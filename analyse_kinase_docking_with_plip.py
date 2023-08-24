import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbqt as pdbqt
from meeko import RDKitMolCreate
from meeko import PDBQTMolecule
from plip.structure.preparation import PDBComplex
import numpy as np
import pandas as pd
from opencadd.databases import klifs
import argparse
import os 
import biotite.database.rcsb as rcsb
from rdkit.Chem import rdFMCS,AllChem
from rdkit import Chem

def find_mcs(crystal_ligand,docked_ligand):
    mcs = rdFMCS.FindMCS([crystal_ligand,docked_ligand])
    return mcs

def get_max_common_substructure_rmsd(crystal_ligand,docked_pose,mcs,poseid):
    mcs_crystal_ligand = Chem.MolFromSmarts(mcs.smartsString)
    match = crystal_ligand.GetSubstructMatch(mcs_crystal_ligand)
    conf = Chem.Conformer(mcs_crystal_ligand.GetNumAtoms())
    for i,mi in enumerate(match):
        conf.SetAtomPosition(i,crystal_ligand.GetConformer().GetAtomPosition(mi))
    mcs_crystal_ligand.AddConformer(conf)
    
    mcs_docked_pose = Chem.MolFromSmarts(mcs.smartsString)
    match = docked_pose.GetSubstructMatch(mcs_docked_pose)
    conf = Chem.Conformer(mcs_docked_pose.GetNumAtoms())
    for i,mi in enumerate(match):
        conf.SetAtomPosition(i,docked_pose.GetConformer(poseid).GetAtomPosition(mi))
    mcs_docked_pose.AddConformer(conf)
    return AllChem.AlignMol(mcs_crystal_ligand,mcs_docked_pose)

def sanitize_pdb(filename='',addHs_pH=7.4,output='',try_renumberResidues=False):
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    fix = PDBFixer(filename=filename)
    fix.findMissingResidues()
    fix.findNonstandardResidues()
    fix.replaceNonstandardResidues()
    fix.findMissingAtoms()
    fix.addMissingAtoms()
    fix.addMissingHydrogens(addHs_pH)
    PDBFile.writeFile(fix.topology, fix.positions, open(filename, 'w'))

def get_active_site_dict(pdbid):
    file = rcsb.fetch(pdbid,"pdb")
    structure = pdb.PDBFile.read(file).get_structure()[0]
    chainA = structure[structure.chain_id==np.unique(structure.chain_id)[0]]
    protein_chainA = chainA[chainA.hetero==False]
    protein_chainA_resid = pd.unique(protein_chainA.res_id)
    protein_chainA_resname = {}
    for i,resid in enumerate(protein_chainA_resid):
        protein_chainA_resname[resid] = protein_chainA[protein_chainA.res_id==resid].res_name[0]
    
    entries = klif_local.structures.by_structure_pdb_id(pdbid)
    binding_site = klif_local.pockets.by_structure_klifs_id(
        entries['structure.klifs_id'][0]
    )
    
    binding_site_dict = {}
    cnt = 0
    for i,val in enumerate(binding_site['residue.id']):
        binding_site_dict[val+protein_chainA_resname[int(val)]+'_bb'] = cnt
        cnt += 1
        binding_site_dict[val+protein_chainA_resname[int(val)]+'_sc'] = cnt
        cnt += 1    

    return binding_site_dict,protein_chainA_resname

def convert_plip_to_ifp(plip_interaction_object,ligand_map,binding_site_dict,protein_chainA_resname,output_prefix):
    # anatomy of an interaction fingerprint[i,j,k]
    # i - interaction type
    # (0: hydrophobic, 1: hbonds_ldon, 2: hbonds_pdon, 
    # 3: pication_laro, 4: halogen_bonds, 5: salt_bridge_lneg,
    # 6: salt_bridge_pneg, 7: pistacking, 8: water bridge) 
    # j - ligand atom number 
    # k - 85 active site residues * 2 (back bone and side chain) 
    ifp_matrix = np.zeros((9,len(plip_interaction_object.ligand.all_atoms),85*2))
    # here all column names have been made 
    # now make row names, ligand atom+pdbatomnumber 
    ligand_row_names = []
    for i in range(len(plip_interaction_object.ligand.all_atoms)):
        ligand_row_names.append(f"{ligand_map[i+1]}")

    
    # hydrophobic interactions 
    for interaction_id, interaction in enumerate(plip_interaction_object.hydrophobic_contacts):
        ligatomid = interaction.ligatom.idx-1
        prot_res_sc_or_bb = str(interaction.resnr)
        prot_res_sc_or_bb += protein_chainA_resname[interaction.resnr]
        prot_res_sc_or_bb += "_sc"
        try:
            ifp_matrix[0,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1
        except KeyError:
            print(f"{prot_res_sc_or_bb} not in active site!")
            continue
     
    # Hydrogen bonds with Ligand donors 
    for interaction_id, interaction in enumerate(plip_interaction_object.all_hbonds_ldon):
        ligatomid = interaction.d.idx-1
        prot_res_sc_or_bb = str(interaction.resnr)
        prot_res_sc_or_bb += protein_chainA_resname[interaction.resnr]
        if interaction.sidechain:
            prot_res_sc_or_bb += "_sc"
        else:
            prot_res_sc_or_bb += "_bb"
        try:
            ifp_matrix[0,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1
        except KeyError:
            print(f"{prot_res_sc_or_bb} not in active site!")
            continue
        ifp_matrix[1,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1
    
    # Hydrogen bonds with protein donors 
    for interaction_id, interaction in enumerate(plip_interaction_object.all_hbonds_pdon):
        ligatomid = interaction.a.idx-1
        prot_res_sc_or_bb = str(interaction.resnr)
        prot_res_sc_or_bb += protein_chainA_resname[interaction.resnr]
        if interaction.sidechain:
            prot_res_sc_or_bb += "_sc"
        else:
            prot_res_sc_or_bb += "_bb"
        try:
            ifp_matrix[0,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1
        except KeyError:
            print(f"{prot_res_sc_or_bb} not in active site!")
            continue
        ifp_matrix[2,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1
    
    # pication with ligand providing pi group 
    for interaction_id, interaction in enumerate(plip_interaction_object.all_pi_cation_laro):
        ligatomidx = []
        for atom in interaction.ring.atoms:
            ligatomidx.append(atom.idx-1)
        prot_res_sc_or_bb = str(interaction.resnr)
        prot_res_sc_or_bb += protein_chainA_resname[interaction.resnr]
        prot_res_sc_or_bb += "_sc"
        try:
            ifp_matrix[0,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1
        except KeyError:
            print(f"{prot_res_sc_or_bb} not in active site!")
            continue
        for ligatomid in ligatomidx:
            ifp_matrix[3,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1
    
    # halogen bonds
    for interaction_id, interaction in enumerate(plip_interaction_object.halogen_bonds):
        ligatomid = interaction.don.x.idx-1
        prot_res_sc_or_bb = str(interaction.resnr)
        prot_res_sc_or_bb += protein_chainA_resname[interaction.resnr]
        if interaction.sidechain:
            prot_res_sc_or_bb += "_sc"
        else:
            prot_res_sc_or_bb += "_bb"
        try:
            ifp_matrix[0,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1
        except KeyError:
            print(f"{prot_res_sc_or_bb} not in active site!")
            continue
        ifp_matrix[4,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1
    
    # salt bridge with negative ligand atom 
    for interaction_id, interaction in enumerate(plip_interaction_object.saltbridge_lneg):
        ligatomidx = []
        for atom in interaction.negative.atoms:
            ligatomidx.append(atom.idx-1)
        prot_res_sc_or_bb = str(interaction.resnr) + protein_chainA_resname[interaction.resnr]+ "_sc"
        try:
            ifp_matrix[0,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1
        except KeyError:
            print(f"{prot_res_sc_or_bb} not in active site!")
            continue
        ifp_matrix[5,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1
    
    # salt bridge with negative protein atom 
    for interaction_id, interaction in enumerate(plip_interaction_object.saltbridge_pneg):
        ligatomidx = []
        for atom in interaction.positive.atoms:
            ligatomidx.append(atom.idx-1)
        prot_res_sc_or_bb = str(interaction.resnr) + protein_chainA_resname[interaction.resnr]+ "_sc"
        try:
            ifp_matrix[0,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1
        except KeyError:
            print(f"{prot_res_sc_or_bb} not in active site!")
            continue
        ifp_matrix[6,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1
    
    # pi stacking
    for interaction_id, interaction in enumerate(plip_interaction_object.pistacking):
        ligatomidx = []
        for atom in interaction.ligandring.atoms:
            ligatomidx.append(atom.idx-1)
        prot_res_sc_or_bb = str(interaction.resnr) + protein_chainA_resname[interaction.resnr]+  '_sc'
        try:
            ifp_matrix[0,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1
        except KeyError:
            print(f"{prot_res_sc_or_bb} not in active site!")
            continue
        ifp_matrix[7,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1
    
    # water bridge 
    for interaction_id, interaction in enumerate(plip_interaction_object.pistacking):
        if interaction.protisdon:
            ligatomidx = interaction.a.idx-1
        else:
            ligatomidx = interaction.d.idx-1
        prot_res_sc_or_bb = str(interaction.resnr) + protein_chainA_resname[interaction.resnr]+ '_sc'
        try:
            ifp_matrix[0,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1
        except KeyError:
            print(f"{prot_res_sc_or_bb} not in active site!")
            continue
        ifp_matrix[8,ligatomid,binding_site_dict[prot_res_sc_or_bb]] = 1   
    
    interaction_type_dict = {
        'hydrophobic_contacts':1,
        'hbonds_ldon':2,
        'hbonds_pdon':3,
        'pication_laro':4,
        'halogen_bonds':5,
        'salt_bridge_lneg':6,
        'salt_bridge_pneg':7,
        'pistacking':8,
        'waterbridge':9
        }
    
    for i in range(9):
        df = pd.DataFrame(ifp_matrix[i,:,:],)
        df.columns = list(binding_site_dict.keys())
        df.index = ligand_row_names
        df.to_csv(f"{output_prefix}/{list(interaction_type_dict.keys())[i]}.csv")
    
    sum_ifp = np.sum(ifp_matrix,axis=0)
    df = pd.DataFrame(sum_ifp)
    df.columns = list(binding_site_dict.keys())
    df.index = ligand_row_names
    df.to_csv(f'{output_prefix}/overall.csv')
    return sum_ifp

#def get_active_site_dict(pdbid):
#    klif_local = klifs.setup_local

def get_ifp_from_pdb(pdbid,binding_site_dict,protein_chainA_resname):
    file = rcsb.fetch(pdbid,"pdb",'.')
    mycomplex = PDBComplex()
    mycomplex.load_pdb(f"{pdbid}.pdb")
    # Make PLIP do it's magic
    mycomplex.analyze()
    interaction_set_name = list(mycomplex.interaction_sets.keys())[0]
    interactions = mycomplex.interaction_sets[interaction_set_name]
    ligand_map = mycomplex.interaction_sets[interaction_set_name].Mapper.ligandmaps[interaction_set_name]
    ifp_matrix = convert_plip_to_ifp(interactions,ligand_map,binding_site_dict,protein_chainA_resname,f"{pdbid}_crystal_structure_interactions")
    return ifp_matrix

def compare_crystal_and_docked_pose_interactions(crystal_ifp,docked_ifp):
    crystal_ifp = np.sum(np.array(crystal_ifp),axis=0)
    docked_ifp = np.sum(np.array(docked_ifp),axis=0)
    score = 0
    for i,val in enumerate(crystal_ifp):
        if val>0:
            if docked_ifp[i]>val:
                score += val
            else:
                score += docked_ifp[i]
    return score,np.sum(crystal_ifp)


    
parser = argparse.ArgumentParser(
    prog="Code to calculate an interaction fingerprint\
        from docked poses of a Kinase with a ligand.",
    description="The code takes a PDB ID as input, the kinase \
        structure from which was used for a docking study. It \
        uses the KLIFS database to get the 85 active site     \
        residues and then uses the PLIP tool to calculate     \
        interactions with the ligand."
)

parser.add_argument('-pdbid',
                    '--PDBID',
                    required=True,
                    type=str)
parser.add_argument('-receptor',
                    '--receptor_pdbqt',
                    required=True,
                    type=str
                    )
parser.add_argument('-poses',
                    '--pose_pdbqt',
                    required=True,
                    type=str)
parser.add_argument('-output',
                    '--output_prefix',
                    required=False,
                    default="docked_complex")
parser.add_argument('-ligand_sdf',
                    '--crystal_ligand_sdf',
                    required=True)
args = parser.parse_args()

# setup local version of KLIFS database
klif_local = klifs.setup_remote()
binding_site_dict,prot_chainA_resname = get_active_site_dict(args.PDBID)

if not os.path.isdir(f"{args.PDBID}_crystal_structure_interactions"):
        os.mkdir(f"{args.PDBID}_crystal_structure_interactions")
crystal_interactions = get_ifp_from_pdb(args.PDBID,binding_site_dict,prot_chainA_resname)

# Read the receptor pdbqt and 
receptor = pdbqt.PDBQTFile.read(args.receptor_pdbqt)
protein_structure = receptor.get_structure()[0]
ligand_poses = pdbqt.PDBQTFile.read(args.pose_pdbqt)
pdbqt_ligand = PDBQTMolecule.from_file(args.pose_pdbqt,skip_typing=True)
rdkitmol_list = RDKitMolCreate.from_pdbqt_mol(pdbqt_ligand)
cocrystal_ligand_sdf = Chem.SDMolSupplier(args.crystal_ligand_sdf)[0]
mcs = find_mcs(cocrystal_ligand_sdf,rdkitmol_list[0])


# create PLIP object to load each docked pose after it is
# converted into a protein ligand complex pdb 
mycomplex = PDBComplex()
interactions = []
docked_pose_score_matrix = np.zeros((len(ligand_poses.get_structure()),4))
for poseid,pose_structure in enumerate(ligand_poses.get_structure()):
    print (f"RUNNING FINGERPRINTING FOR POSE {poseid+1}...")
    output_folder_name = f"{args.output_prefix}_pose_{poseid+1}_interactions"
    if not os.path.isdir(output_folder_name):
        os.mkdir(output_folder_name)
    pose_structure.hetero = [True for i in range(len(pose_structure))]
    combined = protein_structure + pose_structure
    file = pdb.PDBFile()
    file.set_structure(combined)
    filename = output_folder_name+f"/poseid{poseid+1}_complex.pdb"
    file.write(filename)
    # Sanitize this pdb, if you so choose
    # sanitize_pdb(filename)
    # load this complex pdb into PLIP object
    mycomplex.load_pdb(filename)
    # Make PLIP do it's magic
    mycomplex.analyze()

    ifp = convert_plip_to_ifp(
        mycomplex.interaction_sets[f'{pose_structure.res_name[0]}:Z:1'],
        mycomplex.interaction_sets[f'{pose_structure.res_name[0]}:Z:1'].Mapper.ligandmaps[f'{pose_structure.res_name[0]}:Z:1'],
        binding_site_dict,
        prot_chainA_resname,
        output_folder_name
        )
    
    rmsd = get_max_common_substructure_rmsd(cocrystal_ligand_sdf,rdkitmol_list[0],mcs,poseid)
    dock_score,total = compare_crystal_and_docked_pose_interactions(crystal_interactions,ifp)
    docked_pose_score_matrix[poseid,:] = [poseid+1,dock_score,total,rmsd]

sorted_indices = np.lexsort((docked_pose_score_matrix[:,-1],-docked_pose_score_matrix[:,1]))
sorted_matrix = docked_pose_score_matrix[sorted_indices]
np.savetxt("overall_interactions_rmsd_sorting.csv",sorted_matrix,delimiter=',',header="Pose ID,# of Pose Interactions,# of Crystal Interactions,MCS RMSD")