# Kinase inhibitor pose scorer
Codes to calculation interaction fingerprints of Kinase inhibitors from docked poses and MD simulations.

For docked poses, an example use case would be, 
```
python analyse_kinase_docking_with_plip.py -pdbid 6rsh -receptor 6rsh_receptor_hadd.pdbqt -poses out.pdbqt -ligand_sdf 6rsh_ligand.sdf
```
where one has docked a ligand of choice to the protein chain from PDB ID 6RSH (script assumes it's chain A), and docking was run using the receptor file ```"6rsh_receptor_hadd.pdbqt"``` and the docked poses were saved in ```"out.pdbqt"```. The co crystallized ligand is also provided in the ```"6rsh_ligand.sdf"``` file. 

The script outputs folders named "docked_complex_pose_XX_interactions" which contain interactions of pose ID XX with the receptor split into different types of interactions. Each file contains rows equivalent to the number of atoms in the ligand and columns equivalent to the Kinase active site residues (backbone and side chain) which is taken fro m the KLIFS database. These interactions are then summed up in the overall.csv file. The docked pose is also saved in a pdb file within each "docked_complex_pose_XX_interactions" directory. 

To rank the docked poses, the script also outputs another file called "overall_interactions_rmsd_sorting.csv" which contains the ids of the poses sorted according to first, the number of interactions a pose captures compared to the cocrystallized ligand and then according to the RMSD of the maximum common substructure between the docked and cocrystallized ligand. A combination of both of these metrics could be used to decide which one of the top docked poses is close to the actual pose. Keep in mind that with very small maximum common substructures or tanimoto similarity between the docked and cocrystallized ligand, this method of ranking docked poses will lose its effectiveness. 
