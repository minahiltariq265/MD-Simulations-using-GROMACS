# MD-Simulations-using-GROMACS on HPC 
# Extract the ligand coordinates from the protein structure file
grep ligand61 protein.pdb > ligand61.pdb
# Ensure GROMACS and related tools are in your environment path
export PATH=$PATH:/apps/gromacs/gromacs-2020.4/build/bin
export PATH=$PATH:/state/partition1/apps/packages/anaconda3/bin
# Prepare the protein topology
gmx pdb2gmx -f protein.pdb -o protein_processed.gro -ter
# Sort bonds in the ligand and prepare its topology
perl sort_mol2_bonds.pl ligand61.mol2 ligand61_fix.mol2
python run_acpype.py -i ligand61.mol2
# In topol.top, include the ligand topology
# ; Include ligand topology
# include "ligand61_GMX.itp"
# Set up the simulation box
gmx editconf -f complex.gro -o newbox.gro -bt cubic -d 1.0
# solvate
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
# Add ions
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
# Energy Minimization
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
# Position Restraints for Ligand
gmx make_ndx -f ligand61_GMX.gro -o index_ligand61.ndx
# Select 0 & ! a H*  S
gmx genrestr -f ligand61_GMX.gro -n index_ligand61.ndx -o posre_ligand61.itp -fc 1000 1000 1000
# Position Restraints for Ligand
gmx make_ndx -f ligand61_GMX.gro -o index_ligand61.ndx
0 & ! a H*  # Non-hydrogen atoms
gmx genrestr -f ligand61_GMX.gro -n index_ligand61.ndx -o posre_ligand61.itp -fc 1000 1000 1000
# Update Topology File for Position Restraints
# In topol.top ADD:
# ; Ligand position restraints
# ifdef POSRES
# include "posre_ligand61.itp"
# endif
# MOVE TOWARDS NVT and NPT Equilibration:
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -deffnm nvt
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -deffnm npt
# Run the production MD simulation:
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr
nohup gmx mdrun -v -deffnm md_0_100 &
# Analysis
gmx rms -s em.tpr -f md_0_10_center.xtc -n index.ndx -tu ns -o rmsd_ligand61.xvg
gmx rmsf -s em.tpr -f md_0_10_center.xtc -n index.ndx -tu ns -o rmsf.xvg -res
# Always verify the integrity of input files after each step.
# Check GROMACS logs (.log files) for errors or warnings during simulations.
