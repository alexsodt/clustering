# clustering
Configurational clustering of ATP in a protein's binding pocket.

# Extract nucleotide configurations aligned to MgtA_ATP.pdb, using psf and dcd:
# make sure there are < 25k configurations to keep things in memory.
extractNucleotides MgtA_ATP.pdb step5_assembly.psf plus_atp.dcd > all_nucleotides.pdb

# get 5 clusters, aligned on atp_align.pdb ... output to clustering.out
cluster 5 atp_align.pdb all_nucleotides > clustering.out

# generates centers.pdb file.
