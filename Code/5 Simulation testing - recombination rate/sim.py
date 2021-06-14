import msprime as ms
import os
import random
import sys
       
def print_stats(treeseq, filename, Z, seg_sites):
    total_mutations = treeseq.get_num_mutations()
    
    MRCA_time = 0
    for r in treeseq.trees():
        MRCA_time += (r.time(r.root))*8.4
    MRCA_time = MRCA_time/treeseq.num_trees
            
    with open(filename, 'a') as f:
        print(str(Z) + " " + str(treeseq.num_trees) + " " + str(MRCA_time) + " " + str(total_mutations) + " " + str(seg_sites), end='', flush=True, file=f)

# Simulate finite sites data with input recombination rate and mutation map
Z = random.randrange(100000)
with open("mut_map.txt") as f:
    ratemap = [float(i) for i in f.read().splitlines()]
ratemap = [i*0.00002*len(ratemap) for i in ratemap]
pop = [ms.PopulationConfiguration(80, initial_size = 1000000, growth_rate=1.5)]
trees = ms.simulate(length=len(ratemap), random_seed=Z, population_configurations=pop, recombination_rate=float(sys.argv[1]))
model = ms.MatrixMutationModel(alleles=["A","G"], root_distribution=[1.0, 0.0], transition_matrix=[[0.0, 1.0],[1.0, 0.0]])
ratemap = ms.RateMap(position=range(len(ratemap)+1), rate = ratemap)
M = ms.sim_mutations(trees, model=model, rate=ratemap)

# Calculate number of segregating sites
seg_sites = 0
for site in M.sites():
    if(len(site.mutations)>0):
       seg_sites += 1

# Write data to fasta file      
with open("gen_data.fasta", "w") as fa_file:
    M.write_fasta(fa_file)
       
# Add sample stats to results file
print_stats(M, "data_props.txt", Z, seg_sites)
    
# Output a reference sequence (each site in of the ancestral type "A")   
refseq = ["A"]*seg_sites
with open("ref_seq.fasta", 'w') as f:
    print(">Reference", file=f)
    for pos in refseq:
        f.write('%s' % pos)
    print("", file=f)




