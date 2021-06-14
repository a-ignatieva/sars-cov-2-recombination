#!/bin/bash

echo "Time j f Seed Trees MRCA Muts Seg_sites Pmin" > data_props.txt

for j in 0 # This is the number of adjusted sites
do	
	for f in 0 # This is the factor H
	do
		for i in {1..200}
		do
			run_id=$(date -d "today" +"%d%H%M%S")
			echo -n "${run_id} ${j} ${f} " >> data_props.txt
			echo "${run_id} ${i} ... "
			python3 sim.py 0.0 mut_map_adjusted_p_${j}f_${f}.txt
			cat ref_seq.fasta gen_data.fasta > sample_.fasta
			snp-sites -m -o sample.fasta sample_.fasta
			timeout 300 kwarg -Q100 -S1,0.01 -M1.1,0.02 -R-1,1 -C-1,2 -k -n -f < sample.fasta > kwarg_out.txt
			Rscript kwarg_count.R kwarg_out.txt >> data_props.txt
			echo "done"
		done
	done
done

for j in 1 5 10 20 50 100 200 # This is the number of adjusted sites
do	
	for f in 2 5 10 20 50 # This is the factor H
	do
		for i in {1..200}
		do
			run_id=$(date -d "today" +"%d%H%M%S")
			echo -n "${run_id} ${j} ${f} " >> data_props.txt
			echo "${run_id} ${i} ... "
			python3 sim.py 0.0 mut_map_adjusted_p_${j}f_${f}.txt
			cat ref_seq.fasta gen_data.fasta > sample_.fasta
			snp-sites -m -o sample.fasta sample_.fasta
			timeout 300 kwarg -Q200 -S1,0.01 -M1.1,0.02 -R-1,1 -C-1,2 -k -n -f < sample.fasta > kwarg_out.txt
			Rscript kwarg_count.R kwarg_out.txt >> data_props.txt
			echo "done"
		done
	done
done