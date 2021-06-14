#!/bin/bash

echo "Time Rec_rate Seed Trees MRCA Muts Seg_sites Pmin" > data_props.txt

for rec_rate in 0.0000001 0.0000002 0.0000004 0.0000008 0.0000016 0.0000032 0.0000064 0.0000128
do	
	for i in {1..200}
	do
		run_id=$(date -d "today" +"%d%H%M%S")
		echo -n "${run_id} ${rec_rate} " >> data_props.txt
		echo -n "${run_id} ${rec_rate} $i ..."
		python3 sim.py ${rec_rate}
		cat ref_seq.fasta gen_data.fasta > sample_.fasta
		snp-sites -m -o sample.fasta sample_.fasta
		timeout 300 kwarg -Q100 -S1,0.01 -M1.1,0.02 -R-1,1 -C-1,2 -k -n -f < sample.fasta > kwarg_out.txt
		Rscript kwarg_count.R kwarg_out.txt >> data_props.txt
		echo "done"
	done
done


