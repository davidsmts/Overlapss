python3 src/main.py x --opt_get_reads ../data/Pseudoalteromonas_translucida_KMM_520_genomic.fasta --opt_reads_into reads_withreps.fasta --opt_read_len 250 --opt_read_cov 10 # get data and write reads and the full genome into two different files
python3 src/main.py reads_withreps.fasta --remove_repeats 
penguin nuclassemble data/wout_reps/correct.fasta script_assembly.fas tmp --min-contig-len 10000
python3 ../../other/quast-5.3.0/quast.py -r ../../code/Overlapss/long_chromosome.fasta ../../code/Overlapss/assembly_long_chromosome.fas -o quast_script/