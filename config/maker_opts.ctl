genome:/Users/myandell//maker/data/dpp_contig.fasta #genome sequence file.
protein:/Users/myandell/maker/data/dpp_proteins.fasta #protein sequence file.
est:/Users/myandell/maker/data/dpp_transcripts.fasta #EST sequence file.
repeat_protein:/Users/myandell/maker/data/te_proteins.fasta #A transposable Elements Library
clean_up:0 #remove theVoid directory: 1 = yes, 0 = no
rmlib: #An organism specific Repeat Library
use_seq_dir:0 #place files in same directory as sequence file: 1 = yes, 0 = no
split_hit:10000 #length of the splitting of hits
snap_flank:200 #number of bp surrounding SNAP predictions used in attempts to extend gene by MAKER
te_remove:1 #Mask Regions with Excess Similarity to Transposable Element Proteins, 1 = yes, 0 = no 
single_exon:0 #Includes est hits aligning to single exons, 1 = yes, 0 = no
