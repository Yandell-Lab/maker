#-----Genome (Required for De-Novo Annotation)
genome: #genome sequence file in fasta format

#-----Re-annotation Options (Only MAKER derived GFF3)
genome_gff: #re-annotate genome based on this gff3 file
est_pass:0 #use ests in genome_gff: 1 = yes, 0 = no
altest_pass:0 #use alternate organism ests in genome_gff: 1 = yes, 0 = no
protein_pass:0 #use proteins in genome_gff: 1 = yes, 0 = no
rm_pass:0 #use repeats in genome_gff: 1 = yes, 0 = no
model_pass:0 #use gene models in genome_gff: 1 = yes, 0 = no
pred_pass:0 #use ab-initio predictions in genome_gff: 1 = yes, 0 = no
other_pass:0 #passthrough everything else in genome_gff: 1 = yes, 0 = no

#-----EST Evidence (you should provide a value for at least one)
est: #non-redundant set of assembled ESTs in fasta format (classic EST analysis)
est_reads: #unassembled nextgen mRNASeq in fasta format (not fully implemented)
altest: #EST/cDNA sequence file in fasta format from an alternate organism
est_gff: #EST evidence from an external gff3 file
altest_gff: #Alternate organism EST evidence from a seperate gff3 file

#-----Protein Homology Evidence (you should provide a value for at least one)
protein:  #protein sequence file in fasta format
protein_gff:  #protein homology evidence from an external gff3 file

#-----Repeat Masking (leave values blank to skip)
model_org:all=STATIC #model organism for RepBase masking in RepeatMasker
repeat_protein:/home/cholt/usr/local/maker/data/te_proteins.fasta=STATIC #a database of transposable element proteins in fasta format
rmlib: #an organism specific repeat library in fasta format
rm_gff: #repeat elements from an external gff3 file

#-----Gene Prediction Options
organism_type:eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic
predictor: #prediction methods for annotations (seperate multiple values by ',')
unmask:0 #Also run ab-initio methods on unmasked sequence, 1 = yes, 0 = no
snaphmm: #SNAP HMM model
gmhmm: #GeneMark HMM model
augustus_species: #Augustus gene prediction model
fgenesh_par_file: #Fgenesh parameter file
model_gff: #gene models from an external gff3 file (annotation pass-through)
pred_gff: #ab-initio predictions from an external gff3 file

#-----Other Annotation Type Options (features maker doesn't recognize)
other_gff: #features to pass-through to final output from an extenal gff3 file

#-----External Application Specific Options
alt_peptide:C #amino acid used to replace non standard amino acids in blast databases
cpus:1=DISABLED #max number of cpus to use in BLAST and RepeatMasker

#-----MAKER Specific Options
evaluate:0 #run EVALUATOR on all annotations, 1 = yes, 0 = no
max_dna_len:100000 #length for dividing up contigs into chunks (larger values increase memory usage)
min_contig:1 #all contigs from the input genome file below this size will be skipped
min_protein:0 #all gene annotations must produce a protein of at least this many amino acids in length
AED_threshold:1 #Maximum Annotation Edit Distance allowed for annotations (bound by 0 and 1)
softmask:1 #use soft-masked rather than hard-masked seg filtering for wublast
split_hit:10000 #length for the splitting of hits (expected max intron size for evidence alignments)
pred_flank:200 #length of sequence surrounding EST and protein evidence used to extend gene predictions
single_exon:0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length:250 #min length required for single exon ESTs if 'single_exon is enabled'
keep_preds:0 #Add non-overlapping ab-inito gene prediction to final annotation set, 1 = yes, 0 = no
map_forward:0 #try to map names and attributes forward from gff3 annotations, 1 = yes, 0 = no
retry:1 #number of times to retry a contig if there is a failure for some reason
clean_try:0 #removeall data from previous run before retrying, 1 = yes, 0 = no
clean_up:1 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP:=DISABLED #specify a directory other than the system default temporary directory for temporary files

#-----EVALUATOR Control Options
side_thre:5
eva_window_size:70
eva_split_hit:1
eva_hspmax:100
eva_gspmax:100
enable_fathom:0
