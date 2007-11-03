max_dna_len:100000 #the max DNA length used for sequence similarity comparisons.
percov_blastn:0.80 #Blastn Percent Coverage Threhold EST-Genome Alignments
percid_blastn:0.85 #Blastn Percent Identity Threshold EST-Genome Aligments
eval_blastn:1e-10 #Blastn eval cutoff
bit_blastn:40 #Blastn bit cutoff
percov_blastx:0.50 #Blastx Percent Coverage Threhold Protein-Genome Alignments
percid_blastx:0.40 #Blastx Percent Identity Threshold Protein-Genome Aligments
eval_blastx:1e-5 #Blastx eval cutoff
bit_blastx:30 #Blastx bit cutoff
e_perc_cov:50 #Exonerate Percent Coverage Thresshold EST_Genome Alignments
ep_score_limit:10 #Report  alignments scoring at least this percentage of the maximal score exonerate nucleotide
en_score_limit:10 #Report  alignments scoring at least this percentage of the maximal score exonerate protein
model_org:drosophila #Model Organism use for finding repeat sequences by RepeatMasker
snaphmm:/Users/bcantarel/programs/snap/HMM/Dm.hmm #SNAP HMM Model - This can be organism specific by iterating MAKER see HMM BUILDING
