#-----BLAST and Exonerate Statistics Thresholds
blast_type:wublast #set to 'wublast' or 'ncbi'

pcov_blastn:0.8 #Blastn Percent Coverage Threhold EST-Genome Alignments
pid_blastn:0.85 #Blastn Percent Identity Threshold EST-Genome Aligments
eval_blastn:1e-10 #Blastn eval cutoff
bit_blastn:40 #Blastn bit cutoff

pcov_blastx:0.5 #Blastx Percent Coverage Threhold Protein-Genome Alignments
pid_blastx:0.4 #Blastx Percent Identity Threshold Protein-Genome Aligments
eval_blastx:1e-06 #Blastx eval cutoff
bit_blastx:30 #Blastx bit cutoff

pcov_rm_blastx:0.5 #Blastx Percent Coverage Threhold For Transposable Element Masking
pid_rm_blastx:0.4 #Blastx Percent Identity Threshold For Transposbale Element Masking
eval_rm_blastx:1e-06 #Blastx eval cutoff for transposable element masking
bit_rm_blastx:30 #Blastx bit cutoff for transposable element masking

pcov_tblastx:0.8 #tBlastx Percent Coverage Threhold alt-EST-Genome Alignments
pid_tblastx:0.85 #tBlastx Percent Identity Threshold alt-EST-Genome Aligments
eval_tblastx:1e-10 #tBlastx eval cutoff
bit_tblastx:40 #tBlastx bit cutoff

eva_pcov_blastn:0.8 #EVALUATOR Blastn Percent Coverage Threshold EST-Genome Alignments
eva_pid_blastn:0.85 #EVALUATOR Blastn Percent Identity Threshold EST-Genome Alignments
eva_eval_blastn:1e-10 #EVALUATOR Blastn eval cutoff
eva_bit_blastn:40 #EVALUATOR Blastn bit cutoff

ep_score_limit:20 #Exonerate protein percent of maximal score threshold
en_score_limit:20 #Exonerate nucleotide percent of maximal score threshold
