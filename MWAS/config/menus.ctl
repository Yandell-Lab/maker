##Menus from Data::Dumper
$menus = {
           'altest' => {},
           'gmhmm_e' => {
                          'P. ultimum' => '/home/apache/MWS/users/1/pyu.mod',
                          'M. truncatula' => '/home/cholt/usr/local/gmes/HMM/m_truncatula.mod',
                          'C. remanei' => '/home/cholt/usr/local/gmes/HMM/c_remanei.mod',
                          'C. briggsae' => '/home/cholt/usr/local/gmes/HMM/c_briggsae.mod',
                          'O. sativa' => '/home/cholt/usr/local/gmes/HMM/o_sativa.mod',
                          'A. gambiae' => '/home/cholt/usr/local/gmes/HMM/a_gambiae.mod',
                          'A. thaliana' => '/home/cholt/usr/local/gmes/HMM/a_thaliana.mod',
                          'D. melanogaster' => '/home/cholt/usr/local/gmes/HMM/d_melanogaster.mod',
                          'C. reinhardtii' => '/home/cholt/usr/local/gmes/HMM/c_reinhardtii.mod',
                          'C. elegans' => '/home/cholt/usr/local/gmes/HMM/c_elegans.mod',
                          'C. intestinalis' => '/home/cholt/usr/local/gmes/HMM/c_intestinalis.mod'
                        },
           'est' => {
                      'D. melanogaster : example cDNA' => '/data1/cholt/maker/MWAS/bin/../../data/dpp_transcripts.fasta',
	              'De novo/Legacy Annotaion : example ESTs' => '/home/apache/MWS/users/1/pyu-est.fasta'
                    },
           'genome' => {
                         'D. melanogaster : example contig' => '/data1/cholt/maker/MWAS/bin/../../data/dpp_contig.fasta',
	              'De novo Annotation : example contig' => '/home/apache/MWS/users/1/pyu-contig.fasta',
	              'Legacy Annotation : example contig' => '/home/apache/MWS/users/1/legacy-contig.fasta',
	              'Pass-through : example contig' => '/home/apache/MWS/users/1/pass-contig.fasta'
                       },
           'snaphmm' => {
                          'P. ultimum' => '/home/apache/MWS/users/1/pyu.hmm',
                          'A. thaliana' => '/home/cholt/usr/local/snap/HMM/A.thaliana.hmm',
                          'D. melanogaster' => '/home/cholt/usr/local/snap/HMM/D.melanogaster.hmm',
                          'B. malayi' => '/home/cholt/usr/local/snap/HMM/B.malayi.hmm',
                          'C. elegans' => '/home/cholt/usr/local/snap/HMM/C.elegans.hmm',
                          'O. sativa' => '/home/cholt/usr/local/snap/HMM/O.sativa.hmm',
                          'P. marinus' => '/home/cholt/usr/local/snap/HMM/P.marinus.hmm',
                          'A. gambiae' => '/home/cholt/usr/local/snap/HMM/A.gambiae.hmm',
                          'C. intestinalis' => '/home/cholt/usr/local/snap/HMM/C.intestinalis.hmm',
                          'A. mellifera' => '/home/cholt/usr/local/snap/HMM/A.mellifera.hmm',
                        },
           'protein_gff' => {},
           'repeat_gff' => {},
           'model_org' => {
                            'All species' => 'all'
                          },
           'fgenesh_par_file' => {
                                   'Phytophtora' => '/home/cholt/usr/local/fgenesh/Phytophtora',
                                   'Dicots' => '/home/cholt/usr/local/fgenesh/Dicots'
                                 },
           'augustus_species' => {
                                   'Zea mays' => 'maize',
                                   'Saccharomyces cerevisiae' => 'saccharomyces_cerevisiae_rm11-1a_1',
                                   'Homo sapiens' => 'human',
                                   'Kluyveromyces lactis' => 'kluyveromyces_lactis',
                                   'Pichia stipitis' => 'pichia_stipitis',
                                   'Histoplasma capsulatum' => 'histoplasma_capsulatum',
                                   'Botrytis cinerea' => 'botrytis_cinerea',
                                   'Phanerochaete chrysosporium' => 'phanerochaete_chrysosporium',
                                   'Caenorhabditis elegans' => 'caenorhabditis',
                                   'Magnaporthe grisea' => 'magnaporthe_grisea',
                                   'Rhizopus oryzae' => 'rhizopus_oryzae',
                                   'Tribolium castaneum' => 'tribolium',
                                   'Schistosoma mansoni' => 'schistosoma',
                                   'Aspergillus terreus' => 'aspergillus_terreus',
                                   'Laccaria bicolor' => 'laccaria_bicolor',
                                   'Candida tropicalis' => 'candida_tropicalis',
                                   'Brugia malayi' => 'brugia',
                                   'Cryptococcus neoformans gattii' => 'cryptococcus_neoformans_gattii',
                                   'Drosophila melanogaster' => 'fly',
                                   'Lodderomyces elongisporus' => 'lodderomyces_elongisporus',
                                   'Acyrthosiphon pisum' => 'pea_aphid',
                                   'Encephalitozoon cuniculi' => 'encephalitozoon_cuniculi_GB',
                                   'Candida guilliermondii' => 'candida_guilliermondii',
                                   'Cryptococcus neoformans neoformans' => 'cryptococcus_neoformans_neoformans_JEC21',
                                   'Candida albicans' => 'candida_albicans',
                                   'Yarrowia lipolytica' => 'yarrowia_lipolytica',
                                   'Aspergillus oryzae' => 'aspergillus_oryzae',
                                   'Debaryomyces hansenii' => 'debaryomyces_hansenii',
                                   'Chlamydomonas reinhardtii' => 'chlamydomonas',
                                   'Coprinus cinereus' => 'coprinus_cinereus',
                                   'Aspergillus fumigatus' => 'aspergillus_fumigatus',
                                   'Coccidioides immitis' => 'coccidioides_immitis',
                                   'Solanum lycopersicum' => 'tomato',
                                   'Aspergillus nidulans' => 'aspergillus_nidulans',
                                   'Chaetomium globosum' => 'chaetomium_globosum',
                                   'Ustilago maydis' => 'ustilago_maydis',
                                   'Galdieria sulphuraria' => 'galdieria',
                                   'Aedes aegypti' => 'aedes',
                                   'Arabidopsis thaliana' => 'arabidopsis',
                                   'Schizosaccharomyces pombe' => 'schizosaccharomyces_pombe',
                                   'Toxoplasma gondii' => 'toxoplasma',
                                   'Tetrahymena thermophila' => 'tetrahymena',
                                   'Amphimedon queenslandica' => 'amphimedon',
                                   'Neurospora crassa' => 'neurospora_crassa',
                                   'Eremothecium gossypii' => 'eremothecium_gossypii',
                                   'Nasonia vitripennis' => 'nasonia',
                                   'Fusarium graminearum' => 'fusarium_graminearum'
                                 },
           'gmhmm_p' => {},
           'model_gff' => {
			     'Legacy Annotation : example model set 1' => '/home/apache/MWS/users/1/legacy-set1.gff',
			     'Legacy Annotation : example model set 2' => '/home/apache/MWS/users/1/legacy-set2.gff'
	                  },
           'repeat_protein' => {
                                 'RepeatRunner te_proteins' => '/data1/cholt/maker/MWAS/bin/../../data/te_proteins.fasta'
                               },
           'rmlib' => {},
           'pref_gff' => {},
           'protein' => {
                          'D. melanogaster : example proteins' => '/data1/cholt/maker/MWAS/bin/../../data/dpp_proteins.fasta',
		              'De novo/Legacy Annotation : example Protein' => '/home/apache/MWS/users/1/pyu-protein.fasta'
                        },
           'altest_gff' => {},
           'est_gff' => {
    		          'Pass-through : example mRNAseq' => '/home/apache/MWS/users/1/pass-mRNAseq.gff'
	                }
         };
