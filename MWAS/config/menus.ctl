##Menus from Data::Dumper
$menus = {
           'repeat_gff' => {},
           'genome' => {
                         'D. melanogaster : example contig' => '/home/cholt/maker/MWAS/bin/../../data/dpp_contig.fasta'
                       },
           'altest' => {},
           'snaphmm' => {
                          'A. mellifera' => '/home/cholt/usr/local/snap/HMM/A.mellifera.hmm',
                          'C. elegans' => '/home/cholt/usr/local/snap/HMM/C.elegans.hmm',
                          'C. intestinalis' => '/home/cholt/usr/local/snap/HMM/C.intestinalis.hmm',
                          'B. malayi' => '/home/cholt/usr/local/snap/HMM/B.malayi.hmm',
                          'A. gambiae' => '/home/cholt/usr/local/snap/HMM/A.gambiae.hmm',
                          'A. thaliana' => '/home/cholt/usr/local/snap/HMM/A.thaliana.hmm',
                          'D. melanogaster' => '/home/cholt/usr/local/snap/HMM/D.melanogaster.hmm',
                          'O. sativa' => '/home/cholt/usr/local/snap/HMM/O.sativa.hmm',
                          'P. marinus' => '/home/cholt/usr/local/snap/HMM/P.marinus.hmm'
                        },
           'rmlib' => {},
           'protein_gff' => {},
           'model_org' => {
                            'All species' => 'all'
                          },
           'pref_gff' => {},
           'fgenesh_par_file' => {},
           'gmhmm_E' => {
                          'C. elegans' => '/home/cholt/usr/local/gmes/HMM/c_elegans.mod',
                          'C. intestinalis' => '/home/cholt/usr/local/gmes/HMM/c_intestinalis.mod',
                          'C. remanei' => '/home/cholt/usr/local/gmes/HMM/c_remanei.mod',
                          'A. gambiae' => '/home/cholt/usr/local/gmes/HMM/a_gambiae.mod',
                          'M. truncatula' => '/home/cholt/usr/local/gmes/HMM/m_truncatula.mod',
                          'A. thaliana' => '/home/cholt/usr/local/gmes/HMM/a_thaliana.mod',
                          'C. reinhardtii' => '/home/cholt/usr/local/gmes/HMM/c_reinhardtii.mod',
                          'C. briggsae' => '/home/cholt/usr/local/gmes/HMM/c_briggsae.mod',
                          'D. melanogaster' => '/home/cholt/usr/local/gmes/HMM/d_melanogaster.mod',
                          'O. sativa' => '/home/cholt/usr/local/gmes/HMM/o_sativa.mod'
                        },
           'protein' => {
                          'D. melanogaster : example proteins' => '/home/cholt/maker/MWAS/bin/../../data/dpp_proteins.fasta',
                          'UniProt' => '/home/cholt/maker/data/uniprot_swissprot.prot'
                        },
           'est' => {
                      'D. melanogaster : example cDNA' => '/home/cholt/maker/MWAS/bin/../../data/dpp_transcripts.fasta'
                    },
           'gmhmm_P' => {},
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
           'model_gff' => {},
           'altest_gff' => {},
           'repeat_protein' => {
                                 'RepeatRunner te_proteins' => '/home/cholt/maker/MWAS/bin/../../data/te_proteins.fasta'
                               },
           'est_gff' => {}
         };
