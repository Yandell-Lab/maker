#------------------------------------------------------------------------
#----                          maker::sens_spec                      ----
#------------------------------------------------------------------------
package maker::sens_spec;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use FileHandle;
use PostData;
use Exporter;
use PhatHit_utils;
@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub load_anno_hsps {
        my $annotations = shift;

        my @coors;
        foreach my $g (@$annotations){
                foreach my $t_struct (@{$g->{t_structs}}){
                        my $hit = $t_struct->{hit};
                        foreach my $hsp ($hit->hsps()){
                                push(@coors, [$hsp->nB('query'),
                                              $hsp->nE('query'),
                                             ]);
                        }
                }
        }
        return \@coors;
}
#-----------------------------------------------------------------------------
sub load_snap_hps {
        my $snaps = shift;

        my @coors;
        foreach my $hit (@{$snaps}){
                foreach my $hsp ($hit->hsps()){
                        push(@coors, [$hsp->nB('query'),
                                      $hsp->nE('query'),
                                    ]);
                }

        }
        return \@coors;
}
#-----------------------------------------------------------------------------
sub spew_sens_spec_template {
        my $type  = shift;
        my $fasta = shift;
        my $data  = shift;

        my $def   = Fasta::getDef($fasta);
        my $seq   = Fasta::getSeq($fasta);

        my ($seq_id)  = $def =~ /^>(\S+)/;

        my $coors = $type eq 'a' ? load_anno_hsps($data)
                                 : load_snap_hps($data);

        my $masked_seq = Shadower::maskSequence($seq, $coors, 0, '1');

           $$masked_seq =~ s/[^1]/0/g;

        my $m_fasta = Fasta::toFasta($def." Maker sens_spec_template",
                                     $masked_seq,
                                    );


        return $$m_fasta;
}
#------------------------------------------------------------------------
1;


