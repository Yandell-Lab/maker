#<top> ------------------------------------------------------------------------
#     <Name>       Game.pm                    </Name> 
#     <class>      Dumper::
#                  ::XML::Game             </class>       
#     <Author>     M Hadi Islam            </Author>
#     <email>      hadi@genetics.utah.edu  </email>
#     <Does>       Prints xml dump of 
#                  CGL::Annotation Object  </Does>
#</top> ----------------------------------------------------------------------

package Dumper::XML::Game;

use Class::Struct;
use Data::Dumper;
use CGL::Annotation;
use Dumper::XML::Game_Xml;
use strict "vars";
use strict "refs";
use SpaceBase;
use Fasta;
use FastaChunker;
use Iterator::Fasta;
#------------- Requires-------------------------------------------

#<classdef>--------------------------------------------------------------------
struct Dumper::XML::Game =>
{
    Annot          => '$',
    ChaosFile      => '$',
    FastaFile      => '$',
    Iter           => '$',
    Hits           => '$',
    XmlP           => '$',
    AutoAnts       => '$',
    Contig         => '$',
    WriteFile      => '$',
};



sub init{
    my ($self,%args)=@_;
    %{$self}=%args;
    my $gx= new Game_Xml();
  
    $self->XmlP($gx);
    my @a;
    $self->Hits(\@a);
     
    
    return $self;
}

#</classdef>--------------------------------------------------------------------

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#---------------------- SUBS ---------------------------------------------------
#-------------------------------------------------------------------------------


#<sub>--------------------------------------------------------------------------
sub Game{
    my ($self)=@_;
    
    my $game = new Game::game();
    
    my $seq  = $self->Seq();
   
    my $an   = $self->Annotation();
    #print "system out\n";die;
    my $compA=$self->Computational_Analysis();
  

    $game->seq($seq);
    $game->Annotation($an);
  #  print "system out\n";die;
    $game->computational_analysis($compA);
   
   

    $self->XmlP->out($self->{WriteFile});
    $self->XmlP->print_header();
    $self->XmlP->treat($game);
    $self->Write($self->XmlP->sting());
   
}
#</sub>--------------------------------------------------------------------

#<sub>------------------------------------------------------------------------ 
sub Seq{
    my ($self)=@_;
    
    my $annot=$self->{Annot};
    my $ff=$self->{FastaFile};
    my $seq;

    if(defined($annot) && defined($ff)){
	print STDERR "((%): Ambiguous residues\n";
	die "Both CGL and Fasta found\n";
    }
    elsif(defined($annot)){
	$seq=$self->Get_Seq_from_CGL($annot);
    }elsif(defined($ff)){
	$seq=$self->Get_Seq_from_Fasta($ff);
    }else{
	print STDERR "((%): No Residues\n";
	die "Must Provide\nCGL::Annotation || fasta file;\n";
    }
    return $seq;
}
#</sub>-------------------------------------------------------------------- 

#<sub>------------------------------------------------------------------------ 
sub Get_Seq_from_CGL{
    my ($self,$annot)=@_;
    
    my $con=$annot->contig(0);
    my @feats=$annot->features();
    
    my $name=$con->name;
    my $focus ="false";
    my $orga;
    
    if($feats[0]){
	$orga=$feats[0]->organismstr;
    }
    if(defined($con)){
	$focus="true";
    }

    my $res=$con->residues();
    my $type =$con->type();
    my $id = $con->id();

    my $seq_attr=new Game::attr::seq(type  =>$type,
				     id    =>$id,
				     focus =>$focus,
	);
    
    my $seq=new Game::seq(name     => $name,
			  organism => $orga,
			  residues => $res
	);
    
    $seq->attr($seq_attr);
    return $seq;
}
#</sub>-------------------------------------------------------------------- 

#<sub>------------------------------------------------------------------------ 
sub Get_Seq_from_Fasta{
    my ($self,$file)=@_;
    
    my $focus="true";
    my $orga;

    my $fasta_iterator = new Iterator::Fasta($file);
    my $fasta     = $fasta_iterator->nextEntry();
    my $query_def = Fasta::getDef($fasta);
    my ($id)  = $query_def =~ /^>(\S+)/;
    my $res=Fasta::getSeq($fasta);;
    my $type ="FASTA";
    my $name=$id;
    

    my $seq_attr=new Game::attr::seq(type  =>$type,
				     id    =>$id,
				     focus =>$focus,
	);
    
    my $seq=new Game::seq(name     => $name,
			  organism => $orga,
			  residues => $$res
	);
    
    $seq->attr($seq_attr);
    return $seq;
}
#</sub>-------------------------------------------------------------------- 
#<sub>------------------------------------------------------------------------ 
sub Get_Seq_from_Autoant{
    my ($self,$a)=@_;
    
    my @ret;
    my $p_seq_attr=new Game::attr::seq(type  =>'protein',
				       id    =>$a->{t_name},
	);
    
    my $p_seq=new Game::seq(name     => $a->{t_name},
			    residues => $a->{p_seq}
	);
    
    $p_seq->attr($p_seq_attr);

    my $t_seq_attr=new Game::attr::seq(type  =>'transcript',
				       id    =>$a->{t_name},
	);
    
    my $t_seq=new Game::seq(name     => $a->{t_name},
			    residues => $a->{t_seq}
	);
    
    $t_seq->attr($t_seq_attr);

    push @ret,$t_seq;
    push @ret,$p_seq;
    return \@ret;
}
#</sub>-------------------------------------------------------------------- 


#<sub>--------------------------------------------------------------------- 
sub MapPosition{
    my ($self)=@_;
    my $txt="=========>STDERR::GAME Has \'map_position \' as an element
             but could not find something similar in
             CGL::Annotation object\n";
    print STDERR $txt;
    return 0;
}
#</sub>--------------------------------------------------------------------- 

#<sub>---------------------------------------------------------------------- 
sub Annotation{
    my ($self)=@_;
  
    my $CGL=$self->{Annot};
    my $HITS=$self->AutoAnts;
    if(defined($CGL)){
       #) <gene cgl is array..needs> 
       #) <to take care of more gene come along>   
	
	my $gn_arr=$CGL->genes();

	my @ret;
	foreach my $geneCGL(@{$gn_arr}){
	    #my $geneCGL=$CGL->gene(0);
	    my $Contig =$CGL->contig(0);
	    
	    my $tscrCGL =$geneCGL->transcripts();
	    
	    my $geneGM = $self->make_gene($geneCGL);
	    #)<comes from gene>
	    my $dbxrefGM =" ";
	    my $dateGM =" ";
	    
	   
       
	    foreach my $ts(@$tscrCGL){
		my $annotGM = new Game::Annotation;
		
		$self->generalized_stuffing($geneCGL,$annotGM);
		
		my $fst = $self->get_featsets($ts,$Contig);
		
		$annotGM->feature_set($fst);

		$annotGM->gene($geneGM);
		# print Dumper($geneGM);die;
		push @ret,$annotGM;
	    }
	}
	return \@ret;
    }elsif(defined($HITS)){
	my $ant  = $self->Ant();
	#print Dumper($ant->[0]->feature_set);die;
	return $ant;#flag
    }
}
#</sub>--------------------------------------------------------------------- 
#<sub>---------------------------------------------------------------------- 
sub add_hits{
    my ($self,$hit)=@_;
      
    my $hits=$self->Hits;
    push @{$hits},$hit;
    $self->Hits($hits);
    return $hits;
}

#</sub>--------------------------------------------------------------------- 
#<sub>---------------------------------------------------------------------- 
sub add_autoant{
    my ($self,$hit)=@_;
      
    my $hits=$self->AutoAnts;
    push @{$hits},$hit;
    $self->AutoAnts($hits);
    return $hits;
}

#</sub>--------------------------------------------------------------------

#<sub>---------------------------------------------------------------------- 
sub add_contig{
    my ($self,$hit)=@_;
    $self->Contig($hit);
    return $hit;
}

#</sub>--------------------------------------------------------------------- 

#<sub>---------------------------------------------------------------------- 
sub Ant{
    my ($self)=@_;
    
    my $HITS=$self->AutoAnts;
    if (defined($HITS)){
	print STDERR "((%): Proccessing Auto Annotation\n";
    }else{
	print STDERR "((%): No Auto Annotated Hits\n";
	return;
    } 

    my @ants_ret;
    foreach my $PerHit(@$HITS){
	my $annot = new Game::Annotation;
	my $ant= $self->get_ant_fts($PerHit);
	$annot->feature_set($ant);
	push @ants_ret,$annot;	
    }
    
    return \@ants_ret;
}
#</sub>--------------------------------------------------------------------- 
#<sub>---------------------------------------------------------------------- 
sub Computational_Analysis{
    my ($self)=@_;
  
    
    my $HITS=$self->Hits;
    
    if (defined($HITS->[0])){
	print STDERR " ((%): Processing Compuational_Analisis \n";
    }else{
	print STDERR " ((%): No Computational_Analysis Data\n";
	return undef;
    } 

    my @comps_ret;
    foreach my $PerHit(@$HITS){
	my $comp= $self->get_comp_A($PerHit);
	push @comps_ret,$comp;	
    }
    return \@comps_ret;
}
#</sub>--------------------------------------------------------------------- 

#<sub>---------------------------------------------------------------------- 
sub get_comp_A{
    my ($self,$hits)=@_;
    my $compA=new Game::computational_analysis();
    if(defined($hits->[0])){
	$compA->program($hits->[0]->algorithm);
    }
    my @rs;
    foreach my $hit(@$hits){
	my $res_set=$self->get_res_set($hit);
	push @rs,$res_set;
    }
    $compA->result_set(\@rs);
    return $compA;
}

#</sub>--------------------------------------------------------------------- 
#<sub>---------------------------------------------------------------------- 
sub get_ant_fts{
    my ($self,$a)=@_;
    #)add seq_relationship
    my $ft_set= new Game::feature_set();
   
    if(defined($a->{hit})){
	my $hit=$a->{hit};
	my $codon_pos=$a->{t_offset};
	$ft_set->feature_span($self->get_ant_ftspans($hit,$codon_pos));
    }
    $ft_set->seq($self->Get_Seq_from_Autoant($a));
    $ft_set->name($a->{t_name});
    
    return $ft_set;
}

#</sub>--------------------------------------------------------------------- 

#<sub>---------------------------------------------------------------------- 
sub get_res_set{
    my ($self,$hit)=@_;
    #)add seq_relationship
    my $res_set= new Game::result_set();
    $res_set->name($hit->name);
    $res_set->output($self->get_allopts($hit));
    $res_set->result_span($self->get_res_spans($hit));
    return $res_set;
}

#</sub>--------------------------------------------------------------------- 


#<sub>---------------------------------------------------------------------- 
sub get_output{
    my ($self,$tp,$v)=@_;
    
    if(defined($tp) && defined($v)){
	my $oput = new Game::output();
	$oput->type($tp);
	$oput->value($v);
	return $oput;
    }else{
	return undef;
    }
}
#</sub>--------------------------------------------------------------------- 

#<sub>---------------------------------------------------------------------- 
sub get_allopts{
    my ($self,$hsp)=@_;
    my @Oput;
    my $op1;
    my $op2;
    my $op3;
    if(ref($hsp) =~m/PhatHit/){
	$op1=$self->get_output('Percent_Aligned:Query',$hsp->pAq());
        $op2=$self->get_output('Percent_Aligned:Hit',$hsp->pAh());
	$op3=$self->get_output('E/P:Significance',$hsp->significance);

	push @Oput,$op1,$op2,$op3;
	return \@Oput;
    }elsif(ref($hsp)=~m/PhatHsp/){
	$op1=$self->get_output('Fraction_Identical:query',$hsp->frac_identical('query'));
	$op2=$self->get_output('Fraction_Identical:total',$hsp->frac_identical('total'));
	push @Oput,$op1,$op2;
	return \@Oput;
    }
    return undef;
}
#</sub>--------------------------------------------------------------------- 

#<sub>---------------------------------------------------------------------- 
sub get_all_seqrels{
    my ($self,$hsp)=@_;
    my @Oput;
    my $op1;
    my $op2;
    my $op3;
    if(ref($hsp) =~m/PhatHit/){
	$op1=$self->get_output('Percent_Aligned:Query',$hsp->pAq());
        $op2=$self->get_output('Percent_Aligned:Hit',$hsp->pAh());
	$op3=$self->get_output('E/P:Significance',$hsp->significance);
	push @Oput,$op1,$op2,$op3;
	return \@Oput;
    }elsif(ref($hsp)=~m/PhatHsp/){
	$op1=$self->get_output('Fraction_Identical:query',$hsp->frac_identical('query'));
	$op2=$self->get_output('Fraction_Identical:total',$hsp->frac_identical('total'));
	push @Oput,$op1,$op2;
	return \@Oput;
    }
   
}
#</sub>--------------------------------------------------------------------- 

#<sub>---------------------------------------------------------------------- 
sub get_ant_ftspans{
    my ($self,$hit,$codon)=@_;
    
    my $i=0;
    my @ret;
    my $codon_span=$self->get_codon_span("start_codon",$codon,$codon+3);
    push @ret,$codon_span;
    foreach my $hsp($hit->hsps){
	
	my $fspan= new Game::feature_span();

	my $sr2=$self->get_seq_rel($hsp->nB('query'),
				   $hsp->nE('query'),
				   'Exon',
				   $hsp->seq_id,
				   $hsp->hit_string,
	    );
	
	
	$fspan->seq_relationship($sr2);
	push @ret,$fspan;
    }

    return \@ret;
}

#</sub>--------------------------------------------------------------------- 

#<sub>---------------------------------------------------------------------- 
sub get_res_spans{
    my ($self,$hit,$h_id)=@_;
    
    my $i=0;
    my @ret;
    foreach my $hsp($hit->hsps){
	
	my $rspan= new Game::result_span();
	my $attr= new Game::attr::result_span();
	if(defined($h_id)){}else{
	    $h_id=" ";
	}
	$attr->id(++$i);
	$rspan->type($hsp->algorithm);
	$rspan->output($self->get_allopts($hsp));
	$rspan->name($hsp->name);
	my $score_DD=$self->decipher_score($hsp);
	$rspan->score($score_DD);
	my @srs;
	my $sr1=$self->get_seq_rel($hsp->nB('query'),
				   $hsp->nE('query'),
				   'query',
				   $hsp->name,
				   $hsp->query_string,
	    );

	push @srs,$sr1;

	my $sr2=$self->get_seq_rel($hsp->nB('hit'),
				   $hsp->nE('hit'),
				   'subject',
				   $hsp->name,
				   $hsp->hit_string,
	    );
	
	push @srs,$sr2;

	$rspan->seq_relationship(\@srs);
	push @ret,$rspan;
    }
    
    return \@ret;
}

#</sub>--------------------------------------------------------------------- 
sub decipher_score{
    my ($self,$hsp)=@_;
   
    my ($class) = ref($hsp) =~ /.*::(\S+)$/;
    
    my $score;
    if($class eq 'blastx'){
	$score=$hsp->significance();
	
    }
    elsif ($class eq 'blastn'){
	$score=$hsp->significance();
	
    }
    elsif ($class eq 'protein2genome'){
	$score=$hsp->score();
    }
    elsif ($class eq 'est2genome'){
	$score=$hsp->score();
    }
    
    elsif (ref($hsp)  =~ /snap/){
	$score=$hsp->score();
    }
    elsif (ref($hsp)  =~ /repeatmasker/){
	$score=$hsp->score();
    }
    else {
	$score="UNK";
    }
    
    return $score;
}
#<sub>---------------------------------------------------------------------- 

sub make_gene(){
    my($self,$gnCGL)=@_;
    
    my $ret_gn= new Game::gene();
    $ret_gn->name($gnCGL->{name});

    my $gn_attr= new Game::attr::gene();
    $gn_attr->id(new Game_Xml->get_actual($gnCGL->feature_dbxref));
    $gn_attr->association("IS");
    $ret_gn->attr($gn_attr);

    return $ret_gn;


}
#</sub>--------------------------------------------------------------------- 

#<sub>---------------------------------------------------------------------- 
sub generalized_stuffing{
    my($self,$fromft,$toft)=@_;
    
    
    my $prop=$fromft->{properties};
    my $name=$fromft->{name};
    my $type=$fromft->{type};
    my $synonym = $prop->{synonym};
    my $comment =$prop->{note};
    
   
    
    my $to_prop= new Game::property;
    $to_prop->type('locus_tag');
    $to_prop->value($prop->{locus_tag});
    
    if($toft->can("name")){
	$toft->name($name);
    }

    $toft->type($type);
    $toft->synonym($synonym);
    $toft->comment($comment);
    $toft->property($to_prop);

    
}
#</sub>---------------------------------------------------------------------

#<sub>---------------------------------------------------------------------- 
sub get_seq_rel{
    my($self,$start,$end,$type,$seq,$al)=@_;

    my $seq_rel= new Game::seq_relationship;
    my $seq_rel_attr= new Game::attr::seq_relationship;

    my $span = new Game::span;
    $span->start($start);
    $span->end($end);

    $seq_rel_attr->type($type);
    $seq_rel_attr->seq($seq);
        
    $seq_rel->attr($seq_rel_attr);
    $seq_rel->span($span);
    $seq_rel->alignment($al);
    return $seq_rel;
}

#</sub>---------------------------------------------------------------------- 
#<sub>---------------------------------------------------------------------- 
sub get_featspans{
     my($self,$exons)=@_;

     my @f_spans;
    
    foreach my $ex (@$exons){

	my $f_span    =new Game::feature_span;
	my $sp_name   =$ex->name;
	my $sp_type   =$ex->type;
	$f_span->name($sp_name);
	$f_span->type($sp_type);
	
	my $seq_rel= new Game::seq_relationship;
	my $seq_rel_attr= new Game::attr::seq_relationship;

	my $span = new Game::span;
	my ($base1beg,$base1end) = SpaceBase::tobase1base($ex->nbeg,$ex->nend); 
	$span->start($base1beg);
	$span->end($base1end);

	$seq_rel_attr->type('query');

	if($ex->id){
	    $seq_rel_attr->seq($ex->id);
	}elsif($ex->feature_id){
	    $seq_rel_attr->seq($ex->feature_id);
	}

	$seq_rel->attr($seq_rel_attr);
	$seq_rel->span($span);

	$f_span->seq_relationship($seq_rel);
	   
	push @f_spans,$f_span;
    }
     return \@f_spans;
}
#</sub>--------------------------------------------------------------------- 

#<sub>---------------------------------------------------------------------- 
sub get_codon_span{
     my($self,$sp_type,$nbeg,$nend)=@_;

     
     my $f_span=new Game::feature_span;
     $f_span->name($sp_type);
     $f_span->type($sp_type);
     
     my $seq_rel= new Game::seq_relationship;
     my $seq_rel_attr= new Game::attr::seq_relationship;

     my $span = new Game::span;
     my ($base1beg,$base1end) = SpaceBase::tobase1base($nbeg,$nend); 
     $span->start($base1beg);
     $span->end($base1end);
     
     $seq_rel_attr->type('query');
     $seq_rel->attr($seq_rel_attr);
     $seq_rel->span($span);
     
     $f_span->seq_relationship($seq_rel);
     
     return $f_span;
     
}
#</sub>--------------------------------------------------------------------- 
#<sub>---------------------------------------------------------------------- 
sub get_start_codon{
    my($self,$ft,$c)=@_;

    my $start_codon_pos;
    my $p = $ft->translation(0);
    if(defined($p)){
#----	$start_codon_pos = $p->metaPos($ft,0);
    }else{
	print STDERR "No Protein\n";
	return;
    }

#---my $offset=$ft->metaPos($c,$start_codon_pos);
    my $offset=3;#nw
    my $codons=substr($ft->residues,$start_codon_pos,3);
    my $start=$offset;
    my $end=$start+3;
    #print "offset: ",$offset,"\n";
    #print "C O D O N:",$codons;
    #print "Start_codon_start: ",$start,"\n";
    #print "Start_codon_end: ",$end,$c,"\n";
    #die;
#---return $start,$end;
    return 1,2;
}
#</sub>---------------------------------------------------------------------- 
#<sub>---------------------------------------------------------------------- 
sub get_stop_codon{
    my($self,$ft)=@_;
    

    my $stop_codon_pos;
    my $p = $ft->translation(0);
    if(defined($p)){
#-----	$stop_codon_pos = $p->metaPos($ft,length($p->residues));
    }else{
	print STDERR "No Protein\n";
	return;
    }

    my $start=$stop_codon_pos;
    my $end=$start+3;
#----    return $start,$end;
    return 1,2;
}
#</sub>---------------------------------------------------------------------- 

#<sub>---------------------------------------------------------------------- 
sub get_featsets{
    my($self,$ft,$c)=@_;

    my $exons= $ft->exons;
    my $introns= $ft->introns;
    my $trlations=$ft->translations;

    my $f_set_ex=new Game::feature_set;
    my $f_set_in=new Game::feature_set;
    my $f_set_tr=new Game::feature_set;

    $self->generalized_stuffing($ft,$f_set_ex);
    $self->generalized_stuffing($ft,$f_set_in);
    $self->generalized_stuffing($ft,$f_set_tr);

    
    my $f_spans_ex=$self->get_featspans($exons);
    my $f_spans_in=$self->get_featspans($introns);
    my $f_spans_tr=$self->get_featspans($trlations);
   
   # print "system out\n";die;
    my ($start,$end)=$self->get_start_codon($ft,$c);
   
    my $start_codon_span=$self->get_codon_span("start_codon",$start,$end);

    unshift @{$f_spans_ex}, $start_codon_span;
   

    $f_set_ex->feature_span($f_spans_ex);
    $f_set_in->feature_span($f_spans_in);
    $f_set_tr->feature_span($f_spans_tr);

    $f_set_ex->name($ft->{properties}{gene}.":".$ft->{properties}{product});
    $f_set_in->name($ft->{properties}{gene}.":".$ft->{properties}{product});
    $f_set_tr->name($ft->{properties}{gene}.":".$ft->{properties}{product});

    my @ret;

    push @ret,$f_set_ex;
    push @ret,$f_set_in;
    push @ret,$f_set_tr;

    return \@ret;
}
#</sub>--------------------------------------------------------------------- 

sub  Write(){
    my($self)=@_;

    close $self->XmlP->sting(); 
}
#<end>---------------------------------------------------------------------- 



1;
#</end>--------------------------------------------------------------------- 
