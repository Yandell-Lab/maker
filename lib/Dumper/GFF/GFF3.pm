#<top> ------------------------------------------------------------------------
#     <Name>       GFF3.pm                    </Name> 
#     <class>      Dumper::GFF
#                  ::GFF3                     </class>       
#     <Author>     M Hadi Islam               </Author>
#     <email>      hadi@genetics.utah.edu     </email>
#     <Does>       Prints xml dump of 
#                  CGL::Annotation Object     </Does>
#</top> ----------------------------------------------------------------------
package Dumper::GFF::GFF3;

#------------- Requires-------------------------------------------
use Class::Struct;
use Data::Dumper;
use CGL::Annotation;
use Dumper::GFF::GFF3_DEF;
use Dumper::XML::Game_Xml;
use Dumper::XML::Game;
use strict "vars";
use strict "refs";
use SpaceBase;
use Fasta;
use FastaChunker;
use Iterator::Fasta;
#------------- Requires-------------------------------------------

BEGIN {
	$ENV{ZOE} = '/usr/local/SNAP'
}


#<classdef>--------------------------------------------------------------------
struct Dumper::GFF::GFF3 =>
{
    Annot          =>'$',
    ChaosFile      =>'$',
    FastaFile      =>'$',
    Iter           =>'$',
    Hits           =>'$',
    GFF3           =>'$',
    AutoAnts       =>'$',
    Contig         =>'$',
    Version        =>'$',
    BigName        =>'$',
    MRNA           =>'$',
    Game           =>'$',
    CDS            =>'$',
    exon_num       =>'$',
    nbeg           =>'$',
    nend           =>'$',
    strand         =>'$',
    AutoGenes     =>'$',
    WriteFile     =>'$',
    GENTS         =>'$',
    FASTA         =>'$',
};



sub init{
    my ($self,%args)=@_;
    %{$self}=%args;
    $self->Version(3); 
    return $self;
}

#</classdef>--------------------------------------------------------------------



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#---------------------- SUBS ---------------------------------------------------
#-------------------------------------------------------------------------------


#<sub>--------------------------------------------------------------------------
sub GFF3{
    my ($self)=@_;

    my $gff3 = new GFF3::GFF3();
    $self->Seq();
    my $meta   = $self->Meta();
    $self->exon_num(1);
    my $entry    = $self->Entries();
    
    $gff3->Meta($meta);
    $gff3->Entry($entry);
    my $Gee= new GFF3_DEF->print($gff3);
    $self->Write($Gee);
}
#</sub>----------------------------------------------------------------------- 

#<sub>------------------------------------------------------------------------ 
sub Meta{
    my ($self)=@_;
    
    my $annot=$self->{Annot};
    my $HITS=$self->AutoAnts;
      
    if(defined($HITS)){
	my $seq  = $self->Seq();
	my $meta=$self->Get_Meta_from_Auto($seq);
	my $fasta_str=$seq->residues();
	my $form_fasta=Fasta::formatSeq(\$fasta_str,60);
	
	$meta->FASTA("##FASTA\n".">".$seq->name()."\n".$form_fasta."\n");
    	return $meta;
    }
    
    if(defined($annot)){
	my $meta=$self->Get_Meta_from_CGL($annot);
    	return $meta;
    }
}
#</sub>-------------------------------------------------------------------- 


#<sub>------------------------------------------------------------------------ 
sub Get_Meta_from_Auto{
    my ($self,$seq)=@_;
       
    my @on;
    foreach my $cling (@{$self->AutoAnts()}){
	foreach my $cl (@$cling){
	    push @on,$cl;
	}
    }
   
    my $ant  = $self->Ant(\@on);
    $self->order_Auto_mrna($ant);
    my $name=$seq->name();
    $self->BigName($name);  
    
    my $gff_meta=new GFF3::meta(gff_version        =>"##gff-version\t3\n",
				sequence_region    =>"$name\t.\t"."contig"."\t".$self->nbeg."\t".$self->nend()."\t.\t.\t.\tID=$name\n",
				
	);
    
     return $gff_meta;
}
#</sub>-------------------------------------------------------------------- 
#<sub>------------------------------------------------------------------------ 
sub Get_Meta_from_CGL{
    my ($self,$annot)=@_;
    
    my $con=$annot->contig(0);
    my $name=$con->name;
    $self->BigName($name);  
    my $res=$con->residues();
    my $type =$con->type();
    my $id = $con->id();
   
    my $gff_meta=new GFF3::meta(gff_version        =>"##gff-version\t3\n",
			       sequence_region    =>"##sequence-region\t".$name."\t".$con->nbeg."\t".$con->nend()."\n",
			       
	);
    return $gff_meta;
}
#</sub>-------------------------------------------------------------------- 
#<sub>---------------------------------------------------------------------- 
sub Entries{
    my ($self)=@_;
  
    my $CGL=$self->{Annot};
   # print Dumper($self);die;
    my $HITS=$self->AutoAnts;
    my @ret_ar;
    if(defined($HITS)){
	my $i=0;
	my @BMZD;
	#print Dumper($HITS);die;
	foreach my $GEE (@$HITS){
	    my $ant  = $self->Ant($HITS->[$i]);
	    my $seq  = $self->Seq();
	    my $auto_gn=$self->AutoGenes->[$i];
	    my $ret=$self->formation($ant,$seq,$auto_gn,$i);
	    foreach my $stupid (@$ret){
		push @BMZD,$stupid;
	    }
	    print "-------TO FORM---------------------------------------------------------->$i<<<<<\n";
	    #print Dumper(@BMZD);
            #print ref($ret);die;
	    $i++;
	}
	return \@BMZD;
    }elsif(defined($CGL)){
	
	my $game_structs=$self->get_game();
	$self->order_mrna($game_structs);
	$self->exon_num(1); 
	my $geneCGL=$CGL->gene(0);
	my $Contig =$CGL->contig(0);
	my $tscrCGL =$geneCGL->transcripts();
	my $geneGM = $self->make_gene_gff($geneCGL);
	push(@ret_ar,$geneGM);
            
	my $dbxrefGM =" ";
	my $dateGM =" ";
	my @ret;
	my $mrna_ids=1;
	
	

	foreach my $ts(@$tscrCGL){
	    my $mrna=$self->make_mrna($geneCGL,$ts,$mrna_ids++);
	    push(@ret_ar,$mrna);
	}

	foreach my $ts(@$tscrCGL){
	    my $fst = $self->get_featsets($ts,$Contig);
	    foreach my $elm(@$fst){
		foreach my $spn(@{$elm->feature_span}){
		    if( $spn->type =~m/intron/){}else{
			my $gff_spn= $self->make_span($spn);
			if($gff_spn==-1){}else{
			    if($gff_spn->type()=~m/polypeptide/){
				unshift (@ret_ar,$gff_spn);
			    }else{
				push(@ret_ar,$gff_spn);
			    }
			}
		    }
		}
	    }
	}
    }
    
    my $thing=shift @ret_ar;
    push @ret_ar,$thing;
    return \@ret_ar;
}
#</sub>--------------------------------------------------------------------- 
sub order_mrna(){
    my($self,$st_ar)=@_;
    my @A;
    my @CDS;
    my $j=1;
    foreach my $annot(@$st_ar){
	my $fst=$annot->feature_set();
	$j++;
	foreach my $f_spn(@$fst){
	    foreach my $spn(@{$f_spn->feature_span}){
		if($spn->type()=~m/intron/){next;}
		if($spn->type()=~m/start_codon/){next;}
		if($spn->type()=~m/protein/){
		    my $spn_start= $spn->seq_relationship->span->start;
		    my $spn_end= $spn->seq_relationship->span->end;
		    my $cds=$CDS[$spn_start][$spn_end];
		    push @$cds,'mRNA_007'.$j;
		    $CDS[$spn_start][$spn_end]=$cds;
		}
		my $spn_start= $spn->seq_relationship->span->start;
		my $spn_end= $spn->seq_relationship->span->end;
		my $mrna=$A[$spn_start][$spn_end];
		push @$mrna,'mRNA_007'.$j;
		$A[$spn_start][$spn_end]=$mrna;
	    }
	}
    }
    $self->CDS(\@CDS);
    $self->MRNA(\@A);
}
#--------------------------------------------------------------------
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
    my $residue= $$res;
    my $len= length($residue);
    $self->nbeg(1);
    $self->nend($len);

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

#</sub>--------------------------------------------------------------------- 
sub formation(){
    my($self,$st_ar,$seq,$auto_gn,$idx)=@_;
    
    my $mrna_ids=1;
    my @A;
    my @CDS;
    my $j=0;
    my @ret_ar;

    print "-------IN FORM---------------------------------------------------------->$j<<<<<\n";

     #$auto_gn=$self->AutoGenes->[0];
    my $GFFgene=$self->makeAutogene($auto_gn,$idx);

    push(@ret_ar,$GFFgene);

    foreach my $annot(@$st_ar){
	my $mrna=$self->make_mrna_from_auto_annot($annot,$idx.$mrna_ids++,$idx);
	push(@ret_ar,$mrna);
    }
   
    foreach my $annot(@$st_ar){
       

	my $fst=$annot->feature_set->feature_span();

	$j++;
	foreach my $spn(@$fst){
	    if($spn->seq_relationship->attr->type()=~m/intron/){next;}
	    if(defined($spn->type())){
		if($spn->type()=~m/start_codon/){next;}
	    }
	   
	    my $gff_spn= $self->make_span($spn);

	    if($gff_spn==-1){}else{
		if($gff_spn->type()=~m/polypeptide/){
		    unshift (@ret_ar,$gff_spn);
		}else{
		    push(@ret_ar,$gff_spn);
		}
	    }
	}
    }
    return \@ret_ar;
}
#--------------------------------------------------------------------

#</sub>--------------------------------------------------------------------- 
sub order_Auto_mrna(){
    my($self,$st_ar)=@_;
  
    my @A;
    my @CDS;
    my $j=0;
    my $end=1;
    foreach my $annot(@$st_ar){
	my $fst=$annot->feature_set->feature_span();
	$j++;
	foreach my $spn(@$fst){
	    if($spn->seq_relationship->attr->type()=~m/intron/){next;}
	    if(defined($spn->type())){
		if($spn->type()=~m/start_codon/){
		  
		}
	    }
	    if($spn->seq_relationship->attr->type()=~m/protein/){
		my $spn_start= $spn->seq_relationship->span->start;
		my $spn_end= $spn->seq_relationship->span->end;
		if($spn_end > $end){
		    $end=$spn_end;
		}
		my $cds=$CDS[$spn_start][$spn_end];
		push @$cds,'mRNA_008'.$j;
		$CDS[$spn_start][$spn_end]=$cds;
	    }
	    
	    my $spn_start= $spn->seq_relationship->span->start;
	    my $spn_end= $spn->seq_relationship->span->end;
	    if($spn_end > $end){
		$end=$spn_end;
	    }
	    my $mrna=$A[$spn_start][$spn_end];
	    push @$mrna,'mRNA_008'.$j;
	    $A[$spn_start][$spn_end]=$mrna;
	}
    }
  
    
    $self->CDS(\@CDS);
    $self->MRNA(\@A);
}
#--------------------------------------------------------------------
sub get_game(){
    my($self)=@_;

    my $chaos_file=$self->ChaosFile;
    my $CGL= new CGL::Annotation($chaos_file);
    my $IO= new Dumper::CGL::Annotation::XML::Game->init(Annot=>$CGL);
    return $IO->Annotation();

}
#----------------------------------------------------------------------
sub make_struct_del(){
    my($self,$trCGL)=@_;
    
    my $CGL=$self->{Annot};
    my $HITS=$self->AutoAnts;

   if(defined($CGL)){
	#) <gene cgl is array..needs> 
	#) <to take care of more gene come along>    
	my $geneCGL=$CGL->gene(0);
	my $Contig =$CGL->contig(0);
	
	my $tscrCGL =$geneCGL->transcripts();
	
	my $geneGM = $self->make_gene($geneCGL);
	#)<comes from gene>
	my $dbxrefGM =" ";
	my $dateGM =" ";

	my @ret;
	foreach my $ts(@$tscrCGL){
	    my $annotGM = new Game::Annotation;
	    $self->generalized_stuffing($geneCGL,$annotGM);
	    my $fst = $self->get_featsets($ts,$Contig);
	    $annotGM->feature_set($fst);
	    $annotGM->gene($geneGM);
	    push @ret,$annotGM;
	}

	return \@ret;

    }
}



#<sub>---------------------------------------------------------------------- 
sub make_gene_gff(){
    my($self,$gnCGL)=@_;
    
    my $gn_ent= new GFF3::entry();
    my $gn_attr= new GFF3::attr();

    $gn_attr->Name($gnCGL->{name});
    $gn_attr->ID($gnCGL->{type}.new Game_Xml->get_actual($gnCGL->feature_dbxref));
    
    $gn_ent->seqid($self->BigName);
    $gn_ent->start($gnCGL->nbeg());
    $gn_ent->end($gnCGL->nend());
    $gn_ent->type($gnCGL->{type});
   
    
    $gn_ent->attr($gn_attr);
    
    return $gn_ent;


}
#</sub>--------------------------------------------------------------------- 
#<sub>---------------------------------------------------------------------- 
sub makeAutogene(){
    my($self,$an,$j)=@_;
    

    my $g_name     = $an->{g_name};
    my $g_s        = $an->{g_start};
    my $g_e        = $an->{g_end};
    my $g_strand   = $an->{g_strand};
    

    my $gn_ent= new GFF3::entry();
    my $gn_attr= new GFF3::attr();

    $gn_attr->Name($g_name);
    $gn_attr->ID("gene000".$j);
    
    $gn_ent->seqid($self->BigName);
    
    ($g_s,$g_e)=($g_e,$g_s) if $g_s > $g_e;
    $gn_ent->start($g_s);
    $gn_ent->end($g_e);
    $gn_ent->type("gene");
    if($g_strand == 1){
	$gn_ent->strand('+');
	$self->strand('+');
    }elsif($g_strand == 0){
	$gn_ent->strand('-');
	$self->strand('-');
    }
    
    $gn_ent->attr($gn_attr);
    
    return $gn_ent;


}
#</sub>--------------------------------------------------------------------- 


sub make_mrna_from_auto_annot(){
    my($self,$annot,$i,$j)=@_;
    

    my $fst=$annot->feature_set->feature_span();
    
    my $Mend;
    my $Mstart;
    my $M_name=$annot->feature_set->name();
    foreach my $spn(@$fst){
	if($spn->seq_relationship->attr->type()=~m/intron/){next;}
	if(defined($spn->type())){
	    if($spn->type()=~m/start_codon/){
		$Mstart=$spn->seq_relationship->span->start();
		$Mend=$spn->seq_relationship->span->end();
	    }
	}
	my $spn_end= $spn->seq_relationship->span->end;

    }

    my $gn_ent= new GFF3::entry();
    my $gn_attr= new GFF3::attr();
   
    $gn_attr->Name($M_name);
    $gn_attr->ID('mRNA_008'.$i);
    $gn_attr->Parent("gene000$j");
    $gn_ent->seqid($self->BigName);
    #print $Mstart," ". $Mend,"\n";
    ($Mstart,$Mend)=($Mend,$Mstart) if $Mstart > $Mend;
    #print $Mstart," ",$Mend,"\n";die;
    $gn_ent->start($Mstart);
    $gn_ent->end($Mend);
    $gn_ent->type("mRNA");
    $gn_ent->strand($self->strand());
    $gn_ent->attr($gn_attr);
    
    return $gn_ent;
}
#</sub>---------------------------------------------------------------------
#<sub>---------------------------------------------------------------------- 

sub make_mrna(){
    my($self,$gnCGL,$tr,$i)=@_;
    
    my $gn_ent= new GFF3::entry();
    my $gn_attr= new GFF3::attr();

    $gn_attr->Name($tr->{name});
    $gn_attr->ID('mRNA_007'.$i);
    $gn_attr->Parent($gnCGL->{type}.new Game_Xml->get_actual($gnCGL->feature_dbxref));
    $gn_ent->seqid($self->BigName);
    $gn_ent->start($tr->nbeg());
    $gn_ent->end($tr->nend());
    $gn_ent->type($tr->type);
   
    $gn_ent->attr($gn_attr);
    
    return $gn_ent;
}
#</sub>--------------------------------------------------------------------- 
#<sub>---------------------------------------------------------------------- 
sub make_span(){
    my($self,$spn)=@_;
  

    
    my $start=$spn->seq_relationship->span->start;
    my $end=$spn->seq_relationship->span->end;
  
    my $A=$self->MRNA();
   
    my @a;
  
    @a=@{$A};
    my $mrna=$a[$start][$end];
    #print $mrna,$start,$end,"\n";die;  
    if(defined($mrna->[0])){
	my $parent=join(",",@$mrna);
	
	foreach my $some (1..@$mrna){
	    pop @{$mrna};
	}
	$a[$start][$end]=$mrna;

	$self->MRNA(\@a);
	my $gn_ent= new GFF3::entry();
	my $gn_attr= new GFF3::attr();

	$gn_attr->Name($spn->name);
	$gn_ent->seqid($self->BigName);
	my $mmms=$spn->seq_relationship->span->start();
	my $mmme=$spn->seq_relationship->span->end();
	($mmms,$mmme)=($mmms,$mmme) if $mmms > $mmme;
	$gn_ent->start($mmms);
	$gn_ent->end($mmme);
	$gn_ent->type("exon");
   
	if($gn_ent->type=~m/protein/){
	    my $cd=$self->CDS();
	    my @c;
	    @c=@{$cd};
	    my $der=$c[$start][$end];
    
	    $gn_attr->ID("CDS367895".$self->exon_num);
	    $gn_ent->type('polypeptide');
	    $gn_attr->Derives_from($der->[0]);
	}else{
	    $gn_attr->ID("Exon00".$self->exon_num);
	    my $i=$self->exon_num();
	    $self->exon_num($i+1);
	    $gn_attr->Parent($parent);
	}
	$gn_ent->strand($self->strand);
	$gn_ent->attr($gn_attr);
	return $gn_ent;
    }else{
	return -1;
    }
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
	#my ($base1beg,$base1end) = SpaceBase::tobase1base($ex->nbeg,$ex->nend); 
	$span->start($ex->nbeg);
	$span->end($ex->nend);

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
sub add_autoant{
    my ($self,$hit)=@_;
      
    my $hits=$self->AutoAnts;
    push @{$hits},$hit;
    $self->AutoAnts($hits);
    return $hits;
}

#</sub>--------------------------------------------------------------------
#<sub>---------------------------------------------------------------------- 
sub add_auto_gene{
    my ($self,$hit)=@_;
      
    my $hits=$self->AutoGenes;
    push @{$hits},$hit;
    $self->AutoGenes($hits);
    
    return $hits;
}

#</sub>--------------------------------------------------------------------
#<sub>---------------------------------------------------------------------- 
sub add_GENTS{
    my ($self,$hit)=@_;
      
    my $hits=$self->GENTS;
    push @{$hits},$hit;
    $self->GENTS($hits);
    
    return $hits;
}

#</sub>--------------------------------------------------------------------



#<sub>---------------------------------------------------------------------- 
sub Ant{
    my ($self,$WHAT_THE)=@_;
    #print Dumper($WHAT_THE);
    my $HITS=$WHAT_THE;
    if (defined($HITS)){
	print STDERR "((%): Proccessing Auto Annotation\n";
    }else{
	print STDERR "((%): No Auto Annotated Hits\n";
	return;
    } 

#    foreach my $hit (@{$HITS}){
#	print "nB:".$hit->{hit}->nB('query')." nE".$hit->{hit}->nE('query')."\n";
	
#    }
   
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


#<sub>---------------------------------------------------------------------- 
sub get_ant_ftspans{
    my ($self,$hit,$codon)=@_;


   # print ref $hit;die;
    
    my $i=0;
    my @ret;
    my $codon_span=$self->get_codon_span("start_codon",$hit->nB('query'),$hit->nE('query'));
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
	$start_codon_pos = $p->metaPos($ft,0);
    }else{
	print STDERR "No Protein\n";
	return;
    }

    my $offset=$ft->metaPos($c,$start_codon_pos);

    my $codons=substr($ft->residues,$start_codon_pos,3);
    my $start=$offset;
    my $end=$start+3;
    #print "offset: ",$offset,"\n";
    #print "C O D O N:",$codons;
    #print "Start_codon_start: ",$start,"\n";
    #print "Start_codon_end: ",$end,$c,"\n";
    #die;
    return $start,$end;
}
#</sub>---------------------------------------------------------------------- 
#<sub>---------------------------------------------------------------------- 
sub get_stop_codon{
    my($self,$ft)=@_;
    

    my $stop_codon_pos;
    my $p = $ft->translation(0);
    if(defined($p)){
	$stop_codon_pos = $p->metaPos($ft,length($p->residues));
    }else{
	print STDERR "No Protein\n";
	return;
    }

    my $start=$stop_codon_pos;
    my $end=$start+3;
    return $start,$end;
}
#</sub>---------------------------------------------------------------------- 
sub  Write(){
    my($self,$gee)=@_;
    
    my $file=$self->{WriteFile};
    open FFF,">$file";
    print FFF $gee; 
    close FFF;
}


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

#<end>---------------------------------------------------------------------- 
1;
#</end>--------------------------------------------------------------------- 
