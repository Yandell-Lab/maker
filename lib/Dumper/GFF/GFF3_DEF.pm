#<top> ----------------------------------------------------------------------
#     <Name>       GFF3_DEF.pm                            </Name>        
#     <class>      GFF3_DEF                               </class>
#     <Author>     M Hadi Islam                           </Author>
#     <email>      hadi@genetics.utah.edu                 </email>
#     <Does>       Contain game xml definitions 
#                  prints object tree as xml              </Does>
#</top> ---------------------------------------------------------------------

package GFF3_DEF;

use Class::Struct;
use Data::Dumper;


#<classdef>-------------------------------------------------------------------
struct GFF3_DEF =>
{
    Ind     => '$',
        
};
#</classdef>------------------------------------------------------------------

#<sub classdef>---------------------------------------------------------------
struct GFF3::GFF3 =>
{
    Meta      => '$',
    Entry     => '$',
    
};

struct GFF3::meta =>
{
    gff_version        =>'$',
    sequence_region    =>'$',
    feature_ontology   =>'$',
    attribute_ontology =>'$',
    source_ontology    =>'$',
    FASTA              =>'$',
};

struct GFF3::GENTS =>
{
        Gene                =>'$',
	Ants                =>'$',
};


struct GFF3::attr =>
{
    ID 	                   =>'$',
    Name                   =>'$',
    Alias                  =>'$',
    Parent                 =>'$',
    Target                 =>'$',
    Gap                    =>'$',
    Derives_from           =>'$',
    Note                   =>'$',
    Dbxref                 =>'$', 
    Ontology_term          =>'$',
};

struct GFF3::entry =>
{
    seqid       =>'$',
    source      =>'$',
    type        =>'$',
    start       =>'$',
    end         =>'$',
    score       =>'$',
    strand      =>'$',
    phase       =>'$',
    attr        =>'$',
};


my %GFF3Order =
(
 GFF3::attr                   =>[qw(ID Name Alias Parent Target Gap Derives_from Note Dbxref Ontology_term)],
 GFF3::entry                  =>[qw(seqid  source  type start end score strand phase)],
 GFF3::meta                   =>[qw(gff_version sequence_region feature_ontology attribute_ontology source_ontology)],
);


#<constructor>--------------------------------------------------------------
sub init{
    my ($self,%args)=@_;
    
    %{$self}=%args;
    return $self;
}
#</constructor>--------------------------------------------------------------


#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#---------------------- SUBS ------------------------------------------------
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
sub keep_order{
    my ($self,$item)=@_;

    my $ref=ref($item);
    my $order=$GFF3Order{$ref};
 
   return $order;
}


#-----------------------------------------------------------------------
sub print{
    my ($self,$gff3,$elm_nm)=@_;

    my $meta=$gff3->Meta();
    my $st.=$self->print_meta($meta);
  
   # print Dumper($gff3);die;
  
    foreach my $elm (@{$gff3->Entry()}){
	 
	  $st.=$self->print_struct($elm);
	  $st.=$self->print_attr($elm->attr);
	  $st.="\n";
    }
    $st.=$self->print_FASTA_from_meta($meta);
    return $st;
}
#<sub>-----------------------------------------------------------------------
sub print_struct{
    my ($self,$elm)=@_;
    
    my $st;
    foreach my $key(@{$self->keep_order($elm)}){
   
	my $value=$elm->$key;
	if(defined($value)){
	    $st.=$value."\t";
	}else{
	    $st.="."."\t";
	}
    }
   
    return $st;
}
#<sub>-----------------------------------------------------------------------
sub print_attr{
    my ($self,$elm)=@_;
  
    my $st;
    foreach my $key(@{$self->keep_order($elm)}){

	my $value=$elm->$key;
	
	if(defined($value)){
	    $st.=$key."=".$value.";";
	}
    }

    chop $st;
    return $st;
}
#<sub>-----------------------------------------------------------------------
sub print_FASTA_from_meta{
    my ($self,$elm)=@_;

    return $elm->FASTA();;
}


#<sub>-----------------------------------------------------------------------
#<sub>-----------------------------------------------------------------------
sub print_meta{
    my ($self,$elm)=@_;

    my $st;
    foreach my $key(@{$self->keep_order($elm)}){

	my $value=$elm->$key;
	if(defined($value)){
	    $st.=$value;
	}
    }
    return $st;
}


#<sub>-----------------------------------------------------------------------
sub get_ind{
    my ($self)=@_;

    my $i=$self->Ind();
    my $ind="      ";  
    my $x_st = $ind x $i;
    return $x_st;
}


#<sub>-----------------------------------------------------------------------
sub add_level{
    my ($self)=@_;
    
    my $i=$self->Ind();
    $i++;
    $self->setInd($i);
}


#<sub>-----------------------------------------------------------------------
sub less_level{
    my ($self)=@_;
    
    my $i=$self->Ind();
    $i--;
    $self->setInd($i);
}


#<sub>-----------------------------------------------------------------------
sub print_att{
    my ($self,$item,$elm_nm)=@_;
    
    my $x_st=$self->get_ind();
    $x_st.= "<".$elm_nm." ";
    while ( my ($key, $value) = each(%{$item}) ) {
	my $elm_count=$self->get_zero_more($value);
	if($elm_count==0){
	    next;
	}
	#print $key,"--------->\n";
	my $act=$self->get_actual($key);
	#print $act,"--------->\n";
	$x_st.= $act."="."\"".$value."\" ";
    }
    $x_st.=">\n";
    print $x_st;
    return $x_st;
}


#<sub>-----------------------------------------------------------------------
sub print_start{
    my ($self,$item)=@_;
    
    $self->add_level();
    my $elm_name;
    if($self->isObj($item)){
	$elm_name=$self->get_name_from($item);
	if($elm_name eq "val"){
	    #$self->less_level();
	    return;
	}
        if(defined($item->attr())){
	    $self->print_attr($item->attr,$elm_name);
	    #print Dumper($item->attr);die;
	}else{  
            my $x_st=$self->get_ind();
            $x_st.="<".$elm_name.">\n";
            print $x_st;
        }
    }else{
	$elm_name=$item;
	if($elm_name eq "val"){
	   # $self->less_level();
	    return;
	}
	my $x_st=$self->get_ind();
	$x_st.="<".$elm_name.">\n";
	print $x_st;
	return $x_st;
    }
}


#<sub>-----------------------------------------------------------------------
sub print_end{
    my ($self,$item)=@_;

    #$self->less_level();
    
    # check is obj
    my $key=  ref($item);
    my $act;
    if($key=~/^Game/){
	$act=$self->get_actual($key);
    }else{
	$act=$self->get_actual($item);
    }


    if($act eq "val"){
	return;
    }
    
    my $x_st=$self->get_ind();
    $x_st.="</".$act.">\n";
    
    print $x_st;
    return $x_st;
}


#<sub>-----------------------------------------------------------------------
sub setInd{
    my ($self,$i)=@_;
    if(0 <= $i ){
	$self->Ind($i);
    }
}


#<sub>-----------------------------------------------------------------------
sub get_name_from{
    my ($self,$item)=@_;
    my $st=ref($item);
    #print $st,"from actual\n";
    my @a=split(/:/,$st);
    #print Dumper(@a);
    my $ret= pop @a;
    return $ret;
}


#<sub>-----------------------------------------------------------------------
sub get_actual{
    my ($self,$st)=@_;

    my @a=split(/:/,$st);
    #print Dumper(@a);
    my $ret= pop @a;
    return $ret;
}


#<sub>-----------------------------------------------------------------------
sub print_header{
    my ($self)=@_;
    print "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
}


#<sub>-----------------------------------------------------------------------
sub Stuffer{
    my ($self,$from_h,$to_h,$head)=@_;

    my $targ = $to_h->{$head};
    while ( my ($key, $value) = each(%{$from_h}) ) {
        print "$key => $value\n";
	$targ->{$key}=$value;
    }

}



#<end>-----------------------------------------------------------------------
1;

