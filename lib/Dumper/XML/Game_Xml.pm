#<top> ----------------------------------------------------------------------
#     <Name>       Game_Xml.pm                 </Name>        
#     <class>      Game_Xml                    </class>
#     <Author>     M Hadi Islam                </Author>
#     <email>      hadi@genetics.utah.edu      </email>
#     <Does>       Contain game xml definitions 
#                  prints object tree as xml    </Does>
#</top> ----------------------------------------------------------------------

package Game_Xml;

use Class::Struct;
use Data::Dumper;


#<classdef>--------------------------------------------------------------------
struct Game_Xml =>
{
    Ind     => '$',
    sting   => '$',
    out     => '$',
    
};


struct Game::game =>
{
    attr                   =>'$',
    seq                    =>'$',
    map_position           =>'$',
    Annotation             =>'$',
    computational_analysis =>'$',
};


struct Game::seq =>
{
    attr=>'$',
    name=> '$',
    organism=>'$',
    dbxref=>'$',
    description=>'$',
    residues=>'$',
};

struct Game::attr::seq =>
{
    attr=>'$',
    type=>'$',
    focus=> '$',
    md5checksum=>'$',
    version=>'$',
    length=>'$',
    id=>'$',
};

struct Game::dbxref =>
{
    attr  =>'$',
    xref_db  => '$',
    db_xref_id  => '$',
};

struct Game::map_position =>
{
    attr=>'$',
    arm=>'$',
    chromosome=>'$',
    organism=>'$',
    span=> '$',
};

struct Game::attr::map_position=>
{
    attr  =>'$',
    is_elm   =>'$',
    seq      =>'$',
    type     =>'$',
};


struct Game::span =>
{
    attr  =>'$',
    start    =>'$',
    end      => '$',
};



struct Game::feature_set =>
{
    attr            =>'$',
    name            =>'$',
    type            =>'$',
    author          =>'$',
    date            =>'$',
    synonym         =>'$',
    comment         =>'$',
    description     =>'$',
    property        =>'$',
    feature_set     =>'$',
    feature_span    =>'$',
    seq             =>'$',

};

struct Game::result_set =>
{
    attr              =>'$',
    name              =>'$',
    score             =>'$',
    seq_relationship  =>'$',
    output            =>'$',
    result_set        =>'$',
    result_span       =>'$',
};

struct Game::attr::result_set =>
{
    attr              =>'$',
    id                =>'$',
    
};


struct Game::attr::result_span =>
{
    attr              =>'$',
    id                =>'$',
    
};

struct Game::result_span=>
{
    attr             =>'$',
    name             =>'$',
    type             =>'$',
    score            =>'$',
    output           =>'$',
    seq_relationship =>'$',
};

struct Game::output=>
{
    attr             =>'$',
    type             =>'$',
    value            =>'$',
    
};







struct Game::attr::feature_set =>
{
     attr            =>'$',
     id              =>'$',
     problem         =>'$',
     produces_seq    =>'$',
     type            =>'$', 
};



struct Game::date =>
{
    attr            =>'$',
    val             =>'$',
   
};

struct Game::attr::date =>
{
    attr                 =>'$',
    timestamp            =>'$',
    
};

struct Game::comment =>
{
    attr            =>'$',
    text            =>'$',
    person          =>'$',
    date            =>'$',
    internal        =>'$',

};

struct Game::attr::comment =>
{
    attr                =>'$',
    internal            =>'$',
    id                  =>'$',
};


struct Game::property =>
{
    attr                   =>'$',
    type                   =>'$',
    value                  =>'$',
};

struct Game::seq_relationship=>
{
    attr       =>'$',
    span       =>'$',
    alignment  =>'$',
    
}; 

struct Game::computational_analysis=>
{
    attr            =>'$',
    program         =>'$',
    database        =>'$',
    version         =>'$',
    type            =>'$',
    property        =>'$',
    date            =>'$',
    result_set      =>'$',
    result_span     =>'$'
    
};


struct Game::gene=>
{
    attr        =>'$',
    name        =>'$',
    dbxref      =>'$',

};
struct Game::attr::gene=>
{
    attr          =>'$',
    id            =>'$',
    association   =>'$',

};


struct Game::attr::seq_relationship=>
{
    attr               =>'$',
    type               =>'$',
    seq                =>'$',
};

struct Game::feature_span=>
{
    attr             =>'$',
    name             =>'$',
    type             =>'$',
    seq_relationship =>'$',
};

struct Game::attr::feature_span=>
{
    attr           =>'$',
    produces_seq   =>'$',
    type           =>'$',
    id             =>'$'

};

struct Game::Annotation=>
{
    attr         =>'$',
    name         =>'$',
    type         =>'$',
    date         =>'$',
    property     =>'$',
    synonym      =>'$',
    gene         =>'$',
    dbxref       =>'$',
    comment      =>'$',
    description  =>'$',
    feature_set  =>'$',
    feature_span =>'$',
};

my %GameOrder =
(
 Game::game              =>[qw(attr seq map_position Annotation computational_analysis)],
 Game::seq               =>[qw(attr name organism dbxref description residues)],
 Game::attr::seq         =>[qw(attr type focus md5checksum version length id)],
 Game::dbxref            =>[qw(attr xref_db db_xref_id)],
 Game::map_position      =>[qw(attr arm chromosome organism span)],
 Game::span              =>[qw(attr start end)],
 Game::feature_set       =>[qw(attr type name author date synonym comment description property feature_set feature_span seq)],
 Game::attr::feature_set =>[qw(attr id problem produces_seq type)],
 Game::date              =>[qw(attr val)],
 Game::comment           =>[qw(attr text person date internal)],
 Game::attr::comment     =>[qw(attr internal id)],
 Game::property          =>[qw(attr type value)],
 Game::seq_relationship  =>[qw(attr span alignment)],
 Game::gene              =>[qw(attr name dbxref)], 
 Game::attr::gene        =>[qw(attr seq type)],
 
 Game::attr::seq_relationship =>[qw(attr type seq)], 
 Game::feature_span           =>[qw(attr name type seq_relationship)],
 Game::attr::feature_span     =>[qw(attr produces_seq type id)],
 Game::Annotation             =>[qw(attr name type date property synonym gene dbxref comment description feature_set feature_span)],
 Game::computational_analysis =>[qw(attr program database version type property date result_set result_span)],
 Game::attr::feature_set      =>[qw(attr id problem produces_seq type)],
 Game::output                 =>[qw(attr type value)],
 Game::result_span            =>[qw(attr name type score output seq_relationship)],
 Game::attr::result_span      =>[qw(attr id)],
 Game::attr::result_set       =>[qw(attr id)],
 Game::result_set             =>[qw(attr name score seq_relationship output result_set result_span)],
 

);



sub init{
    my ($self,%args)=@_;

  
    %{$self}=%args;
    
    $self->Ind(0);
    return $self;
}



#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#---------------------- SUBS ------------------------------------------------
#----------------------------------------------------------------------------


sub isObj{
    my($self,$obj)=@_;
    my $ref=ref($obj);
    if( $ref=~/^Game/){
	return 1;
    }else{
	return 0;
    }
}

#----------------------------------------------------------------------------
sub keep_order{
    my ($self,$item)=@_;

    my $ref=ref($item);
    my $order=$GameOrder{$ref};
 
   return $order;
}


#------------------my baby-----------------------------------------------------
sub treat{
    my ($self,$item,$elm_nm)=@_;
    

    if($self->isObj($item)){                                                              
        $self->print_start($item);                               
	
	foreach my $key(@{$self->keep_order($item)}){

	    my $value=$item->$key;
	    my $elm_count=$self->get_zero_more($value);
            
	    if($elm_count == 0){
                next;
	    }elsif($elm_count == 1){
		my $nm=$self->get_actual($key);
		$self->treat($value,$nm);
		$self->print_end($key,$value); 
		$self->less_level();
	    }else{
		foreach my $more (@{$value}){
		    my $nm=$self->get_actual($key);
		    $self->treat($more,$nm);
		    $self->print_end($key,$more); 
		    $self->less_level();  
		}
	    }
        }
            
        if($self->Ind() eq 1){
            $self->print_end($item); 
        }   
    }else{
	$self->print_start($elm_nm);
	$self->print_item($item,$nm);
	return;
    }
}


#-----------------------------------------------------------------------
sub get_zero_more{
    my ($self,$item)=@_;    
    
    if(defined($item)){
        if(ref($item) =~/ARRAY/){
            return 2;
        }elsif(ref($item)=~/attr/){
            return 0;
        }
    return 1;
    }else{
        return 0;
    }
    
}


#<sub>-----------------------------------------------------------------------
sub print_item{
    my ($self,$item)=@_;
    my $st;
    $self->add_level();
    my $ind=$self->get_ind();
    $st= $ind."\t".$item."\n";
    $self->add2print($st);
    $self->less_level();
    return $st;
}

#<sub>-----------------------------------------------------------------------
sub add2print{
    my ($self,$item)=@_;
    my $thing=$self->sting();
    print $thing $item;
    return;
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
sub print_attr{
    my ($self,$item,$elm_nm)=@_;
    
    my $x_st=$self->get_ind();
    $x_st.= "<".$elm_nm." ";
    while ( my ($key, $value) = each(%{$item}) ) {
	my $elm_count=$self->get_zero_more($value);
	if($elm_count==0){
	    next;
	}
	my $act=$self->get_actual($key);
	$x_st.= $act."="."\"".$value."\" ";
    }
    $x_st.=">\n";
    $self->add2print($x_st);
    return $x_st;
}


#<sub>-----------------------------------------------------------------------
sub print_start{
    my ($self,$item)=@_;
    
    my $x_st;
    
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
            $x_st=$self->get_ind();
            $x_st.="<".$elm_name.">\n";
            #print $x_st;
	    $self->add2print($x_st);
        }
    }else{
	$elm_name=$item;
	if($elm_name eq "val"){
	   # $self->less_level();
	    return;
	}
	$x_st.=$self->get_ind();
	$x_st.="<".$elm_name.">\n";
	$self->add2print($x_st);
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
    
    #print $x_st;
    $self->add2print($x_st);
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
     my $file=$self->out();
    
    if(defined($file)){
        open FFF,">$file";
	$self->sting(*FFF{IO});
    }else{
	open FFF,STDOUT;
	$self->sting(*STDOUT{IO});
    }
#    print Dumper(FFF);die;
#    print Dumper($self);die;
    
    my $st=  "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
    print FFF $st;
    #close $self->sting;
    #$self->add2print($st);
    #$self->sting(*FFF{IO});
    return *FFF{IO};
}


#<sub>-----------------------------------------------------------------------
sub Stuffer{
    my ($self,$from_h,$to_h,$head)=@_;

    my $targ = $to_h->{$head};
    while ( my ($key, $value) = each(%{$from_h}) ) {
        #print "$key => $value\n";
	$targ->{$key}=$value;
    }

}



#<end>-----------------------------------------------------------------------
1;

