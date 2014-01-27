#-------------------------------------------------------------------------------
#------                            Fasta                               ---------
#-------------------------------------------------------------------------------
package Fasta;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use URI::Escape;
use Carp;

@ISA = qw(
          );

#-------------------------------------------------------------------------------
#------------------------------- Methods ---------------------------------------
#-------------------------------------------------------------------------------
sub new {
	my $class = shift;
	my $fasta = shift;

	my $self = {};

	bless $self, $class;

	$self->fasta($fasta);

	return $self;
}
#-------------------------------------------------------------------------------
sub identifier {
	my $self = shift;
	my $id   = shift;

        if (defined($id)){
                $self->{identifier} = $id;
        }
        elsif (defined($self->{identifier})){
                return $self->{identifier};
        }
	else {
		my ($id) = $self->def() =~ /^>(\S+)/;
		$self->{identifier} = $id;
		return $self->{identifier};
	}
}
#-------------------------------------------------------------------------------
sub fasta {
	my $self  = shift;
	my $fasta = shift;

	if (defined($fasta)){
		$self->{fasta} = $fasta;
	}
	else {
		return $self->{fasta};
	}
}
#-------------------------------------------------------------------------------
sub seq {
        my $self = shift;
        my $seq  = shift;

        if  (defined($seq) && $seq eq 'revComp' && !defined($self->{seq})){
                my @a = split(/\n/, $self->fasta());
                shift(@a);
                $self->{seq} = join('', @a);;
                return revComp($self->{seq});
        }
        if (defined($seq) && $seq eq 'revComp'){
                return revComp($self->{seq});
        }
        elsif (defined($seq)){
                $self->{seq} = $seq;
        }
        elsif (defined($self->{seq})){
                return $self->{seq};
        }
        else {
                my @a = split(/\n/, $self->fasta());
		shift(@a);
                $self->{seq} = join('', @a);;
                return $self->{seq};
        }
}
#-------------------------------------------------------------------------------
sub def {
	my $self = shift;
	my $def  = shift;

        if (defined($def)){
                $self->{def} = $def;
        }
	elsif (defined($self->{def})){
		return $self->{def};
	}
        else {
		my @a = split(/\n/, $self->fasta());
		$self->{def} = shift(@a);
		return $self->{def};
        }
}
#-------------------------------------------------------------------------------
#------------------------------- SUBS ------------------------------------------
#-------------------------------------------------------------------------------
sub getWantedFromMulti {
        my $multiFasta = shift;
        my $wanted     = shift;

        my $fh = new FileHandle;
           $fh->open("$multiFasta");

	my @fastas;
	{
	    local $/ = "\n>";
	    while(my $line = <$fh>){
                $line =~ s/>//;
                $line = ">".$line;
		if (!defined($wanted) || $line =~ /$wanted/){
		    push(@fastas, $line);
                }
	    }
	}
        $fh->close;

	return \@fastas;
}
#-----------------------------------------------------------------------------
sub getDef {
        my $fasta = shift;

	#always work with references
	$fasta = $$fasta while(ref($fasta) eq 'REF');
	my $fasta_ref = (ref($fasta) eq '') ? \$fasta : $fasta;

	my ($def) = $$fasta_ref =~ /(>[^\n\cM]+)/;
	
	return $def;
}
#-----------------------------------------------------------------------------
sub getSeqID {
        my $fasta = shift;

	#always work with references
	$fasta = $$fasta while(ref($fasta) eq 'REF');
	my $fasta_ref = (ref($fasta) eq '') ? \$fasta : $fasta;

	my $def = Fasta::getDef($fasta_ref);
	my $seq_id = def2SeqID($def);

	return $seq_id;
}
#-----------------------------------------------------------------------------
sub def2SeqID {
        my $def = shift;

	my ($seq_id)  = $def =~ /^>(\S+)/;

	return $seq_id;
}
#-----------------------------------------------------------------------------
sub getSafeID {
        my $fasta = shift;

	#always work with references
	$fasta = $$fasta while(ref($fasta) eq 'REF');
	my $fasta_ref = (ref($fasta) eq '') ? \$fasta : $fasta;

	my $seq_id = getSeqID($fasta_ref);
	my $safe_id = SeqID2SafeID($seq_id);

	return $safe_id;
}
#-----------------------------------------------------------------------------
sub seqID2SafeID {
    my $seq_id = shift;
    
    my $safe_id = uri_escape($seq_id,
			     '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\.\+'
			     );
    
    return $safe_id;
}
#-----------------------------------------------------------------------------
sub safeID2SeqID {
    my $seq_id = shift;

    my $safe_id = uri_unescape($seq_id);
    
    return $safe_id;
}
#-------------------------------------------------------------------------------
sub getSeq {
    my $fasta = shift;

    return ${getSeqRef(\$fasta)};
}
#-------------------------------------------------------------------------------
sub getSeqRef {
    my $fasta = shift;

    #always work with references
    $fasta = $$fasta while(ref($fasta) eq 'REF');
    my $fasta_ref = (ref($fasta) eq '') ? \$fasta : $fasta;

    my @fasta = split(/\n/, $$fasta_ref);

    my $seq = '';
    while(my $l = shift(@fasta)){
	chomp($l);
	next if $l =~ /^>/;
	#remove contaminating whitespace
	$l =~ s/\s+//g;	
	$seq .= $l;
    }

    return \$seq;
}
#-------------------------------------------------------------------------------
sub getSeqLength {
    my $fasta = shift;

    #always work with references
    $fasta = $$fasta while(ref($fasta) eq 'REF');
    my $fasta_ref = (ref($fasta) eq '') ? \$fasta : $fasta;

    my $def = getDef($fasta);
    $fasta_ref = fasta2seqRef($fasta_ref);
    my $length = length($$fasta_ref);
    $fasta_ref = seq2fastaRef($def, $fasta_ref);

    return $length;
}

#-------------------------------------------------------------------------------
sub fasta2seqRef {
    my $fasta = shift;

    #always work with references
    $fasta = $$fasta while(ref($fasta) eq 'REF');
    my $fasta_ref = (ref($fasta) eq '') ? \$fasta : $fasta;

    $$fasta_ref =~ s/(>[^\n\cM]+)//; #remove header
    my $length = length($$fasta_ref);
    my $limit = 10000000;
    if($length < $limit){ #weird large regular expression limit
	$$fasta_ref =~ s/[^A-Za-z]//g; #remove whitespace and binary characters
    }
    else{
	my $new= '';
	my $offset = 0;	
	while($offset < $length){
	    my $l = substr($$fasta_ref, $offset, $limit);
	    $l =~ s/[^A-Za-z]//g; #remove whitespace and binary characters
            $new .= $l;
	    $offset += $limit;
        }
        $fasta_ref = \$new;
    }

    return $fasta_ref;
}
#-------------------------------------------------------------------------------
sub fasta2seq {
    my $fasta = shift;
    return ${fasta2seqRef(\$fasta)};
}
#-------------------------------------------------------------------------------
sub seq2fastaRef {
    my $def = shift;
    my $seq = shift;

    #always work with references
    $seq = $$seq while(ref($seq) eq 'REF');
    my $fasta_ref = (ref($seq) eq '') ? \$seq : $seq;

    $def = ">$def" unless($def =~ /^>/);
    my $width = 60;
    my $limit = int(10000000/$width)*$width;
    my $length = length($$fasta_ref);
    if($length < $limit){ #weird large regular expression limit
	$$fasta_ref =~ s/(.{1,$width})/$1\n/g;
    }
    else{
	my $new= '';
	my $offset = 0;	
	while($offset < $length){
	    my $l = substr($$fasta_ref, $offset, $limit);
	    $l =~ s/(.{1,$width})/$1\n/g;
            $new .= $l;
	    $offset += $limit;
        }
        $fasta_ref = \$new;
    }
    $$fasta_ref =~ s/^(.)/$def\n$1/;

    return $fasta_ref;
}
#-------------------------------------------------------------------------------
sub seq2fasta {
    my $def = shift;
    my $seq = shift;
    
    return ${seq2fastaRef($def, \$seq)};
}
#-------------------------------------------------------------------------------
sub ucFasta{
    my $fasta = shift;

    return ${ucFastaRef(\$fasta)};
}
#-------------------------------------------------------------------------------
sub ucFastaRef{
        my $fasta = shift;

	#always work with references
	$fasta = $$fasta while(ref($fasta) eq 'REF');
	my $fasta_ref = (ref($fasta) eq '') ? \$fasta : $fasta;	

	my $def = getDef($fasta_ref);
	$fasta_ref = fasta2seqRef($fasta_ref);
	$$fasta_ref =~ tr/a-z/A-Z/;
	$fasta_ref = seq2fastaRef($def, $fasta_ref);

        return $fasta_ref;
}
#-------------------------------------------------------------------------------
sub revComp {
    my $seq = shift;

    return ${revCompRef(\$seq)};
}
#-------------------------------------------------------------------------------
sub revCompRef {
	my $seq = shift;

	#always work with references
	$seq = $$seq while(ref($seq) eq 'REF');
	my $seq_ref = (ref($seq) eq '') ? \$seq : $seq;	

        my($rc);

        $rc = $$seq_ref;
        $rc =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX
                 /tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
        $rc = reverse scalar ($rc);

        return \$rc;
}
#-------------------------------------------------------------------------------
sub toFasta {
    my $def = shift;
    my $seq = shift;

    return ${toFastaRef($def, \$seq)};
}
#-------------------------------------------------------------------------------
sub toFastaRef {
        my $def = shift;
        my $seq = shift;

	#always work with references
	$seq = $$seq while(ref($seq) eq 'REF');
	my $seq_ref = (ref($seq) eq '') ? \$seq : $seq;

        my $fasta = $def."\n";

	$fasta .=  ${_formatSeq($seq_ref, 60)};

        return \$fasta;
}
#-------------------------------------------------------------------------------
sub _formatSeq {
    my $seq = shift;
    my $l   = shift;
    
    #always work with references
    $seq = $$seq while(ref($seq) eq 'REF');
    my $seq_ref = (ref($seq) eq '') ? \$seq : $seq;
    
    confess "ERROR: Not a SCALAR reference\n" if(ref $seq_ref ne 'SCALAR');
    
    my $f_seq = '';
    for (my $i=0; $i< length($$seq_ref);$i+=$l){
	$f_seq .= substr($$seq_ref, $i, $l)."\n";
    }
    
    return \$f_seq;
}
#-------------------------------------------------------------------------------
sub toBpos {
        my $def = shift;
        my $seq = shift;

	#always work with references
	$seq = $$seq while(ref($seq) eq 'REF');
	my $seq_ref = (ref($seq) eq '') ? \$seq : $seq;

        my $fasta = $def."\n";

	my $bpos = "0 " x length($$seq_ref);

        for (my $i=0; $i< length($bpos);$i+=60){
                $fasta .= substr($bpos, $i, 60)."\n";
        }

        return $fasta;
}

#-------------------------------------------------------------------------------
sub toQual {
	my $def = shift;
	my $seq = shift;

	#always work with references
	$seq = $$seq while(ref($seq) eq 'REF');
	my $seq_ref = (ref($seq) eq '') ? \$seq : $seq;

	$$seq_ref  =~ s/^\s+//;

	my @values = split(/\s+/, $$seq_ref);
        my $fasta  = $def."\n";

	my $j = 0;
        for (my $i=0; $i< @values ;$i++){
		if ($j < 20){
			my $v = $values[$i];

			$fasta .= length($values[$i]) == 1 ? $values[$i]."  " 
			                                   : $values[$i]." ";
			$j++;
		}
		elsif ($j == 20){
			$fasta .= $values[$i]."\n";
			$j = 0;
		}
		else {
			confess "dead in toQual\n";
		}
        }
	$fasta .= "\n" unless $fasta =~ /\n$/;

        return $fasta;
}
#-------------------------------------------------------------------------------
sub shift2 (\@) {
	return splice(@{$_[0]}, 0, 2);
}
#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        print STDERR "Fasta::AutoLoader called for: ",
              "\$self->$call","()\n";
        print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if (defined($arg)){
                $self->{$call} = $arg;
        }
        else {
                return $self->{$call};
        }
}
#----------------------------------------------------------------------------

1;
