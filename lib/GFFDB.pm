#----------------------------------------------------------------------------
#----                               GFFDB                                ---- 
#----------------------------------------------------------------------------
package GFFDB;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use DBI;
use DBD::SQLite;
use File::Temp;
use URI::Escape;
use Bio::Search::Hit::PhatHit::gff3;
use Bio::Search::HSP::PhatHSP::gff3;
use Cwd;
use File::NFSLock;
use FastaSeq;

@ISA = qw(
       );

#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
sub new {
        my ($class, @args) = @_;

        my $self = {};
        bless ($self, $class);

	my $in = shift @args;
	my $dbfile = 'dbfile.db';
	if(ref $in eq 'HASH'){
	    my $CTL_OPT = $in;
	    $dbfile = "$CTL_OPT->{out_base}/$CTL_OPT->{out_name}.db";
	    $CTL_OPT->{_dbfile} = $dbfile;

	    #rebuild database from scratch on force
	    unlink($dbfile) if($CTL_OPT->{force} && !$CTL_OPT->{_multi_chpc});

	    $self->{dbfile} = $dbfile;
	    $self->{last_build} = undef;
	    $self->{next_build} = $CTL_OPT->{out_name}."::1.00";
	    $self->{go_gffdb} = $CTL_OPT->{go_gffdb};

	    my $lock;
	    while(! $lock || ! $lock->maintain(30)){
		$lock = new File::NFSLock($self->{dbfile}, 'EX', 1200, 60);
	    }	    

	    $self->initiate();
	    return $self unless($self->{go_gffdb});
	    
	    $self->add_maker($CTL_OPT->{maker_gff},$CTL_OPT);
	    $self->add_repeat($CTL_OPT->{rm_gff});
	    $self->add_est($CTL_OPT->{est_gff});
	    $self->add_altest($CTL_OPT->{altest_gff});
	    $self->add_protein($CTL_OPT->{protein_gff});
	    $self->add_pred($CTL_OPT->{pred_gff});
	    $self->add_model($CTL_OPT->{model_gff});
	    $self->add_other($CTL_OPT->{other_gff});
	    $self->do_indexing();
	    $lock->unlock;
	}
	elsif(defined $in){
	    $dbfile = $in;
	    $self->{dbfile} = $dbfile;
	    $self->{in_memory} = shift @args;
	    $self->{last_build} = undef;
	    $self->{next_build} = "Build::1.00";
	    $self->{go_gffdb} = 1;

	    my $lock;
	    while(! $lock || ! $lock->maintain(30)){
		$lock = new File::NFSLock($self->{dbfile}, 'EX', 1200, 60);
	    }
	    $self->initiate();
	    $lock->unlock;
	}

	return $self;
}
#-------------------------------------------------------------------------------
sub initiate {
    my $self = shift;

    my $dbfile = $self->{dbfile};

    my $val;
    if($self->{go_gffdb}){
	my $dbh;
	if($self->{in_memory}){
	    $dbh = DBI->connect("dbi:SQLite::memory:", "", "", {sqlite_use_immediate_transaction => 1});
	    $self->{DBH} = $dbh;
	}
	else{
	    $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","",{AutoCommit => 0});
	    $dbh->do(qq{PRAGMA default_synchronous = OFF}); #improve performance
	    $dbh->do(qq{PRAGMA default_cache_size = 10000}); #improve performance
	}
	
	my $tables = $dbh->selectcol_arrayref(qq{SELECT name FROM sqlite_master WHERE type = 'table'});
	if (! grep( /^sources$/, @{$tables})){
	    $dbh->do(qq{CREATE TABLE sources (name TEXT, source TEXT)});
	}
	
	$dbh->commit unless($self->{in_memory});
	$dbh->disconnect unless($self->{in_memory});
    }
    elsif($dbfile && -e $dbfile){
       unlink($dbfile);
    }
 }
#-------------------------------------------------------------------------------
sub do_indexing {
    my $self = shift;

    return unless($self->{go_gffdb});
    my $dbfile = $self->{dbfile};

    my $dbh;
    if($self->{in_memory}){
	$dbh = $self->{DBH};
    }
    else{
	$dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","",{AutoCommit => 0});
	$dbh->do(qq{PRAGMA default_synchronous = OFF}); #improve performance
	$dbh->do(qq{PRAGMA default_cache_size = 10000}); #improve performance
    }
    
    my $tables = $dbh->selectcol_arrayref(qq{SELECT name FROM sqlite_master WHERE type = 'table'});
    my $indices = $dbh->selectcol_arrayref(qq{SELECT name FROM sqlite_master WHERE type = 'index'});
    
    foreach my $table (@{$tables}){
	next if($table eq 'sources');
	unless (grep(/^$table\_inx$/, @{$indices})){
	    $dbh->do("CREATE INDEX $table\_inx ON $table(seqid)");
	}
    }
	
    $dbh->commit unless($self->{in_memory});
    $dbh->disconnect unless($self->{in_memory});
}
#-------------------------------------------------------------------------------
sub add_maker {
    my $self = shift @_;
    my $gff_file = shift @_;
    my %codes = %{shift @_};

    return unless($self->{go_gffdb});

    my @types;
    push(@types, 'repeat_maker')  if($codes{rm_pass});
    push(@types, 'est_maker')     if($codes{est_pass});
    push(@types, 'altest_maker')  if($codes{altest_pass});
    push(@types, 'protein_maker') if($codes{protein_pass});
    push(@types, 'pred_maker')    if($codes{pred_pass});
    push(@types, 'model_maker')   if($codes{model_pass});
    push(@types, 'other_maker')   if($codes{other_pass});    

    my $dbfile = $self->{dbfile};

    my $dbh;
    if($self->{in_memory}){
	$dbh = $self->{DBH};
    }
    else{
	$dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","",{AutoCommit => 0});
	$dbh->do(qq{PRAGMA default_synchronous = OFF}); #improve performance
	$dbh->do(qq{PRAGMA default_cache_size = 10000}); #improve performance
    }
    
    #check to see if tables need to be created, erased, or skipped
    my $tables = $dbh->selectcol_arrayref(qq{SELECT name FROM sqlite_master WHERE type = 'table'});
    @$tables = grep {/\_maker$/} @$tables; #filter out the source table
    
    my %all;
    foreach my $t (@types, @$tables){
	$all{$t}++;
    }
    
    my %skip;
    foreach my $table (keys %all){
	my $source = (grep {/^$table$/} @types) ? $gff_file : 'empty';
	
	if (grep {/^$table$/} @{$tables}){
	    my ($o_source) = $dbh->selectrow_array(qq{SELECT source FROM sources WHERE name = '$table'});
	    
	    if($source ne $o_source){
		$dbh->do(qq{DROP TABLE $table});
		$dbh->do(qq{CREATE TABLE $table (seqid TEXT, source TEXT, parent TEXT, start INT, end INT, line TEXT)});
		$dbh->do(qq{UPDATE sources SET source = '$source' WHERE name = '$table'});
		}
	    else{
		$skip{$table}++;
	    }
	}
	else{
	    $dbh->do(qq{CREATE TABLE $table (seqid TEXT, source TEXT, parent TEXT, start INT, end INT, line TEXT)});
	    $dbh->do(qq{INSERT INTO sources (name, source) VALUES ('$table', '$source')});
	}
    }
    
    #parse gff3
    if(scalar(keys %skip) < @types && $gff_file){
	my @files = split(/\,/, $gff_file);
	open (my $IN, "< $gff_file") or die "ERROR: Could not open file: $gff_file\n" if(@files == 1);
	open ($IN, "$FindBin::Bin/gff3_merge -l -s -n ".join(' ', @files)." |") if(@files > 1);
	my $count = 0;
	my $line;
	while(defined($line = <$IN>)){
	    chomp($line);
	    if($line =~ /^\#\#genome-build maker ([^\s\n\t]+)/){
		my $build = $1;
		my ($id, $count) = split("::", $build);
		next unless($count =~ /^\d+\.\d\d$/);
		
		$self->{last_build} = $build;
		$self->{next_build} = sprintf $id.'::%.2f', $count;
	    }
	    last if ($line =~ /^\#\#FASTA/);
	    next if ($line =~ /^\s*$/);
	    next if ($line =~ /^\#/);
	    
	    #for line with multiple parents
	    my ($parent) = $line =~ /Parent=([^\;\n]+)/;
	    my @parents = split(",", $parent);
	    
	    my @lines;
	    foreach my $p (@parents){
	       $line =~ s/Parent=[^\;\n]+/Parent=$p/;
	       push(@lines, $line);
	    }
	    push(@lines, $line) if(!@lines);

	    foreach my $ln (@lines){
	       my $l = $self->_parse_line(\$ln);
	       my $table;
	       if($l->{source} =~ /^repeatmasker|^blastx\:repeat|^repeatrunner|^repeat_gff\:/i){
		  next if (! $codes{rm_pass});
		  next if ($skip{repeat_maker});
		  $table = 'repeat_maker';
	       }
	       elsif($l->{source} =~ /^blastn|^est2genome|^est_gff\:/i){
		  next if (! $codes{est_pass});
		  next if ($skip{est_maker});
		  $table = 'est_maker';
	       }
	       elsif($l->{source} =~ /^tblastx|^altest_gff\:/i){
		  next if (! $codes{altest_pass});
		  next if ($skip{altest_maker});
		  $table = 'altest_maker';
	       }
	       elsif($l->{source} =~ /^blastx|^protein2genome|^protein_gff\:/i){
		  next if (! $codes{protein_pass});
		  next if ($skip{protein_maker});
		  $table = 'protein_maker';
	       }
	       elsif($l->{source} =~ /^snap\_?|^augustus\_?|^fgenesh\_*?|^genemark\_?|^pred_gff\:/i){
		  next if (! $codes{pred_pass});
		  next if ($skip{pred_maker});
		  $table = 'pred_maker';
	       }
	       elsif($l->{source} =~ /^maker|^model_gff\:/i){
		  next if (! $codes{model_pass});
		  next if ($skip{model_maker});		  
		  $table = 'model_maker';
	       }
	       elsif($l->{source} =~/^\.$/){
		  next;  #this is just the contig line
	       }
	       else{
		  next if (! $codes{other_pass});
		  next if ($skip{other_maker});
		  $table = 'other_maker';
	       }
	       
	       $self->_add_to_db($dbh, $table, $l);
	       if($count == 10000){ #commit every 10000 entries
		  $dbh->commit unless($self->{in_memory});
		  $count = 0;
	       }
	       else{
		  $count++;
	       }
	    }
	}
	close($IN);
    }
    
    #commit changes
    $dbh->commit unless($self->{in_memory});
    $dbh->disconnect unless($self->{in_memory});
}
#-------------------------------------------------------------------------------
sub add_repeat {
    my $self = shift;
    my $gff_file = shift;
    my $table = 'repeat_gff';

    $self->_add_type($gff_file, $table);
}
#-------------------------------------------------------------------------------
sub add_est {
    my $self = shift;
    my $gff_file = shift;
    my $table = 'est_gff';

    $self->_add_type($gff_file, $table);
}
#-------------------------------------------------------------------------------
sub add_altest {
    my $self = shift;
    my $gff_file = shift;
    my $table = 'altest_gff';

    $self->_add_type($gff_file, $table);
}
#-------------------------------------------------------------------------------
sub add_protein {
    my $self = shift;
    my $gff_file = shift;
    my $table = 'protein_gff';

    $self->_add_type($gff_file, $table);
}
#-------------------------------------------------------------------------------
sub add_pred {
    my $self = shift;
    my $gff_file = shift;
    my $table = 'pred_gff';

    $self->_add_type($gff_file, $table);
}

#-------------------------------------------------------------------------------
sub add_model {
    my $self = shift;
    my $gff_file = shift;
    my $table = 'model_gff';

    $self->_add_type($gff_file, $table);
}
#-------------------------------------------------------------------------------
sub add_other {
    my $self = shift;
    my $gff_file = shift;
    my $table = 'other_gff';

    $self->_add_type($gff_file, $table);
}
#-------------------------------------------------------------------------------
sub _add_type {
    my $self = shift;
    my $gff_file = shift;
    my $table = shift;

    return unless($self->{go_gffdb});

    my $dbfile = $self->{dbfile};

    my $dbh;
    if($self->{in_memory}){
        $dbh = $self->{DBH};
    }
    else{
        $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","",{AutoCommit => 0});
        $dbh->do(qq{PRAGMA default_synchronous = OFF}); #improve performance
        $dbh->do(qq{PRAGMA default_cache_size = 10000}); #improve performance
    }
    
    #see if table needs to be created, erased or skipped
    my $source = ($gff_file) ? $gff_file : 'empty';
    my $tables = $dbh->selectcol_arrayref(qq{SELECT name FROM sqlite_master WHERE type = 'table'});
    my $skip = 0;
    if (grep(/^$table$/, @{$tables})){
	my ($o_source) = $dbh->selectrow_array(qq{SELECT source FROM sources WHERE name = '$table'});
	
	if($source ne $o_source){
	    my ($index) = $dbh->selectrow_array(qq{SELECT name FROM sqlite_master WHERE name = '$table\_inx'});
	    $dbh->do(qq{DROP TABLE $table});
	    $dbh->do(qq{CREATE TABLE $table (seqid TEXT, source TEXT, parent TEXT, start INT, end INT, line TEXT)});
	    $dbh->do(qq{UPDATE sources SET source = '$source' WHERE name = '$table'});
	}
	else{
	    $skip = 1;
	}
    }
    else{
	$dbh->do(qq{CREATE TABLE $table (seqid TEXT, source TEXT, parent TEXT, start INT, end INT, line TEXT)});
	$dbh->do(qq{INSERT INTO sources (name, source) VALUES ('$table', '$source')});
    }
    
    #parse gff3
    if(! $skip && $gff_file){
	my @files = split(/\,/, $gff_file);
	open (my $IN, "< $gff_file") or die "ERROR: Could not open file: $gff_file\n" if(@files == 1);
	open ($IN, "$FindBin::Bin/gff3_merge -l -s -n ".join(' ', @files)." |") if(@files > 1);
	my $count = 0;
	while(defined(my $line = <$IN>)){
	    chomp($line);
	    last if ($line =~ /^\#\#FASTA/);
	    next if ($line =~ /^\s*$/);
	    next if ($line =~ /^\#/);
	    
	    #for line with multiple parents
	    my ($parent) = $line =~ /Parent=([^\;\n]+)/;
	    my @parents = split(",", $parent);
	    
	    my @lines;
	    foreach my $p (@parents){
	       $line =~ s/Parent=[^\;\n]+/Parent=$p/;
	       push(@lines, $line);
	    }
	    push(@lines, $line) if(!@lines);

	    foreach my $ln (@lines){
	       my $l = $self->_parse_line(\$ln, $table);
	       $self->_add_to_db($dbh, $table, $l) unless($l->{type} eq 'contig');
	    
	       if($count == 10000){ #commit every 10000 entries
		  $dbh->commit unless($self->{in_memory});
		  $count = 0;
	       }
	       else{
		  $count++;
	       }
	    }
	}
	close($IN);
    }
    
    #commit changes
    $dbh->commit unless($self->{in_memory});
    $dbh->disconnect unless($self->{in_memory});
}
#-------------------------------------------------------------------------------
sub _parse_line{
    my $self = shift;
    my $line = shift;
    my $tag = shift;

    chomp $$line;
    my @data = split(/\t/, $$line);

    foreach my $d (@data){
	$d = uri_escape($d,'\'\"\%');
    }
    
    #add tag to source
    if($tag && $data[1] !~ /^$tag\:/){
	$data[1] = "$tag:$data[1]";
    }
    elsif($data[1] eq 'maker' && !$tag){
	$data[1] = ($self->{last_build}) ?
	    "model_gff:maker_".$self->{last_build} : 'model_gff:maker';
    }

    $data[6] = '+' if($data[6] =~ /^\.$|^0$/); #fixes some repeat entries
    my ($parent) = $data[8] =~ /Parent=([^\n\;]+)/;

    my %l = (seqid  => $data[0],
	     source => $data[1],
	     parent => ($parent || '.'),
	     type   => $data[2], 
	     start  => $data[3], 
	     end    => $data[4], 
	     line   => join("\t", @data)
	    );

    return (\%l);
}
#-------------------------------------------------------------------------------
sub _add_to_db {
    my $self = shift;
    my $dbh = shift;
    my $table = shift;
    my $l = shift;

    $dbh->do("INSERT INTO $table (seqid, source, parent, start, end, line) ".
	     "VALUES (\'".$l->{seqid}."\', \'".$l->{source}."\', \'".$l->{parent}."\', ".
	     $l->{start}.", ".$l->{end}.", \'".$l->{line}."\')"
	    );
}
#-------------------------------------------------------------------------------
sub phathits_on_chunk {
    my $self = shift;
    my $chunk = shift;
    my $seq = shift;
    my $h_type = shift;
    my $seq_len = shift || length_o($seq); #sometimes slow step

    return [] unless($self->{go_gffdb});
    
    my $dbfile = $self->{dbfile};

    my $c_start = $chunk->start;
    my $c_end = $chunk->end;
    my $seqid = $chunk->seqid;

    my $ref1 = [];
    my $ref2 = [];

    my $dbh;
    if($self->{in_memory}){
        $dbh = $self->{DBH};
    }
    else{
        $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","",{AutoCommit => 0});
        $dbh->do(qq{PRAGMA default_synchronous = OFF}); #improve performance
        $dbh->do(qq{PRAGMA default_cache_size = 10000}); #improve performance
    }
    
    my $tables = $dbh->selectcol_arrayref(qq{SELECT name FROM sqlite_master WHERE type = 'table'});
    
    #get gff annotations
    if (grep(/^$h_type\_gff$/, @{$tables})){
	$ref1 = $dbh->selectall_arrayref(qq{SELECT line FROM $h_type\_gff }.
					 qq{WHERE seqid = '$seqid' }.
					 qq{AND start BETWEEN $c_start AND $c_end }.
					 qq{AND parent = '.' });
	my $safe = quotemeta($seqid);
	@$ref1 = grep {$_->[0] !~ /[\t\;]ID=$safe[\;\n]/} @$ref1; #removes contig/chromosome
	
	my %IDs;
	foreach my $row (@$ref1){
	   my $line = $row->[0];
	   if($line =~ /ID=([^\;\n]+)/){
	      $IDs{$1}++;
	   }
	}
	my @check = keys %IDs;

	while(@check){
	   my @subset;
	   if(@check > 300){
	       @subset = @check[0..299];
	       @check = @check[300..$#check];
	   }
	   else{
	       @subset = @check;
	       undef @check;
	   }

	   my $dsn = "parent = '".join("' OR parent = '", @subset)."'";
	   my $ref = $dbh->selectall_arrayref(qq{SELECT line FROM $h_type\_gff }.
					      qq{WHERE seqid = '$seqid' }.
					      qq{AND ( $dsn )});
	   
	   push(@$ref1, @$ref);

	   %IDs = ();
	   foreach my $row (@$ref){
	      my $line = $row->[0];
	      if($line =~ /ID=([^\;\n]+)/){
		 $IDs{$1}++;
	      }
	   }
	   push(@check, keys %IDs);
	}
     }
    
    #get maker annotations
    if (grep(/^$h_type\_maker$/, @{$tables})){
	$ref2 = $dbh->selectall_arrayref(qq{SELECT line FROM $h_type\_maker }.
					 qq{WHERE seqid = '$seqid' }.
					 qq{AND start BETWEEN $c_start AND $c_end }.
					 qq{AND parent = '.' });

	my $safe = quotemeta($seqid);
	@$ref2 = grep {$_->[0] !~ /[\t\;]ID=$safe[\;\n]/} @$ref2; #removes contig/chromosome

	my %IDs;
	foreach my $row (@$ref2){
	   my $line = $row->[0];
	   if($line =~ /ID=([^\;\n]+)/){
	      $IDs{$1}++;
	   }
	}
	my @check = keys %IDs;

	while(@check){
	   my $dsn = "parent = '".join("' OR parent = '", @check)."'";
	   my $ref = $dbh->selectall_arrayref(qq{SELECT line FROM $h_type\_maker }.
					      qq{WHERE seqid = '$seqid' }.
					      qq{AND ( $dsn )});
	   
	   push(@$ref2, @$ref);

	   %IDs = ();
	   foreach my $row (@$ref){
	      my $line = $row->[0];
	      if($line =~ /ID=([^\;\n]+)/){
		 $IDs{$1}++;
	      }
	   }
	   @check = keys %IDs;
	}
    }
    
    $dbh->disconnect unless($self->{in_memory});
    
    my $features = _ary_to_features($ref1, $ref2);
    
    my $structs;
    if($h_type eq 'model'){
	$structs = _get_genes($features, $seq, $seq_len);
    }
    elsif($h_type eq 'repeat'){
	$structs = _get_structs($features, $seq, $seq_len);
    }
    elsif($h_type eq 'est'){
	$structs = _get_structs($features, $seq, $seq_len);
    }
    elsif($h_type eq 'altest'){
	$structs = _get_structs($features, $seq, $seq_len);
    }
    elsif($h_type eq 'protein'){
	$structs = _get_structs($features, $seq, $seq_len);
    }
    elsif($h_type eq 'pred'){
	$structs = _get_genes($features, $seq, $seq_len);
	my $structs2 = _get_structs($features, $seq, $seq_len);
	push(@$structs, @$structs2); 
    }
    elsif($h_type eq 'other'){
	die "ERROR: Can not build phathits for type: \'other\'\n";
    }
    else{
	die "ERROR: no recognized type in GFFDB::phathits_on_chunk\n";
    }
    
    my @phat_hits;    
    foreach my $s (@{$structs}){
	next unless ($c_start <= $s->{start} && $s->{start} <= $c_end);
	push(@phat_hits, @{_load_hits($s, $seq, $seq_len)});
    }

    return \@phat_hits;
}
#-------------------------------------------------------------------------------
sub lines_for_chunk {
    my $self = shift;
    my $chunk = shift;
    my $h_type = shift;

    return [] unless($self->{go_gffdb});
    
    my $dbfile = $self->{dbfile};

    my $c_start = $chunk->start;
    my $c_end = $chunk->end;
    my $seqid = $chunk->seqid;

    my $dbh;
    if($self->{in_memory}){
        $dbh = $self->{DBH};
    }
    else{
        $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","",{AutoCommit => 0});
        $dbh->do(qq{PRAGMA default_synchronous = OFF}); #improve performance      
        $dbh->do(qq{PRAGMA default_cache_size = 10000}); #improve performance     
    }
    
    my $tables = $dbh->selectcol_arrayref(qq{SELECT name FROM sqlite_master WHERE type = 'table'});
    
    #get gff annotations
    my @lines;
    if (grep(/^$h_type\_gff$/, @{$tables})){
	my $ref1 = $dbh->selectall_arrayref(qq{SELECT line FROM $h_type\_gff }.
					    qq{WHERE seqid = '$seqid' }.
					    qq{AND start BETWEEN $c_start AND $c_end }.
					    qq{AND parent = '.' });
	my $safe = quotemeta($seqid);
	@$ref1 = grep {$_->[0] !~ /[\t\;]ID=$safe[\;\n]/} @$ref1; #removes contig/chromosome
	
	my %IDs;
	foreach my $row (@$ref1){
	   my $line = $row->[0];
	   push(@lines, $line);
	   if($line =~ /ID=([^\;\n]+)/){
	      $IDs{$1}++;
	   }
	}
	my @check = keys %IDs;

	while(@check){
	   my $dsn = "parent = '".join("' OR parent = '", @check)."'";
	   my $ref = $dbh->selectall_arrayref(qq{SELECT line FROM $h_type\_gff }.
					      qq{WHERE seqid = '$seqid' }.
					      qq{AND ( $dsn )});
	   
	   %IDs = ();
	   foreach my $row (@$ref){
	      my $line = $row->[0];
	      push(@lines, $line);
	      if($line =~ /ID=([^\;\n]+)/){
		 $IDs{$1}++;
	      }
	   }
	   @check = keys %IDs;
	}
     }
    
    #get maker annotations
    if (grep(/^$h_type\_maker$/, @{$tables})){
	my $ref2 = $dbh->selectall_arrayref(qq{SELECT line FROM $h_type\_maker }.
					    qq{WHERE seqid = '$seqid' }.
					    qq{AND start BETWEEN $c_start AND $c_end }.
					    qq{AND parent = '.' });

	my $safe = quotemeta($seqid);
	@$ref2 = grep {$_->[0] !~ /[\t\;]ID=$safe[\;\n]/} @$ref2; #removes contig/chromosome

	my %IDs;
	foreach my $row (@$ref2){
	   my $line = $row->[0];
	   push(@lines, $line);
	   if($line =~ /ID=([^\;\n]+)/){
	      $IDs{$1}++;
	   }
	}
	my @check = keys %IDs;

	while(@check){
	   my $dsn = "parent = '".join("' OR parent = '", @check)."'";
	   my $ref = $dbh->selectall_arrayref(qq{SELECT line FROM $h_type\_maker }.
					      qq{WHERE seqid = '$seqid' }.
					      qq{AND ( $dsn )});	   
	   %IDs = ();
	   foreach my $row (@$ref){
	      my $line = $row->[0];
	      push(@lines, $line);
	      if($line =~ /ID=([^\;\n]+)/){
		 $IDs{$1}++ if($line =~ /ID=([^\;\n]+)/);
	      }
	   }
	   @check = keys %IDs;
	}
    }
    
    $dbh->disconnect unless($self->{in_memory});

    return \@lines;
}
#-------------------------------------------------------------------------------
sub phathits_on_contig {
    my $self = shift;
    my $seqid = shift;
    my $seq = shift;
    my $h_type = shift;
    my $seq_len = shift || length_o($seq); #sometimes slow step

    return [] unless($self->{go_gffdb});

    my $ref1 = [];
    my $ref2 = [];

    my $dbfile = $self->{dbfile};
    my $dbh;
    if($self->{in_memory}){
        $dbh = $self->{DBH};
    }
    else{
        $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","",{AutoCommit => 0});
        $dbh->do(qq{PRAGMA default_synchronous = OFF}); #improve performance      
        $dbh->do(qq{PRAGMA default_cache_size = 10000}); #improve performance     
    }
    
    my $tables = $dbh->selectcol_arrayref(qq{SELECT name FROM sqlite_master WHERE type = 'table'});
    
    #get gff annotations
    if (grep(/^$h_type\_gff$/, @{$tables})){
	$ref1 = $dbh->selectall_arrayref(qq{SELECT line FROM $h_type\_gff WHERE seqid = '$seqid'});
    }
    
    #get maker annotations
    if (grep(/^$h_type\_maker$/, @{$tables})){
	$ref2 = $dbh->selectall_arrayref(qq{SELECT line FROM $h_type\_maker WHERE seqid = '$seqid'});
    }
    
    $dbh->disconnect unless($self->{in_memory});
    
    my $features = _ary_to_features($ref1, $ref2);
    
    my $structs;
    if($h_type eq 'model'){
	$structs = _get_genes($features, $seq, $seq_len);
    }
    elsif($h_type eq 'repeat'){
	$structs = _get_structs($features, $seq, $seq_len);
    }
    elsif($h_type eq 'est'){
	$structs = _get_structs($features, $seq, $seq_len);
    }
    elsif($h_type eq 'altest'){
	$structs = _get_structs($features, $seq, $seq_len);
    }
    elsif($h_type eq 'protein'){
	$structs = _get_structs($features, $seq, $seq_len);
    }
    elsif($h_type eq 'pred'){
	$structs = _get_structs($features, $seq, $seq_len);
	my $structs2 = _get_genes($features, $seq, $seq_len);
	push(@{$structs}, @{$structs2});
    }
    elsif($h_type eq 'other'){
	die "ERROR: Can not build phathits for type: \'other\'\n";
    }
    else{
	die "ERROR: no recognized type in GFFDB::phathits_on_contig\n";
    }
    
    my @phat_hits;    
    foreach my $s (@{$structs}){
	push(@phat_hits, @{_load_hits($s, $seq, $seq_len)});
    }
    
    return \@phat_hits;
}
#-------------------------------------------------------------------------------
sub get_existing_gene_names {
    my $self = shift;
    my $seqid = shift;

    return {} unless($self->{go_gffdb});

    my %names;

    my $dbfile = $self->{dbfile};
    my $dbh;
    if($self->{in_memory}){
        $dbh = $self->{DBH};
    }
    else{
        $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","",{AutoCommit => 0});
        $dbh->do(qq{PRAGMA default_synchronous = OFF}); #improve performance      
        $dbh->do(qq{PRAGMA default_cache_size = 10000}); #improve performance     
    }
    
    my $tables = $dbh->selectcol_arrayref(qq{SELECT name FROM sqlite_master WHERE type = 'table'});
    
    #---get names for genes
    my $h_type = 'model';
    
    #get gff annotations
    my $ref1 = [];
    if (grep(/^$h_type\_gff$/, @{$tables})){
	$ref1 = $dbh->selectall_arrayref(qq{SELECT line FROM $h_type\_gff WHERE seqid = '$seqid'});
    }
    
    #get maker annotations
    my $ref2 = [];
    if (grep(/^$h_type\_maker$/, @{$tables})){
	$ref2 = $dbh->selectall_arrayref(qq{SELECT line FROM $h_type\_maker WHERE seqid = '$seqid'});
    }
    
    my $safe = quotemeta($seqid);
    my $features = _ary_to_features($ref1, $ref2);
    foreach my $f (@$features){
	my $tag = $f->primary_tag();
	
	if ($tag eq 'gene') {    
	    my $id   = _get_annotation($f,'ID');
	    my $name = _get_annotation($f,'Name');
	    
	    $name = $id if($name eq '');
	    
	    #get old names
	    ($name) = $name =~ /^([^\s\t\n]+)/;
	    $names{$name}++;
	    
	    #get old maker cluster ids
	    my ($c_id) = $name =~ /$safe\-[\-]+\-gene\-(\d+\.*\d*)/;
	    $names{$c_id}++ if(defined $c_id);
	}
    }
    
    #---get names for preds
    $h_type = 'pred';
    
    #get gff annotations
    $ref1 = [];
    if (grep(/^$h_type\_gff$/, @{$tables})){
	$ref1 = $dbh->selectall_arrayref(qq{SELECT line FROM $h_type\_gff WHERE seqid = '$seqid'});
    }
    
    #get maker annotations
    $ref2 = [];
    if (grep(/^$h_type\_maker$/, @{$tables})){
	$ref2 = $dbh->selectall_arrayref(qq{SELECT line FROM $h_type\_maker WHERE seqid = '$seqid'});
    }
    
    $features = _ary_to_features($ref1, $ref2);
    foreach my $f (@$features){
	my $tag = $f->primary_tag();
	
	if ($tag eq 'match' || $tag eq 'gene') { 
	    my $id   = _get_annotation($f,'ID');   
	    my $name = _get_annotation($f,'Name');
	    
	    $name = $id if($name eq '');
	    
	    #get old names
	    ($name) = $name =~ /^([^\s\t\n]+)/;
	    $name =~ s/\-mRNA\-\d+$//;
	    $names{$name}++;
	    
	    #get old maker cluster ids
	    my ($c_id) = $name =~ /$safe\-[\-]+\-gene\-(\d+\.*\d*)/;
	    $names{$c_id}++ if(defined $c_id);
	}
    }
    
    $dbh->disconnect unless($self->{in_memory});

    return \%names;
}
#-------------------------------------------------------------------------------
sub last_build {
    my $self = shift;

    return $self->{last_build} || undef;
}
#-------------------------------------------------------------------------------
sub next_build {
    my $self = shift;

    return $self->{next_build} || undef;
}
#-------------------------------------------------------------------------------
sub dbfile {
    my $self = shift;

    return $self->{dbfile};
}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
#based on Bio::Tools::GFF code
sub _ary_to_features{
   my @features;

   while(my $ary_ref = shift){
      foreach my $row (@{$ary_ref}){
	 my $string = $row->[0];
	 my $feat = Bio::SeqFeature::Generic->new();
	 
	 my ($seqname, $source, $primary, $start, $end, 
	     $score, $strand, $frame, $groups) = split(/\t/, $string);
	 
	 if ( ! defined $frame ) {
	    die "ERROR: [$string] does not look like GFF3\n";
	 }
	 
	 $feat->seq_id($seqname);
	 $feat->source_tag($source);
	 $feat->primary_tag($primary);
	 $feat->start($start);
	 $feat->end($end);
	 $feat->frame($frame);
	 $feat->score($score) unless ( $score eq '.' );
	 if ( $strand eq '-' ) { $feat->strand(-1); }
	 if ( $strand eq '+' ) { $feat->strand(1); }
	 if ( $strand eq '.' ) { $feat->strand(0); }
	 my @groups = split(/\s*;\s*/, $groups);
	 
	 for my $group (@groups) {
	    my ($tag,$value) = split (/=/,$group);
	    $tag             = uri_unescape($tag);
	    my @values       = map {uri_unescape($_)} split(/,/,$value);
	    for my $v ( @values ) {  $feat->add_tag_value($tag,$v); }
	 }
	 push (@features, $feat);
      }
   }

   return \@features;
}
#-------------------------------------------------------------------------------
sub _load_hits {
    my $g   = shift;
    my $seq = shift;
    my $seq_len = shift || length_o($seq);

    my $gene_id   = $g->{id};
    my $gene_name = $g->{name};

    #strip off unrecognized gene attributes for storage
    my @anns;
    if(@anns = $g->{f}->annotation->get_all_annotation_keys()){ #sometimes bioperl does this, version?
	@anns = grep {!/^ID$|^Name$|^Target$|^Parent$|^_AED$|^_eAED$|^_QI$/} @anns;
	foreach my $ann (@anns){
	    my @list = $g->{f}->annotation->get_Annotations();
	    @list = map {$_->value()} @list;
	    $ann = $ann.'='.join(',', @list);
	}
    }
    elsif(@anns = $g->{f}->get_all_tags){ #sometimes bioperl does this, version?
	@anns = grep {!/^ID$|^Name$|^Target$|^Parent$|^_AED$|^_eAED$|^_QI$/} @anns;
        foreach my $ann (@anns){
            my @list = $g->{f}->get_tag_values($ann);
            $ann = $ann.'='.join(',', @list);
        }
    }
    my $gene_attrib = join(';', @anns);
        
    my @phat_hits;
    foreach my $t (@{$g->{mRNAs}}){
	my $tran_id   = $t->{id};
	my $tran_name = $t->{name};
	   ($tran_name) = $tran_name =~ /(^[^\s]+)/;
	    $tran_name =~ s/(mRNA\-\d+)\-AED\:.*/$1/; #for backwards compatability to older maker

	my $description = "g_name=$gene_name;g_id=$gene_id;t_name=$tran_name;t_id=$tran_id" unless(! $gene_name);
	
	my $f = new Bio::Search::Hit::PhatHit::gff3('-name'         => $tran_name,
						    '-description'  => $description,
						    '-algorithm'    => $t->{f}->source_tag,
						    '-length'       => length($t->{seq}),
						    '-score'        => $t->{f}->score || '.',
						    );

	$f->{gene_id} = $gene_id unless(! $gene_id);
	$f->{gene_name} = $gene_name unless(! $gene_name);
	$f->{_tran_name} = $tran_name unless(! $gene_name);
	$f->{_tran_id}   = $tran_id unless(! $tran_id);
	$f->{maker_qi}  = $t->{maker_qi} unless(! $t->{maker_qi});
	$f->{gene_attrib} = $gene_attrib unless(! $gene_attrib);

	#strip off unrecognized hit attributes for storage
	my @anns;
	if(@anns = $t->{f}->annotation->get_all_annotation_keys()){ #sometimes bioperl does this, version?
	    @anns = grep {!/^ID$|^Name$|^Target$|^Parent$|^_AED$|^_eAED$|^_QI$/} @anns;
	    foreach my $ann (@anns){
		my @list = $t->{f}->annotation->get_Annotations();
		@list = map {$_->value()} @list;
		$ann = $ann.'='.join(',', @list);
	    }
	}
	elsif(@anns = $t->{f}->get_all_tags){ #sometimes bioperl does this, version?
	    @anns = grep {!/^ID$|^Name$|^Target$|^Parent$|^_AED$|^_eAED$|^_QI$/} @anns;
	    foreach my $ann (@anns){
		my @list = $t->{f}->get_tag_values($ann);
		$ann = $ann.'='.join(',', @list);
	    }
	}
	$f->{-attrib} = join(';', @anns) if(@anns);


	my $type = $t->{f}->primary_tag;
	$f->{transcript_type}=$type;

	$f->queryLength(length($t->{seq}));

	#this are CDS entries but no exons!!!
	if((!$t->{exons} || !@{$t->{exons}}) && $t->{cdss} && @{$t->{cdss}}){
	    $t->{exons} = $t->{cdss};
	    $t->{seq} = $t->{cds_seq};
	}

	if(defined $t->{cdss} && @{$t->{cdss}}){
	    my ($t_offset, $t_end) = _get_t_offset_and_end($t);

	    if($t_offset != -1){ #only happens on bad CDS entries
		$f->{translation_offset} = $t_offset;
		$f->{translation_end}    = $t_end;
		my $cdss = _load_cdss($t, $seq, $seq_len);
		$f->{cdss} = $cdss;
	    }
	    else{
		warn "WARNING: Problem cause by bad CDS entries in GFF3 file for ".
		    $t->{name}."\n".
		    "Maker will just figure out a new CDS entry internally\n\n";
	    }
	}

	#$f->{seq} = $t->{seq};

	my $hsps = _load_hsps($t, $seq, $seq_len);

	foreach my $hsp (@{$hsps}){
	    $f->add_hsp($hsp);
	}

	push(@phat_hits, $f);
    }

    return \@phat_hits;
}
#-------------------------------------------------------------------------------
sub _get_t_offset_and_end {
    my $t  = shift;
                
    my $t_seq = $t->{seq};
    my $c_seq = $t->{cds_seq};
                
    my $t_offset = index($t_seq, $c_seq);
                
    my $t_end = $t_offset + length($c_seq) +1;
                
    warn "WARNING: Problem in GFFDB::_get_t_offset_and_end\n" if $t_offset == -1;

    return ($t_offset, $t_end);
}
#-------------------------------------------------------------------------------
sub _load_cdss {
    my $t   = shift;
    my $seq = shift;
    my $seq_len = shift || length_o($seq);

    my @hsps;
    my $hit_start = 1;
    my $hit_strand = 1;
    my $hit_name = $t->{name};
    ($hit_name) = $hit_name =~ /(^[^\s]+)/;
    $hit_name =~ s/(mRNA\-\d+)\-AED\:.*/$1/; #for backward compatability to old maker
    
    foreach my $e (@{$t->{cdss}}){
	my @args;
	my $hit_end = $e->{f}->end - $e->{f}->start + $hit_start;

	my $value = _get_annotation($e->{f}, 'Target');
	if($value ne ''){
	    my @dats = split(/\s/, $value);

	    $hit_name  = $dats[0];
	    $hit_start = $dats[1];
	    $hit_end   = $dats[2];
	    $hit_strand = (defined ($dats[3]) && $dats[3] eq '-') ? -1 : 1;
	}

	push(@args, '-query_start');
	push(@args, $e->{f}->start);

	push(@args, '-query_seq');
	push(@args, $e->{seq});

	push(@args, '-score');
	push(@args, $e->{f}->score);
	
	push(@args, '-homology_seq');
	push(@args, $e->{seq});
	
	push(@args, '-hit_start');
	push(@args, $hit_start);
	
	push(@args, '-hit_seq');
	push(@args, $e->{seq});
	
	push(@args, '-hsp_length');
	push(@args, $e->{f}->end - $e->{f}->start + 1);
	
	push(@args, '-identical');
	push(@args, $e->{f}->end - $e->{f}->start + 1);
	
	push(@args, '-hit_length');
	push(@args, $e->{f}->end - $e->{f}->start + 1);
	
	push(@args, '-query_name');
	push(@args, $e->{f}->seq_id);
	
	push(@args, '-algorithm');
	push(@args, $e->{f}->source_tag);
	
	push(@args, '-bits');
	push(@args, 2*($e->{f}->end - $e->{f}->start + 1));
	
	push(@args, '-evalue');
	push(@args, 0.0);
	
	push(@args, '-pvalue');
	push(@args, 0.0);
	
	push(@args, '-query_length');
	push(@args, $seq_len);

	push(@args, '-query_end');
	push(@args, $e->{f}->end);

	push(@args, '-conserved');
	push(@args, length($e->{seq}));

	push(@args, '-hit_name');
	push(@args, $hit_name);

	push(@args, '-hit_end');
	push(@args, $hit_end);

	push(@args, '-query_gaps');
	push(@args, 0);

	push(@args, '-hit_gaps');
	push(@args, 0);

	my $hsp = new Bio::Search::HSP::PhatHSP::gff3(@args);
	   $hsp->queryName($e->{f}->seq_id);
	#-------------------------------------------------
	# setting strand because bioperl is all messed up!
	#------------------------------------------------
	if ($e->{f}->strand == 1 ){
	    $hsp->{_strand_hack}->{query} = 1;
	    $hsp->{_strand_hack}->{hit}   = $hit_strand;
	}
	else {
	    $hsp->{_strand_hack}->{query} = -1;
	    $hsp->{_strand_hack}->{hit}   = $hit_strand;
	}

	$hit_start += $e->{f}->end - $e->{f}->start + 1;

	push(@hsps, $hsp);
    }

    return \@hsps;
}
#-------------------------------------------------------------------------------
sub _load_hsps {
    my $t   = shift;
    my $seq = shift;
    my $seq_len = shift || length_o($seq);

    my @hsps;
    my $hit_start = 1;
    my $hit_strand = 1;
    my $hit_name = $t->{name};
    ($hit_name) = $hit_name =~ /(^[^\s]+)/;
    $hit_name =~ s/(mRNA\-\d+)\-AED\:.*/$1/; #for backward compatability to old maker

    #catch error caused by malformed GFF3
    if(! $t->{exons}->[0] || ! $t->{exons}->[0]->{f}){
	die "ERROR: Failed on $hit_name\n".
	    "Check your input GFF3 file for errors! (from GFFDB)\n";
    }

    #added 3/19/2009
    #check for single and double base pair overhangs
    if($t->{exons}->[0]->{f}->source_tag =~ /^snap|^augustus|^fgenesh|^genemark/){
	my $features = $t->{exons};
	@{$features} = sort {$a->{f}->start <=> $b->{f}->start} @{$features};
	my $length = 0;
	foreach my $e (@{$features}){
	    $length += abs($e->{f}->end - $e->{f}->start) + 1;
	}
	
	my $overhang = $length % 3;
	if($overhang != 0){
	    if($features->[0]->{f}->strand == 1){
		my $last = $features->[-1];
		my $l_length = abs($last->{f}->end - $last->{f}->start) + 1;
		
		while($l_length <= $overhang){
		    pop(@{$features});
		    $overhang -= $l_length;
		    $last = $features->[-1];
		    $l_length = abs($last->{f}->end - $last->{f}->start) + 1;
		}
		
		$last->{f}->end($last->{f}->end - $overhang);
	    }
	    elsif($features->[0]->{f}->strand == -1){
		my $last = $features->[0];
		my $l_length = abs($last->{f}->end - $last->{f}->start) + 1;
		
		while($l_length <= $overhang){
		    shift(@{$features});
		    $overhang -= $l_length;
		    $last = $features->[0];
		    $l_length = abs($last->{f}->end - $last->{f}->start) + 1;
		}
		
		$last->{f}->start($last->{f}->start + $overhang);
	    }
	    else{
		die "FATAL: No exon strand in Widget::snap\n";
	    }
	}
    }

    #build hsps
    foreach my $e (@{$t->{exons}}){
	my @args;
	my $hit_end = $e->{f}->end - $e->{f}->start + $hit_start;

	my $value = _get_annotation($e->{f}, 'Target');
        if($value ne ''){
            my @dats = split(/\s/, $value);

            $hit_name  = $dats[0];
            $hit_start = $dats[1];
            $hit_end   = $dats[2];
            $hit_strand = (defined ($dats[3]) && $dats[3] eq '-') ? -1 : 1;
        }

	push(@args, '-query_start');
	push(@args, $e->{f}->start);

	push(@args, '-query_seq');
	push(@args, $e->{seq});

	push(@args, '-score');
	push(@args, $e->{f}->score);

	push(@args, '-homology_seq');
	push(@args, $e->{seq});

	push(@args, '-hit_start');
	push(@args, $hit_start);

	push(@args, '-hit_seq');
	push(@args, $e->{seq});

	push(@args, '-hsp_length');
	push(@args, $e->{f}->end - $e->{f}->start + 1);

	push(@args, '-identical');
	push(@args, $e->{f}->end - $e->{f}->start + 1);

	push(@args, '-hit_length');
	push(@args, $e->{f}->end - $e->{f}->start + 1);

	push(@args, '-query_name');
	push(@args, $e->{f}->seq_id);

	push(@args, '-algorithm');
	push(@args, $e->{f}->source_tag);

	push(@args, '-bits');
	push(@args, 2*($e->{f}->end - $e->{f}->start + 1));
	
	push(@args, '-evalue');
	push(@args, 0.0);

	push(@args, '-pvalue');
	push(@args, 0.0);
	
	push(@args, '-query_length');
	push(@args, $seq_len);

	push(@args, '-query_end');
	push(@args, $e->{f}->end);

	push(@args, '-conserved');
	push(@args, length($e->{seq}));
	
	push(@args, '-hit_name');
	push(@args, $hit_name);
	
	push(@args, '-hit_end');
	push(@args, $hit_end);
	
	push(@args, '-query_gaps');
	push(@args, 0);

	push(@args, '-hit_gaps');
	push(@args, 0);

	#strip off unrecognized hit attributes for storage
	my @anns;
	if(@anns = $e->{f}->annotation->get_all_annotation_keys()){ #sometimes bioperl does this, version?
	    @anns = grep {!/^ID$|^Name$|^Target$|^Parent$|^_AED$|^_eAED$|^_QI$/} @anns;
	    foreach my $ann (@anns){ 
		my @list = $e->{f}->annotation->get_Annotations();
		@list = map {$_->value()} @list;
		$ann = $ann.'='.join(',', @list);
	    }
	}
	elsif(@anns = $e->{f}->get_all_tags){ #sometimes bioperl does this, version?
	    @anns = grep {!/^ID$|^Name$|^Target$|^Parent$|^_AED$|^_eAED$|^_QI$/} @anns;
	    foreach my $ann (@anns){
		my @list = $e->{f}->get_tag_values($ann);
		$ann = $ann.'='.join(',', @list);
	    }
	}
	my $attrib = join(';', @anns) if(@anns);
	push(@args, '-attrib');
	push(@args, $attrib);

	my $hsp = new Bio::Search::HSP::PhatHSP::gff3(@args);
	   $hsp->queryName($e->{f}->seq_id);
	#-------------------------------------------------
	# setting strand because bioperl is all messed up!
	#------------------------------------------------
	if ($e->{f}->strand == 1 ){
	    $hsp->{_strand_hack}->{query} = 1;
	    $hsp->{_strand_hack}->{hit}   = $hit_strand;
	}
	else {
	    $hsp->{_strand_hack}->{query} = -1;
	    $hsp->{_strand_hack}->{hit}   =  $hit_strand;
	}

	$hit_start += $e->{f}->end - $e->{f}->start + 1;

	push(@hsps, $hsp);
    }

    return \@hsps;
}

#-------------------------------------------------------------------------------
sub _get_genes {
    my $features = shift;
    my $seq    = shift;
    my $seq_len = shift || length_o($seq);

    my $exons = _grab(['exon'], $features);
    my $cdss  = _grab(['CDS'],  $features);
    my $mRNAs = _grab(['mRNA'], $features);
    my $UTRs = _grab(['five_prime_UTR', 'three_prime_UTR'], $features);

    foreach my $p_id (keys %{$mRNAs}){
	for (my $i = 0; $ i < @{$mRNAs->{$p_id}}; $i++) {
	    my $f  = $mRNAs->{$p_id}->[$i]->{f};
	    my $id = $mRNAs->{$p_id}->[$i]->{id}; 

	    $mRNAs->{$p_id}->[$i]->{exons} = $exons->{$id};
	    $mRNAs->{$p_id}->[$i]->{cdss}  = $cdss->{$id};
	    $mRNAs->{$p_id}->[$i]->{maker_qi} = _get_maker_qi($mRNAs->{$p_id}->[$i]);
	}
    }

    my @genes;
    foreach my $f (@{$features}){
	my $tag = $f->primary_tag();

	if ($tag eq 'gene') {	    
	    my $id=_get_annotation($f,'ID');
	    my $name=_get_annotation($f,'Name');

	    $name = $id if ($name eq '');

	    #take care of wormbases incorrect parentage
	    if(exists $exons->{$id}){
		foreach my $mRNA (@{$mRNAs->{$id}}){
		    my $t_id = $mRNA->{id};
		    $mRNA->{exons} = _fix_wormbase($exons->{$id}, $cdss->{$t_id}, $UTRs->{$t_id});
		}
	    }
	    
	    push(@genes, {'f'       => $f,
			  'mRNAs'   => $mRNAs->{$id},
			  'part_of' => [],
			  'id'      => $id,
			  'name'    => $name,
			  'start'   => $f->start,
			  'end'     => $f->end,
		         }
		);
	}
    }

    my @valid_genes;
    for my $gene (@genes) {
	push @valid_genes, $gene if _validate_gene($gene);
    }

    _load_seqs(\@valid_genes, $seq, $seq_len);

    return (\@valid_genes);
}
#-------------------------------------------------------------------------------
#try and discover the exon parantage for wormbase entries
sub _fix_wormbase {
    my $exons = shift || [];
    my $cdss = shift || [];
    my $UTRs = shift || [];

    my @keepers;
    foreach my $exon (@$exons){
	my $e = $exon->{f};
	my $eB = $e->start;
	my $eE = $e->end;

	my $ok = 0;
	my $okB = 0;
	my $okE = 0;

	#check if exon goes with this CDS
	foreach my $piece (@$cdss){
	    my $p = $piece->{f};
	    my $pB = $p->start;
	    my $pE = $p->end;

	    if($pB == $eB && $pE == $eE){
		$ok = 1;
	    }
	    elsif($pB == $eB && $pE != $eE){
		$okB = 1;
	    }
	    elsif($pB != $eB && $pE == $eE){
		$okE = 1;
	    }

	    last if($ok);
	}

	if($ok){
	    push(@keepers, $exon);
	    next;
	}

	#check if exon goes with or is completed by this UTR
	foreach my $piece (@$UTRs){
            my $p = $piece->{f};
            my $pB = $p->start;
            my $pE = $p->end;

            if($pB == $eB && $pE == $eE){
                $ok = 1;
            }
            elsif($pB == $eB && $pE != $eE){
		$okB = 1;
                $ok = 1 if($okE);
            }
            elsif($pB != $eB && $pE == $eE){
		$okE = 1;
                $ok = 1 if($okB);
            }

	    last if($ok);
        }

	if($ok){
            push(@keepers, $exon);
            next;
        }
    }

    return \@keepers;
}
#-------------------------------------------------------------------------------
#try too build a sructure similar to that producd by _get_genes
#but for non gene hits in the gff3
sub _get_structs {
    my $features = shift;
    my $seq    = shift;
    my $seq_len = shift || length_o($seq);

    my @bases;
    my %index;
    foreach my $f (@{$features}){
	my $tag_t = $f->primary_tag();

	next if($tag_t =~ /^gene$|^mRNA$|^exon$|^CDS$/);

	my $id = _get_annotation($f,'ID');
	my $name = _get_annotation($f, 'Name');	
	my $p_ids = _get_p_ids($f);

	$name = $id if($name eq '');

	my $struct = { f        => $f,
		       id       => $id,
		       name     => $name,
		       part_of  => $p_ids
		     };

	if (! @{$p_ids}){
	    push(@bases, $struct);
	}
	else{
	    foreach my $p_id (@{$p_ids}){
		push(@{$index{$p_id}}, $struct);
	    }
	}
    }

    foreach my $b (@bases){
	if (exists $index{$b->{id}}){
	    $b->{exons} = $index{$b->{id}};
	}
	else{
	    push(@{$b->{exons}}, { f       => $b->{f},
				   id      => $b->{id},
				   name    => $b->{name},
				   part_of => [$b->{id}]
				 }
		);
	}
    }
    
    #keep compatible with old code, gene based hits
    my @genes;
    foreach my $b (@bases){
	push(@genes, { f        => $b->{f},
		       mRNAs    => [$b],
		       part_of  => [],
		       id       => '',
		       name     => '',
		       start    => $b->{f}->start,
		       end      => $b->{f}->end
		     }
	     );
    }

    _load_seqs(\@genes, $seq, $seq_len);

    return (\@genes);
}
#-------------------------------------------------------------------------------
sub _grab {
    my $types    = shift;
    my $features = shift;

    my %booty;

    for my $type (@{$types}) {
	foreach my $f (@{$features}){
	    my $tag_t = $f->primary_tag();

	    if ($tag_t eq $type) {
		my $id = _get_annotation($f,'ID');
		my $name = _get_annotation($f, 'Name');
		my $p_ids = _get_p_ids($f);

		$name = $id if($name eq '');

		foreach my $p_id (@{$p_ids}){
		    push(@{$booty{$p_id}}, {f        => $f,
					    id       => $id,
					    name     => $name,,
					    part_of  => $p_ids,
					});
		}
	    }
	}
    }
    return \%booty;
}
#-------------------------------------------------------------------------------
sub _get_annotation {
    my ($f,$type)=@_;;
    my $annotation_collection = $f->annotation;
    my ($annotation) = $annotation_collection->get_Annotations($type);

    my $value;
    if (defined $annotation){
	$value = $annotation->value();
    }#fix for weird Bioperl error on Linux vs Mac OS
    elsif( grep {/^$type$/} $f->get_all_tags ){
	($value) = $f->get_tag_values($type);
    }#end fix
    else{
	$value = '';
    }

    return $value;
}
#-------------------------------------------------------------------------------
sub _get_p_ids {
    my $f = shift;

    my @parents = $f->annotation->get_Annotations('Parent');
    my @p_ids;

    foreach my $p (@parents){
	push(@p_ids, $p->{value});
    }

    #fix for weird Bioperl error on Linux vs Mac OS
    if (! @p_ids && grep {/^Parent$/} $f->get_all_tags ){
        @p_ids  = $f->get_tag_values('Parent');
    }
    #end fix

    return \@p_ids;
}
#-------------------------------------------------------------------------------
sub _validate_gene {
    my $gene = shift;
                
    my @strands;
        
    #Check strand and fail if we don't get a valid strand value.
    my $g_strand = $gene->{f}{_location}{_strand};
    if ($g_strand != 1 && $g_strand != -1) {
	warn "Invalid strand in gene caught at " .
	    "FlyBase::validate_gene\n";
	return undef;
    }
        
    #Push the strand onto an array so that we can check all strands for
    #internal consistancy at the end.
    push @strands, $g_strand;
        
    for my $transcript (@{$gene->{transcripts}}) {
	#Check strand and fail if we don't get a valid strand value.
	my $t_strand = $transcript->{f}{_location}{_strand};
	if ($t_strand != 1 && $t_strand != -1) {
	    warn "Invalid strand in transcript caught " .
		"at FlyBase::validate_gene\n";
	    return undef;
	}
	#Push the strand onto an array so that we can check 
	#all strands for internal consistancy at the end.
	push @strands, $t_strand;
                
	for my $cds (@{$transcript->{cdss}}) {
	    #Check strand and fail if we don't get a valid strand value.
	    my $c_strand = $transcript->{f}{_location}{_strand};
	    if ($c_strand != 1 && $c_strand != -1) {
		warn "Invalid strand in cds caught " .
		    "at FlyBase::validate_gene\n";
		return undef;
	    }
	    #Push the strand onto an array so that we can check
	    #all strands for internal consistancy at the end.
	    push @strands, $c_strand;
	}
	
	for my $exon (@{$transcript->{exons}}) {
	    #Check strand and fail if we don't get a valid strand value.
	    my $e_strand = $transcript->{f}{_location}{_strand};
	    if ($e_strand != 1 && $e_strand != -1) {
		warn "Invalid strand in exon caught " .
		    "at FlyBase::validate_gene\n";
		return undef;
	    }
	    #Push the strand onto an array so that we can check
	    #all strands for internal consistancy at the end.
	    push @strands, $e_strand;
	}
    }
    #Check that all strands are the same, and fail if they are not.
    my $first = shift @strands;
    return undef if grep {$first != $_} @strands;
    
    #No failures, so return success.
    return 1;
}
#--------------------------------------------------------------------------------
sub _get_gene_seq {
    my $g   = shift;
    my $seq = shift;
    my $seq_len = shift || length_o($seq);

    my $g_b = $g->{f}->start();
    my $g_e = $g->{f}->end();

    my ($src_s, $src_e);
    $src_s = (($g_b - 1) < 0) ? 1 : $g_b;
    $src_e = ($g_e > $seq_len) ? $seq_len : $g_e; 

    my $g_seq = substr_o($seq, $src_s-1, abs($src_e-$src_s)+1);

    return ($g_seq, $src_s, $src_e);
}
#-------------------------------------------------------------------------------
sub _load_seqs {
    my $genes = shift;
    my $seq   = shift;
    my $seq_len = shift || length_o($seq);

    foreach my $g (@{$genes}){
	my ($g_seg_seq, $src_start, $src_end) = _get_gene_seq($g, $seq, $seq_len);
	
	$g->{seq}     = $g_seg_seq;
	$g->{src_s}   = $src_start;
	$g->{src_e}   = $src_end;
	$g->{i_start} = 1;
	$g->{i_end}   = length($g_seg_seq);

        foreach my $t (@{$g->{mRNAs}}){
	    foreach my $e (@{$t->{exons}}){
		my $e_seq = _get_exon_seq($e, $seq); 
		
		$e->{seq} = $e_seq;
                $e->{i_start} = $e->{f}->start() - $g->{src_s} + 1;

		$e->{i_end}   = $e->{f}->end() - $g->{src_s} + 1;  
	    }
	    foreach my $c (@{$t->{cdss}}){
		my $c_seq = _get_exon_seq($c, $seq);

		$c->{seq} = $c_seq;
		$c->{i_start} = $c->{f}->start() - $g->{src_s} + 1;

		$c->{i_end} = $c->{f}->end() - $g->{src_s} + 1;
	    }

	    my $t_seq   = _get_mRNA_seq($t, $seq); 
	    my $cds_seq = _get_cds_seq($t, $seq);

	    #debug
	    my $index = index($t_seq, $cds_seq) if($cds_seq);

	    $t->{seq} = $t_seq;
	    $t->{cds_seq} = $cds_seq;

	    $t->{i_start} = $t->{f}->start() - $g->{src_s} + 1;
	    $t->{i_end} = $t->{f}->end()   - $g->{src_s} + 1; 

	}
    }
}
#-------------------------------------------------------------------------------
sub _get_cds_seq {
    my $t   = shift;
    my $seq = shift;

    my $sorted = _sort_cdss($t);

    my $cds_seq;
    foreach my $c (@{$sorted}){
	my $c_seq = $c->{seq} || _get_exon_seq($c, $seq);
	$cds_seq .= $c_seq;

    }
    return $cds_seq;
}
#-------------------------------------------------------------------------------
sub _get_mRNA_seq {
    my $t   = shift;
    my $seq = shift;
    
    my $sorted = _sort_exons($t);

    my $transcript;
    foreach my $e (@{$sorted}){
	my $exon_seq = $e->{seq} || _get_exon_seq($e, $seq);
	$transcript .= $exon_seq;
	
    }
    return $transcript;
}
#-------------------------------------------------------------------------------
sub _get_exon_seq {
    my $e   = shift;
    my $seq = shift;

    my $e_b = $e->{f}->start();
    my $e_e = $e->{f}->end();

    my $exon_seq = substr_o($seq, $e_b-1, abs($e_e-$e_b)+1);

    $exon_seq = Fasta::revComp($exon_seq) if $e->{f}->strand() == -1;

    return $exon_seq;
}
#-------------------------------------------------------------------------------
sub _sort_exons {
    my $t = shift;

    my @sorted;

    if ($t->{f}->strand() ==  1){
	@sorted = 
	    sort {$a->{f}->start <=> $b->{f}->start} @{$t->{exons}}; 
    }
    elsif ($t->{f}->strand() == -1){
	@sorted = 
	    sort {$b->{f}->start <=> $a->{f}->start} @{$t->{exons}}; 
    }
    else {
	die "ERROR: Unknown strand in GFFDB::_sort_exons!\n";
    }
    return \@sorted;
}
#-------------------------------------------------------------------------------
sub _sort_cdss {
    my $t = shift;

    my @sorted;
    if    ($t->{f}->strand() ==  1){
                @sorted =
		    sort {$a->{f}->start <=> $b->{f}->start} @{$t->{cdss}};
	    }
    elsif ($t->{f}->strand() == -1){
                @sorted =
		    sort {$b->{f}->start <=> $a->{f}->start} @{$t->{cdss}};
	    }
    else {
	die "Error: unknown strand in GFFDB::_sort_cdss!\n";
    }
    return \@sorted;
}
#-------------------------------------------------------------------------------
sub _get_maker_qi {
    my $mRNA = shift;

    my $name = $mRNA->{name};

    return 'NA' unless defined $name;

    my ($maker_qi) = $name =~ /QI:([\d\|\.-]+)/;
    return $maker_qi;
}
#-------------------------------------------------------------------------------

1;


