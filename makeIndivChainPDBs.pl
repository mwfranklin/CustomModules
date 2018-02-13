use strict;
use Carp qw(croak);
use File::Copy;
#use List::MoreUtils qw(firstidx);
{
    my ($pdb_id, $one_chain)= @ARGV;
    
    #Step 2: Obtaining chain ids for all the chains in protein complex
    my ($chain_ids_ref, $structure)=obtain_chain_ids_present_in_protein($pdb_id);
    my @chain_ids = @$chain_ids_ref;
	#print @chain_ids;
    my $ext = "_chain_";
    if ($structure eq "NMR"){
        eliminate_extra_models($pdb_id);
    }
    
	if ($one_chain eq "True"){
		make_chain_specific_protein_file($pdb_id, @chain_ids[0])
	}
    #Step 3: Making chain specific protein files out of protein complex
    else{
		foreach(@chain_ids)
		{
	        #print $_, "\n";
	        make_chain_specific_protein_file($pdb_id, $_);
	    }
    }
	
}

sub obtain_chain_ids_present_in_protein {
    my ($pdb) = @_;
    my $file = "$pdb.pdb";
    my @chain_ids;
    my @chain_ids_uniq;
    my $structure_type;
    open IN, "<", $file or croak "can't open $file\n";
    while (my $line= <IN>){
        if ($line =~ /^ATOM/){
            substr($line , 21, 1) =~ /(\S+)/;
            my $chain_id = $1;
            push (@chain_ids, $chain_id);
        }
        if ($line =~ /^EXPDTA/){
            if ($line =~ m/SOLUTION NMR/){
                $structure_type = "NMR";
            }
            else {
                $structure_type = "Other";
            }
        }
    }
    my @chain_ids_uniq = do { my %seen; grep{ !$seen{$_}++} @chain_ids}; # to get chain ids for all the chains in a protein
    close IN;
    return (\@chain_ids_uniq, $structure_type);
}
sub eliminate_extra_models {
    my ($pdb) =@_;
    my $file = "$pdb.pdb";
    my $ext = "model";
    my $count = 0;
    open IN, "<", $file or croak "can't open $file\n";
    open OUTFILE, ">$pdb$ext.pdb";
    while (my $line = <IN>){
        if ($line =~/^ENDMDL/ && $count ne 0){
            last;
        }
        if (($line =~/^ATOM/) || ($line =~/^HETATM/)) {#Includes all ligands except HOH in protein with same chain
            if (substr ($line, 17,3) ne "HOH" && substr($line, 16,1) ne "B"){
                print OUTFILE $line;
                $count ++;
            }
        }
    }
    close IN;
    close OUTFILE;
    my $cmd1 = "rm $pdb.pdb";
    my $cmd2 = "cp $pdb$ext.pdb $pdb.pdb";
    my $cmd3 = "rm $pdb$ext.pdb";
    system ($cmd1);
    system ($cmd2);
    system ($cmd3);
}

sub make_chain_specific_protein_file {
    my ($pdb, $chain_id) =@_;
    my $file = "$pdb.pdb";
    my $ext = "_chain_";
    my $count = 0;
    open IN, "<", $file or croak "can't open $file\n";
    open OUTFILE, ">$pdb$ext$chain_id.pdb";
    while (my $line = <IN>){
        if ($line =~/^END/ && $count ne 0){
            last;
        }
        if (($line =~/^ATOM/) || ($line =~/^HETATM/)) {#Includes all ligands except HOH in protein with same chain
            if (substr ($line, 21,1) eq $chain_id && substr ($line, 17,3) ne "HOH" && substr($line, 16,1) ne "B"){
                print OUTFILE $line;
                $count ++;
            }
        }
    }
    close IN;
    close OUTFILE;
}