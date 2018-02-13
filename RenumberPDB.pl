use strict;
use Carp qw(croak);
use File::Copy;
use Cwd;

{
	my ($PDB) = @ARGV;
	my $ext = "_renumbered";
	my $dir = getcwd;
	#print $dir;
	my $oldresnu;

	open(FILE, '<', "$PDB.pdb") or die "Cannot open file: $!";
	while (my $line=<FILE>){
	    if ($line =~ /^ATOM/){
	        substr($line, 22, 5) =~ /(\S+)/;
        
	        $oldresnu =$1  ;
        
	        #print $oldresnu ,"\n\n\n";;
	        last;
	    }
    
	}
	close FILE;

	my $new_resnu;
	my $t1;
	my $t2;
	my $count = 1;

	open (OUTFILE, ">>$PDB$ext.pdb");

	open(FILE1, '<', "$PDB.pdb") or die "Cannot open file: $!";
	while (my $line1=<FILE1>){
    
	    if ($line1 =~ /^ATOM/){
	        #substr($line1, 21,1) = "A";
	        substr($line1, 22, 5) =~ /(\S+)/;
	        $new_resnu =$1  ;
	        if ($new_resnu eq $oldresnu){
	            $t1 = sprintf "%4d", $count;
	            substr($line1, 22, 4) = $t1;
	            print OUTFILE $line1;
            
	        }
	        else{
	            $count++;
	            $t2 =sprintf "%4d", $count;
	            substr($line1, 22, 4) = $t2;
	            print OUTFILE $line1;
	            $oldresnu = $new_resnu;
	        }
	    }
	}
	close OUTFILE;
	close FILE1;
	my $cmd1 = "rm $PDB.pdb";
	my $cmd2 = "cp $PDB$ext.pdb $PDB.pdb";
	my $cmd3 = "rm $PDB$ext.pdb";
	system ($cmd1);
	system ($cmd2);
	system ($cmd3);
}

    