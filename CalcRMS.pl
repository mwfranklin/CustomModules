#!/usr/bin/perl 
use strict;
use Carp qw(croak);

my ($PDB1,$PDB2,$number_of_residues_1,$number_of_residues_2)=@ARGV;
my $file = "$PDB1.pdb";
my $file2 = "$PDB2.pdb";

my ($line,$res_start, $line2, $res_start2, @array_residues, @array_residues2, $in, $in2, $resname2, $resname, $atomname, $atomname2, @sum);
my ($i, $i2);
#my $count=0;
#my $count2=0;
open FILE, '<', $file or croak "can't open $file\n";
open FILE2, '<', $file2 or croak "can't open $file2\n";
my $n;
my $rms;
my @rmsd;
if ($number_of_residues_1 == $number_of_residues_2){
    $n = $number_of_residues_1-1;#because looping is from 0 to n
}
else{
    if ($number_of_residues_1<$number_of_residues_2){
        $n = $number_of_residues_1-1;
    }
    else{
        $n = $number_of_residues_2-1;
    }
}

for $i (0..$n){
    my ($x1, $y1, $z1, @x1, @y1, @z1);
    my ($x2,$y2,$z2,@x2,@y2,@z2);
    my $count=0;
	my $count2=0;
    my $atomname_start;
    my $atomname_start2;
    $line = <FILE>;
    
    #print $line , "\n";
    $res_start = substr($line,17,3);
    #print $res_start,"TAB\n\n";
    substr($line, 12, 4) =~ /(\S+)/;
    $atomname_start =$1;
    #print $line;
    
    while (substr($line,17,3)eq $res_start){
        #print $line;
        substr($line, 30, 8) =~ /(\S+)/;
        $x1 = $1;
        
        substr($line, 38, 8) =~ /(\S+)/;
        $y1 = $1;
        
        substr($line, 46, 8) =~ /(\S+)/;
        $z1 = $1;

        #if ($atomname_start eq "CA"){

        push (@x1, $x1);
        push (@y1, $y1);
        push (@z1, $z1);
        #print $line, "\n";
        #print $x1, "\t", $y1,"\t", $z1, "\n";
        #}
        
        $line = <FILE>;
        $count++;
        substr($line, 12, 4) =~ /(\S+)/;
        $atomname_start =$1;
        if ($atomname_start eq "N"){
            last;
        }
    }
    #print "\n\n";
    seek(FILE, -length($line), 1);#unread a line
    #print @x1, "\t", @y1,"\t", @z1, "\n";
    #print $count, "\n\n";
    
    push (@array_residues, $count);
    
    $line2 = <FILE2>;
    #print $line2, "\n";
    $res_start2 = substr($line2,17,3);
    #print $res_start2,"\tWAY\n\n";
    substr($line2, 12, 4) =~ /(\S+)/;
    $atomname_start2 =$1;
    #print $line;
    while (substr($line2,17,3)eq $res_start2){
        #print $line2;
        substr($line2, 30, 8) =~ /(\S+)/;
        $x2 = $1;
        
        substr($line2, 38, 8) =~ /(\S+)/;
        $y2 = $1;
        
        substr($line2, 46, 8) =~ /(\S+)/;
        $z2 = $1;

        #if ($atomname_start2 eq "CA"){

        push (@x2, $x2);
        push (@y2, $y2);
        push (@z2, $z2);
        #print $line2, "\n";
        #print $x2, "\t", $y2,"\t", $z2, "\n";
        #}

        $line2 = <FILE2>;
        $count2++;
        substr($line2, 12, 4) =~ /(\S+)/;
        $atomname_start2 =$1;
        if ($atomname_start2 eq "N"){
            last;
        }
    }
   # print "$count2\n\n\n";
    seek(FILE2, -length($line2), 1);
    #my $k = 0;
    my $c = 0;
if ($count2 == $count){
    $c = $count;
}else{
    if ($count2 > $count){
    	$c = $count;
    }else{
	$c = $count2;
    }
}
my $sum2;
    for my $k (0..($c-1)){
    #print "($x1[$k]-$x2[$k])**2 + ($y1[$k]-$y2[$k])**2 + ($z1[$k]-$z2[$k])**2\n";
        $sum2 += ($x1[$k]-$x2[$k])**2 + ($y1[$k]-$y2[$k])**2 + ($z1[$k]-$z2[$k])**2;
    #print $sum2, "\n\n";
    }
    $rms = sqrt($sum2/$c);
   print $rms, "\t", $res_start,"\t", $res_start2, "\n";
    push (@rmsd, $rms);
}
close (FILE);
close (FILE2);
#print (\@rmsd);