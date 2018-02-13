#!/usr/bin/perl

################################################################################
#                                                                              #
#                                  pdb2seq1h.pl                                #
#                                                                              #
# Prints horizontally on screen the amino acid sequence of a protein chain     #
# given its ATOM records in PDB format.                                        #
#                                                                              # 
# ARGUMENT:                                                                    #
# ---> PDBFILE: path to the PDB file containing ONLY the ATOM records of the   #
#               chain. The last record does have a line terminator.            #
#                                                                              #
# OUTPUT:                                                                      # 
#  The amino acid sequence of the chain, printed with one-letter codes. No     #
#  separators are put between the codes. The sequence starts at the start of   #
#  the first line of output and is followed by a line terminator.              #
#                                                                              #
################################################################################

%AMINO1 = ("ALA" => "A",
           "VAL" => "V",
           "PHE" => "F",
           "PRO" => "P",
           "MET" => "M",
           "ILE" => "I",
           "LEU" => "L",
           "ASP" => "D",
           "GLU" => "E",
           "LYS" => "K",
           "ARG" => "R",
           "SER" => "S",
           "THR" => "T",
           "TYR" => "Y",
           "HIS" => "H",
           "CYS" => "C",
           "ASN" => "N",
           "GLN" => "Q",
           "TRP" => "W",
           "GLY" => "G");

open(PDB, "<$ARGV[0]");

if($line=<PDB>) {
    substr($line, 12, 4) =~ /(\S+)/;
    $anam= $1;
    if ($line =~/^ATOM/ && $anam eq "CA"){
  $oldnam=substr($line, 17, 3);
  $oldidx=substr($line, 22, 4);
  $oldcod=substr($line, 26, 1);
  print $AMINO1{$oldnam};
    }
}

while($line=<PDB>) {
    substr($line, 12, 4) =~ /(\S+)/;
    $anam= $1;
    if ($line =~/^ATOM/ && $anam eq "CA"){
  $newnam=substr($line, 17, 3);
  $newidx=substr($line, 22, 4);
  $newcod=substr($line, 26, 1);

  if(($newidx ne $oldidx) || ($newcod ne $oldcod)) {
    print $AMINO1{$newnam};
  }
  
  $oldnam=$newnam;
  $oldidx=$newidx;
  $oldcod=$newcod;
    }
}

print "\n";

close(PDB);
