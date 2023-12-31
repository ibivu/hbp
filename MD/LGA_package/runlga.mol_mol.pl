#!/usr/bin/perl
#
# runlga.mol_mol.pl
# script to evaluate molecule pairs (model.target) 
# using LGA program and collect_PDB.pl script
#
# Author: Adam Zemla
# Email: adamz@llnl.gov
#
# version: 08/19/2002 
# updated by Lucy Forrest 10/16/2004 (email: lrf2103@columbia.edu)
#
# usage: runlga.mol_mol.pl model target selected_parameters
#

$|=1;

# subdirectory with executables
#$bin='bin';
#$bin='.';
$bin='.';

# subdirectory for results
$dirres = 'RESULTS/';
$pdbs = 'MOL2/';
$tmp = 'TMP/';

# --------------------------------------------------------------------- #

# usage 
if (@ARGV < 2) { die "\nUsage:\n\t$0 mol1 mol2 params\n\nSee http://as2ts.llnl.gov/ for help\n\n"; }

$par="@ARGV";
$res='.res';
$pdb='.pdb';
$lga='.lga';
$d='.';

print "\n";
if (!-d $dirres) { mkdir($dirres); print "Making $dirres for output (*$res) files\n"; }
# print "LGA results (*$res data) are stored in ./$dirres\n";

# read input parameters and filenames
$entrmol1=$ARGV[0];
$entrmol2=$ARGV[1];
$ARGV[0]="";
$ARGV[1]="";
# if (!-e $entrmol1) { die "\nError: $entrmol1 doesn't exist\n\n"; }
# if (!-e $entrmol2) { die "\nError: $entrmol2 doesn't exist\n\n"; }

@LINE=split(/\//,$entrmol1);
@R=reverse(@LINE);
$mol1=$R[0];
@LINE=split(/\//,$entrmol2);
@R=reverse(@LINE);
$mol2=$R[0];

$par="@ARGV";

# subdirectory for input. LGA takes input from MOL2 directory
if (!-d $pdbs) { mkdir($pdbs); print "Making $pdbs for input files\n"; }
# print "\nPutting PDB input for LGA processing in ./$pdbs\n";
sleep 1;

# create input file for LGA
$model="$mol1$d$mol2";
$mol1_new="";
$mol2_new="";
while ($_="@ARGV") {
  shift;
  ($xx)=/(\S+)/;
  if($xx =~ /-mol1:/) {
    ($tmp1,$mol1_new) = split /\:/,$xx;
    $model="$mol1_new";
  }
  if($xx =~ /-mol2:/) {
    ($tmp2,$mol2_new) = split /\:/,$xx;
    $model="$mol2_new";
  }
}
if($mol1_new ne "" && $mol2_new ne "") {
  $model="$mol1_new$d$mol2_new";
}

# print "Processing structures in: $pdbs$model\n";
$model=~s/\*//g;
system "rm -rf $tmp$model$pdb";
system "echo 'MOLECULE $mol1' > $pdbs$model ";
system "$bin/collect_PDB.pl $entrmol1 | grep -v '^MOLECULE ' | grep -v '^END' >> $pdbs$model ";
system "echo 'END' >> $pdbs$model ";
system "echo 'MOLECULE $mol2' >> $pdbs$model ";
system "$bin/collect_PDB.pl $entrmol2 | grep -v '^MOLECULE ' | grep -v '^END' >> $pdbs$model ";
system "echo 'END' >> $pdbs$model ";

# temp directory with outputs from LGA. LGA puts calculation results to TMP directory
if (!-d $tmp) { mkdir($tmp); print "Making $tmp for output files\n"; }
sleep 1;

# run LGA - it knows to look in MOL2 for the input...
print "Running: lga $model $par\n";
system "$bin/lga $model $par > $dirres$model$res";
if (-z "$dirres$model$res") { die "\nError: Problem running LGA program. Check settings (see INFO).\n"; }
sleep 1;

# tidy up
if (!-e "$tmp$model$pdb") { 
  print "Warning: LGA didn\'t create file: $tmp$model$pdb \n"; 
}
else {
  system "cat $tmp$model$pdb >> $dirres$model$res";
}
sleep 1;

# removing input and output files after the processing is done ...
system "rm -rf $pdbs$model $tmp$model$lga $tmp$model$pdb";

print "Done! \n";

exit;
