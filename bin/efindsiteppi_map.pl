#!/usr/bin/perl -w

#===============================================================================
#        ______ _           _  _____ _ _       _____  _____ _____ 
#       |  ____(_)         | |/ ____(_) |     |  __ \|  __ \_   _|
#    ___| |__   _ _ __   __| | (___  _| |_ ___| |__) | |__) || |  
#   / _ \  __| | | '_ \ / _` |\___ \| | __/ _ \  ___/|  ___/ | |  
#  |  __/ |    | | | | | (_| |____) | | ||  __/ |    | |    _| |_ 
#   \___|_|    |_|_| |_|\__,_|_____/|_|\__\___|_|    |_|   |_____|
#
#
#   eFindSitePPI - prediction of protein binding sites from meta-threading
#
#   Computational Systems Biology Group
#   Department of Biological Sciences
#   Center for Computation & Technology
#   Louisiana State University
#   407 Choppin Hall, Baton Rouge, LA 70803, USA
#
#   http://www.brylinski.org
#
#   Report bugs and issues to smahes2@lsu.edu michal@brylinski.org
#
#   Copyright 2013 Michal Brylinski Surabhi Maheshwari
#
#   This file is part of eFindSitePPI package.
#
#   eFindSitePPI is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   eFindSitePPI is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with eFindSitePPI. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

 use strict;
 use File::Path;
 use File::Copy;
 use Benchmark;
 use Cwd;
 use File::Slurp;
 use File::Temp qw/ tempfile tempdir /;
 
 local $| = 1;
 
 print "------------------------------------------------------------\n";
 print "                      efindsiteppi_map\n";
 print "                        version 1.0\n";
 print "                 Threading library mapping\n\n";
 print "       report bugs and issues to michal\@brylinski.org\n";
 print "------------------------------------------------------------\n\n";
 
 if ($#ARGV < 5)
 {
  print "efindsite_map -t <threading library in FASTA format>\n";
  print "              -p <eFindSitePPI library in FASTA format>\n";
  print "              -o <output filename>\n";
  print "              -a <number of processors to use, default 1>\n";
  die "\n";
 }
 
 die "Could not find formatdb\n" if ( !( `which formatdb` ) );
 die "Could not find blastall\n" if ( !( `which blastall` ) );
 
 my $ftpl1 = '';
 my $fpdb1 = '';
 my $fout1 = '';
 my $fcpu1 = 1;
 
 for ( my $i = 0; $i <= $#ARGV; $i++ )
 {
  $ftpl1 = $ARGV[$i+1] if ( $ARGV[$i] eq '-t' );
  $fpdb1 = $ARGV[$i+1] if ( $ARGV[$i] eq '-p' );
  $fout1 = $ARGV[$i+1] if ( $ARGV[$i] eq '-o' );
  $fcpu1 = $ARGV[$i+1] if ( $ARGV[$i] eq '-a' );
 }
 
 die "Provide template library in FASTA format\n" if ( !( -e $ftpl1 ) or !length($ftpl1) );
 die "Provide PDB library in FASTA format\n"      if ( !( -e $fpdb1 ) or !length($fpdb1) );
 die "Provide output filename\n"                  if ( !length($fout1) );
 
 my $bt0 = Benchmark->new;
 
 my @pdb1 = read_file($fpdb1); chomp(@pdb1);
 
 my @tpl1 = read_file($ftpl1); chomp(@tpl1);
 
 my $dir1 = getcwd();
 
 my $dir2 = tempdir( CLEANUP => 1 );
 
 my ($fh1, $tmpfil1) = tempfile( DIR => $dir2, UNLINK => 1);
 my ($fh2, $tmpfil2) = tempfile( DIR => $dir2, UNLINK => 1);
 
 copy($ftpl1, "$tmpfil1");
 copy($fpdb1, "$tmpfil2");
 
 printf("Tempdir created: %s\n\n", $dir2);
 
 chdir($dir2) or die "Cannot chdir to $dir2 $!";
 
 open (FOR, "formatdb -i $tmpfil1 -p T -o T 2>&1 |") || die "Cannot run formatdb -i $tmpfil1 -p T -o T\n";
  my @for1=<FOR>;
  chomp(@for1);
 close (FOR);
 
 print "Mapping ... ";
 
 open (BLA, "blastall -p blastp -d $tmpfil1 -i $tmpfil2 -m 8 -v 1 -b 1 -a $fcpu1 2>&1 |") || die "Cannot run blastall -p blastp -d $tmpfil1 -i $tmpfil2 -m 8 -v 1 -b 1 -a $fcpu1\n";
  my @bla1=<BLA>;
  chomp(@bla1);
 close (BLA);
 
 chdir($dir1) or die "Cannot chdir to $dir1 $!";
 
 my %bla3 = ();
 my %bla4 = ();
 
 foreach my $wbla1 (@bla1)
 {
  while ( $wbla1 =~ /\t/ ) { $wbla1 =~ s/\t/\ /g; }
  while ( $wbla1 =~ /\ \ / ) { $wbla1 =~ s/\ \ /\ /g; }
  
  substr($wbla1,  0, 1) = '' if ( substr($wbla1,  0, 1) eq ' ' );
  substr($wbla1, -1, 1) = '' if ( substr($wbla1, -1, 1) eq ' ' );
  
  my @bla2 = split(/\ /, $wbla1);
  
  my $nbla2 = @bla2;
  
  if ( $nbla2 == 12 )
  {
   $bla2[10] = '1'.$bla2[10] if ( substr($bla2[10],0,1) eq 'e' ); $bla2[10] *= 1.0;
   
   if ( exists $bla3{$bla2[0]} )
   {
    if ( $bla2[10] < $bla4{$bla2[0]} )
    {
     $bla3{$bla2[0]} = $bla2[1];
     $bla4{$bla2[0]} = $bla2[10];
    }
   }
   else
   {
    $bla3{$bla2[0]} = $bla2[1];
    $bla4{$bla2[0]} = $bla2[10];
   }
  }
 }
 
 my %tpl2 = ();
 
 foreach my $wbla3 ( keys %bla3 )
 {
  if ( exists $tpl2{$bla3{$wbla3}} )
  {
   $tpl2{$bla3{$wbla3}} .= '&'.$wbla3;
  }
  else
  {
   $tpl2{$bla3{$wbla3}} = $wbla3;
  }
 }
 
 my @out1 = ();
 
 my $n1 = 0;
 my $n2 = 0;
 
 foreach my $wtpl2 ( sort { $a cmp $b } keys %tpl2 )
 {
  my @out2 = split(/\&/, $tpl2{$wtpl2});
  
  my $nout2 = @out2;
  
  if ( $nout2 )
  {
   my $out3 = sprintf("%15s%5d", $wtpl2, $nout2);
   
   foreach my $wout2 (@out2)
   {
    $out3 .= " $wout2";
   }
   
   push(@out1, "$out3\n");
   
   $n1 += $nout2;
   
   $n2++;
  }
 }
 
 print "done: $n1 -> $n2\n\n";
 
 if ( @out1 )
 {
  write_file($fout1, @out1);
  
  print "Mapping file written to $fout1\n\n";
 }
 
 my $bt1 = Benchmark->new;
 
 printf("------------------------------------------------------------\n");
 printf("Walltime: %s\n", timestr(timediff($bt1, $bt0)));
 printf("------------------------------------------------------------\n");
 
 exit(0);
 
