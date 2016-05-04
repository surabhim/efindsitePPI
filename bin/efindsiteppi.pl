#!/project/michal/apps/perl/bin/perl
##!/usr/bin/perl -w

#===============================================================================
#        ______ _           _  _____ _ _       _____  _____ _____ 
#       |  ____(_)         | |/ ____(_) |     |  __ \|  __ \_   _|
#    ___| |__   _ _ __   __| | (___  _| |_ ___| |__) | |__) || |  
#   / _ \  __| | | '_ \ / _` |\___ \| | __/ _ \  ___/|  ___/ | |  
#  |  __/ |    | | | | | (_| |____) | | ||  __/ |    | |    _| |_ 
#   \___|_|    |_|_| |_|\__,_|_____/|_|\__\___|_|    |_|   |_____|
#                                                  
#   eFindSitePPI - prediction of protein binding sites from meta-threading
#
#   This software is distributed WITHOUT ANY WARRANTY (but with best wishes)
#
#   Report bugs and issues to smahes2@tigers.lsu.edu  michal@brylinski.org
#
#   Computational Systems Biology Group
#   Department of Biological Sciences
#   Center for Computation & Technology
#   Louisiana State University
#   407 Choppin Hall, Baton Rouge, LA 70803, USA
#
#   http://www.brylinski.org
#
#===============================================================================

use strict;
use warnings ;
use Benchmark;
use File::Slurp;
use Cwd;
use File::Copy;
use File::Temp qw/ tempfile tempdir /;
use Compress::Zlib;
use Algorithm::NeedlemanWunsch;
use List::Util qw( min max );
use AI::NaiveBayes1;
use YAML;

 local $| = 1;

 print "------------------------------------------------------------\n";
 print "                        efindsiteppi\n";
 print "                        version 1.0\n";
 print "              protein binding site prediction\n\n";
 print "       report bugs and issues to smahes2\@tigers.lsu.edu\n";
 print "                                 michal\@brylinski.org\n";
 print "------------------------------------------------------------\n\n";
 
 die "PPI_FRTMALIGN is not set\n" if !( $ENV{'PPI_FRTMALIGN'} );
 die "PPI_NACCESS is not set\n" if !( $ENV{'PPI_NACCESS'} );
 die "PPI_LIBSVM is not set\n" if !( $ENV{'PPI_LIBSVM'} );
 die "PPI_MODELS is not set\n" if !( $ENV{'PPI_MODELS'} );
 die "PPI_LIBRARY is not set\n" if !( $ENV{'PPI_LIBRARY'} );
 die "PPI_ETHREADMAP is not set\n" if !( $ENV{'PPI_ETHREADMAP'} );
 
 if ($#ARGV < 7)
 {
  print "ethread -s <target structure in PDB format>\n";
  print "        -t <output from eThread>\n";
  print "        -e <sequence profile>\n";
  print "        -o <output filename>\n\n";
  print "additional options:\n\n";
  print "        -b <sequence identity threshold for benchmarks, default 1.0>\n";
  print "        -m <TMscore threshold, default 0.4>\n";
  print "        -x <max number of templates, default 1000>\n";
  die "\n";
 }

 my $env_frtmalign = $ENV{'PPI_FRTMALIGN'};
 my $env_naccess = $ENV{'PPI_NACCESS'};
 my $env_libsvm = $ENV{'PPI_LIBSVM'};
 my $env_tempdataset = $ENV{'PPI_LIBRARY'};
 my $env_svm = ($ENV{'PPI_MODELS'}.'/residueSVM.model');
 my $env_nb = ($ENV{'PPI_MODELS'}.'/residueNBC.model');
 my $mfile = $ENV{'PPI_ETHREADMAP'};

 die "Could not find frtmalign\n" if ( !( -e $env_frtmalign ) );
 die "Could not find naccess\n" if ( !( -e $env_naccess ) );
 die "Could not find svm-predict\n" if ( !( -e $env_libsvm ) );
 die "Could not find template library\n" if ( !( -d $env_tempdataset ) or !( -e $env_tempdataset.'/pdb.ver' ) );
 die "Could not find SVM model\n" if ( !( -e $env_svm ) );
 die "Could not find NBC model\n" if ( !( -e $env_nb ) );
 die "Could not find mapping file\n" if ( !( -e $mfile ) );
 
 my $thresh_tmscore = 0.4 ;
 my $thresh_seqid = 1.00;
 my $thresh_temp = 1000;
 
 my $i; my $efile=''; my $sfile=''; my $target_id='';my $ffile='';my $pfile='';
 for ($i = 0; $i <= $#ARGV; $i++)
 {
     if ($ARGV[$i] eq '-t') {$efile = $ARGV[$i+1] ;} 
  elsif ($ARGV[$i] eq '-s') {$sfile = $ARGV[$i+1] ;}
  elsif ($ARGV[$i] eq '-o') {$target_id = $ARGV[$i+1] ;}
  elsif ($ARGV[$i] eq '-e') {$pfile = $ARGV[$i+1] ;}
  elsif ($ARGV[$i] eq '-b') {$thresh_seqid = $ARGV[$i+1] ;}
  elsif ($ARGV[$i] eq '-m') {$thresh_tmscore = $ARGV[$i+1] ;}
  elsif ($ARGV[$i] eq '-x') {$thresh_temp = $ARGV[$i+1] ;}
 }
 
 die "Provide target structure in PDB format\n" if ( !( -e $sfile ) );
 die "Provide output from eThread\n" if ( !( -e $efile ) );
 die "Provide sequence profile\n" if ( !( -e $pfile ) );
 die "Provide output filename\n" if ( !( length($target_id) ) );
 
 printf("!!! Benchmarking mode activated with max sid of %.2f !!!\n\n", $thresh_seqid) if ( $thresh_seqid < 1.0 );
 
 my @version1 = read_file($env_tempdataset.'/pdb.ver'); chomp(@version1);
 
 my $dir1 = getcwd();
 my $dir2 = tempdir(CLEANUP => 1);
 
 print "Tempdir created: $dir2\n\n";
 
 print "eFindSitePPI library path: $env_tempdataset\n";
 print "eFindSitePPI library version: $version1[0]\n";
 print "eFindSitePPI mapping: $mfile\n\n";
 
#NOTE:Target pdb should have only 1 chain , the "receptor"
##section1#################################################################################################
## Reading input files 
 
 my $t0 = Benchmark->new; 

 print "Reading input files ... ";

 my @target_pdb1 = read_file($sfile);
 my @target_pdb = grep(/^ATOM/,@target_pdb1);
 
 # fix chain id
 
 foreach my $wtarget_pdb (@target_pdb)
 {
  substr($wtarget_pdb, 21, 1) = 'A' if ( substr($wtarget_pdb, 21, 1) eq ' ' );
 }
 
 my @ethread = read_file($efile);
 my @prf = read_file($pfile);
 my @map_file_1 = read_file($mfile);

 my %h0 = (); my $key0 = ();
 foreach (@map_file_1)
  {
    my @split_line = split (/ +/, $_);
    $key0 = $split_line[1];
    my $values = substr($_, 21);
    my @array = split (/ /,$values);
    $h0{$key0} = [@array];
  }

 print "done\n\n";

##section2.1###############################################################################################
## Making a function for finding TM-score
## usage > tmscore(template.pdb,target.pdb);

 sub tmscore 
 {
   #my $tm_score = '/project/michal/apps/frtmalign/frtmalign';
   #my $tm_score = '/usr/local/frtmalign/frtmalign';
   my $tm_score =$env_frtmalign;
   open (FH, "$tm_score -m 1 $_[0] $_[1] 2>&1 |");
     my @tmscore_out = <FH>;
   close FH ; 
   my @grep_tms1 = grep(/TM-score=/,@tmscore_out);
   my @split_tms1 = split(/ +/,$grep_tms1[0]);
   my $tmscore1 = $split_tms1[5];   
   $tmscore1 =~ s/,// ; $tmscore1 =~ s/TM-score=//;
   my $aln_length1= $split_tms1[2];$aln_length1 =~ s/,// ;
   my $rmsd1 = $split_tms1[4]; $rmsd1 =~ s/,// ;
   my $seq_id1 = $split_tms1[6]; chomp $seq_id1;$seq_id1 =~ s/,// ;$seq_id1 =~ s/ID=//; 
   my $aln_temp1 = $tmscore_out[16];chomp $aln_temp1;
   my $aln_dots1 = $tmscore_out[17];chomp $aln_dots1;
   my $aln_targ1 = $tmscore_out[18];
   my @mat = read_file("trf.mat") if ( -e "trf.mat" );
   my @grep_tms2 = grep(/Chain 1/,@tmscore_out);
   my @split_tms2 = split(/ +/,$grep_tms2[0]);
   # my $length_temp1 = $split_tms2[3];chomp $length_temp1;
   my $length_temp1=substr($grep_tms2[0],25,4)*1;
   my @grep_tms3 = grep(/Chain 2/,@tmscore_out);
   my @split_tms3 = split(/ +/,$grep_tms3[0]);
   # my $length_targ1 = $split_tms3[3];chomp $length_targ1;
    my $length_targ1=substr($grep_tms3[0],25,4)*1;
   return ($length_targ1,$length_temp1,$tmscore1,$rmsd1,$aln_length1,$seq_id1,$aln_temp1,$aln_dots1,$aln_targ1,@mat);      
  }

##section2.2###############################################################################################
## Translate PDB files to sequence string file
## usage > pdb2seq(@input_pdb_array);

 sub pdb2seq
  {
  my %aa=qw(
	ALA 	A
	CYS 	C
	ASP 	D
	GLU 	E 
	PHE 	F
	GLY	G
	HIS	H
	ILE	I
	LYS	K
	LEU	L
	MET	M
	ASN	N
	PRO	P
	GLN	Q
	ARG	R
	SER	S
	THR	T
	VAL	V
	TRP	W
	TYR	Y);

  my @pdb_in = @_;my $residue ; my $resno;
  my $oldresno=-1;my $seq=(); #print "@pdb_in";
  foreach (@pdb_in)
   {
    if (/^ATOM/)
    {
     my $type = substr($_,13,2); 
     if ($type eq "CA")	
      {	
       my $res = substr($_, 17, 3); 
       chomp($residue=$aa{$res});$residue=~ s/^\s+//;$residue=~s/\s+$//;
       my $resno=substr($_, 22, 4);
       if ($resno>$oldresno)
        {
         $seq=$seq.$residue ; 
         $oldresno=$resno;
         }
       }
     }
    }
   return $seq;
   }

##section2.3###############################################################################################
# Gives the global sequence identity between two sequence strings.
# usage > get_identity(seq_string1, seq_string2);

 use vars qw( %nwmat4 @nwseq3 @nwseq4 $nwseq3o $nwseq4o );
 my %nwmat1 = (); 
 my @nwmat3 = qw(A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X);
 
 while ( my $wdat1 = <DATA> )
 {
  chomp $wdat1;
  if ( length($wdat1) == 70 and $wdat1 ne '   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X' )
  {
   my $nwr1 = substr($wdat1, 0, 1);   
   for ( my $xg = 0; $xg < 23; $xg++ )
   {
    $nwmat1{$nwr1.$nwmat3[$xg]} = substr($wdat1, 1 + $xg * 3, 3) * 1;
   }
  }
 }
 
 my @nwseq1 = ();my @nwseq2 = ();
 my $nwseq1o = '';my $nwseq2o = '';
 
 sub blosum 
 { 
  my ($anw, $bnw) = @_;  
  my $snw = 0;  
  $snw = $nwmat1{$anw.$bnw} if ( exists $nwmat1{$anw.$bnw} );  
  return ($snw);
 }
 
 my $matcher1 = Algorithm::NeedlemanWunsch->new(\&blosum);
 $matcher1->gap_open_penalty(-10);
 $matcher1->gap_extend_penalty(-2);
 
 sub prepend_align1 
 { 
  my ($i, $j) = @_;
  $nwseq1o = $nwseq1[$i] . $nwseq1o;
  $nwseq2o = $nwseq2[$j] . $nwseq2o;
 }
 
 sub prepend_first_only1 
 { 
  my $i = shift;  
  $nwseq1o = $nwseq1[$i] . $nwseq1o;
  $nwseq2o = "-$nwseq2o";
 }
 
 sub prepend_second_only1 
 {
  my $j = shift;
  $nwseq1o = "-$nwseq1o";
  $nwseq2o = $nwseq2[$j] . $nwseq2o;
 }
 
 sub get_identity 
 {
  my ($iseq1, $iseq2) = @_;  
  my $iss = 0.0; 
  @nwseq1 = split(//, $iseq1);
  @nwseq2 = split(//, $iseq2);
  
  $nwseq1o = '';$nwseq2o = '';
 
  my $score = $matcher1->align(\@nwseq1, \@nwseq2, { align => \&prepend_align1, shift_a => \&prepend_first_only1, shift_b => \&prepend_second_only1, });
  
  my @nwseq1a = split(//, $nwseq1o);
  my @nwseq2a = split(//, $nwseq2o);
  
  my $niseq1 = @nwseq1a;my $niseq2 = @nwseq2a;
  my $isid1 = 0; my $isid2 = 0;
  
  if ( $niseq1 == $niseq2 )
  {
   for ( my $ixa = 0; $ixa < $niseq1; $ixa++ )
   {
    $isid1++ if ( $nwseq1a[$ixa] ne '-' and $nwseq2a[$ixa] ne '-' and $nwseq1a[$ixa] eq $nwseq2a[$ixa] );
    $isid2++ if ( $nwseq1a[$ixa] ne '-' and $nwseq2a[$ixa] ne '-' );
   }
  }
  
  $iss = $isid1 / $isid2 if ( $isid2 );  
  return $iss;
 }

##section2.4############################################################################################### 

 sub map_alignments
   {
    my $aln_target = $_[0];
    my $aln_template = $_[1];
    my %hash1 = %{$_[2]};
    my %hash2 = %{$_[3]};
    my $target_length = $_[4];
    my $length =  length ($aln_target);
    my $m1 = 0 ; my $m2 =0;
    for (my $i=0; $i<= $length; $i++)
     {
      my $count_m1=0; my $count_m2=0;
      if ( substr($aln_target,$i,1) ne '-') { $m1++; $count_m1 = 1;}
      if ( substr($aln_template,$i,1) ne '-') { $m2++; $count_m2 = 1;}
      if ( exists $hash1{$m2} && $count_m1 == 1 && $count_m2 == 1 && $m1 <= $target_length)
       {
         if (exists $hash2{$m1})
          {
           $hash2{$m1}=$hash2{$m1}.":".$hash1{$m2};
           }
         else 
          {
           $hash2{$m1}=$hash1{$m2}; 
           } 
        }
       elsif ( !exists $hash1{$m2} && $count_m1 == 1 && $count_m2 == 1 && $m1 <= $target_length)
        {
         if (exists $hash2{$m1})
          {
           $hash2{$m1}=$hash2{$m1}.":"."0";
           }
         else 
          {
           $hash2{$m1}="0"; 
           }
         }
       }
    return \%hash2 ; 
   }

  sub find_fraction
   {
    my %hash3 = %{$_[0]};
    my %hash4 = () ;
    foreach my $key (keys %hash3)
      {
       my @split = split(/:/,$hash3{$key});
       my $denominator=0;
       foreach (@split)
         {
          if ($_ == 0){$denominator = $denominator+1 ;}
          else {$denominator = $denominator + $_ ;}        
          }
       my $sum=0; my $fraction=0;
       if ($denominator > 0)
        {
         foreach (@split)
         {
         $sum= $sum + $_;
         }
        $fraction = ($sum/$denominator) ;
        $hash4{$key}=$fraction; 
        }
      }
    return \%hash4 ;
    }

##section2.5############################################################################################### 

 sub scale 
  {
  my $in = $_[0] ;
  my $lb = $_[1]  ;
  my $ub = $_[2] ;
  my $out = ($in - $lb)/($ub-$lb)* 2.0 - 1.0;
  $out = -1.0 if ( $out < -1.0 );
  $out = 1.0 if ( $out > 1.0 );
  return $out;
  }

##section3.1###############################################################################################
## Finding the TM-score of all the templates identified in the ethread file (mapped on ethread mapping file)

print "Calculating structure alignments .";

 my $target_seq = pdb2seq(@target_pdb);
 my @ethread_temp1 = grep(/>P1;/,@ethread);
 my %ethread_temp2; my @count_template ; my $count_temp;
 foreach ( @ethread_temp1)
 {
  my $ethread_temp = substr($_,4,5); # finds the pdbid of the template
  my $tmscore_temp = substr($_,10,6); 
  $ethread_temp2{$ethread_temp}=$tmscore_temp ; 
  #push(@ethread_temp2,"$ethread_temp \n"); 
 }

 my ($fh0, $tmpfil0) = tempfile( DIR => $dir2,SUFFIX => ".pdb", UNLINK => 1);
 write_file( $tmpfil0, @target_pdb ) ;
 
 chdir $dir2; 

 my $target_length ; my $k1 = 0; 
 my (%h1, $value1, $key1,%h2,$value2,%h3,$value3,%h4,%h20);
 my (%hbond, %salt,%contacts);
 foreach my $template (reverse sort { $ethread_temp2{$a} <=> $ethread_temp2{$b} } keys(%ethread_temp2) )
 { 
  if ($k1 < $thresh_temp)
  { #print "$template $ethread_temp2{$template} $k1 \n";
   $value1=();$key1=();
   $template=~ s/^\s+|\s+$//g;
   foreach my $key_chk0 ( keys %h0)
   {
    if ($template eq $key_chk0)
      {
        foreach my $TEMPLATE (@{$h0{$key_chk0}})
          {
           chomp $TEMPLATE;
            my $file1 = $env_tempdataset.'/data/'.substr($TEMPLATE, 1, 2).'/'.$TEMPLATE; 
            ##my $file1 = $env_tempdataset."/$TEMPLATE";  
            if (-e $file1)
            {
             my @file2 = ();
             
             my $gz_tpl = gzopen($file1, "rb")  or die "Cannot open $file1 for reading: $gzerrno\n";
             while ( $gz_tpl->gzreadline($_) > 0) { push(@file2, $_); }
             die "Error reading from $file1: $gzerrno\n" if $gzerrno != Z_STREAM_END ;
             $gz_tpl->gzclose();
             
             my @grep_1 = grep(/^PRT ATOM/, @file2); foreach (@grep_1) {$_ =~ s/PRT //g ;}  
             
             my ($fh1, $tmpfil1) = tempfile( DIR => $dir2, UNLINK => 1);
             write_file( $tmpfil1, @grep_1 ) ;
             my $template_seq = pdb2seq(@grep_1); 
             my $seq_global = get_identity($target_seq,$template_seq);
             
             if ($seq_global <= $thresh_seqid)
             {
              print '.';
             
              my @grep_2 = grep (/^IAL PRT/,@file2);       
              my ($length_targ2,$length_temp2,$tmscore2,$rmsd2,$aln_length2,$seq_id2,$aln_temp2,$aln_dots2,$aln_targ2, @mat2) = &tmscore( $tmpfil1,$tmpfil0);
              $target_length = $length_targ2 ;
              if ($tmscore2 >= $thresh_tmscore)
                {
                 $k1++;
                 $count_temp = scalar (@grep_2);
                 push (@count_template, $count_temp);
                 my @split_1 = split (/\s+/,$mat2[2]);
                 my @split_2 = split (/\s+/,$mat2[3]);
                 my @split_3 = split (/\s+/,$mat2[4]);

                my ($m10,$m11,$m12,$m13) = ($split_1[2],$split_1[3],$split_1[4],$split_1[5]);
                my ($m20,$m21,$m22,$m23) = ($split_2[2],$split_2[3],$split_2[4],$split_2[5]);
                my ($m30,$m31,$m32,$m33) = ($split_3[2],$split_3[3],$split_3[4],$split_3[5]);

                $value1 = sprintf("TEMPLTE %s %4d %6.3f %6.2f %4d %6.3f %6.3f",$TEMPLATE,$length_temp2,$tmscore2,$rmsd2,$aln_length2,$seq_global,$seq_id2);
                $value2 = "ALIGNMENT $aln_temp2 $aln_targ2";
                #$value3 = sprintf("ROTMTRX %-6s %-16s %-16s %-16s %-16s %-16s %-16s %-16s %-16s %-16s %-16s %-16s %-16s",$TEMPLATE,$m10,$m11,$m12,$m13,$m20,$m21,$m22,$m23,$m30,$m31,$m32,$m33);
                ####$value3 = sprintf("ROTMTRX %-6s %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f",$TEMPLATE,$m10,$m11,$m12,$m13,$m20,$m21,$m22,$m23,$m30,$m31,$m32,$m33);
                $value3 = "ROTMTRX $TEMPLATE $m10 $m11 $m12 $m13 $m20 $m21 $m22 $m23 $m30 $m31 $m32 $m33";
                $key1 = $TEMPLATE;
                $h1{$key1} = $value1 ; 
                $h2{$key1} = $value2 ;$h20{$key1}=$aln_dots2;
                $h3{$key1} = $value3 ;
                $h4{$key1} = [@grep_2];

                 my @grep_3 = grep (/^HBN Donor:/,@file2);foreach (@grep_3) {$_ =~ s/HBN //g ;}
                 my @grep_4 = grep (/^SLT Cation/,@file2);foreach (@grep_4) {$_ =~ s/SLT //g ;}
                 my @grep_5 = grep (/^CON/,@file2);foreach (@grep_5) {$_ =~ s/CON //g ;}
                 $hbond{$key1}=[@grep_3] if [@grep_3] ;
                 $salt{$key1}=[@grep_4] if [@grep_4];
                 $contacts{$key1}=[@grep_5] if [@grep_4];
                }
               }
              unlink($fh1, $tmpfil1) if ( -e $tmpfil1 );
             }
          }
       } 
    }
   }
 }
 
 my @out_file_1 ;
 foreach my $y(keys %h2)
  {
   my @split_line1 = split(/\s+/,$h2{$y});
   my @split_line2 = split(/\s+/,$h1{$y});
   my $line1=">$y $split_line2[2] $split_line2[5] $split_line2[3] $split_line2[4] $split_line2[6]";
   my $line2=$split_line1[1];
   my $line4=$split_line1[2];
   my $line3=$h20{$y};
   my $line5="*";
   push (@out_file_1,"$line1\n$line2\n$line3\n$line4\n$line5\n");
  } 
 
 my $t1 = Benchmark->new;
 my $td01 = timediff($t1, $t0);
 
 if ( $k1 > 1 )
 {
  print " $k1 templates found\n\n";
 }
 elsif ( $k1 == 1 )
 {
  print " $k1 template found\n\n";
 }
 else
 {
  print " no templates found\n";
  
  my $t4 = Benchmark->new;
  
  printf("\n------------------------------------------------------------\n");
  printf("Walltime: %s\n", timestr(timediff($t4, $t0)));
  printf("------------------------------------------------------------\n");
  
  exit(0);
 }
 
##section3.2###############################################################################################
# This section identifies all the ASA residues of the input target residues
# Assigns residue interface propensity to all surface residues.

 my $y2 = substr($tmpfil0, 0, index($tmpfil0, '.')); 
 my $y3 = $y2.".rsa";
   sub naccess 
 {
   my $run_naccess =$env_naccess;
   open (FH, "$run_naccess $_[0] |");
     my @naccess_out = <FH>;
   close FH ; 
   my @rsa = read_file("$y3") if ( -e "$y3" );
   return (@rsa);
  }

 my $remark = 'REM';
 my (@rsa1)=naccess($tmpfil0);  
 my @rsa2 = grep( !/^$remark/,@rsa1); 

 my %res_name ;
 my %naccess; my $key_naccess; my $start = 'RES';
 my %rip ; 
 foreach my $y7 (@rsa2)
  {
   if ( substr($y7, 0, 3) eq "RES" )
    {
      my @split_naccess = split (/\s+/,$y7);
      my $key_naccess = $split_naccess[3];
      my $res_type = $split_naccess[1];
      $res_name{$key_naccess} = $res_type ;
      if ($split_naccess[5]>5)
       {
        my $ripx = 0 ; 
        my $nacess_scale = scale($split_naccess[5],5,125);
        $naccess{$key_naccess} = $nacess_scale;
         if ($res_type eq "TRP"){$ripx= 0.83;}
         elsif ($res_type eq "PHE"){$ripx= 0.82;}
         elsif ($res_type eq "MET"){$ripx= 0.66;}
         elsif ($res_type eq "TYR"){$ripx= 0.66;}
         elsif ($res_type eq "ILE"){$ripx= 0.44;}
         elsif ($res_type eq "CYS"){$ripx= 0.43;}
         elsif ($res_type eq "HIS"){$ripx= 0.41;}
         elsif ($res_type eq "LEU"){$ripx= 0.40;}
         elsif ($res_type eq "ARG"){$ripx= 0.27;}
         elsif ($res_type eq "VAL"){$ripx= 0.27;}
         elsif ($res_type eq "ASN"){$ripx= 0.12;}
         elsif ($res_type eq "ALA"){$ripx= -0.17;}
         elsif ($res_type eq "ASP"){$ripx= -0.38;}
         elsif ($res_type eq "GLN"){$ripx= -0.11;}
         elsif ($res_type eq "GLU"){$ripx= -0.13;}
         elsif ($res_type eq "GLY"){$ripx= -0.07;}
         elsif ($res_type eq "LYS"){$ripx= -0.36;}
         elsif ($res_type eq "PRO"){$ripx= -0.25;}
         elsif ($res_type eq "SER"){$ripx= -0.33;}
         elsif ($res_type eq "THR"){$ripx= -0.18;}
         my $ripx_scale = scale($ripx,-0.38,0.83);
         $rip{$key_naccess} = $ripx_scale;
       }
    }
  }  

##section3.3###############################################################################################
# This section  parses .prf file and finds the entropy of all surface residues of the input target residues 

 my @value_prf; my $key_prf=(); my %prf1;my %entropy;
 my $last_key_prf ;
 foreach my $y1 (@prf)
  {
    my $string1 = substr($y1,0,6);$string1 =~ s/^\s+//;$string1 =~ s/\s+$//;
    my @split_3 = split (/\s+/,$string1);
    $y1 =~ s/^\s+//;
    if (exists $split_3[1] && $split_3[1] =~ /[a-zA-Z]/)
     {
      $prf1{$key_prf} = [@value_prf] if $key_prf;
      my @split_4 = split(/\s+/,$y1);
      $key_prf = $split_4[0];
      @value_prf = ();
      foreach (3 .. $#split_4)
       {
         push(@value_prf,$split_4[$_]);
        }
      }
     else
      {
       my @split_5 = split(/\s+/,$y1);
       foreach (0 .. $#split_5)
        {
         push(@value_prf,$split_5[$_]);
         $last_key_prf = $key_prf;
         }      
       }
     }
 
 $prf1{$last_key_prf} = [@value_prf] ;

 foreach my $y2 (sort {$a<=>$b} keys %prf1)
   {
    my $E = 0 ;
    foreach (@{$prf1{$y2}})
     {
      if ($_ != 0)
       {
        my $E1 = ($_)*(log($_)/log(2));
        $E = $E + $E1 ;
       }
     }
    my $E_scale=scale(-$E,0,4.24870939634279);
    $entropy{$y2}=$E_scale;
    }

##section3.4###############################################################################################

 my %psip ;
 my @propensity = (0.43,0.66,0.82,0.44,0.40,0.27,0.83,0.66,-0.17,-0.07,-0.18,-0.33,-0.11,0.12,-0.13,-0.38,0.41,0.27,-0.36,-0.25);
 foreach my $y3 (keys %prf1)
   {
    my $PSIP = 0 ;
    for (0 ..19)
     {
      if (${$prf1{$y3}}[$_]!= 0)
       {
        my $E1 = ${$prf1{$y3}}[$_]*$propensity[$_];
        $PSIP = $PSIP + $E1 ;
       }
     }
    my $PSIP_scale = scale($PSIP,-0.38,0.83);
    $psip{$y3}=$PSIP_scale;
    }
 
##section4#################################################################################################

 my (%h5,%h6);
 
 foreach my $key_chk1 (keys %h1)
  {
      my @split_4 = split (/\s+/,$h1{$key_chk1});
      my $template_1 = $split_4[1];     
      my @temp_interface =  @{$h4{$key_chk1}};
      my @temp_inter_resi=(); my @temp_inter_resno=();
      foreach (@temp_interface)
       {
        my @split_5 = split(/\s+/,$_);
        my $inter_resi = $split_5[4];
        my $inter_res_no = substr($_,23,5);$inter_res_no=~ s/^\s+//;$inter_res_no=~s/\s+$//;
        push(@temp_inter_resi,$inter_resi);
        push(@temp_inter_resno,$inter_res_no);
       }
     $h5{$key_chk1} = [@temp_inter_resi];
     $h6{$key_chk1} = [@temp_inter_resno]; 
  }

##section5#################################################################################################

 my (%h7); my (%p1); my ($m1,$m2);  
 foreach my $key_chk2 (keys %h6)
  {
   my $q = $h2{$key_chk2};
   my @split_5 = split (/\s+/,$q);
   my $aln_temp3 = $split_5[1];
   my $aln_targ3 = $split_5[2];
   my @inter_resi = @{$h6{$key_chk2}};
   my $length = length ($aln_targ3);
   undef %p1;
   foreach (@inter_resi)
     {
      if (exists $p1{$_}){$p1{$_}=$p1{$_}+1;}
      else {$p1{$_}=1;}
      } 
    my $out = map_alignments($aln_targ3,$aln_temp3,\%p1, \%h7, $target_length);
    %h7 = %{$out};
    }

  my $h8_out = find_fraction(\%h7);
  my %h8 = %{$h8_out};my %h81;

  foreach (keys %naccess)
   { 
    if (!exists $h8{$_}) {$h81{$_}=-1;}
    else {$h81{$_}=scale($h8{$_},0,1)}  
   }
 
##section6#################################################################################################
# making input for lib-svm

 my @value9; my %h9;
 my @value10; my %h10;
 foreach my $y (keys %naccess)
  {
   @value9 = ();  @value10 = ();
   push (@value9,$naccess{$y});push (@value10,"1:".$naccess{$y});
   push (@value9,$rip{$y})  ;push (@value10,"2:".$rip{$y});
   push (@value9,$entropy{$y});push (@value10,"3:".$entropy{$y});
   push (@value9,$psip{$y});push (@value10,"4:".$psip{$y});
   push (@value9,$h81{$y}) if ($h81{$y}); push (@value10,"5:".$h81{$y}) if ($h81{$y});
   push (@value9,"-1") if !($h81{$y}); push (@value10,"5:-1") if !($h81{$y});
   $h9{$y} = [@value9]; $h10{$y} = [@value10];
  }



##section7.1###############################################################################################

print "Applying SVM ... ";

## Prepare data for svm_predict##

 my @svm_input;
 foreach (sort {$a<=>$b} keys %h10)
   {
    my $i = 0;
    foreach (@{$h10{$_}})
     {
      $i++;
      if ($i == 1) {push(@svm_input,"0 $_ ");}
      else {push(@svm_input,"$_ ") ;}
     }
    push(@svm_input,"\n") ;
   }

 
write_file("svm_input_scale",@svm_input );

## svm_predict ##

 my $svm_predict = $env_libsvm;  
 open (FH, "$svm_predict -b 1 svm_input_scale $env_svm svm.prediction 2>&1|");
   my @svm_out2 = <FH>;
 close FH ;  

 my @svm_prediction;
 if ( -e "svm.prediction" ) {@svm_prediction = read_file("svm.prediction") ;}
 else {die "svm.prediction failed"};

## parse svm.prediction ##
 my @svm_prediction1 = grep (!/labels/,@svm_prediction);
 my @svm_in_keys = sort {$a<=>$b} keys %h10 ;
 if ($#svm_in_keys != $#svm_prediction1){die "could not map keys of svm prediction";}
 
 my %sv_prediction ; 
 foreach (0..$#svm_prediction1)
  {
   my $key = $svm_in_keys[$_];
   my @split = split (/\s+/,$svm_prediction1[$_]);
   $sv_prediction{$key}=$split[2];
   }

 chdir $dir1;
 print "done\n\n";

##section7.2###############################################################################################

 print "Applying NBC ... ";

 my $nb_model = AI::NaiveBayes1->import_from_YAML_file("$env_nb");

 my %nb_prediction ;
foreach (sort {$a<=>$b} keys %h9)
 {
  my $p1 = 0.0;
  my $p2 = 0.0;
  my @array = @{$h9{$_}};
  my $p3 = $nb_model->predict(attributes=>{RSA=>$array[0],RIP=>$array[1],ENTROPY=>$array[2],PSIP=>$array[3],FRACTION=>$array[4]});
  $p1 = $p3->{'TRUE'};
  $p2 = $p3->{'FALSE'};
  $nb_prediction{$_}=$p1;
 } 
 my $t2 = Benchmark->new;
 my $td12 = timediff($t2, $t1);
 
 print "done\n\n";
 
##section8#################################################################################################
## This section predictshydrogen bonds, hydrophobic bcontacts , aromatic contacts an dsalt bridges.


 my (%hbond2, %salt2,%hydrophobic2,%aromatic2);
 my (%hbond3, %salt3,%hydrophobic3,%aromatic3); 
 my ($salt1,$key_s1);

 ## HBOND ##  
 foreach my $hkey1(keys %hbond)
  {
   my %hbond1 = (); 
   my $template_id = $hkey1 ; 
   my $template_chn = substr($template_id,4,1);
   my @template_hbond = @{$hbond{$hkey1}};

   my $q2 = $h2{$hkey1};
   my @split_h0 = split (/\s+/,$q2);
   my $aln_temp4 = $split_h0[1];
   my $aln_targ4 = $split_h0[2];
   my $length_aln = length ($aln_targ4);
    
   foreach my $h1(@template_hbond)
    {
     $h1=~ s/^\s+|\s+$//g;
     my @split_h1 = split(/\s+/,$h1);
     my $chn1_h1 = $split_h1[2]; my $resno1_h1 = $split_h1[3];
     my $chn2_h1 = $split_h1[7]; my $resno2_h1 = $split_h1[8];

     if ($chn1_h1 eq $template_chn)
       {$hbond1{$resno1_h1} = 1;}
     elsif ($chn2_h1 eq $template_chn)
       {$hbond1{$resno2_h1}= 1;}
     } 
    my $out1 = map_alignments($aln_targ4,$aln_temp4,\%hbond1, \%hbond2, $target_length);
    %hbond2 = %{$out1};
   }
  my $hbond3_out = find_fraction(\%hbond2);
  %hbond3 = %{$hbond3_out};
 
 ## SALT BRIDGES ##
 foreach my $hkey2(keys %salt)
  {
   my $template_id = $hkey2 ; 
   my $template_chn = substr($template_id,4,1);
   my @template_salt = @{$salt{$hkey2}};
   my $q2 = $h2{$hkey2};
   my @split_h0 = split (/\s+/,$q2);
   my $aln_temp4 = $split_h0[1];my $aln_targ4 = $split_h0[2];
   my $length_aln = length ($aln_targ4);
   my %salt1 = (); 
   foreach my $s1(@template_salt)
    {
     $s1=~ s/^\s+|\s+$//g;
     my @split_s1 = split(/\s+/,$s1);
     my $chn1_s1 = $split_s1[2]; my $resno1_s1 = $split_s1[3];
     my $chn2_s1 = $split_s1[7]; my $resno2_s1 = $split_s1[8];

     if ($chn1_s1 eq $template_chn)
       {$salt1{$resno1_s1} = 1;}
     elsif ($chn2_s1 eq $template_chn)
       {$salt1{$resno2_s1}= 1;}
     } 
    my $out2 = map_alignments($aln_targ4,$aln_temp4,\%salt1, \%salt2, $target_length);
    %salt2 = %{$out2};
   }
 my $salt3_out = find_fraction(\%salt2);
 %salt3 = %{$salt3_out}; 
 
  my @aromatic=qw/TYR PHE TRP HIS/;
  my @apolar=qw/ALA ILE LEU PHE VAL PRO GLY/;
  my @salt=qw/ARG LYS GLU ASP/;

 foreach my $hkey3(keys %contacts)
  {
   my $template_id = $hkey3 ; 
   my $template_chn = substr($template_id,4,1);
   my @template_contacts = @{$contacts{$hkey3}};

   my $q2 = $h2{$hkey3};
   my @split_h0 = split (/\s+/,$q2);
   my $aln_temp4 = $split_h0[1];my $aln_targ4 = $split_h0[2];
   my $length_aln = length ($aln_targ4);
   my %aromatic1 = () ; my %hydrophobic1 = () ;

   ## AROMATIC CONTACTS ## 
   foreach my $a1 (@template_contacts)
    {
     $a1=~ s/^\s+|\s+$//g;
     my @split_line_a1 = split(/\s+/,$a1);
     my $chn1=$split_line_a1[2];
     my $chn2=$split_line_a1[6];
     my $i_res1=$split_line_a1[0];
     my $i_res2=$split_line_a1[4];
     my $i_type=$split_line_a1[9];
     #my @hphobic=qw/TYR PHE TRP HIS/;
     #my @apolar=qw/ALA ILE LEU PHE VAL PRO GLY/;
     if($chn1 ne $chn2 && $i_type=~m/(S-S)/i && /$i_res1/i ~~ @aromatic && /$i_res2/i ~~ @aromatic )
      {
       my $chn1_resno=$split_line_a1[1] ;
       my $chn2_resno=$split_line_a1[5] ;
       if ($chn1 eq $template_chn)
        {
         $aromatic1{$chn1_resno} = 1;
         }
       elsif ($chn2 eq $template_chn)
        {
         $aromatic1{$chn2_resno}=1;
        }
      }
     if($chn1 ne $chn2 && $i_type=~m/(S-S)/i && /$i_res1/i ~~ @apolar && /$i_res2/i ~~ @apolar )
      {
       my $chn1_resno=$split_line_a1[1] ;
       my $chn2_resno=$split_line_a1[5] ;
       if ($chn1 eq $template_chn)
        {
         $hydrophobic1{$chn1_resno} = 1;
         }
       elsif ($chn2 eq $template_chn)
        {
         $hydrophobic1{$chn2_resno}=1;
        }
      }
    }
    my $out3 = map_alignments($aln_targ4,$aln_temp4,\%aromatic1, \%aromatic2, $target_length);
    %aromatic2 = %{$out3}; 
    my $out4 = map_alignments($aln_targ4,$aln_temp4,\%aromatic1, \%aromatic2, $target_length);
    %hydrophobic2 = %{$out4};
 }

  my $aromatic3_out = find_fraction(\%aromatic2);
  %aromatic3 = %{$aromatic3_out};

 my $hydrophobic3_out = find_fraction(\%hydrophobic2);
 %hydrophobic3 = %{$hydrophobic3_out};

 my @out_file_2;my $j1 = 0; my $j2 = 0 ; my %star;
 my $confidence_sum=0;my $confidence =0;  
 foreach ( sort{$a<=>$b} keys %naccess)
 {
   $j1++;
  if ($sv_prediction{$_}>0.202 && $nb_prediction{$_}>0.178 && exists $h8{$_} )
   {
    push (@out_file_2,sprintf("RESIDUE %4s %4d %s %4d %8.5f %8.5f %9.6f\n","<*>",$j1,$res_name{$_},$_,$h8{$_},$sv_prediction{$_},$nb_prediction{$_}));
    $star{$_}=1;
    $j2++;
    $confidence_sum = $confidence_sum + $sv_prediction{$_}*$nb_prediction{$_};
   }
  elsif ( exists $h8{$_} ) 
   {
    push (@out_file_2,sprintf("RESIDUE %4s %4d %s %4d %8.5f %8.5f %9.6f\n"," ",$j1,$res_name{$_},$_,$h8{$_},$sv_prediction{$_},$nb_prediction{$_}));
   }
 } 
 
 $confidence = $confidence_sum/$j2 if ($j2>0) ;

 my $confidence_out = 0;
 if ($confidence>=0.5){$confidence_out="HIGH";}
 elsif ($confidence>=0.25 && $confidence<0.5){$confidence_out="MEDIUM";}
 else {$confidence_out="LOW";}
 unshift (@out_file_2, "CONFDNC $confidence_out\n");
 unshift (@out_file_2, "LIBRARY $version1[0]\n");
 
 print "Prediction confidence: $confidence_out\n"; 

 foreach (sort{$a cmp $b} keys %h1){push (@out_file_2,"$h1{$_}\n");}
 foreach (sort{$a cmp $b} keys %h3){push (@out_file_2,"$h3{$_}\n");}

 foreach (sort{$a <=> $b} keys %hbond3){
   if ($hbond3{$_} > 0 && exists $naccess{$_}){
    if (exists $star{$_} && $hbond3{$_} > 0.041){ push (@out_file_2, sprintf("INTRCTN HBND %-4s%-6s %-6s %9.6f\n","<*>",$res_name{$_},$_,$hbond3{$_}));}
    else {push (@out_file_2, sprintf("INTRCTN HBND %-4s%-6s %-6s %9.6f\n"," ",$res_name{$_},$_,$hbond3{$_}));}
   }}

 foreach (sort{$a <=> $b} keys %salt3){
  if ($salt3{$_} > 0 && exists $naccess{$_}  && /$res_name{$_}/i ~~ @salt ) {
   if (exists $star{$_} && $salt3{$_} > 0.006 ) { push (@out_file_2, sprintf("INTRCTN SALT %-4s%-6s %-6s %9.6f\n","<*>",$res_name{$_},$_,$salt3{$_}));}
   else { push (@out_file_2, sprintf("INTRCTN SALT %-4s%-6s %-6s %9.6f\n"," ",$res_name{$_},$_,$salt3{$_}));}
  }}

 foreach (sort{$a <=> $b} keys %hbond3){
  if ($hydrophobic3{$_}>0 && exists $naccess{$_} && /$res_name{$_}/i ~~ @apolar){
   if (exists $star{$_} && $hydrophobic3{$_}>0.005) { push (@out_file_2, sprintf("INTRCTN HYFB %-4s%-6s %-6s %9.6f\n","<*>",$res_name{$_},$_,$hydrophobic3{$_}));}
    else { push (@out_file_2, sprintf("INTRCTN HYFB %-4s%-6s %-6s %9.6f\n"," ",$res_name{$_},$_,$hydrophobic3{$_}));}
   }}

 foreach (sort{$a <=> $b} keys %hbond3){
  if ($aromatic3{$_} > 0 && exists $naccess{$_} && /$res_name{$_}/i ~~ @aromatic){ 
   if  (exists $star{$_} && $aromatic3{$_} > 0.005 ) {push (@out_file_2, sprintf("INTRCTN AROM %-4s%-6s %-6s %9.6f\n","<*>",$res_name{$_},$_,$aromatic3{$_}));}
    else {push (@out_file_2, sprintf("INTRCTN AROM %-4s%-6s %-6s %9.6f\n"," ",$res_name{$_},$_,$aromatic3{$_}));}
   }} 

 my $out_file_1_name;my $out_file_2_name;
 $out_file_1_name= $target_id.".alignments.dat";
 write_file($out_file_1_name,@out_file_1);
 $out_file_2_name= $target_id.".sites.dat";
 write_file($out_file_2_name,@out_file_2);
##section9#################################################################################################

 my $t3 = Benchmark->new;
 
 printf("\n------------------------------------------------------------\n");
 printf("Walltime: %s\n", timestr(timediff($t3, $t0)));
 printf("------------------------------------------------------------\n");
 
 exit(0);

__DATA__
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
A  5 -2 -1 -2 -1 -1 -1  0 -2 -1 -2 -1 -1 -3 -1  1  0 -3 -2  0 -2 -1 -1
R -2  7 -1 -2 -4  1  0 -3  0 -4 -3  3 -2 -3 -3 -1 -1 -3 -1 -3 -1  0 -1
N -1 -1  7  2 -2  0  0  0  1 -3 -4  0 -2 -4 -2  1  0 -4 -2 -3  4  0 -1
D -2 -2  2  8 -4  0  2 -1 -1 -4 -4 -1 -4 -5 -1  0 -1 -5 -3 -4  5  1 -1
C -1 -4 -2 -4 13 -3 -3 -3 -3 -2 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1 -3 -3 -2
Q -1  1  0  0 -3  7  2 -2  1 -3 -2  2  0 -4 -1  0 -1 -1 -1 -3  0  4 -1
E -1  0  0  2 -3  2  6 -3  0 -4 -3  1 -2 -3 -1 -1 -1 -3 -2 -3  1  5 -1
G  0 -3  0 -1 -3 -2 -3  8 -2 -4 -4 -2 -3 -4 -2  0 -2 -3 -3 -4 -1 -2 -2
H -2  0  1 -1 -3  1  0 -2 10 -4 -3  0 -1 -1 -2 -1 -2 -3  2 -4  0  0 -1
I -1 -4 -3 -4 -2 -3 -4 -4 -4  5  2 -3  2  0 -3 -3 -1 -3 -1  4 -4 -3 -1
L -2 -3 -4 -4 -2 -2 -3 -4 -3  2  5 -3  3  1 -4 -3 -1 -2 -1  1 -4 -3 -1
K -1  3  0 -1 -3  2  1 -2  0 -3 -3  6 -2 -4 -1  0 -1 -3 -2 -3  0  1 -1
M -1 -2 -2 -4 -2  0 -2 -3 -1  2  3 -2  7  0 -3 -2 -1 -1  0  1 -3 -1 -1
F -3 -3 -4 -5 -2 -4 -3 -4 -1  0  1 -4  0  8 -4 -3 -2  1  4 -1 -4 -4 -2
P -1 -3 -2 -1 -4 -1 -1 -2 -2 -3 -4 -1 -3 -4 10 -1 -1 -4 -3 -3 -2 -1 -2
S  1 -1  1  0 -1  0 -1  0 -1 -3 -3  0 -2 -3 -1  5  2 -4 -2 -2  0  0 -1
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  2  5 -3 -2  0  0 -1  0
W -3 -3 -4 -5 -5 -1 -3 -3 -3 -3 -2 -3 -1  1 -4 -4 -3 15  2 -3 -5 -2 -3
Y -2 -1 -2 -3 -3 -1 -2 -3  2 -1 -1 -2  0  4 -3 -2 -2  2  8 -1 -3 -2 -1
V  0 -3 -3 -4 -1 -3 -3 -4 -4  4  1 -3  1 -1 -3 -2  0 -3 -1  5 -4 -3 -1
B -2 -1  4  5 -3  0  1 -1  0 -4 -4  0 -3 -4 -2  0  0 -5 -3 -4  5  2 -1
Z -1  0  0  1 -3  4  5 -2  0 -3 -3  1 -1 -4 -1  0 -1 -2 -2 -3  2  5 -1
X -1 -1 -1 -1 -2 -1 -1 -2 -1 -1 -1 -1 -1 -2 -2 -1  0 -3 -1 -1 -1 -1 -1
