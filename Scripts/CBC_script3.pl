#! usr/local/bin/perl

#THIS SCRIPT IS RELATED TO THE FOLLOWING STUDY PUBLISHED ON SCIENTIFIC REPORTS:

#Armanhi, J. S. L.; de Souza, R. S. C.; de Araujo, L. M.; Okura, V. K.; Mieczkowski, P.; Imperial, J.; Arruda, P. Multiplex amplicon sequencing for microbe identification in community-based culture collections. Sci. Rep. 6, 29543; doi:10.1038/srep29543 (2016) (www.nature.com/articles/srep29543/).

#If you use this script or part of it please cite us.

use strict;

my $usearch = '$uv8';
my @name_array = split ("/",$ARGV[0]);
my $name = pop @name_array;

if (!(scalar(@ARGV) == 3)){
   die ("\n\nUSAGE: perl script.pl <dataset> <16S_database> <output_dir>\n\n
(1) This scripts uses a .fasta file as input.
(2) The first step is to align these sequences against a database (usearch_global) with identity 0.97.
(3) Sequences without hit are then aligned to the original .fasta file (usearch_global) with identity 0.97.
(4) Now sequences without hit against the original dataset are discarded.
(5) The output is a file (.filtered) which comprises original dataset of sequences except those discarded.\n\n\n");
}

#Create folders
my @folders = ($ARGV[2]);
sub makeFolder {
   my @params = @_;
   my @f = @{$params[0]};
   mkdir "$params[1]";
   foreach my $var (@f){
      mkdir "$params[1]/$var";
   }
}
makeFolder(\@folders, $ARGV[2]);

#1.Alignment
my $results=system ("$usearch -usearch_global $ARGV[0] -db $ARGV[1] -strand both -id 0.97 -blast6out $ARGV[2]/$name.db_m8 -alnout $ARGV[2]/$name.db_aln -threads 120");
if ($results != 0){
   die ("\nERROR: dataset vs database alignment");
}

#2.1.Creating file of IDs of CCSs that have match
open (my $file3, ">", "$ARGV[2]/$ARGV[0].db_hit.id");
   open (FILE1, "$ARGV[2]/$name.db_m8");
      my $ids;
      while (<FILE1>){
         my $linhaAtual1 = $_;
         my @pad = split ("\t", $linhaAtual1);
         print {$file3} ">"."$pad[0]\n";	
         $ids = "$ids"."/".">"."$pad[0]";
      }
      $ids = "$ids"."/";
   close (FILE1, "$ARGV[2]/$name.db_m8");
close ($file3);

#2.2.Looking for .fasta of IDs that have no match
open (my $file4, ">", "$ARGV[2]/$name.db_nohit.fasta");
   open (DATASET, $ARGV[0]);
      my $escrevi = 0;
      while (<DATASET>){
         chomp;
         my $linhaAtual2 = $_;
         if ($linhaAtual2 =~ m/>/){
            if (!($ids =~ m/$linhaAtual2/)){
               print {$file4} "$linhaAtual2\n";
               $escrevi = 1;
            }
         }else{
            if ($escrevi == 1){
               print {$file4} "$linhaAtual2\n";
               $escrevi = 0;
            }
         }
      }
   close (DATASET, $ARGV[0]);
close ($file4);

#2.3 Creating file .db_hit.fasta
open (my $file41, ">", "$ARGV[2]/$name.db_hit.fasta");
   open (DATASET5, $ARGV[0]);
      my $escrevi4 = 0;
      while (<DATASET5>){
         chomp;
         my $linhaAtual3 = $_;
         if ($linhaAtual3 =~ m/>/){
            if (($ids =~ m/$linhaAtual3/)){
               print {$file41} "$linhaAtual3\n";
               $escrevi4 = 1;
            }
         }else{
            if ($escrevi4 == 1){
               print {$file41} "$linhaAtual3\n";
               $escrevi4 = 0;
            }
         }
      }
   close (DATASET5, $ARGV[0]);
close ($file41);

#3.Alignment against the dataset
my $results = system ("$usearch -usearch_global $ARGV[2]/$name.db_nohit.fasta -db $ARGV[0] -strand both -blast6out $ARGV[2]/$name.db_nohit.fasta.dataset_m8 -alnout $ARGV[2]/$name.db_nohit.fasta.dataset_aln -id 0.97 -top_hits_only -self");
if ($results != 0){
   die ("\nERROR: dataset vs dataset alignment");
}

#4.1.Creating file with IDs that have match
open (my $file5, ">", "$ARGV[2]/$name.db_nohit.fasta.dataset_hit.id");
   open (FILE2, "$ARGV[2]/$name.db_nohit.fasta.dataset_m8");
      my $ids2;
      while (<FILE2>){
         my $linhaAtual4 = $_;
         my @pad = split ("\t", $linhaAtual4);
         print {$file5} ">"."$pad[0]\n";
         $ids2 = "$ids2"."/".">"."$pad[0]";
      }
      $ids2 = "$ids2"."/";
   close (FILE2, "$ARGV[2]/$name.db_nohit.fasta.dataset_m8");
close ($file5);

#4.2.Looking for IDs that have no match
open (my $file6, ">", "$ARGV[2]/$name.db_nohit.fasta.dataset_nohit.fasta");
   open (DATASET2, "$ARGV[2]/$name.db_nohit.fasta");
      my $escrevi2 = 0;
      while (<DATASET2>){
         chomp;
         my $linhaAtual5 = $_;
         if ($linhaAtual5 =~ m/>/){
            if (!($ids2 =~ m/$linhaAtual5/)){
               print {$file6} "$linhaAtual5\n";
               $escrevi2 = 1;
            }
         }else{
            if ($escrevi2 == 1){
               print {$file6} "$linhaAtual5\n";
               $escrevi2 = 0;
            }
         }
      }
   close (DATASET2, "$ARGV[2]/$name.db_nohit.fasta");
close ($file6);

#5.1 Annotating IDs that have no match against both databases
open (FILE3, "$ARGV[2]/$name.db_nohit.fasta.dataset_nohit.fasta");
   my $ids3;
   while (<FILE3>){
      chomp;
      my $linhaAtual6 = $_;
      if ($linhaAtual6 =~ m/>/){
         $ids3 = "$ids3"."/"."$linhaAtual6";
      }
   }
   $ids3 = "$ids3"."/";
close (FILE3, "$ARGV[2]/$name.db_nohit.fasta.dataset_nohit.fasta");

#5.2 Reading original seq file and looking for those seqs with no match against both databases
open (my $file7, ">", "$ARGV[2]/$name.filtered");
   open (DATASET3, $ARGV[0]);
      my $escrevi3 = 0;
      while (<DATASET3>){
         chomp;
         my $linhaAtual7 = $_;
         if ($linhaAtual7 =~ m/>/){
            if (!($ids3 =~ m/$linhaAtual7/)){
               print {$file7} "$linhaAtual7\n";
               $escrevi3 = 1;
            }
         }else{
            if ($escrevi3 == 1){
               print {$file7} "$linhaAtual7\n";
               $escrevi3 = 0;
            }
         }
      }
   close (DATASET3, $ARGV[0]);
close ($file7);
