#! usr/local/bin/perl

#THIS SCRIPT IS RELATED TO THE FOLLOWING STUDY PUBLISHED ON SCIENTIFIC REPORTS:

#Armanhi, J. S. L.; de Souza, R. S. C.; de Araujo, L. M.; Okura, V. K.; Mieczkowski, P.; Imperial, J.; Arruda, P. Multiplex amplicon sequencing for microbe identification in community-based culture collections. Sci. Rep. 6, 29543; doi:10.1038/srep29543 (2016) (www.nature.com/articles/srep29543/).

#If you use this script or part of it please cite us.

use strict;

if (!(scalar (@ARGV) == 8)){
   die ("\n\nUSAGE: perl <program> -i <input_directory> -o <output_root_directory> -ch <chimera_filtering_database> -db <database_for_16S_filter>\n\n
(1) The first step of this script is to create folders which will harbor files from its different steps.
(2) The script recognizes a directory that contains .demult files and concatenates them into a .cat file.
(3) UCHIME algorithm is run to the .cat file and chimeras are discarded.
(4) Sequences are alignment against a database (0.75 id) and no-hit-sequences (nonspecific amplification) are discarded.
(5) Sequences are then redemultiplexed (based on its head, from demultiplexing script).
(6) A report is created to summarize steps during the process.\n\n\n");
}

my $usearch = '$uv8';
my @folders = ('01_CAT_DEMULTIPLEXED','02_CHIMERA_FILTERING','03_ALIGNMENT_FOR_16S_FILTER', '04_16S_FILTERED');#,'05_REDEMULTIPLEXED'); #Folders to be created
my $inputDir; #Input directory from -i parameter
my $outputDir; #Output root directory from -o parameter
my $chimeraDb; #Database for chimera filtering from -d parameter
my $identity; #Identity value from -id parameter for OTU clustering
my $db16S; #Database directory/file to database for 16S filtering
my $demultPattern = "*.demult"; #Pattern searched of file that contains demultiplexed well (use bc.._....demult for demultiplexed files and .demult for cat)
my $nonchimerasPattern = ".demult.cat.nonchimeras"; #Pattern searched of file that contains nonchimeric sequences from demultiplexed well
my $filterPattern = ".redemult"; #".demult.cat.nonchimeras.16Sfiltered";

sub makeFolder {
   my @params = @_;
   my @f = @{$params[0]};
   mkdir "$params[1]";
   foreach my $var (@f){
      mkdir "$params[1]/$var";
   }
}

my $out;
sub catchARGV {
   my @params = @_;
   my $cont = 0;
   my $in;
   my $ch;
   my $ident;
   my $db;
   while ($cont < (scalar @params)){
      if ($params[$cont] eq '-i'){
         $in = $params[$cont+1];
      }
      if ($params[$cont] eq '-o'){
         $out = $params[$cont+1];
      }
      if ($params[$cont] eq '-ch'){
         $ch = $params[$cont+1];
      }
      if ($params[$cont] eq '-id'){
         $ident = $params[$cont+1];
      }
      if ($params[$cont] eq '-db'){
         $db = $params[$cont+1];
      }
      $cont++;
   }
   return ($in, $out, $ch, $ident, $db);
}

sub chimeraFilter {
   my @params = @_;
   print "\n\nCommand line:\n\n$usearch -uchime_ref $params[0]/$params[1] -db $params[3] -strand plus -nonchimeras $params[2]/$params[1].nonchimeras -chimeras $params[2]/$params[1].chimeras -uchimeout $params[2]/$params[1].uchimeout -uchimealns $params[2]/$params[1].uchimealns -threads 120\n\n";
   my $results = system ("$usearch -uchime_ref $params[0]/$params[1] -db $params[3] -strand plus -nonchimeras $params[2]/$params[1].nonchimeras -chimeras $params[2]/$params[1].chimeras -uchimeout $params[2]/$params[1].uchimeout -uchimealns $params[2]/$params[1].uchimealns -threads 120");
   if ($results != 0){\
      die "\n\nERROR: Chimera filter error.\n";
      }
}

my %barcodes;
sub idCount {
   my @params = @_;
   open (my $input, $params[0]);
      while (<$input>){
         chomp;
         my $linhaAtual = $_;
         if ($linhaAtual =~ m/>/){
            my @barcodeLabel = split (";", "$linhaAtual");
            my @Label = split ("barcodelabel=", $barcodeLabel[1]);
            $barcodes{$Label[1]}[$params[1]]=$barcodes{$Label[1]}[$params[1]]+1;
         }
      }
   close ($input, $params[0]);
}

sub filter {
   my @params = @_;
   opendir (my $dirFilter, $params[0]);
      while (my $inputFile = readdir ($dirFilter)){
         if ($inputFile =~ m/$nonchimerasPattern/){
            print "\n\nCommand line:\n\n$usearch -usearch_global $params[0]/$inputFile -db $params[3] -strand both -id 0.75 -blast6out $params[1]/$inputFile.m8 -threads 120\n\n";
            my $results = system "$usearch -usearch_global $params[0]/$inputFile -db $params[3] -strand both -id 0.75 -blast6out $params[1]/$inputFile.m8 -threads 120";
            if ($results != 0){
               die "\n\nERROR: Usearch_global for 16S filter error.\n";
            }
         }
      }
   close ($dirFilter);
   my %hashFiltering;
   opendir (my $dirFiltering, $params[1]);
      while (my $inputFile = readdir ($dirFiltering)){
         open (my $file, "$params[1]/$inputFile");
            while (<$file>){
               my $linhaAtual = $_;
               my @filterColumn = split ("\t", $linhaAtual);
               $hashFiltering{$filterColumn[0]}=1;
            }
         close ($file, "$params[1]/$inputFile");
      }
   close ($dirFiltering);
   opendir (my $dirRead, $params[0]);
      while (my $inputRead = readdir ($dirRead)){
         if ($inputRead =~ m/$nonchimerasPattern/){
            open (my $fileRead, "$params[0]/$inputRead");
               open (my $file, ">", "$params[2]/$inputRead.16Sfiltered");
                  my $fileTemp;
                  my $idPrinted = 0;
                  while (<$fileRead>){
                     chomp;
                     my $linhaAtualRead = $_;
                     if ($linhaAtualRead =~ m/>/){
                        my $idRead = substr($linhaAtualRead, 1);
                        if (exists($hashFiltering{$idRead})){
                           $fileTemp = "$fileTemp"."$linhaAtualRead\n";
                           $idPrinted = 1;
                        }else{
                           $idPrinted = 0;
                        }
                     }elsif ($idPrinted == 1){
                        $fileTemp = "$fileTemp"."$linhaAtualRead\n";
                     }
                  }
                  print {$file} $fileTemp;
               close ($file);
            close ($fileRead);
         }
      }
   close ($dirRead);
}

sub redemult {
   my @params = @_;
   open (my $inputFile, $params[0]);
      my $fileName;
      my $escrevi;
      my $ondeEscrevi;
      my %hashRedemult;
      while (<$inputFile>){
         chomp;
         my @id;
         my $linhaAtual = $_;
         if ($linhaAtual =~ m/>/){
            $escrevi = 0;
            @id = split ("barcodelabel=", $linhaAtual);
            $fileName = "$id[1].redemult";
            if (!(exists $hashRedemult{$fileName})){
               open (my $arq, ">", "$params[1]/$fileName");
                  $hashRedemult{$fileName} = $arq;
               #close ($arq);
            }
            print {$hashRedemult{$fileName}} "$linhaAtual\n";
            $escrevi = 1;
            $ondeEscrevi = $fileName;
         }elsif ($escrevi == 1){
            print {$hashRedemult{$ondeEscrevi}} "$linhaAtual\n";
         }
      }
   close ($inputFile);
}

#Catch parameters
print "\n***SCRIPT STARTING";
print "\n\n***Catching arguments";
($inputDir, $outputDir, $chimeraDb, $identity, $db16S) = catchARGV (@ARGV);

#Create folder in output directory
print "\n\n***Creating folders";
makeFolder(\@folders, $outputDir);

#Cat for .demult files
print "\n\n***Cat demultiplexed files";
my @concatenado = split ("/", $out);
my $name = pop @concatenado;
my $cat = "$name".".demult.cat";
opendir (my $dirChimera, $inputDir);
   my $results = system ("cat $inputDir/$demultPattern > $outputDir/01_CAT_DEMULTIPLEXED/$cat");
   if ($results != 0){
      die "\n\nERROR: Cat .demult error.\n";
   }
closedir ($dirChimera);

#Create chimera filtering
print "\n\n***Removing chimeras";
opendir (my $dirChimera, "$outputDir/01_CAT_DEMULTIPLEXED");
   while (my $inputFile = readdir ($dirChimera)){	
      if ($inputFile =~ m/$cat/){
         chimeraFilter("$outputDir/01_CAT_DEMULTIPLEXED", $inputFile, "$outputDir/02_CHIMERA_FILTERING", $chimeraDb);
         idCount("$outputDir/01_CAT_DEMULTIPLEXED/$inputFile", "0");
      }
   }	
closedir ($dirChimera);

#Reading chimera.filtered.file to count
opendir (my $countChime, "$outputDir/02_CHIMERA_FILTERING");
   while (my $inputFile = readdir ($countChime)){
      if ($inputFile =~ m/$nonchimerasPattern/){
      idCount("$outputDir/02_CHIMERA_FILTERING/$inputFile","1");
      }
   }
closedir ($countChime);

#Runing 16S filter
print "\n\n***Filtering 16S";
filter("$outputDir/02_CHIMERA_FILTERING", "$outputDir/03_ALIGNMENT_FOR_16S_FILTER", "$outputDir/04_16S_FILTERED", $db16S);

#Reading 16S.filtered.file to count
opendir (my $count16S, "$outputDir/04_16S_FILTERED");
   while (my $inputFile = readdir ($count16S)){
      if ($inputFile =~ m/$nonchimerasPattern/){
         idCount ("$outputDir/04_16S_FILTERED/$inputFile", "2");
      }
   }
closedir ($count16S);

#printar hash que contabiliza quantos arquivos sobram depois de demultiplex, chimera, e 16S filtering.
open (my $report, ">>", "$out/report.report");
   print {$report} "SHORT REPORT\n";
   print {$report} "WELL\t\tDEMULT\tNONCHIM\t16S\n";
   my @barcodesArray = sort keys %barcodes;
   foreach my $var (@barcodesArray){
      print {$report} "$var\t$barcodes{$var}[0]\t$barcodes{$var}[1]\t$barcodes{$var}[2]\n";
   }
close ($report, "$out/report.report");
