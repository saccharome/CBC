#! usr/local/bin/perl

#THIS SCRIPT IS RELATED TO THE FOLLOWING STUDY PUBLISHED ON SCIENTIFIC REPORTS:

#Armanhi, J. S. L.; de Souza, R. S. C.; de Araujo, L. M.; Okura, V. K.; Mieczkowski, P.; Imperial, J.; Arruda, P. Multiplex amplicon sequencing for microbe identification in community-based culture collections. Sci. Rep. 6, 29543; doi:10.1038/srep29543 (2016) (www.nature.com/articles/srep29543/).

#If you use this script or part of it please cite us.

use strict;

if (!(scalar (@ARGV) == 6)){
        die ("\n\nUSAGE:\nperl <program> -i <input_directory> -o <output_root_directory> -id <uparse_ID>\n\n
(1) This script creates folders in the input directory.
(2) It will perform a pipeline based on UPARSE (see ref. on github page) when dereplication, sorting and clustering.
(3) OTUs obtained on the fist step of clustering (withing wells) will be re-clustered.
(4) Then, a final OTU table will be obtained from OTUs of the second-step of clustering.
\n\n");
}

my $usearch = '$uv8';
my @folders = ('06_PER_WELL_DEREPLICATION', '07_PER_WELL_SORTING', '08_PER_WELL_CLUSTERING','09_RELABEL_OTU_WELL', '09.1_MAPPING_READS', '11_DEREPLICATION_CAT_OTU_WELL', '12_SORTING_CAT_OTU_WELL', '13_CAT_OTU_WELL_CLUSTERING', '10_CAT_OTU_WELL', '15_OTU_COLLECTION_MAP', '16_OTU_COLLECTION_TABLE', '14_OTU_COLLECTION_RELABEL'); #Folders to be created
my $inputDir; #Input directory from -i parameter
my $outputDir; #Output root directory from -o parameter
my $chimeraDb; #Database for chimera filtering from -d parameter
my $identity; #Identity value from -id parameter for OTU clustering
my $db16S; #Database directory/file to database for 16S filtering
my $demultPattern = "*.demult"; #Pattern searched of file that contains demultiplexed well (use bc.._....demult for demultiplexed files and .demult for cat)
my $nonchimerasPattern = ".demult.cat.nonchimeras"; #Pattern searched of file that contains nonchimeric sequences from demultiplexed well
my $redemult = ".redemult"; #".demult.cat.nonchimeras.16Sfiltered";
my $derepPattern = ".redemult.derep"; #".demult.cat.nonchimeras.16Sfiltered.derep"; #Pattern searched of file that cintains derep sequences from nonchimeric wells
my $sortPattern = ".redemult.derep.sorted"; #".demult.cat.nonchimeras.16Sfiltered.derep.sorted"; #Pattern searched of file that contains sorted sequences frm derep ones
my $relabelPattern = "otu.relabeled";
my $relabelCatPattern = "*otu.relabeled";
my $otusWellPattern = "otuswell.cat";
my $derepCollPattern = "otuswell.cat.derep";
my $sortCollPattern = "otuswell.cat.derep";
my $otuCollPattern = "otuswell.cat.derep.sorted.....otu";

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
      if ($params[$cont] eq '-id'){
         $ident = $params[$cont+1];
      }
      $cont++;
   }
   return ($in, $out, $ident);
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
   close ($input);
}

sub uparsePipelineDerep {
   my @params = @_;
   print "\n\nCommand line:\n\n$usearch -derep_fulllength $params[0]/$params[1] -fastaout $params[2]/$params[1].derep -sizeout -threads 120\n\n";
   my $results = system ("$usearch -derep_fulllength $params[0]/$params[1] -fastaout $params[2]/$params[1].derep -sizeout -threads 120");
   if ($results != 0){
      die "\n\nERROR: UparseDerep error.\n";
   }
}

sub uparsePipelineSort {
   my @params = @_;
   print "\n\nCommand line:\n\n$usearch -sortbysize $params[0]/$params[1] -fastaout $params[2]/$params[1].sorted\n\n";
   my $results = system ("$usearch -sortbysize $params[0]/$params[1] -fastaout $params[2]/$params[1].sorted");
   if ($results != 0){
      die "\n\nERROR: UparseSort error.\n";
   }
}

sub uparsePipelineClustering {
   my @params = @_;
   print "\n\nCommand line:\n\n$usearch -cluster_otus $params[0]/$params[1] -id $params[3] -otus $params[2]/$params[1].$params[3]otu -strand both\n\n";
   my $results = system ("$usearch -cluster_otus $params[0]/$params[1] -id $params[3] -otus $params[2]/$params[1].$params[3]otu -strand both");
   if ($results != 0){
      die "\n\nERROR: UparseCluster error.\n";
   }
}

sub relabelOtuWell {
   my @params = @_;
   open (my $file, $params[0]);
      open (my $fileRelabel, ">", $params[1]);
         my $otuNumber = 1;
         while (<$file>){
            chomp;
            my $linhaAtual = $_;
            if ($linhaAtual =~ m/>/){
               my @otuWellId = split (";", $linhaAtual);
               my $idSubstr = substr ($otuWellId[0], 1);
               my $relabelId = ">OTUwell_"."$otuNumber"."_"."$idSubstr".";"."$otuWellId[1]".";";
               print {$fileRelabel} "$relabelId\n";
               $otuNumber++;
            }else{
               print {$fileRelabel} "$linhaAtual\n";
            }
         }
      close ($fileRelabel);
   close ($file);
}

sub relabelOTUcollection {
   my @params = @_;
   print "\n\nCommand line: python /home/jaderson.armanhi/programs/src/drive5_py/fasta_number.py $params[0] OTU_ > $params[1]";
   my $results = system ("python /home/jaderson.armanhi/programs/src/drive5_py/fasta_number.py $params[0] OTU_ > $params[1]");
}

sub uparseMap {
   my @params = @_;
   print ("\n\nCommand line:$usearch -usearch_global $params[0] -db $params[1] -strand both -id 0.97 -uc $params[2]");
   my $results = system ("$usearch -usearch_global $params[0] -db $params[1] -strand both -id 0.97 -uc $params[2]");
   if ($results != 0){
      die ("\nUparse map reads error");
   }
}

sub uparseTable {
   my @params = @_;
   print "\n\nCommand line: python /home/jaderson.armanhi/programs/src/drive5_py/uc2otutab.py $params[0] > $params[1]";
   my $results = system ("python /home/jaderson.armanhi/programs/src/drive5_py/uc2otutab.py $params[0] > $params[1]");
   if ($results != 0){
      die ("\nUparse otu_table error");
   }
}

#Catch parameters
print "\n***SCRIPT STARTING";
print "\n\n***Catching arguments";
($inputDir, $outputDir, $identity) = catchARGV (@ARGV);

#Create folder in output directory
print "\n\n***Creating folders";
makeFolder(\@folders, $outputDir);

#Runing well Uparse dereplication
opendir (my $dirDerep, $inputDir);
   while (my $inputFile = readdir ($dirDerep)){
      if ($inputFile =~ m/$redemult/){
         my $inputDirFile = "$inputDir/$inputFile";
         if (-s $inputDirFile){
            print "\n\n***Uparse dereplication";
            uparsePipelineDerep($inputDir, $inputFile, "$outputDir/06_PER_WELL_DEREPLICATION");
         }
      }
   }
closedir ($dirDerep);

#Runing well Uparse sorting
opendir (my $dirSort, "$outputDir/06_PER_WELL_DEREPLICATION");
   while (my $inputFile = readdir ($dirSort)){
      if ($inputFile =~ m/$derepPattern/){
         print "\n\n***Uparse sorting";
         uparsePipelineSort("$outputDir/06_PER_WELL_DEREPLICATION", $inputFile, "$outputDir/07_PER_WELL_SORTING");
      }
   }
closedir ($dirSort);

#Runing well Uparse clustering
opendir (my $dirClust, "$outputDir/07_PER_WELL_SORTING");
   while (my $inputFile = readdir ($dirClust)){
      if ($inputFile =~ m/$sortPattern/){
         print "\n\n***Uparse clustering";
         uparsePipelineClustering("$outputDir/07_PER_WELL_SORTING", $inputFile, "$outputDir/08_PER_WELL_CLUSTERING", $identity);
         print "\n\n***Relabeling";
         relabelOtuWell("$outputDir/08_PER_WELL_CLUSTERING/$inputFile.$identity"."otu", "$outputDir/09_RELABEL_OTU_WELL/$inputFile.$identity"."otu".".relabeled");
      }
   }
closedir ($dirClust);

#Creating file .uc for each well
my $ucPattern = ".derep.sorted.0.97otu.relabeled";
opendir (my $dirFasta, "$outputDir/05_REDEMULTIPLEXED");
   opendir (my $dirOTU, "/$outputDir/09_RELABEL_OTU_WELL");
      while (my $inputFile = readdir ($dirFasta)){
         if ($inputFile =~ m/$redemult/){
            my $inputFile2 = "$inputFile"."$ucPattern";
            print "\n\n***Uparse uc files for wells";
            print ("\n\nCommand line:$usearch -usearch_global $outputDir/05_REDEMULTIPLEXED/$inputFile -db $outputDir/09_RELABEL_OTU_WELL/$inputFile2 -strand both -id 0.97 -uc $outputDir/09.1_/$inputFile.well_uc");
            my $results = system ("$usearch -usearch_global $outputDir/05_REDEMULTIPLEXED/$inputFile -db $outputDir/09_RELABEL_OTU_WELL/$inputFile2 -strand both -id 0.97 -uc $outputDir/09.1_/$inputFile.well_uc");
            if ($results != 0){
               die ("\nUparse map reads error");
            }
	 }
      }
   closedir ($dirOTU);
closedir ($dirFasta);

#Reading OTUwells to count
opendir (my $countOTUs, "$outputDir/09_RELABEL_OTU_WELL");
   while (my $inputFile = readdir ($countOTUs)){
      if ($inputFile =~ m/$relabelPattern/){
         idCount ("$outputDir/09_RELABEL_OTU_WELL/$inputFile", "3");
      }
   }
closedir ($countOTUs);

#COLLECTION OTUS
#Cat for .redemult files
print "\n\n***Cat redemultiplexed files";
my @concatenado = split ("/", $out);
my $name = pop @concatenado;
my $cat = "$name".".otuswell.cat";
opendir (my $dirRedemult, "$outputDir/09_RELABEL_OTU_WELL");
   my $results = system ("cat $outputDir/09_RELABEL_OTU_WELL/$relabelCatPattern > $outputDir/10_CAT_OTU_WELL/$cat");
   if ($results != 0){
      die "\n\nERROR: Cat .redemult error.\n";
   }
closedir ($dirRedemult);

#Runing collection dereplication
opendir (my $dirDerep, "$outputDir/10_CAT_OTU_WELL");
   while (my $inputFile = readdir ($dirDerep)){
      if ($inputFile =~ m/$otusWellPattern/){
         my $inputDirFile = "$outputDir/10_CAT_OTU_WELL/$inputFile";
         if (-s $inputDirFile){
            print "\n\n***Uparse COLLECTION dereplication";
            uparsePipelineDerep("$outputDir/10_CAT_OTU_WELL", $inputFile, "$outputDir/11_DEREPLICATION_CAT_OTU_WELL");
         }
      }
   }
closedir ($dirDerep);

#Runing collection Uparse sorting
opendir (my $dirSort, "$outputDir/11_DEREPLICATION_CAT_OTU_WELL");
   while (my $inputFile = readdir ($dirSort)){
      if ($inputFile =~ m/$derepCollPattern/){
         print "\n\n***Uparse COLLECTION sorting";
         uparsePipelineSort("$outputDir/11_DEREPLICATION_CAT_OTU_WELL", $inputFile, "$outputDir/12_SORTING_CAT_OTU_WELL");
      }
   }
closedir ($dirSort);

#Runing collection Uparse clustering
opendir (my $dirClust, "$outputDir/12_SORTING_CAT_OTU_WELL");
   while (my $inputFile = readdir ($dirClust)){
      if ($inputFile =~ m/$sortCollPattern/){
         print "\n\n***Uparse COLLECTION clustering";
         uparsePipelineClustering("$outputDir/12_SORTING_CAT_OTU_WELL", $inputFile, "$outputDir/13_CAT_OTU_WELL_CLUSTERING", $identity);
      }
   }
closedir ($dirClust);

#Runung relabel otu collection
print "\n\n***Uparse COLLECTION relabel";
relabelOTUcollection("$outputDir/13_CAT_OTU_WELL_CLUSTERING/$name.otuswell.cat.derep.sorted*", "$outputDir/14_OTU_COLLECTION_RELABEL/$name.otuswell.cat.derep.sorted.relabel");

#Runing read_map
print "\n\n***Uparse COLLECTION read map";
uparseMap("$outputDir/10_CAT_OTU_WELL/$name.otuswell.cat", "$outputDir/14_OTU_COLLECTION_RELABEL/$name.otuswell.cat.derep.sorted.relabel", "$outputDir/15_OTU_COLLECTION_MAP/$name.otuswell.cat.derep.sorted.map.uc");

#Runing otu_table
print "\n\n***Uparse COLLECTION otu_table";
uparseTable("$outputDir/15_OTU_COLLECTION_MAP/$name.otuswell.cat.derep.sorted.map.uc","$outputDir/16_OTU_COLLECTION_TABLE/$name.otuswell.cat.derep.sorted.map.uc.otu_table.txt");

#printar hash que contabiliza quantos arquivos sobram depois de demultiplex, chimera, e 16S filtering.
print "\n\nSHORT REPORT\n";
print "WELL\tOTUS\n";
my @barcodesArray = sort keys %barcodes;
foreach my $var (@barcodesArray){
   print "$var\t$barcodes{$var}[3]\n";
}
