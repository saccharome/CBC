#! usr/local/bin/perl

#THIS SCRIPT IS RELATED TO THE FOLLOWING STUDY PUBLISHED ON SCIENTIFIC REPORTS:

#Armanhi, J. S. L.; de Souza, R. S. C.; de Araujo, L. M.; Okura, V. K.; Mieczkowski, P.; Imperial, J.; Arruda, P. Multiplex amplicon sequencing for microbe identification in community-based culture collections. Sci. Rep. 6, 29543; doi:10.1038/srep29543 (2016) (www.nature.com/articles/srep29543/).

#If you use this script or part of it please cite us.

use strict;

my $usearch = '$uv8';

if (!(scalar(@ARGV) == 10)){
   die ("\n\nUSAGE: perl script.pl <input_fa> <plate_bc.fa> <row_bc.fa> <column_bc.fa> <fwd_nex.fa> <rev_nex.fa> <rev_primer.fa> <error_margin> <max_length> <output_dir>\n\n
(0) This script has a .FASTA FILE as input and creates a directory as OUTPUT, given by ARGV parameter.
(1) This script first identifies the length of Illumina adapters and reverse primer.
(2) Filter sequences from .fasta file based on maximum length accepted for each sequence. Sequences out of this filter are discarded.
(3) Sequences filtered are aligned (ublast) against nextera transposase forward, reverse and reverse primer.
(4) Based on its alignments it is estimated the position of the barcodes to be identified.
(5) Barcodes are aligned with no mismatch allowed (search_oligodb) and only those with correct predicted positions are accepted.
(6) A margin of error can be given (it is the number of nucleotides permited right and left based on predicted positions).
(7) Sequences are then barcode trimmed and outputted in plus-strand based on the reverse primer alignment.
(8) Sequences that lack at least one of the three barcodes are considered undemultiplexed. Sequences with barcores identified but not in the expected order.
(9) A report is criated to organized the generated informations during the process.\n\n\n");
}

#(0).CREATING folders
my @folders = ($ARGV[9]);

sub makeFolder {
   my @params = @_;
   my @f = @{$params[0]};
   mkdir "$params[1]";
   foreach my $var (@f){
      mkdir "$params[1]/$var";
   }
}
makeFolder(\@folders, $ARGV[9]);

#(1)Recognition of Illumina adaptors length step
open (FLENGTH, @ARGV[4]);
   my $lengthF;
   while (<FLENGTH>){
      chomp;
      my $linhaAtualFL = $_;
      if (!($linhaAtualFL =~ />/)){
         $lengthF = length $linhaAtualFL;
      }
   }
close (FLENGTH, @ARGV[4]);

open (RLENGTH, @ARGV[5]);
   my $lengthR;
   while (<RLENGTH>){
      chomp;
      my $linhaAtualRL = $_;
      if (!($linhaAtualRL =~ />/)){
         $lengthR = length $linhaAtualRL;
      }
   }
close (RLENGTH, @ARGV[5]);

open (PRLENGTH, @ARGV[6]);
   my $lengthP;
   while (<PRLENGTH>){
      chomp;
      my $linhaAtualPL = $_;
      if (!($linhaAtualPL =~ />/)){
         $lengthP = length $linhaAtualPL;
      }
   }
close (PRLENGTH, @ARGV[6]);

#Concatenating file step if needed
open (FASTA, @ARGV[0]);
my $escrevi = 0;
while (<FASTA>){
   chomp;
   my $linhaAtual = $_;
   open (my $file3, ">>", "$ARGV[9]/concatenated.fasta");
      if ($linhaAtual =~ />/){
         if ($escrevi == 1){
            print {$file3} "\n";
         }
         print {$file3} "$linhaAtual"."\n";
         $escrevi = 1;
      }else{
         print {$file3} "$linhaAtual";
      }
   close ($file3);
}
close (FASTA, @ARGV[0]);

#(1.1)Size filtering of CCSs (1 to 1600 bp)
my $fastaEntrada = 0;
my $fastaSaida = 0;
open (FASTA, "@ARGV[9]/concatenated.fasta");
   my $linhaAnterior;
   my $filtro_menor = 1;
   my $filtro_maior = $ARGV[8];
   open (my $fastafiltred, ">>", "$ARGV[9]/FASTAFILTERED.output");
      open (my $outoflength, ">>", "$ARGV[9]/outoflength.output");
         while (<FASTA>){
            chomp;
            my $linhaAtual = $_;
            if ($linhaAtual =~ />/){
               $linhaAnterior = $linhaAtual;
               $fastaEntrada = $fastaEntrada+1;
            }else{
               if ((length $linhaAtual >= $filtro_menor) and (length $linhaAtual <= $filtro_maior)){
                  print {$fastafiltred} "$linhaAnterior\n$linhaAtual\n";
                  $fastaSaida = $fastaSaida+1;
               }else{
                  print {$outoflength} "$linhaAnterior\n$linhaAtual\n";
               }
            }
         }
         my $outoffilter = $fastaEntrada - $fastaSaida;
      close ("$ARGV[9]/outoflength.output");
   close ("$ARGV[9]/FASTAFILTERED.output");
close (FASTA, "@ARGV[9]/concatenated.fasta");

#(2)Step of alignment of .fasta with adaptors/primers
my $result = system ("$usearch -ublast $ARGV[9]/FASTAFILTERED.output -db @ARGV[4] -evalue 10 -strand both -blast6out $ARGV[9]/FNEXT.output");
if ($result == 1){
   die ("wrong adaptadorF.fasta file");
}

my $result = system ("$usearch -ublast $ARGV[9]/FASTAFILTERED.output -db @ARGV[5] -evalue 10 -strand both -blast6out $ARGV[9]/RNEXT.output");
if ($result == 1){
   die ("wrong adaptadorR.fasta file");
}

my $result = system ("$usearch -ublast $ARGV[9]/FASTAFILTERED.output -db @ARGV[6] -evalue 10 -strand both -blast6out $ARGV[9]/RPRIM.output");
if ($result == 1){
   die ("wrong primerreve.fasta file");
}

#(3)Step for marking position of adaptors and estimating position of barcodes
my $mismatches = 3;
my %posicoesbc;

open (FNEXTADAPTRED, "$ARGV[9]/FNEXT.output");
   my $bcRstart;
   my $bcRend;
   my $start16S;
   my $end16S;
   while (<FNEXTADAPTRED>){
      chomp;
      my $linhaAtualFA = $_;
      my @FA = split ("\t",$linhaAtualFA);
      if (($lengthF-$FA[3]+$FA[4]) <= $mismatches){
         if (!((defined $posicoesbc{$FA[0]}[2]) and (defined $posicoesbc{$FA[0]}[3]))){
            if (@FA[8]<$FA[9]){
               $bcRstart = ($FA[6]-($FA[8]-1))-8-@ARGV[7];
               $bcRend =   ($FA[6]-($FA[8]-1))-1+@ARGV[7];
               $posicoesbc{$FA[0]}[2]=$bcRstart;
               $posicoesbc{$FA[0]}[3]=$bcRend;
               $start16S = ($FA[7]+(14-$FA[9]))+19+9+2+1;
               $posicoesbc{$FA[0]}[6]=$start16S;
            }else{
               $bcRstart = ($FA[7]+($FA[9]-1))+1-@ARGV[7];
               $bcRend =   ($FA[7]+($FA[9]-1))+8+@ARGV[7];
               $posicoesbc{$FA[0]}[2]=$bcRstart;
               $posicoesbc{$FA[0]}[3]=$bcRend;
               $end16S = ($FA[6]-(14-$FA[8]))-19-9-2-1;
               $posicoesbc{$FA[0]}[7]=$end16S;
            }
         }
      }
   }
close (FNEXTADAPTRED, "$ARGV[9]/FNEXT.output");

open (RNEXTADAPTRED, "$ARGV[9]/RNEXT.output");
   my $bcCstart;
   my $bcCend;
   while (<RNEXTADAPTRED>){
      chomp;
      my $linhaAtualRA = $_;
      my @RA = split ("\t",$linhaAtualRA);
      if ($lengthR-$RA[3]+$RA[4] <= $mismatches){
         if (!((defined $posicoesbc{$RA[0]}[4]) and (defined $posicoesbc{$RA[0]}[5]))){
            if ($RA[8]>$RA[9]){
               $bcCstart = ($RA[7]+($RA[9]-1))+1-@ARGV[7];
               $bcCend =   ($RA[7]+($RA[9]-1))+8+@ARGV[7];
               $posicoesbc{$RA[0]}[4]=$bcCstart;
               $posicoesbc{$RA[0]}[5]=$bcCend;
            }
            if ($RA[8]<$RA[9]){
               $bcCstart = ($RA[6]-($RA[8]-1))-8-@ARGV[7];
               $bcCend =   ($RA[6]-($RA[8]-1))-1+@ARGV[7];
               $posicoesbc{$RA[0]}[4]=$bcCstart;
               $posicoesbc{$RA[0]}[5]=$bcCend;
            }
         }
      }
   }
close (RNEXTADAPTRED, "$ARGV[9]/RNEXT.output");

open (PRIMERREVERSO, "$ARGV[9]/RPRIM.output");
   my $bcPstart;
   my $bcPend;
   my $end16S;
   my $start16S;
   while (<PRIMERREVERSO>){
      chomp;
      my $linhaAtualPR = $_;
      my @PR = split ("\t",$linhaAtualPR);
      if ($lengthP-@PR[3]+@PR[4] <= $mismatches){
         if (!((defined $posicoesbc{$PR[0]}[0]) and (defined $posicoesbc{$PR[0]}[1]))){
            if ($PR[8]>$PR[9]){
               $bcPstart = ($PR[7]+($PR[9]-1))+2+1-@ARGV[7];
               $bcPend =   ($PR[7]+($PR[9]-1))+2+9+@ARGV[7];
               $posicoesbc{$PR[0]}[0]=$bcPstart;
               $posicoesbc{$PR[0]}[1]=$bcPend;
               $end16S = ($PR[7]+($PR[9]-1));
               $posicoesbc{$PR[0]}[7]=$end16S;
            }
            if ($PR[8]<$PR[9]){
               $bcPstart = ($PR[6]-($PR[8]-1))-2-9-@ARGV[7];
               $bcPend =   ($PR[6]-($PR[8]-1))-2-1+@ARGV[7];
               $posicoesbc{$PR[0]}[0]=$bcPstart;
               $posicoesbc{$PR[0]}[1]=$bcPend;
               $start16S = ($PR[6]-($PR[8]-1));
               $posicoesbc{$PR[0]}[6]=$start16S;
            }
         }
      }
   }
close (PRIMERREVERSO, "$ARGV[9]/RPRIM.output");

#(5)Step for search_olibodb (USEARCH) run for barcodes
my $result = system ("$usearch -search_oligodb $ARGV[9]/FASTAFILTERED.output -db @ARGV[1] -strand both -maxdiffs 0 -blast6out $ARGV[9]/PLA.output");
if ($result == 1){
   die ("wrong plate.fasta file");
}

my $result = system ("$usearch -search_oligodb $ARGV[9]/FASTAFILTERED.output -db @ARGV[2] -strand both -maxdiffs 0 -blast6out $ARGV[9]/ROW.output");
if ($result == 1){
   die ("wrong row.fasta file");
}

my $result = system ("$usearch -search_oligodb $ARGV[9]/FASTAFILTERED.output -db @ARGV[3] -strand both -maxdiffs 0 -blast6out $ARGV[9]/COL.output");
if ($result == 1){
   die ("wrong column.fasta file");
}

#(6)Reading files PLA.OUTPUT, ROW. and COL.
my %demultiplex;
open (PLA_OUTPUT, "$ARGV[9]/PLA.output");
   while (<PLA_OUTPUT>){
      my $linhaAtualP = $_;
      my @PLA = split ("\t", $linhaAtualP);
      if ((($PLA[6] >= $posicoesbc{$PLA[0]}[0]) and ($PLA[7] <= $posicoesbc{$PLA[0]}[1])) and (!(defined $demultiplex{$PLA[0]}[0]))){
         $demultiplex{$PLA[0]}[0]=$PLA[1];
      }
   }
close (PLA_OUTPUT, "$ARGV[9]/PLA.output");

open (ROW_OUTPUT, "$ARGV[9]/ROW.output");
   while (<ROW_OUTPUT>){
      my $linhaAtualR = $_;
      my @ROW = split ("\t", $linhaAtualR);
      if ((($ROW[6] >= $posicoesbc{$ROW[0]}[2]) and ($ROW[7] <= $posicoesbc{$ROW[0]}[3])) and (!(defined $demultiplex{$ROW[0]}[1]))){
         $demultiplex{$ROW[0]}[1]=$ROW[1];
      }
   }
close (ROW_OUTPUT, "$ARGV[9]/ROW.output");

open (COL_OUTPUT, "$ARGV[9]/COL.output");
   while (<COL_OUTPUT>){
      my $linhaAtualC = $_;
      my @COL = split ("\t", $linhaAtualC);
      if ((($COL[6] >= $posicoesbc{$COL[0]}[4]) and ($COL[7] <= $posicoesbc{$COL[0]}[5])) and (!(defined $demultiplex{$COL[0]}[2]))){
         $demultiplex{$COL[0]}[2]=$COL[1];
      }
   }
close (COL_OUTPUT, "$ARGV[9]/COL.output");

#(7)Printing nt of begining and endind of each barcode
my @posicoes = keys %posicoesbc;
my $cont = 0;
while ($cont < scalar @posicoes){
   if ((($posicoesbc{$posicoes[$cont]}[2] < $posicoesbc{$posicoes[$cont]}[0]) and ($posicoesbc{$posicoes[$cont]}[0] < $posicoesbc{$posicoes[$cont]}[4]) and ($posicoesbc{$posicoes[$cont]}[2] < $posicoesbc{$posicoes[$cont]}[4]))){
   $posicoesbc{$posicoes[$cont]}[8] = "+";
   }
   if (($posicoesbc{$posicoes[$cont]}[4] < $posicoesbc{$posicoes[$cont]}[0]) and ($posicoesbc{$posicoes[$cont]}[0] < $posicoesbc{$posicoes[$cont]}[2]) and ($posicoesbc{$posicoes[$cont]}[4] < $posicoesbc{$posicoes[$cont]}[2])){
   $posicoesbc{$posicoes[$cont]}[8] = "-";
   }
   $cont++;
}

#(8)Printing seq of wells
open (FASTA, "$ARGV[9]/FASTAFILTERED.output");
   my %hashtemp;
   my $escrevi;
   my $ondeEscrevi;
   my @soID;
   my $file;
   my $semBarcode = 0;
   my $wellDemultiplexado = 0;
   my $CCSdemultiplexado = 0;
   my $quimera = 0;
   while (<FASTA>){
      chomp;
      my $linhaAtual = $_;
      open (my $erro, ">>", "$ARGV[9]/undemultiplexed.undemult");
      if ($linhaAtual =~ />/){
         @soID = split (">", $linhaAtual);
         if ((!(defined $demultiplex{$soID[1]}[0])) or (!(defined $demultiplex{$soID[1]}[1])) or (!(defined $demultiplex{$soID[1]}[2]))){
            print {$erro} ">"."$soID[1]\t$demultiplex{$soID[1]}[0]\t$demultiplex{$soID[1]}[1]\t$demultiplex{$soID[1]}[2]\n";
            $semBarcode = $semBarcode+1;
         }else{
            if ((($posicoesbc{$soID[1]}[2] < $posicoesbc{$soID[1]}[0]) and ($posicoesbc{$soID[1]}[0] < $posicoesbc{$soID[1]}[4]) and ($posicoesbc{$soID[1]}[2] < $posicoesbc{$soID[1]}[4])) or (($posicoesbc{$soID[1]}[4] < $posicoesbc{$soID[1]}[0]) and ($posicoesbc{$soID[1]}[0] < $posicoesbc{$soID[1]}[2]) and ($posicoesbc{$soID[1]}[4] < $posicoesbc{$soID[1]}[2]))){
            $file = "$demultiplex{$soID[1]}[0]"."_"."$demultiplex{$soID[1]}[1]"."$demultiplex{$soID[1]}[2]".".demult";
               if (!(exists $hashtemp{$file})){
                  open (my $arq, ">", "$ARGV[9]/$file");
                  $wellDemultiplexado = $wellDemultiplexado+1;
                  $hashtemp{$file} = $arq;
               }
            print {$hashtemp{$file}} "$linhaAtual".";barcodelabel="."$demultiplex{$soID[1]}[0]"."_"."$demultiplex{$soID[1]}[1]"."$demultiplex{$soID[1]}[2]"."\n";
            $CCSdemultiplexado = $CCSdemultiplexado+1;
            $escrevi = 1;
            $ondeEscrevi = $file;
            }else{
               print {$erro} ">"."$soID[1]\t"."barcodes desordenados: possível quimera\n";
               $quimera = $quimera+1;
            }
         }
      }elsif ($escrevi == 1){
         my $corteInicial = $posicoesbc{$soID[1]}[6]-1;
         my $tamanho16S   = $posicoesbc{$soID[1]}[7]-$corteInicial;
         my $substr16S = substr ($linhaAtual, $corteInicial, $tamanho16S);
         if ($posicoesbc{$soID[1]}[8] eq "+"){
            print {$hashtemp{$ondeEscrevi}} "$substr16S\n";
            $escrevi = 0;
         }else{
            $substr16S =~ tr/ATCGatcg/TAGCtagc/;
            print {$hashtemp{$ondeEscrevi}} (scalar reverse $substr16S);
            print {$hashtemp{$ondeEscrevi}}"\n";
            $escrevi = 0;
         }
      }
      close ($erro, "$ARGV[9]/undemultiplexed.undemult");
   }
close (FASTA, "/$ARGV[9]/FASTAFILTERED.output");

open (my $report, ">>", "$ARGV[9]/report.report");
print {$report} "Input file____________________\t@ARGV[0]
Output dir____________________\t@ARGV[9]
Number of CCS sequences_______\t$fastaEntrada
Length filter_________________\t$filtro_menor to $filtro_maior
Filtered sequences____________\t$fastaSaida
Out of length filter__________\t$outoffilter
Adaptor mismatches accepted___\t$mismatches
Barcode margin of error_______\t@ARGV[7]
Maximum length of sequence____\t@ARGV[8]
Undemult. CCS (>=1 unkn. bc)__\t$semBarcode
Demultiplexed wells___________\t$wellDemultiplexado
Demultiplexed CCS_____________\t$CCSdemultiplexado
Putative chimera (bc desord.)_\t$quimera";
close ("report.report");
