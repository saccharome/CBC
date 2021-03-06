# Community-based culture collection (CBC) scripts

## _Background_
Microbial collections based exclusively on axenic cultures are time consuming, costly and may result in the loss of relevant biological information, since there are microorganisms that might depend on microbe-microbe interactions for their growth. As an alternative, we explore the concept of community-based culture collections (CBC) that is based on the isolation of colonies containing single or multiple microorganisms. By using this approach, we created a community-based culture collection with isolates from sugarcane plant organs. We developed a method to annotate the community-based isolates through a [multiplex strategy](https://www.nature.com/articles/srep29543) of pooling and sequencing near-full-length 16S rRNA gene amplicons. This method provides a high-throughput strategy for identifying isolated microbes in our sugarcane community-based culture collection and allows fast and accurate cross-referencing with the [sugarcane microbiome](https://www.nature.com/articles/srep28774).

# Documentation

## 1. Introduction
"CBC" is a program that combine four independent major scripts. Scripts were designed to optimize [demultiplexing](https://github.com/saccharome/CBC/blob/master/README.md#11-script-1-demultiplexing), [quality filtering](https://github.com/saccharome/CBC/blob/master/README.md#12-script-2-filtering-chimeras-and-non-16s-sequences) and [reliability](https://github.com/saccharome/CBC/blob/master/README.md#13-script-3-filtering-for-reliability) of sequences, and [clustering](https://github.com/saccharome/CBC/blob/master/README.md#14-script-4-two-step-clustering). The scripts were tested in a dataset of the __Sugarcane Community-Based Culture Collection (CBC)__ and described in this [paper](https://www.nature.com/articles/srep29543).

### 1.1. Script #1: Demultiplexing

* 1.1.1. This script aims to demultiplex large datasets of FASTA sequences of a multiplexing 16S rRNA gene amplicon sequencing based on two-step PCR amplifications with tagged primers for plates, rows, and columns, as illustrated below. See [Fig. 1](https://www.nature.com/articles/srep29543) for image details.

![figure 1](https://cloud.githubusercontent.com/assets/20187988/25688998/226b3fec-305a-11e7-9b47-246a3671b964.png)

* 1.1.2. This script uses circular consensus sequences (CCSs) in FASTA format as input, given by PacBio RS II sequencing, and creates a directory as output. This script is exclusively designed for demultiplexing sequences that contain the same amplicon structure as illustrated below. See [Supp. Fig. 1b](https://www.nature.com/articles/srep29543) for image details.

![supplementary figure 1](https://cloud.githubusercontent.com/assets/20187988/25688996/223e299e-305a-11e7-9f6c-cd12f6933572.png)

* 1.1.3. Files containing sequences of _barcodes of plates_, _columns_ and _rows_, _forward_ and _reverse Nextera_ sequences, and _reverse primer_ sequence must be also inputted in FASTA format.

* 1.1.4. The script first aligns _forward_ and _reverse Nextera_ sequences and _reverse primer_ sequence to identify putative positions of neighboring barcodes. The expected position of each barcode can be expanded by considering more nucleotides to right and left side of barcodes, using the parameter _margin of error_.

* 1.1.5. Barcode-trimmed sequences are outputted in plus-strand based on the reverse primer alignment.

* 1.1.6. Demultiplexed sequences header are modified to include the label of demultiplexing. Headers receive _;barcodelabel=_ followed by a given label. The lack of at least one of the three barcodes leads a sequence to be considered undemultiplexed.

* 1.1.7. Parameters must be inputted in the following order:

  * __input_fa__ _Input file containing CCS in FASTA format_
  * __plate_bc.fa__ _FASTA containing barcodes of plates_
  * __row_bc.fa__ _FASTA containing barcodes of rows_
  * __column_bc.fa__ _FASTA containing barcodes of columns_
  * __fwd_nex.fa__ _FASTA containing Nextera forward sequence_
  * __rev_nex.fa__ _FASTA containing Nextera reverse sequence_
  * __rev_primer.fa__ _FASTA containing reverse primer sequence_
  * __error_margin__ _Number of margin of error for barcode prediction; e.g. 3 nt_
  * __max_length__ _Maximum length accepted of sequences; e.g. 1600 nt_
  * __output_dir__ _Directory of outputted files_

* 1.1.8. Example of command line:

```bash
perl CBC_script1.pl ccs_example.fa plate_bc.fa rows_bc.fa column_bc.fa fwd_nex.fa rev_nex.fa rev_primer.fa 3 1600 output_directory
```

* 1.1.9. Output files:

	* __concatenated.fasta__ _FASTA concatenating sequences if broken in lines_
	* __FASTAFILTERED.output__ _FASTA containing CCSs of expected size_
	* __outoflength.output__ _FASTA containing out-of-size CCSs_
	* __PLA.output__ _Tab-separated output of alignment of plates barcodes_
	* __ROW.output__ _Tab-separated output of alignment of rows barcodes_
	* __COL.output__ _Tab-separated output of alignment of columns barcodes_
	* __FNEX.output__ _Tab-separated output of alignment of fwd_nex sequence_
	* __RNEX.output__ _Tab-separated output of alignment of rev_nex sequence_
	* __RPRIM.output__ _Tab-separated output of alignment of rev_primer sequence_
	* __.demult__ _Files containing demultiplexed CCSs_
	* __undemultiplexed.undemult__ _IDs of undemultiplex CCSs_
	* __report.report__ _Numbers obtained and summary of CBC_script1.pl run_

### 1.2. Script #2: Filtering chimeras and non-16S sequences

* 1.2.1. This script aims to filter demultiplexed 16S rRNA gene sequences from chimeras and non-16S sequences based on chimera-free and a large 16S rRNA databases, respectively. The script is based on [UCHIME](https://doi.org/10.1093/bioinformatics/btr381) algorithm.

* 1.2.2. This script uses demultiplexed sequences, outputted from _CBC_script1.pl_.

* 1.2.3. In order to optimize processing, a _cat_ command can be used to join all barcode-labeled sequences in a single file. Files extension must be _.demult_.

* 1.2.4. After filtering and chimera/non-16S sequences discarded, sequences are redemultiplexed.

* 1.2.5. Parameters:

  * __-i__ _Input directory containing .demult file_
  * __-o__ _Output root directory_
  * __-ch__	_Chimera-free database_
  * __-db__	_16S rRNA database_

* 1.2.6. Example of command line:

```bash
perl CBC_script2.pl -i demultiplexed.demult -o output_directory -ch chimera_free.fa -db 16S_db.fa
```

* 1.2.7. Output folders and file:
  * __01_CAT_DEMULTIPLEXED__ _Folder of concatenated .demult files_
  * __02_CHIMERA_FILTERING__ _Folder of files containing chimeras, nonchimeras and alignment output of chimera filtering_
  * __03_ALIGNMENT_FOR_16S_FILTER__ _Folder of tab-separated output of alignment with 16S database_
  * __04_16S_FILTERED__ _Folder of FASTA sequences filtered from chimera and non-16S sequences_
  * __report.report__ _Number of sequences in each step for each well: demultiplexed, nonchimeras and 16S_

* 1.2.8. FASTA sequences are broken in lines and must be concatenated before inputted in _CBC_script3.pl_.

### 1.3. Script #3: Filtering for reliability

* 1.3.1. This script aims to filter sequences based on their reliability. A 16S rRNA gene sequence is considered reliable if it is above a threshold of similarity to a sequence in a 16S database and/or to any other sequence within the dataset.

* 1.3.2. This script uses a FASTA file as input and creates a directory as output.

* 1.3.3. Alignments were performed using _usearch_global_ algorithm from [USEARCH](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btq461) package considering threshold of 0.97.

* 1.3.4. Sequences maintained after filtering are outputted in a _.filtered_ file.

* 1.3.5. Parameters must be inputted in the following order:

  * __dataset__ _Dataset containing sequences that must be filtered_
  * __16S_database__ _16S rRNA gene database_
  * __output_dir__ _Output directory_
  
* 1.3.6. Example of command line:

```bash
perl CBC_script3.pl seqs.fa 16S_db.fa output_directory
```

* 1.3.7. Output files:

  * __.db_m8__ _Tab-separated file of alignment of dataset with 16S database_
  * __.db_aln__ _Alignment file of dataset with 16S database_
  * __.db_hit.fasta__ _FASTA of sequences with hit against 16S database_
  * __.db_nohit.fasta__ _FASTA of sequences without hit against 16S database_
  * __.db_nohit.fasta.dataset_m8__ _Tab-separated file of alignment of dataset with itself_
  * __.db_nohit.fasta.dataset_aln__ _Alignment file of dataset with itself_
  * __.db_nohit.fasta.dataset_hit.id__ _ID of sequences with hit against dataset_
  * __.db_nohit.fasta.dataset_nohit.fasta__ _FASTA of sequences without hit against dataset_
  * __.filtered__ _FASTA of reliable sequences_

* 1.3.8. FASTA sequences must be redemultiplexed in each well before inputted in _CBC_script4.pl_.

### 1.4. Script #4: Two-step clustering

* 1.4.1. This script aims to cluster sequences within wells to obtain OTUs (operational taxonomic units) named well-OTUs (wOTUs) in each single well of 96-well plates. Since might be a redundancy of microbes among wells, a second clustering is performed to obtain unique collection-OTUs (cOTUs). The script is based on [UPARSE](http://www.nature.com/nmeth/journal/v10/n10/full/nmeth.2604.html) algorithm.

* 1.4.2. Files extension must be _.demult_ in inputted folder.

* 1.4.3. Parameters:

  * __-i__ _Input directory containing “.redemult” files_
  * __-o__ _Output directory_
  * __-id__ _Identity of clustering (e.g. 0.97)_

* 1.4.4. Example of command line:

```bash
perl CBC_script4.pl -i input_directory -o output_root_directory -id uparse_ID
```

* 1.4.5. Output folders:

  * __06_PER_WELL_DEREPLICATION__ _Dereplication of sequences per well .redemult_
  * __07_PER_WELL_SORTING__ _Sorting of sequences per well .redemult_
  * __08_PER_WELL_CLUSTERING__ _Clustering of sequences within wells (wOTUs)_
  * __09.1_MAPPING_READS__ _Mapping reads of clustering on above-mentioned step_
  * __09_RELABEL_OTU_WELL__ _Relabeling wOTU names_
  * __10_CAT_OTU_WELL__ _Concatenating wOTUs in a single file_
  * __11_DEREPLICATION_CAT_OTU_WELL__ _Dereplication of sequences of wOTUs_
  * __12_SORTING_CAT_OTU_WELL__ _Sorting of sequences of wOTUs_
  * __13_CAT_OTU_WELL_CLUSTERING__ _Clustering of sequences wOTUs_
  * __14_OTU_COLLECTION_RELABEL__ _Relabeling unique OTUs (cOTUs)_
  * __15_OTU_COLLECTION_MAP__ _Mapping reads of clustering on above-mentioned step_
  * __16_OTU_COLLECTION_TABLE__ _OTU_table construction_

## 2. Pipeline overview

Raw reads are demultiplexed and filtered for reliability and quality of CCSs. The first clustering step allows determination of which and how many microorganisms are present in each well, while the second one provides information on the redundancy of OTUs. See [Supp. Fig. 2](https://www.nature.com/articles/srep29543) for image details.

![supplementary figure 2](https://cloud.githubusercontent.com/assets/20187988/25689108/30b7808c-305b-11e7-88d4-5208b6981d05.png)

## 3. Installation and dependencies

The scripts were written in Perl language and run in Linux environment via command line.

### 3.1. Softwares

#### 3.1.1 Linux

* 3.1.1.1. We used the following version of Linux:

```bash
Linux version 2.6.32-642.el6.x86_64 (mockbuild@x86-033.build.eng.bos.redhat.com)
(gcc version 4.4.7 20120313 (Red Hat 4.4.7-17) (GCC))
```

#### 3.1.2. USEARCH package

* 3.1.2.1. We used the following version of [USEARCH](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btq461):

```bash
usearch v8.0.1517_i86linux64
```

* 3.1.2.2. An environmental variable must be created to USEARCH in PATH, e.g. _$uv9_. Change the variable for USEARCH inside each script as required.

### 3.2. Databases

The following databases must be downloaded before proceed scripts run:

#### 3.2.1. Greengenes

* 3.2.1.1 Greengenes 16S rRNA gene database was downloaded from the [webpage](http://greengenes.secondgenome.com/downloads).

* 3.2.1.2. A release must be chosen, then chose _*.fasta.gz_ file.

* 3.2.1.3. In our work, we used the version “The Greengenes Database files from May, 2013”. Any other release should work, we have not tested though.

#### 3.2.2. UCHIME chimera-free 16S rRNA gene database

* 3.2.2.1. [UCHIME](https://doi.org/10.1093/bioinformatics/btr381) is an algorithm installed with USEARCH and requires a _gold_ database of chimera-free 16S rRNA gene sequences. The [database](http://drive5.com/uchime/uchime_download.html) was downloaded directly from [drive5.com](http://www.drive5.com/) website.

## 4. Example

We uploaded a small dataset of CCSs as an [example](https://github.com/saccharome/CBC/tree/master/example/ccs_dataset_example), as well as the files containing the [expected reports](https://github.com/saccharome/CBC/tree/master/example/expected_reports) of each CBC_script using this file.

## 5. Bugs

Please submit any problems, suggestions or requests [here](https://github.com/saccharome/CBC/issues).

## 6. Citation
Armanhi, J. S. L.; de Souza, R. S. C.; de Araújo, L. M.; Okura, V. K.; Mieczkowski, P.; Imperial, J.; Arruda, P. Multiplex amplicon sequencing for microbe identification in community-based culture collections. Sci. Rep. 6, 29543, 2016. [doi:10.1038/srep29543](http://www.nature.com/articles/srep29543/).
