# PURC: Pipeline for Untangling Reticulate Complexes
# v2.0
#### Last updated 2021 June 22

## Major changes from v1 ##
* PURC updated for compatibility with Python3 and macOS 10+
* Usearch replaced with [vsearch](https://github.com/torognes/vsearch) for Linux and macOS
* Amplicon sequence variant (ASVs) identification using [DADA2](https://benjjneb.github.io/dada2/index.html) introduced
* Linux only: PacBio's [lima](https://github.com/pacificbiosciences/barcoding/) replaces BLAST methods for demultiplexing
* Input file path handling. Files no longer have to be in the working directory (output will still be written to working directory)
* Input sequence file can be gzip compressed
* Primer order in config now must match locus order
* Multiple reference sequences for a locus must be in the same orientation
* Dependencies expected to be in PATH by default

## **Overview** ##
PURC is a pipeline for inferring the underlying biological sequences (alleles, paralogs, or homeologs) from amplicon sequencing data (PacBio, Illumina, etc), de-multiplexing them (labeling each sequence with its locus and source sample), and cleaning them (removing PCR errors, sequencing errors, and chimeras). It is geared toward analyzing polyploid species complexes but is also effective for other applications; the final output of a full run includes an alignment for each locus with each homeolog or allele sequence in the amplicon data labeled with the source sample information and amount of coverage.

PURC allows users to extract and analyze all the copies of a given locus present in an amplicon pool, so is particularly useful in cases such as the study of polyploid complexes or large gene families, which were historically tractable only via time- and expense-intensive cloning approaches. Within a single run users can analyze as many accessions as they like, limited only by the desired coverage/allele and the number of sequences that can be uniquely linked back to their source samples. PURC can perform this linking using barcode sequences, locus identity, or phylogenetic information. For example, within a single run an individual barcode can be used multiple times if the locus is different, or if the accessions can be distinguished phylogenetically (i.e., of different genera or other clades that can be distinguished using BLAST). While it is most useful in cases where multiple copies of fairly long sequences have been amplified from individual accessions, it also (when used on PacBio sequencing data) provides considerable cost- and time savings for sequencing plastid or mitochondrial markers.


### Rough outline of PURC's workflow: ###

* Optional: scan for concatemers and split them into their component sequences (more on concatemers [here](https://github.com/PacificBiosciences/cDNA_primer/wiki/Artificial-concatemers,-PCR-chimeras,-and-fusion-genes))
* Identify and remove barcode sequences
* Trim primers and other adapter sequences
* Assign each read to its source accession based on the barcode and a user-specified list of reference sequences
* Cluster sequences and remove chimeric sequences, iteratively
* Compute consensus sequence for each cluster
* Produce sequence alignments ready for downstream phylogenetic analyses

PURC should run on most recent versions of macOS and Linux. Windows support is deprecated, instead on PCs we recommend running an Ubuntu virtual machine (see tutorial [here](https://itsfoss.com/install-linux-in-virtualbox/)).

## **Quick Start** ##
### Step 0: Clone this repository ###
To download the program, enter your terminal, move to the directory where you want to download it, and clone from the repo. For example:
```bash
cd ~
git clone https://peter_schafran@bitbucket.org/peter_schafran/purc.git && cd purc
```
### Step 1: Setup ###
PURC consists of purc.py (and another variation--purc_recluster.py--that we describe below) and relies on a number of dependencies. We recommend using the [Miniconda](https://conda.io/en/latest/miniconda.html) package manager for installing dependencies. Once installed (and the terminal rebooted), you should be able to run on of these commands to install dependencies, depending on your operating system:
```bash
# macOS
conda env create -n purc --file purc_macos.yaml && conda activate purc

# Linux
conda env create -n purc --file purc_linux.yaml && conda activate purc
```
If PURC fails to run it is likely due to dependency issues. See advanced installation at the bottom of this page. In addition to the main log file created in the working directory, R logs are created in each locus/sample directory in the output folder.

### Step 2: Get files ready ###
PURC requires the following files:

* **Configuration file** - This contains the file names and parameters that PURC needs. There are several examples provided with the PURC distribution (open and edit them in TextWrangler or similar text editor).

* **Barcode sequence file** - In fasta format. PURC uses this file to identify barcodes in the sequence reads. The sequence name is the barcode name. E.g.:

        >BC01
        ACTACATATGAGATG
        >BC02
        TCATGAGTCGACACTA
        >BC03
        TATCTATCGTATACGC


	Note that if you are using the dual barcode function, the names have to be BCF1, BCF2, ..., BCF9, and BCR1, BCR2, ..., BCR9, where         the F and R denote the barcodes on forward and reverse primers respectively. This is used to check whether there are invalid F-F or R-R situation; only F-R or R-F are kept.

* **Reference sequence file** - In fasta format. PURC uses these references to assign reads to loci (and optionally to phylogenetic group). The sequences themselves must be free of gaps ("-") or else Blast will choke when it tries to make the Blast database. Each reference seq name must specify the locus ('locus=') that this sequence represents; you can optionally note where the sequence came from ('ref_taxon='). The ref_taxon information can contain anything (as long as no spaces are special characters are included) and is not used by PURC--it is only useful in allowing users to keep track of which sequences are included in their reference file more easily. Each designation is separated by a backslash ("/"). For example:

        >locus=ApP/ref_taxon=Cystopteris_bulbifera
        TGCCACACTGGTGAGTATTATCTTACTCTTTTTCTTGGTGAGAAAGGGTAGTGTGCATGGCATATTCATGTCAACCCCCGCTTGGGACCGGGGGCGACAAGGATGTTACTGTTGGGTGATACCTGTGATGCCCAATTGGAGCAAGAGTAAAATCAACTTTGTAAATATCATCTATTTGAAGGATTAACAGACATGGTATTTAAATTCCTCTCACATTTAAAACAGGGTGGTTGCAGAACTGGTATGGCCAAAGTAACGAACGCTTACGATTTGCCTGCAAGGTAAAAGTTGAACAATGCTCAAGGTGGGGCTAGTTCTTTTGTCACTTAAGCAAGGATCTTCAAGCATGTAAAATTATTCTCCCTCAACTTTGCTTTACAAAAGAAATTTAACATATTGACTACTTTATCAATGGAATTTGAGCAGCTATCACATGTCGATGTTTATTTTGAGGCAAGGGGTTCTTGTTTGCATGCGGTTGTAGAAATGTCTTATCACATTTCTATTTGGGTTTTTTGACATTGCTATTTTTGCATAAATGCTAAGTTACAAATTAGAATTGTTTACTTGTTTGTTTGGAGGAAATCTCAAACGACTGTCTTTTGCTCTTGTATGCTCATTTGATGATTGCATGCGTACACCTTTATGTTCATTTAGGCTATGTTTTGTCAGCTCACAAATTTTTGATGTTAAACCTAACATGACAGGAAAGTTATACATACTGTTGGTCCAAGATATGCTGTTAAATATCATACAGCTGCAGAAAATGCTCTAAGTCATTGTTACCGATCTTGTTTAGAGGCTTTGATTGACTTAGGCCTTCAAAGGTACCAGCTGCTTGTTTAAACAGCTCAAAATTAAAGGAGAGTTTTTTCCTTTTGGATTAAAGTTTATCTCCCTTGTAATTCTTGCAGCATTGCCCTGGGGTGTATTTACACAGAGTCTAAAGGCTAT
        >locus=GAP/ref_taxon=Cystopteris_bulbifera
        AAAGTAATTTGATGCGTTTTTTTGGTGAAGTCAGTGCATTTTATTCTTTGAAATGGTGAGATCTTATGTTTGCAATCCACTGTTTACCAGGTAGTCCACCAGGAGTTTGGAATTGCCGAGGCCCTCATGACAACTGTGCATGCTACCACTGGTACCGTATTGTCTCCAGCCTCTTTGATTCTGATGTGCTGTAGTTTATACAGAGGCTGGAATTGAATACCCTTGTCATGACTAGCTCTTATATGTCGGTATTTGTGCTTCCAGCTACTCAAAAGACTGTGGATGGTCCCTCTGGCAAAGATTGGCGGGGTGGCCGAGGTGCTGGGCAGAATATCATTCCGAGCTCAACCGGTGCTGCCAAGGTACACTGCTAACTTAACCAAAACCTAAAACATCCTACTAATGTGGTGGCTATTTTATTGCCTTGTAAACATGTTGGTCTTTGTGGTAAACTGCAGGCAGTAGGAAAGGTGCTTCCAGAGTTAAATGGAAAGCTCACTGGAATGGCATTCCGAGTGCCTACCCCCAATGTCTCAGTTGTGGACTTGACTTGCCGCCTAGAGAAAGGGGCATCATACGATACCATCAAAGCAGCTGTGAAGTATGAACATTACATGCTCTAGTCATTGAAACAATTTTTTTTGTCTAAATTAATTTATTCTACTGATGCTGTCTAATTTTTAGCCTACCTGAAAGCTACATCTAATGTGATTTATTACACTTTTTGTAGGGCTGCATCTGAAGGTCCCATGAAAGGAATTCTGGGATACACTGAAGATGATGTTGTATCCACGGACTTTGTCGGTGATGAAAGGTAACTTCCACTTTAGAATTGCAGAAGGTAGTTTCAGTCATGCTCGCTCAAATTCCTCTAACGGTGGACTGGGTTGTTTGCAGATCAAGCATATTTGATGCTAAAGCAGGCATTGCCCTCAGTGATCGGTTTGTAAAGCTTGTTTCGT


	When sequences from multiple specimens are pooled together and share the same barcode, PURC can demultiplex them using user-provided "group" information. A reference sequence must be provided for each group, and the groups must be sufficiently phylogenetically distinct that BLAST can match each sequences to a reference sequences from the appropriate group. You need to supply the group information in the reference sequences ('group='). In this example, group "A" is the genus *Cystopteris*, group "B" is the genus *Acystopteris*, and group "C" is *Gymnocarpium*:

        >locus=ApP/group=A/ref_taxon=Cystopteris_fragilis_Utah
        TGCCACACTGGTGAGTATTGTCTTACTTTTTGTTATCCTTTTTCTTGGTGAGAAAGGGTACTGTGCATGGCATATTCACGTCAGAATCCAAGACCCCCGCTTGGGGCCGAGGGTGACAAGGATGTTTCTGTTGGGTGATACCTGTGATGCCAGTTGGAGCAAGAGTAAAATCAACTTTGTAAACATCATCTATTTGAAGGATTAACAGACATGGTATTTAAATTCCTCTCACATTCAAAACAGGGTGGTTGCAGAACTGGTATGGCCAAAGTAACGAATGCTTACGATTTGCCTGCAAGGTAAAAGTTGCACAATGCTCAAGGTGGGGCTAGTTCTTTTGTCACTTAAGCAAGGATCTTCAAGCATGTAAAATTATTCTCCCTCAACTTTGCTTTACAAAAGAAATTTAATATATTGACTACTTCATGCATGGAATTCGAGCAGCTATCACATGTTGATGTTTTTTTTTGAGGCGAGGGGTTCTTTGCATGTGGTTGTAGAAATGTTTTATCACATTTCTATGTGCTATTTTTGCATAAATGCTACGTTACAAATTAGAATTGTTTACTTGTTTGTTTGTAGGAAATCTCAAACGACTGTCTTTTGCTCTTGTATGCTTAGTTGATGATTGCATGCGTACACCTTTATGTTCATTTCAGGCTATGTTTTGTCAGCTCACAAGTTTTTGATGTTTAACCTAACATGACAGGAAAGTTATACATACTGTTGGTCCAAGATATGCTGTAAAATATCATACAGCTGCAGAAAATGCTCTAAGTCATTGTTACCGATCTTGTTTAGAGGCTTTGATTGACTTAGGCCTTCAAAGGTACCAGCTGCTTGTTTAAACAGCTCAAAATTAAAGGAGAGTGTATTCCTTTTGGATTAAAGTTTATCTCCCTTGTAATTCTTGCAGCATTGCCCTGGGGTGTATTTACACAGAGTCTAAAGGCTAT
        >locus=ApP/group=B/ref_taxon=Acystopteris_japonica_Taiwan
        TGCCACACTGGTGAGTATTGTCTTACTTTTTGTTATCCTTTTTCTTGGTGAGAAAGGGTACTGTGTATGGCATATTCACGTCATAATCCAAGACCCCCGCTTGGGGCTGGGGGGTGACAAGGATGTTTCTGTTGGGTGATACCTGTGATGCCAGTTGGAGCAAGAGTAAAATCAACTTTGTAAACATCATCTATTTGAAGGATTAACAGACATGGTATTTAAATTCCTCTCACATTCAAAACAGGGTGGTTGCAGAACTGGTATGGCCAAAGTAACGAATGCTTACGATTTGCCTGCAAGGTAAAAGTTGCACAATGCTCAAGGTGGGGCTAGTTCTTTTGTCACTTAAGCAAGGATCTTCAAGCATGTAAAATTATTCTCCCTCAACTTTGCTTTACAAAAGAAATTTAATATATTGACTACTTCATGCATGGAATTTGAGCAGCTATCACATGTTGATGTTTTTTTTTCAGGCGAGGGGTTCTTTGCATGTGGTTGTAGAAATGTTTTATCACATTTCTATGTGGGTTTTTTGACATGGCTATTTTTGCATAAATGCTACGTTACAAATTAGAATTGTTTACTTGTTTGTTTGTAGGAAATCTCAAACGACTGTCTTTTGCTCTTGTATGCTTAGTTGATGATTGCATGCGTACACCTTTATGCTCATTTCAGGCTATGTTTTGTCAGCTCACAAGTTTTTGATGTTTAACCTAACATGACAGGAAAGTTATACATACTGTTGGTCCAAGATATGCTGTAAAATATCATACAGCTGCAGAAAATGCTCTAAGTCATTGTTACCGATCTTGTTTAGAGGCTTTGATTGACTTAGGCCTTCAAAGGTACCAGCTGCTTGTTTAAACAGCTCAAAATTAAAGGAGAGTTTATTTCTTTTGGATTAAAGTTTATCTCCCTTGTATTTCTTGCAGCATTGCCCTGGGGTGTATTTACACAGAGTCTAAAGGCTAT
        >locus=ApP/group=C/ref_taxon=Gymnocarpium_dryopteris_MA
        TGCCACACTGGTGAGTATTGTCTTACTTTTTGTTATCCTTTTTCTTGGTGAGAAAGGGTACTGTGTATGGCATATTCACGTCATAATCCAAGACCCCCGCTTGGGGCTGGGGGGTGACAAGGATGTTTCTGTTGGGTGATACCTGTGATGCCAGTTGGAGCAAGAGTAAAATCAACTTTGTAAACATCATCTATTTGAAGGATTAACAGACATGGTATTTAAATTCCTCTCACATTCAAAACAGGGTGGTTGCAGAACTGGTATGGCCAAAGTAACGAATGCTTACGATTTGCCTGCAAGGTAAAAGTTGCACAATGCTCAAGGTGGGGCTAGTTCTTTTGTCACTTAAGCAAGGATCTTCAAGCATGTAAAATTATTCTCCCTCAACTTTGCTTTACAAAAGAAATTTAATATATTGACTACTTCATGCATGGAATTTGAGCAGCTATCACATGTTGATGTTTTTTTTTCAGGCGAGGGGTTCTTTGCATGTGGTTGTAGAAATGTTTTATCACATTTCTATGTGGGTTTTTTGACATGGCTATTTTTGCATAAATGCTACGTTACAAATTAGAATTGTTTACTTGTTTGTTTGTAGGAAATCTCAAACGACTGTCTTTTGCTCTTGTATGCTTAGTTGATGATTGCATGCGTACACCTTTATGCTCATTTCAGGCTATGTTTTGTCAGCTCACAAGTTTTTGATGTTTAACCTAACATGACAGGAAAGTTATACATACTGTTGGTCCAAGATATGCTGTAAAATATCATACAGCTGCAGAAAATGCTCTAAGTCATTGTTACCGATCTTGTTTAGAGGCTTTGATTGACTTAGGCCTTCAAAGGTACCAGCTGCTTGTTTAAACAGCTCAAAATTAAAGGAGAGTTTATTTCTTTTGGATTAAAGTTTATCTCCCTTGTATTTCTTGCAGCATTGCCCTGGGGTGTATTTACACAGAGTCTAAAGGCTAT

* **Map files** - One tab-delimited text file for each locus. These "maps" allow PURC to match sequences to the source specimens, based on the barcode (and the group, if used). The first column is the barcode (the barcode names need to be identical to those listed in the barcode sequence file), and the second is the specimen name. For example:

        BC01	Cystopteris_fragilis_Utah
        BC02	Cystopteris_fragilis_Arizona
        BC03	Cystopteris_fragilis_Taiwan

	When one barcode is used for multiple specimens, you need to add the group designation in the second column. In this case, one *Acystopteris* and one *Cystopteris* specimen were labeled with barcode 1, but they are divergent enough that PURC can pull them apart (based on the reference sequence file):

        BC01	A	Acystopteris_japonica_Taiwan
        BC01	B	Cystopteris_fragilis_Utah
        BC02	A	Acystopteris_japonica_Japan

    Again, in this case PURC will go through each sequence, find out what barcode it has, and what "group" it matches to (based on the reference sequence file described above). And it will then use that barcode and group information to find the corresponding taxon name using the map, above.

    If two barcodes are used (one on each primer), then the first two columns in the map files are the barcodes and the third the specimen names:

    	BCF1	BCR1	Cystopteris_fragilis_Utah
    	BCF2	BCR1	Cystopteris_fragilis_Arizona
    	BCF1	BCR2	Cystopteris_fragilis_Taiwan
    	BCF2	BCR2	Cystopteris_fragilis_AMagicPlace


    Again that if you are using the dual barcode function, the names have to be BCF1, BCF2, ..., BCF9, and BCR1, BCR2, ..., BCR9, where         the F and R denote the barcodes on forward and reverse primers respectively. This is used to check whether there are invalid F-F or R-R situation; only F-R or R-F are kept.

### Step 3: Run ###
PURC can be run by navigating to the directory that contains the configuration, barcode, reference sequence, and map files, and calling the program from there:
```bash
/Users/fayweili/Programs/purc/purc.py purc_configuration.txt
```
This assumes that the purc script (and Dependencies directory) is in /Users/fayweili/Programs/purc. **DO NOT** copy purc.py to your working directory; instead, call purc.py from there. Alternative, you can add the purc directory into your PATH, and in this case, you can run by:
```bash
purc.py purc_configuration.txt
```

After PURC finishes successfully, you should find the final clustered sequences in ```[prefix]_4_[locus]_clustered_reconsensus.fa```, and the aligned sequences in ```[prefix]_4_[locus]_clustered_reconsensus.aligned.fa```

### Step 4: recluster ###
If you want to adjust clustering parameters, instead of re-running the whole thing, you can start from the annotated fasta file by using ```purc_recluster.py.```

Usage:
```
./purc_recluster.py -f annotated_seq_file -o output_folder -c clustID1 clustID2 clustID3 clustID4 -s sizeThreshold1 sizeThreshold2 -a abundance_skew --clean
```

Example:
```
./purc_recluster.py -f purc_3_annotated.fa -o recluster -c 0.997 0.995 0.99 0.997 -s 1 4 --clean
```
```
required arguments:
  -f, --annotated_file
                        The annotated fasta file, generated by purc
  -o, --output_folder
                        The output folder
  -c, --clustering_identities
                        The similarity criterion for the first, second, third
                        and forth USEARCH clustering
  -s, --size_threshold
                        The min. number of sequences/cluster necessary for
                        that cluster to be retained (set to 2 to remove
                        singletons, 3 to remove singletons and doubles, etc)
optional arguments:
  -h, --help            show this help message and exit
  -a, --abundance_skew
                        The parameter to control chimera-killing; the default
                        is 1.9
  --clean               Remove the intermediate files

```
### Examples ###
The PURC package comes with three examples (see below). To run, simply:

```
#!shell

cd /Users/fayweili/Programs/purc/Example_1
/Users/fayweili/Programs/purc/purc.py purc_configuration.txt
```


### Example 1 - PacBio ###
In this dataset, four loci were amplified from 30 *Cystopteris* specimens. Each specimen was labeled with a unique barcode, and all the PCR reactions were pooled together and sequenced in one PacBio SMRT cell. Because each barcode corresponds to a single specimen, in the configuration file we specify:
```
Multiplex_per_barcode	= 0
```


### Example 2 - PacBio ###
In this dataset, four loci were amplified from 28 *Cystopteris*, four *Acystopteris* and 18 *Gymnocarpium* specimens. To save on barcode expenses, we applied some barcodes more than once (e.g., we used the same barcode for a *Cystopteris* and an *Acystopteris* accession). For example, G_dry_7000, A_ten_4225 and C_mil_6761 were all labeled with BC03, but were assigned as different "groups" (A, B, C) in the map files (and the appropriate reference sequences for each group were added to the reference sequence file). In the configuration file we need to specify:
```
Multiplex_per_barcode	= 1
```

And in the map files we need to add the group designation in the second column:
```
BC03	A	G_dry_7000
BC03	B	A_ten_4225
BC03	C	C_mil_6761
...
```
Also for each group of each locus, we need to include the reference sequences, so that PURC can tell them apart:
```
>locus=ApP/group=A/ref_taxon=yourFavGymnocarpium
TGCCACACTGGTGAGTATTGTCTTACTTTTTGTTATCCTTTTTCTTGGTGAGAAAGGGTACTGTGCATGGCATATTCACGTCAGAATCCAAGACCCCCGCTTGGGGCCGAGGGTGACAAGGATGTTTCTGTTGGGTGATACCTGTGATGCCAGTTGGAGCAAGAGTAAAATCAACTTTGTAAACATCATCTATTTGAAGGATTAACAGACATGGTATTTAAATTCCTCTCACATTCAAAACAGGGTGGTTGCAGAACTGGTATGGCCAAAGTAACGAATGCTTACGATTTGCCTGCAAGGTAAAAGTTGCACAATGCTCAAGGTGGGGCTAGTTCTTTTGTCACTTAAGCAAGGATCTTCAAGCATGTAAAATTATTCTCCCTCAACTTTGCTTTACAAAAGAAATTTAATATATTGACTACTTCATGCATGGAATTCGAGCAGCTATCACATGTTGATGTTTTTTTTTGAGGCGAGGGGTTCTTTGCATGTGGTTGTAGAAATGTTTTATCACATTTCTATGTGCTATTTTTGCATAAATGCTACGTTACAAATTAGAATTGTTTACTTGTTTGTTTGTAGGAAATCTCAAACGACTGTCTTTTGCTCTTGTATGCTTAGTTGATGATTGCATGCGTACACCTTTATGTTCATTTCAGGCTATGTTTTGTCAGCTCACAAGTTTTTGATGTTTAACCTAACATGACAGGAAAGTTATACATACTGTTGGTCCAAGATATGCTGTAAAATATCATACAGCTGCAGAAAATGCTCTAAGTCATTGTTACCGATCTTGTTTAGAGGCTTTGATTGACTTAGGCCTTCAAAGGTACCAGCTGCTTGTTTAAACAGCTCAAAATTAAAGGAGAGTGTATTCCTTTTGGATTAAAGTTTATCTCCCTTGTAATTCTTGCAGCATTGCCCTGGGGTGTATTTACACAGAGTCTAAAGGCTAT-------------------------------------------------------------------------------------
>locus=ApP/group=B/ref_taxon=tokenAcystopteris
TGCCACACTGGTGAGTATTGTCTTACTCTTTTTCTTGGTGAGAAAGGGTAAAGTGCATGGCATATTCACGTCAACCCCCGCTTGGGGCCGGGAGTGACAAGGATGTTACTGTTGGGTGATACCTGTGATGCCCAGTTGGAGCAGGAGTAAAATCGACTTTGTAAATATCATTTATTTGAGGGATTAACAGACATGGTATTTAAATTCCTCACATTCAAAACAGGGTGGTTGCAGAACTGGTATGGCCAAAGTAACGAACGCTTACGATTTGCCTGCAAGGTAAAAGTTGCACAATGCTCAAGGTGGGGCTAGTTCTTGTCACGTAAGCAAGGATCATCAAGCATGTAAAATTATTCTCCCTCAACTTTGCTTTACAAAAGAAATTTAACATATTGACCACTTCATCCATGGAATTTGAGCAGCTATCACATGTCGATGTTTTTTTTTGAGGCGAGGGGTTCTTGTTTGCATGTGGTTGTAGAAATGTTTTATCACATTTCTATTCGGATTTTTTGACATGGCTATTTTTGCATAAATGCTAAGTTCCAAATTAGAATTGTTTACTTGTTTGTTTGTAGGAAATCTCAAACGACTGTCTTTTGCTCTTGTATGCTCATTTGATGATTGCATGCGTACACCTTTATGTTCATTTCAGGCTATGTTTTGTCAGCTCACAAGTTTTTGATGTTAAACCTAACATGACAGGAAAGTTATACATACTGTTGGTCCAAGATATGCTGTAAAATATCATACAGCTGCAGAAAATGCTCTAAGTCATTGTTACCGATCTTGTTTAGAGACTTTGATTGACTTAGGCCTTCAAAGGTACCAGCTGCTTGTTTAAACAGATCAAAATTAAAGAGAGTTTTTTCCTTTTGGATTAAAGTTTATCTCCCTTGTAATTCTTGCAGCATTGCCCTGGGGTGTATTTACACAGAGTCTAAAGGCTAT------------------------------------------------------------------------------------------
>locus=ApP/group=C/ref_taxon=CystopterisBoss
TGCCACACTGGTGAGTATTATCTTACTCTTTTTCTTGGTGAGAAAGGGTAGTGTGCATGGCATATTCATGTCAACCCCCGCTTGGCACCGGGGGTGATAAGGATGTTACTGTTGGGTGATACCTGTGATGCCCAATTGGAGCAAATAGTAAAATCAACTTCGTAAATATCATTTATTTGAAGGATTAATAGACATGGTATTTAAATTCCTCTCACATTCAAAACAGGGTGGTTGCAGAACTGGTATGGCCAAAGTAACGAACGCTTACGATTTGCCTGCAAGGTAAAAGTTGAACAATGCTCAAGGTGGGGCTACTCCTTTTGTCACTTAAGCAAGGATCTTCAAGCATGTAAAATTATTCTCCCTCAACTTTGCTTTACAAAAGAAATTTAACATATTGACTACTTCATCCATGGAATTTGAGCAGCTACACATGTAGATGTTTTTTTTTTGAGGCAAGGGGTTCTTGTTTGCATGTGGTTGTAGAAATGTTTATCACATTTCTATTTGGGTCTTTTGACATGGCTATTTTTGCATAACTGCCAAGTTACAAATTAGAATTGTTCACTTGCTGTTTGTACGAAATCTCAATTGACTGTCTTTTGCTCTTGTATGCTCATTTGATGATTGCATGCGTACACCTTTATGTTCATTTCAGGCTATGTTTTGTCAGCTCACAAGATTTTGATGCTAAACCTAACGTGACAGGAAAGTTATACATACTGTTGGTCCAAGATATGCTGTAAAATATCATACAGCTGCAGAAAATGCTCTAAGTCGTTGCTACCGATCTTGTTTAGAGGCTTTGATTGACTTAGGCCTTCAAAGGTACCAGCTGCTTGTTTAACAGCTCAAAATTAAAGGAGAGTTTTTTCCTTTTGGATTAAAGTTTATCTCCCATGTAATTCTTGCAGCATTGCCCTGGGGTGTATTTACACAGAGTCTAAAGGCTAT---------------------------------------------------------------------------------------
...
```

### Example 3 - Illumina ###
PURC can work with Illumina data too. You'd need to merge the pair-end reads first, for example, using [PEAR](https://github.com/xflouris/PEAR). In this dataset, each sequence is labeled by two barcodes, one on the forward and the other on the reverse primer. In the map file, we tell PURC which barcode combination corresponds to which specimen. We invoke dual-barcode mode by:
```
Dual_barcode		= 1
```
And in the map file:
```
BCF1	BCR1	species1
BCF2	BCR1	species2
BCF3	BCR1	species3
...
```

Again that if you are using the dual barcode function, the names have to be BCF1, BCF2, ..., BCF9, and BCR1, BCR2, ..., BCR9, where         the F and R denote the barcodes on forward and reverse primers respectively. This is used to check whether there are invalid F-F or R-R situation; only F-R or R-F are kept.

### Example 4 - PacBio with ASV (FASTQ only) ###
Newer PacBio CCS (HiFi) sequencing data is produced with quality scores in FASTQ format. This type of data can be processed by the 'DADA2' R package to reduce reads into amplicon sequence variants (ASVs), which tend to more accurately represent the source compared to OTU clustering. To turn on ASV inference, change to `Clustering_method = 0` in your config file. Some settings, such as minimum/maximum length of reads to keep, and maximum number of expected errors to allow, can be set with `minLen`, `maxLen`, and `maxEE`, respectively.



### Citation ###
This script relies heavily on VSEARCH, MUSCLE, DADA2, and BLAST.
If this script assisted with a publication, please cite the following papers
(or updated citations, depending on the versions of VSEARCH, etc., used).

PURC:  
Rothfels, C.J., K.M. Pryer, and F-W. Li. 2016. Next-generation polyploid
phylogenetics: rapid resolution of hybrid polyploid complexes using PacBio
single-molecule sequencing. New Phytologist

VSEARCH:  
Rognes T, T. Flouri, B. Nichols, C. Quince, and F. Mahé. 2016. VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584 [https://doi.org/10.7717/peerj.2584](https://doi.org/10.7717/peerj.2584)

Cutadapt:  
Martin, M. 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads.
EMBnet.journal 17:10-12.

MUSCLE:  
Edgar, R.C. 2004. MUSCLE: Multiple sequence alignment with high accuracy and high throughput.
Nucleic Acids Research 32:1792-1797.

BLAST:  
Camacho, C., G. Coulouris, V. Avagyan, N. Ma, J. Papadopoulos, et al. 2009.
BLAST+: Architecture and applications. BMC Bioinformatics 10: 421.

DADA2:  
Callahan, B., McMurdie, P., Rosen, M. et al. 2016. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583. [https://doi.org/10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)

**Deprecated dependencies**  
USEARCH/UCLUST:  
Edgar, R.C. 2010. Search and clustering orders of magnitude faster than BLAST.
Bioinformatics 26(19), 2460-2461.

UCHIME:  
Edgar, R.C., B.J. Haas, J.C. Clemente, C. Quince, R. Knight. 2011.
UCHIME improves sensitivity and speed of chimera detection, Bioinformatics 27(16), 2194-2200.


### FAQ ###
* What does "ERRmidBC" flag mean?

     It indicates that PURC identified a barcode sequence in the middle of the read (not at the ends). The reason could be PCR/sequencing artifacts, or a stretch of reads that resemble one of the barcodes. If latter, then use

        Barcode_detection = 1

     in configuration file to restrict barcode identification to the ends of sequences.

* What are the output files produced (e.g., BC16_SlC1_0.997dCh1Ss1C2_0.995dCh2Ss2C3_0.99dCh3Ss3.fa, the .uc ones, etc)?

     PURC produces a number of intermediate files during the clustering/chimera-killing steps. The '.uc' and '.uchime' files are from USEARCH and UCHIME, respectively, and contain detailed clustering and chimera detection results. The fasta files from each clustering/chimera-killing steps are saved as well. 'Sl' means fasta sorted by length, 'C1_0.997' means the first clustering with identity of 0.997, 'dCh1' means the first chimera detection, 'Ss' means fasta sorted by cluster size. If you don't want to keep all these intermediate files (to save space for example), in configuration file use

        Remove_intermediates = 1


* Why do some of my sequences come out of the pipeline much shorter than they went in?

     Barcodes and primer sequences will be trimmed in the pipeline, so the final sequences should be somewhat shorter
than the ones that enter the pipeline. However, if sequences are much shorter, that suggests that the primer-detection
settings are too lax, and thus primers are being "found" in the middle of the sequences. If this happens, then that region and everything downstream will be erased. To fix this problem, you have to manual change the cutadapt settings in purc.py. There are instructions on how to do this in purc.py itself, around line 528.


### Who do I talk to? ###
Peter Schafran ([ps997@cornell.edu](mailto:ps997@cornell.edu))

Fay-Wei Li ([fl329@cornell.edu](mailto:fl329@cornell.edu))

Carl Rothfels ([crothfels@berkeley.edu](mailto:crothfels@berkeley.edu))

### Advanced Installation ###
If PURC is not functioning correctly when installed with the YAML files, you can try manual installation. Run these commands depending on your OS:
```bash
# macOS
conda create -n purc -c bioconda -c conda-forge cutadapt blast muscle vsearch r-base=4.1 r-essentials bioconductor-dada2 r-ggplot2 r-reshape2 r-gridextra r-rcolorbrewer python">=3.7"

# Linux
conda create -n purc -c bioconda -c conda-forge cutadapt blast muscle vsearch r-base=4.1 r-essentials bioconductor-dada2 r-ggplot2 r-reshape2 r-gridextra r-rcolorbrewer python">=3.7" lima

# macOS or Linux
conda activate purc && pip install BioPython
```
If R is already installed, it is best not to install multiple instances. To install without R, use these commands:
```bash
# macOS
conda create -n purc -c bioconda -c conda-forge cutadapt blast muscle vsearch python">=3.7"

# Linux
conda create -n purc -c bioconda -c conda-forge cutadapt blast muscle vsearch python">=3.7" lima

# macOS or Linux
conda activate purc && pip install BioPython
```
Make sure `R` and `Rscript` are in your PATH (type the command and it runs from anywhere). E.g. if installed with the installer on macOS, you may need to add `/Library/Frameworks/R.framework/Versions/4.0/Resources/`. Note this will need be redone each time your open a new Terminal:
```bash
PATH="/Library/Frameworks/R.framework/Versions/4.0/Resources/:$PATH"
```
Once R is installed and accessible from the command line, open `R` and run these commands to install required packages. If install fails, try changing BiocManager version based on your R version.
```R
# R version = BiocManager Version
# 4.0.2+    = 3.12
# 4.0       = 3.11
# 3.6       = 3.10
# 3.5       = 3.8
# 3.4       = 3.6
# 3.3       = 3.4
# 3.2       = 3.2

### Install ###
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", quiet = TRUE)
BiocManager::install(version = '3.13', ask = FALSE)
if (!require("dada2", quietly = TRUE)) BiocManager::install("dada2", ask = FALSE)
if (!require("gridExtra", quietly = TRUE)) install.packages("gridExtra", quiet = TRUE)
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2", quiet = TRUE)
if (!require("reshape2", quietly = TRUE)) install.packages("reshape2", quiet = TRUE)
if (!require("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer", quiet = TRUE)
```
If you have another R installation, it may have a different library path that can be prioritized over the conda installed path (even in a conda environment). If R errors appear in the DADA2 logs, such as `Error: package or namespace load failed for ‘Rcpp’: package ‘Rcpp’ was installed before R 4.0.0: please re-install it` try entering an interactive R session and type:
```R
.libPaths()
```
You may see output like this:
```R
> .libPaths()
[1] "/home/ps997/R_libs"                            
[2] "/home/ps997/miniconda3/envs/purc/lib/R/library"
```
showing that a previous R library is prioritized over the conda installed one. You can reorder the list by editing your `.Renviron` file located in your home (~) directory. **Backup the original file**, and change the `R_LIBS` variable to the purc install path.
```bash
R_LIBS=/home/ps997/miniconda3/envs/purc/lib/R/library
```
Now when R is started, the `.libPaths()` command should return:
```R
> .libPaths()
[1] "/home/ps997/miniconda3/envs/purc/lib/R/library"
```

If you get an error message like "ERROR: Cython is not installed", install/update [Cython](http://docs.cython.org/src/quickstart/install.html) and try again.

### Known Bugs ###
* Inconsistent number of sequences reported during barcode removal (does not always equal total number of sequences)
* CCS files without forward-slashes around ZMW code in read name will break lima
* Fix vsearch output in log
