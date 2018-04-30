# shRNAscreen
A computational pipeline for analyzing shRNA-Seq generated from in vivo shRNA screen. All 19-mers were firstly extracted from short reads of shRNA-Seq and collapsed, then searching against annotation of shRNAs sequence using R package "stringdist" with mismatches allowed.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

## Prerequisites

### Operating Systems

#### Supported Unix distributions

- Ubuntu
- CentOS
- Red Hat Enterprise Linux (please use the CentOS packages and instructions)

#### Job scheduler

- Univa Grid Engine
- TORQUE Resource Manager

#### Tools or packages
- python >= 3.5.1 (for snakemake)
- R >= 3.1.0
- [pigz](https://zlib.net/pigz/) >= 2.3.1
- [mawk](http://invisible-island.net/mawk/) >= 1.3.4
- [snakemake](http://snakemake.readthedocs.io/en/stable/index.html) >= 3.12.0

#### Libraries or modules

##### R
- [stringdist](https://cran.r-project.org/web/packages/stringdist/)

can be installed by following commands in R:

```{r message = FALSE}
install.packages("stringdist");
```
## Installing

### Install snakemake

Install snakemake into a virtual environment 

```
git clone https://bitbucket.org/snakemake/snakemake.git
cd snakemake
virtualenv -p python3 snakemake
source snakemake/bin/activate
python setup.py install
```
Or you can install it using bioconda

```
conda install snakemake
```

### Download scripts and configuration files from github and add directory of scripts into PATH variable 

```
git clone https://github.com/shenyang1981/shRNAscreen.git
cd shRNAscreen/; export PATH="${PWD}/scripts:$PATH"
```
You may consider put 'export PATH=${PWD}/scripts:$PATH' <replace ${PWD}/scripts with real path to shRNAscreen scripts> into your .bashrc file.
  
### Prepare annotation file for shRNA library.

* hairpin_ERWOOD_good.txt -- shRNA library annotation

The annotation file includes three columns: shRNA ID, sequence and corresponding gene name (**no header**).

ERWOOD_1 |GAGAAGATCCTCTTCATCA|Gas2
:--------|:------------------|:-------
ERWOOD_2 |AGCTTTGACCAGCTTCTTC|Crtc3
ERWOOD_3 |GCCAGCGAATGCAGTACAT|Vmn1r181
ERWOOD_4 |AGGACACATCTGCCAGCAT|Msl3
ERWOOD_5 |AGCAGTACAGGCTGGTACA|Gabrr1
ERWOOD_6 |AGGCTCATACTCTCCTTCT|Dcpp3
ERWOOD_7 |GCCAGGATGTGACTCAGAT|Tmem171
ERWOOD_8 |CATGGAGAAGTACAACATA|Dynll2
ERWOOD_9 |GCCTTCATCATTGGTGCAG|Kcnj2
ERWOOD_10|ACAGAAACATTAGAATTAC|Mtpap



### Prepare input files and sample information

* sampleList.txt -- information of each sequenced library, including species (**Species**), library ID (**LibID**), sequencing batch (**SeqBatch**), Analytic ID(**AnalyticID**). Samples belonged to **the same AnalyticID** would be selected for searching against same shRNA annotation. 

The format of sampleList.txt is like:

Species|LibID |SeqBatch|AnalyticID
:------|:-----|:-------|:---------
Mouse  |AML025_rep1|AML025  |batch3
Mouse  |AML025_rep2|AML025  |batch3


** Note: LibID should be unique as the corresponding sequence file should be named as {LibID}.fastq.gz.

* input reads files -- Reads are single-end. Name of each file should be {LibID}.fastq.gz (LibID should be the same as in sampleList.txt). All of reads files from the same sequencing batch should be put into one folder named by {SeqBatch} as indiciated in the sampleList.txt. For example, reads files "AML025_rep1.fastq.gz", "AML025_rep2.fastq.gz" can be put into folder "input/AML025/"

```
tree input/
input/
├── AML025
│   ├── AML025_rep1.fastq.gz
│   └── AML025_rep2.fastq.gz
└── sampleList.txt
```

### generate config file 

To generate a configuration file for snakemake, several variables need to be defined:
- SHRNASCREEN_ANADIR: path to root folder of analysis
- SAMPLEINFO: path to the sampleList.txt file
- HAIRPIN: path to the file of shRNA library annotation
- HPID: shRNA library ID of screen
- ANAID: analytic ID indicating which libraries should be selected


Configuration file can be generated using script **generateConfigureFile.sh**

```
SHRNASCREEN_ANADIR=. SAMPLEINFO=./input/sampleList.txt HAIRPIN=annotation/hairpin_ERWOOD_good.txt HPID=ERWOOD ANAID=batch3 generateConfigureFile.sh pipeline/configTemplate/conf-shRNAScreen.json > config-batch3.json
```
**config-batch3.json**

Now, files should be organized like:

```
.
├── annotation
│   └── hairpin_ERWOOD_good.txt
├── config-batch3.json
├── input
│   ├── AML025
│   │   ├── AML025_rep1.fastq.gz
│   │   └── AML025_rep2.fastq.gz
│   └── sampleList.txt
├── LICENSE
├── pipeline
│   ├── configTemplate
│   │   └── conf-shRNAScreen.json
│   └── counthpreads.sk
├── README.md
└── scripts
    ├── extractMatureSimple.sh
    ├── generateConfigureFile.sh
    └── summaryShCount.R
```
## Running Pipeline

### Local Mode
The pipeline can be simply run in local mode with the configuration file.

```
source {$pathtosnakemake}/snakemake/bin/activate
snakemake -s pipeline/counthpreads.sk --configfile config-batch3.json -j 32
```

### Submit snakemake jobs to Cluster
Or run it by submitting to the job scheduler

```
runsnake.sh pipeline/parcek.sk conf.batch1.json testjob 24 24
```

### Results

After running pipeline, results would be stored in "result/" and "hpreads" folder.

- **AML025_rep1.mature.count.txt.gz** -- collapsed reads.
- **AML025_rep1_vs_ERWOOD_hp.result.count.txt** -- shRNA hits.

```
.
hpreads/
└── AML025
    ├── AML025_rep1.mature.count.txt.gz
    └── AML025_rep2.mature.count.txt.gz
result/
└── batch3
    └── AML025
        ├── AML025_rep1_vs_ERWOOD_hp.result.count.txt
        ├── AML025_rep1_vs_ERWOOD_hp.result.count.txt.all
        ├── AML025_rep1_vs_ERWOOD_hp.result.count.txt.allmature
        ├── AML025_rep2_vs_ERWOOD_hp.result.count.txt
        ├── AML025_rep2_vs_ERWOOD_hp.result.count.txt.all
        └── AML025_rep2_vs_ERWOOD_hp.result.count.txt.allmature
```

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Yang Shen** - Develope the pipeline

## Contact

Please contact us if you find bugs, have suggestions, need help etc. You can either use our mailing list or send us an email:

* [Yang Shen](mailto:sheny@gis.a-star.edu.sg)

shRNAscreen is developed in the Genome Institute of Singapore

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
