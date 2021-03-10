# De novo Assemblers
This is a comprehensive document that aims to record the installation/configuration and usage of De novo assembly tools. In case if you are using conda package, it's always desirable to make a separate environment.

1. **Allpaths-LG:** Source Code: [Here](http://software.broadinstitute.org/allpaths-lg/blog/?page_id=12)
ALLPATHS‐LG requires a minimum of 2 paired‐end libraries – one short and one long. The short library average separation size must be slightly less than twice the read size, such that the reads from a pair will likely overlap – for example, for 100 base reads the insert size should be 180 bases.

**Steps:**

```bash
##Installing GMP
sudo apt install -y libgmp-dev
## Installing MPFR
sudo apt install libmpfr-dev
# MPC
sudo apt install libmpc-dev
# Picard
sudo apt-get install -y picard
# graphviz 
sudo apt install -y graphviz
# g++ (v4.7.0 or higher according to AllPaths-LG manual but higher version makes ALLpaths to fail). In Ubuntu 16 and 20, installing g++4.8 is difficult
# IMPORTANT:
# Ubuntu 18.04 comes with gcc v7.3.0
# use v7.3.0 causes "make" to fail with these error messages
# 
# Makefile:2654: recipe for target 'CleanEfasta.o' failed
# make[1]: *** [CleanEfasta.o] Error 1
# make[1]: Leaving directory '/usr/local/allpathslg-52488/src'
# Makefile:284: recipe for target 'all-recursive' failed
# make: *** [all-recursive] Error 1
#
# install g++ v4.8 and use this version in the "./configure" step (see below; CXX/CXXPP options)
sudo apt install -y g++-4.8
# libieee
# IMPORTANT:
# Ubuntu 18.04 uses glibc v2.27
# which does not support libieee and causes "make" to fail with these error messages
#
# /usr/bin/ld: cannot find -lieee
# collect2: error: ld returned 1 exit status
# Makefile:2096: recipe for target 'RemodelGaps' failed
# make[1]: *** [RemodelGaps] Error 1
# make[1]: Leaving directory '/usr/local/allpathslg-52188/src'
# Makefile:284: recipe for target 'all-recursive' failed
# make: *** [all-recursive] Error 1
#
# install libieee1284-3 and create a symbolic link
sudo apt install -y libieee1284-3
sudo ln -s /usr/lib/x86_64-linux-gnu/libieee1284.so.3 /usr/lib/libieee.so
#
# start AllPaths-LG installation
sudo wget ftp://ftp.broadinstitute.org/pub/crd/ALLPATHS/Release-LG/latest_source_code/LATEST_VERSION.tar.gz
sudo tar xzf LATEST_VERSION.tar.gz
sudo chown -R root:root allpathslg-52488
sudo chmod 755 allpathslg-52488
cd allpathslg-52488
sudo ./configure --prefix=/usr/local CXX=/usr/bin/g++-4.8 CXXPP=/usr/bin/cpp-4.8
sudo make
sudo make install
```
**Other means:**
Docker Image: [here](https://github.com/remiolsen/allpathslg)

Conda: [here](https://anaconda.org/biobuilds/allpathslg)

How to use: [Thread 1](https://www.biostars.org/p/320596/), [Thread 2](http://software.broadinstitute.org/allpaths-lg/blog/?page_id=336), [Thread 3](http://evomics.org/wp-content/uploads/2012/01/Allpaths_exercises.pdf), [Batch Script for cluster](https://www.msi.umn.edu/sw/allpaths-lg)

2. **Velvet**

**Steps**
```bash
wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz
tar zxvf velvet_1.2.10.tgz
sudo chown -R root:root velvet_1.2.10
sudo chmod 755 velvet_1.2.10
cd velvet_1.2.10
##compilation settings for human genome
#choose the appropriate Kmer length according to your data. Appropriate Kmer length can be found by Velvet Advisor:http://dna.med.monash.edu.au/~torsten/velvet_advisor/
sudo make 'CATEGORIES=10' 'MAXKMERLENGTH=301' 'LONGSEQUENCES=1' 'OPENMP=1'
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc

#Velvet Optimiser: VelvetOptimiser is a multi-threaded Perl script for automatically optimising the 
#three primary parameter options (K, -exp_cov, -cov_cutoff) for the Velvet de novo sequence assembler.
#install cpanm & bioperl. This will install BioSeqIO too required by the tool
sudo apt install -y cpanminus
sudo apt-get install -y bioperl
#Velvet Optimizer
wget https://github.com/tseemann/VelvetOptimiser/archive/master.zip
unzip master.zip
cd VelvetOptimiser-master/
dir=$(pwd)
export PATH="$PATH:$dir"

```
**Other means**
Velvet Conda: [here](https://anaconda.org/bioconda/velvet)
Velvet Optimizer: [here](https://anaconda.org/bioconda/perl-velvetoptimiser)
How to use: [K-mer size](http://seqanswers.com/forums/showthread.php?t=67844), [tutorial](https://bpa-csiro-workshops.github.io/btp-manuals-md/modules/btp-module-velvet/velvet/), PBS [Script](https://wiki.gacrc.uga.edu/wiki/VelvetOptimiser-Sapelo2)

3. **IDBA**

**Steps**

```bash
git clone https://github.com/loneknightpy/idba.git
cd idba
./build.sh
cd bin
```

**Other means**
Conda: [here](https://anaconda.org/bioconda/idba)
How to use: [Source 1](https://dzif-metagenomics-workshop.readthedocs.io/en/latest/assembly/idba_ud.html)

4. **Abyss**

**Steps**
```bash

sudo apt-get install -y zsh
sudo apt-get install -y abyss

```
5. **SOAPDenovo2**

**Steps**

```bash
sudo apt-get install -y soapdenovo2
#WARNING: In case if you are using conda beware it will downgrade your samtools version.
```
usage: SOAPdenovo-127mer or SOAPdenovo-63mer
##Warning 

**Other means** Conda: [here](https://anaconda.org/bioconda/soapdenovo2)

6. **SGA**
**Steps**
```bash
conda create -n sga
conda activate sga
conda install mamba -c conda-forge
mamba install -c bioconda sga
```
7. **Celera Assembler**

```bash
wget https://master.dl.sourceforge.net/project/wgs-assembler/wgs-assembler/wgs-8.3/wgs-8.3rc2.tar.bz2
tar -jxvf wgs-8.3rc2.tar.bz2
cd wgs-8.3rc2/
cd kmer && make install && cd ..
cd src && make && cd ..
cd Linux-amd64/bin
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```
8. Masurca

```bash
git clone https://github.com/alekseyzimin/masurca
git clone https://github.com/swig/swig.git
cd swig
./autogen.sh
./configure
sudo apt-get install bison flex
make
sudo make install
sudo apt-get install -y yaggo
cd ..
cd masurca
git submodule init
git submodule update
make
tar zxf MaSuRCA-3.4.2.tar.gz
cd MaSuRCA-3.4.2/
./install.sh
cd bin
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```
9. Ray

```bash
wget https://master.dl.sourceforge.net/project/denovoassembler/Ray-v2.0.0.zip
sudo apt-get install -y mpi
unzip Ray-v2.0.0.zip
cd Ray-v2.0.0/
make PREFIX=ray-build
sudo make install
ls ray-build
#testing it
#mpiexec -n 1 ray-build/Ray -o test -p test_1.fastq test_2.fastq -k 31 
cd ray-build/
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```
10. **ssake**
```bash
wget https://github.com/bcgsc/SSAKE/releases/download/v4.0/ssake_v4-0.tar.gz
tar zxf ssake_v4-0.tar.gz
cd ssake/
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```

11. Edena

```bash
wget http://www.genomic.ch/edena/EdenaV3.131028.tar.gz
tar zxf EdenaV3.131028.tar.gz
cd EdenaV3.131028/
make
cd bin
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```

12. VCAKE

```bash
wget https://master.dl.sourceforge.net/project/vcake/hybrid%20assembly%20pipeline/AssemblyPipeline_v1.1/AssemblyPipeline_1.1.tar.gz
tar zxf AssemblyPipeline_1.1.tar.gz
cd AssemblyPipeline/
sudo apt-get install ncbi-blast+
gcc -o soapsorttext soapsorttext.c
gcc -o soapsortnum soapsortnum.c
chmod +x *
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
##In case if tool doesnot work because of NEWBLER Assembler, abondon this tool because NEWBLER assembler of Roche is no longer available.
```
13. 
```bash
wget http://sharcgs.molgen.mpg.de/software/1.2.11/sharcgs.pl
chmod +x sharcgs.pl
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```
14. SAGE2
```bash
conda create -n sage2
conda activate sage2
conda install mamba -c conda-forge
mamba install -c mmolnar racer
mamba install -c mmolnar sage2
```
15. JR-Assembler
```bash
wget http://jr-assembler.iis.sinica.edu.tw/packages/JR-Assembler_v1.0.4.tar.gz
tar zxf JR-Assembler_v1.0.4.tar.gz
cd JR-Assembler_v1.0.4/
wget ftp://155.52.47.33/pub/bio/tgi/software/seqclean/mdust.tar.gz
cd mdust
make
chmod +x *
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
##incomplete
```
16. Taipan
```bash
mkdir taipan
cd taipan
wget https://master.dl.sourceforge.net/project/taipan/taipan/taipan%201.0/taipan.tar
tar xf taipan.tar
make
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```
17. SOPRA
```bash
wget http://www.physics.rutgers.edu/~anirvans/SOPRA/SOPRA_v1.4.6.zip
unzip SOPRA_v1.4.6.zip
cd source_codes_v1.4.6/
```
18. Arachane 
```bash
wget ftp://ftp.broadinstitute.org/pub/crd/ARACHNE/latest_source_code/arachne-46233.tar.gz
cd arachne-46233/
##incomplete
```
19. DISCOVER denovo

```bash
wget ftp://69.173.70.223/pub/crd/DiscovarDeNovo/latest_source_code/LATEST_VERSION.tar.gz
tar zxf LATEST_VERSION.tar.gz
cd discovardenovo-52488/
./configure --prefix=/usr/local CXX=/usr/bin/g++-4.8 CXXPP=/usr/bin/cpp-4.8
make all
sudo make install
##DiscovarDeNovo
```
20. Light Assembler
```bash
git clone https://github.com/SaraEl-Metwally/LightAssembler.git
cd LightAssembler/
##Default k <= 31. Use make k=kmersize for k > 31, e.g. make k=49
make 
```
21. GAML
```bash
git clone https://github.com/usamec/GAML.git
cd GAML/
cmake .
make
```
22. Minia
```bash
wget https://github.com/GATB/minia/releases/download/v3.2.4/minia-v3.2.4-bin-Linux.tar.gz
cd minia-v3.2.4-bin-Linux
cd bin
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```
23. ngopt
```bash
wget https://excellmedia.dl.sourceforge.net/project/ngopt/a5_miseq_linux_20160825.tar.gz
tar zxf a5_miseq_linux_20160825.tar.gz
cd a5_miseq_linux_20160825/bin/
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```
24. SWAP2
```bash
wget https://master.dl.sourceforge.net/project/swapassembler/SWAP2/SWAP2.tar.bz2
bunzip2 -d SWAP2.tar.bz2
tar xf SWAP2.tar
cd SWAP2/
make
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```
25. Stride
```bash
mkdir stride
cd stride
wget https://raw.githubusercontent.com/ythuang0522/StriDe/master/StriDe_Linux64bit_v1.0.tar.gz
tar zxf StriDe_Linux64bit_v1.0.tar.gz
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```
26. Meraculous-2D
```bash
conda create -n meraculous
conda activate meraculous
conda install mamba -c conda-forge
mamba install -c bioconda meraculous
##Usage
meraculousTh_56mer
meraculousTh_56mer
```
27. ScaffoldScaffolder
```bash
wget http://bioresearch.byu.edu/scaffoldscaffolder/scaffoldscaffolder-0.1.tar.gz
tar zxf scaffoldscaffolder-0.1.tar.gz
cd ScaffoldScaffolder/
sudo make compile
##usage: java -Xmx2g -cp class::lib/* scaffoldscaffolder.ScaffoldScaffolder
```
28. EPGA2
```bash
git clone https://github.com/bioinfomaticsCSU/EPGA2.git
cd EPGA2/
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
##usage: perl EPGA.pl
```
29. Pandaseq
```bash
git clone https://github.com/neufeld/pandaseq.git
cd pandaseq/
sudo apt-get install build-essential libtool automake zlib1g-dev libbz2-dev pkg-config
./autogen.sh && ./configure && make && sudo make install
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```
30. Readjoiner
```bash
sudo apt install genometools
#usage: gt readjoiner -help
```
31. Fermi

```bash
wget https://github.com/downloads/lh3/fermi/fermi-1.1.tar.bz2
bunzip2 -d fermi-1.1.tar.bz2
tar xf fermi-1.1.tar
cd fermi-1.1/
make
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```
32. Oases

```bash
git clone --recursive https://github.com/dzerbino/oases
cd oases
make
./oases --version
```
 33. PhrapUmd2
 
 ```bash
 wget http://terpconnect.umd.edu/~ALEKSEYZ/PhrapUMDV2/PhrapUMDV2.tgz
 tar zxf PhrapUMDV2.tgz
 cd PhrapUmd2/
 ./install.sh
 ```

34. curtain
```bash
sudo apt install maq
sudo apt install velvet
wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/curtain/curtain-exec_0.2.3-BETA-dist.tar.gz
tar zxf curtain-exec_0.2.3-BETA-dist.tar.gz
cd curtain-exec_0.2.3-BETA/bin
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```
35. Price
```bash
wget http://derisilab.ucsf.edu/software/price/PriceSource140408.tar.gz
tar zxf PriceSource140408.tar.gz
cd PriceSource140408/
make
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```
36. Swalo
```bash
wget https://atifrahman.github.io/SWALO/swalo-0.9.7-beta.tar.gz
tar zxf swalo-0.9.7-beta.tar.gz
cd swalo
make
```
For Long Reads

 37. Canu
```bash
wget https://github.com/marbl/canu/releases/download/v2.1/canu-2.1.Linux-amd64.tar.xz
tar -xJf canu-2.1.Linux-amd64.tar.xz
cd canu-2.1/bin
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```
38. PCAP
```bash
wget http://seq.cs.iastate.edu/PCAP/pcap.rep.linux.xeon64.tar
tar xf pcap.rep.linux.xeon64.tar
cd pcap.rep.linux.xeon64/
dir=$(pwd)
export PATH="$PATH:$dir"
source ~/.bashrc
```
39. HapCol
```bash
wget https://github.com/AlgoLab/HapCol/archive/v1.1.1.tar.gz
tar zxf v1.1.1.tar.gz
cd HapCol-1.1.1/
mkdir -p build
cd build
cmake ../src
make
```
40. pb-assembly
```bash
conda create -n pbassembly
conda activate pbassembly
conda install pb-assembly
```
41. smartdenovo
```bash
git clone https://github.com/ruanjue/smartdenovo.git
cd smartdenovo/
#make minor modification in make file and wtlay.h as mentioned here: https://github.com/ruanjue/smartdenovo/issues/18
#1.  Changes in the  `Makefile`  as mentioned here
#2.  Removing lines 95 to 97 and 832 to 834 in  `wtlay.h`  as mentioned here
#3.  Add  `static`  before  **all**  five occurrences of  `inline`  (lines 491 to 518) in  `wtlay.h`
make
```
42. NextDenovo
```bash
pip install psutil
pip install drmaa
wget https://github.com/Nextomics/NextDenovo/releases/download/v2.3.1/NextDenovo.tgz
tar zxf NextDenovo.tgz
cd NextDenovo/
