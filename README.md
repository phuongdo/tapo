# TAPO

## 1. Downloading TAPO
### 1.1 System Requirements

* Java Development Kit (JDK) at least JDK 1.6
* Memory: No minimum requirement
* Disk: Approximately 100MB
* Operating System: Linux OS, Ubuntu 12, CentOS 6.0 (x32 or x64). Window 7 (non-tested)

### 1.2 Dependencies

* 3DCOMB (http://ttic.uchicago.edu/~jinbo/software.htm )
* CD-Hit (http://weizhong-lab.ucsd.edu/cd-hit/download.php)
* DSSP (http://swift.cmbi.ru.nl/gv/dssp/DSSP_5.html )
* MUSTANG(http://www.csse.monash.edu.au/~karun/Site/mustang.html)
* Wget ( yum install wget in Centos)

```
$ cd $YOUR_INSTALL_DIR
$ wget https://github.com/phuongdo/tapo/releases/download/v1.1.2/bioapps.tar.gz
```

### 1.3 Source Code
```
$ cd $YOUR_INSTALL_DIR
$ git clone https://github.com/phuongdo/tapo.git
```

How to build
```
$ mvn clean package -DskipTests

```


## 2. TAPO How to Install 
The installation of TAPO is a simple process of extracting the archive and modifying configuration file. In BioApps, we included all necessary files and dependencies.

```
$ tar -zxvf bioapps.tar.gz
$ cd BioApps
```
Ensure that you have the correct 32bit/64bit version for your hardware and that the dependencies are made executable:

```
$ cd $YOUR_INSTALL_DIR
$ cd BioApps
$ chmod +x <YOUR DEPENCENCY FILE>
```

## 3.  TAPO Configuration 

In the HOME_DIR/BioApps/ tapo-v1.1.2/conf directory. The red lines are the vital parameters.
You must modify these parameters to fit your system OS. Others are not necessary.

Properties file for TAPO v1.1.2

```
#### PROGRAMME  ##### 

# <HOME_DIR> TAPO directory
# <XXX_TMP_DIR> TAPO temporary file directory
HOME_DIR= /home/mchateau/save/BioApps/tapo
DB_TMP_DIR = /home/mchateau/work/tmp
TMP_DIR = /home/mchateau/work/tmp

#### TAPO PARAMETERS   ##### 
#To change the score parameters by which a repeating units selected
#(BETWEEN 0.0 and 1.0)
TM_THRES = 0.5
#To change the score parameters by which a repeating units extented
#(BETWEEN 0.0 and 1.0)
TM_THRES_EXTEND = 0.5
#To change the score parameters to verify quality of one repeat
#(BETWEEN 0.0 and 1.0)
QA_THRES = 0.0
#To change the score parameters by which a repeating vectors selected
#(BETWEEN 0.0 and 1.0)
VECTOR_THRES = 0.5
# turn off QA score yielding performance
CALS_QA = on

#### TAPO DEPENDENCIES ##### 
#Dali program
# Holm L, Rosenström P (2010) Dali server: conservation mapping in 3D. Nucl. Acids Res. 38, W545-549.
# see more http://ekhidna.biocenter.helsinki.fi/dali_lite/start
# Deprecated ( no longer used, only for the developer mode)
DALI_DIR = /home/mchateau/Downloads/DaliLite_3.3
DALI_TMP_DIR = /home/mchateau/work/tmp
#DSSP program
# Deprecated ( no longer used, only for the developer mode)
DSSP_EXECUTABLE = /home/mchateau/save/BioApps/dsspcmbi/dssp/dsspcmbi
# The DSSP program was designed by Wolfgang Kabsch and Chris Sander to standardize secondary structure assignment.
# DSSP is a database of secondary structure assignments (and much more) for all protein entries in the Protein Data Bank (PDB).
# DSSP is also the program that calculates DSSP entries from PDB entries.
# see more http://swift.cmbi.ru.nl/gv/dssp/
DSSP_EXECUTABLE_v2 = /home/mchateau/save/BioApps/dsspcmbi/dssp-2.0
DSSP_TMP_DIR = /home/mchateau/work/tmp

# Structure Alphabet program
# based on 3D blast, save the file in the type of
# proCode.sa . eg 1BPO.sa
# J.-M. Yang and C.-H. Tung "Protein structure database search and evolutionary classification,"
# Nucleic Acids Research, vol. 34, pp. 3646-3659, 2006
# see more http://3d-blast.life.nctu.edu.tw/
# Deprecated ( no longer used, only for the developer mode)
SADB_EXECUTABLE = /home/mchateau/BioApps/3d-blast/3d-blast
SADB_TMP_DIR = /home/mchateau/work/tmp

#CD-HIT
# Weizhong Li and Adam Godzik. Cd-hit: a fast program for clustering and comparing large sets of protein
# or nucleotide sequences. Bioinformatics, 2006(22): 1658-1659.
# see more http://weizhongli-lab.org/cd-hit/
CDHIT_EXECUTABLE =  /home/mchateau/save/BioApps/cd-hit/cd-hit
CDHIT_TMP_DIR = /home/mchateau/work/tmp

#   MUSTANG (v3.2.1): A MUltiple STuructural AligNment alGorithm.
# see more http://www.ncbi.nlm.nih.gov/pubmed/16736488
MUSTANG_EXECUTABLE = /home/mchateau/save/BioApps/MUSTANG_v3.2.2/bin/mustang-3.2.2
MUSTANG_TMP_DIR = /home/mchateau/work/tmp
#   3DCOMB : Multiple Protein Structure Alignment
# downloaded at  http://ttic.uchicago.edu/~jinbo/software.htm.
# see more http://www.ncbi.nlm.nih.gov/pubmed/21791532
DCOMB_EXECUTABLE = /home/mchateau/save/BioApps/3DCOMB/3DCOMB_linux
DCOMB_TMP_DIR = /home/mchateau/work/tmp

#CMO : Contact Map Optimization
# Deprecated ( no longer used, only for the developer mode)
CMO_DIR = /home/mchateau/ABC
# CM_OUTPUT = HOME_DIR &/ CM_OUTPUT_DIR
CM_OUTPUT_DIR = /home/mchateau/work/tmp

# Deprecated ( no longer used, only for the developer mode)
JMOL_WEB_DATA= C:\\Program Files (x86)\\Zend\\Apache2\\htdocs\\3DRepeat\\data
#JMOL_WEB_DATA= /home/mchateau/save/public_html/3DRepeat/data

# Deprecated ( no longer used, only for the developer mode)
#TM-Align structure alignment program (Zhang and Skolnick 2003)
#http://zhanglab.ccmb.med.umich.edu/TM-align/
TMALIGN_EXECUTABLE = /home/pdoviet/BioApps/TMalignc/TMalign.exe

#### LARGE ANALYSIS  ##### 

####---------------NOTE----------------------

# By-default Leave the folder empty, TAPO can automatically download and install most of the data files that it needs.
# Those downloads will happen only once. Future requests for the data file will re-use the local copy.
####--------------END NOTE----------------------


# The directory <PDB_LOCAL> is the entry directory for the PDB database,
# If you already have a local PDB installation
# see more http://www.rcsb.org/pdb/static.do?p=download/ftp/index.html
PDB_LOCAL = /home/mchateau/work/bank/local_pdb/
# The directory <FASTA_LOCAL> is the entry directory for the FASTA database
# By-default: Leave this folder empty
FASTA_LOCAL = /home/mchateau/work/bank/fasta
# The directory <DSSP_LOCAL> is the entry directory for the DSSP database
# http://swift.cmbi.ru.nl/gv/dssp/
# By-default: Leave this folder empty
DSSP_LOCAL = /home/mchateau/work/bank/dssp
# The directory <SADB_LOCAL> is the entry directory for the SADB database
# see more http://3d-blast.life.nctu.edu.tw/
# By-default: Leave this folder empty
SADB_LOCAL = /home/mchateau/work/bank/sadb
# The directory <SADB_LOCAL> is the entry directory for the SCOPe(extended) database
# see more http://3d-blast.life.nctu.edu.tw/
# By-default: Leave this folder empty
SCOPE_LOCAL = /home/mchateau/work/bank/scopedb
# The directory <CATH_LOCAL> is the entry directory for the CATH database
# see more http://cathdb.info
# download lastest release of CATH (http://release.cathdb.info/)
# By-default: Leave this folder empty
CATH_LOCAL = /home/mchateau/work/bank/cathdb
# temporary files
CLUSTER_WORKING_DIR = /home/mchateau/work

#### GENO-CLUSTER OPTIONS ##### 
#### CLUSTER PARALLEL SMP ##### 
# the number of core
PARALLEL_SMP = 6
USER_NAME = mchateau
# LOCAL OR CLUSTER
CLUSTER_MODE = CLUSTER
#JAVA_OPT = java -version:1.7 -Xmx2042m -Djava.awt.headless=true -Dlog4j.skipJansi=true -cp "target/classes;dependencies/bugs/vecmath-1.3.1.jar;dependencies/*"
JAVA_OPT = java -Xmx2042m -Djava.awt.headless=true -Dlog4j.skipJansi=true -cp  "target/*:dependencies/bugs/vecmath-1.3.1.jar:dependencies/*"
### WEBSERVICE MODE #### 

# port TAPO server is using
WEB_PORT:9999

###########RUN BASH SCRIPT####################
#	     FOR WINDOW USER
##############################################

# put CYGWIN_PATH = empty if you run on Linux  environment
#by getting your current path in your cygwin application, type "cygpath -w" in windown
#for further information, type cygpath -help
#CYGWIN_PATH =
CYGWIN_PATH =
BASH_EXECUTABLE = /bin/bash
```

## 4 Running TAPO


### 4.1 Local Run

Available options as follow:
* -c        protein chain id eg. A
* -help     help
* -nCore    number of cores you want to use. Default is 1 core
* -o        ouput file.
* -p        protein code eg. 1bpo
* -f        file mode ( in PDF formmat)


Run PDF code
```
$ cd $YOUR_INSTALL_DIR
$ cd BioApps/tapo-v1.1.2
$ sh bin/run apps.TaPo "-p 2gsc -c A -nCore 1 -o /home/mchateau/work/2gscA.o"
```
Run PDF file mode
```
$ cd $YOUR_INSTALL_DIR
$ cd BioApps/tapo-v1.1.2
$ sh bin/run apps.TaPo "-f /home/phuongdv/our_dir/afile.pdb -p 2gsc -c A -nCore 1 -o /home/mchateau/work/2gscA.o"
```

### 4.2 Large-Scale Analysis

You must segment your file. This command show how to segment your file into a certain number of proteins per file. Available options as follow:

* <input/pdbList.in> : input file
* <10> : the number of protein structures per file.


```
$ cd $YOUR_INSTALL_DIR
$ cd BioApps/ tapo-v1.1.2
$ sh bin/run apps.GenerateJobInput "input/pdbList.in 10"
```

Before run
*  <input/pdbList.in> : input file


```
$ cd $YOUR_INSTALL_DIR
$ cd BioApps/ tapo-v1.1.2
$ sh bin/run apps.UpdatePDB "input/pdbList.in"
```

How to submit your job to GenoTul.

```
$ cd $YOUR_INSTALL_DIR
$ cd BioApps/ tapo-v1.1.2
$ sh bin/submit_jobarrs
```

How to start TAPO webserver

```
$ cd $YOUR_INSTALL_DIR
$ cd BioApps/tapo-v1.1.2
$ screen -R tapo
$ sh bin/start-webserver
```




4.3 TAPO Output

```

*  COLUMNS         DATA TYPE           CONTENTS
* ---------------------------------------------- 
*  $1              String              PDB_ID
*  $2              String              Finder Name
*  $3              String              Cluster name(selected mark is the best prediction)
*  $4              Integer             Number of repeat units
*  $5              Real                Avg. length of repeat
*  $6              String              TRs regions
*  $7              String              Repeat Units regions
*  $8              Real                QA score
*  $9              Real                R Score
*/
```


```

>1v3w_A|TM-Score=0.736;Psim-Score=0.803;CE-Score=0.593;V-Score=0.637;L-Score=0.332;CA-Score=0.000;S-Score=0.000;RA-Score=0.406|SVM-Score=1.000
1v3w_A	TRUST-23A	cl1_selected	11	6.000	36-131	36-41;51-56;57-62;73-78;79-84;90-95;102-107;108-113;114-119;120-125;126-131	0.895	0.596
1v3w_A	CESymm	cl1	7	19.000	3-136	3-24;25-42;43-63;64-85;86-102;104-120;121-136	0.446	0.264
1v3w_A	RMSD	cl1	4	30.000	16-135	16-47;48-87;88-110;111-135	0.564	0.171
```

Convert TAPO format to sequence residue number format.

Available options as follow:

* <fileTapoArrayResNumber.in>: input file
* <fileTapoSeqResNumber.ResNumber.out>: output file.

```
$ cd $YOUR_INSTALL_DIR
$ cd BioApps/tapo-v1.1.2
$ sh bin/run  org.cnrs.crbm.nextgen.ConvertTaPoFormatToSeqResNumber " fileTapoArrayResNumber.in fileTapoSeqResNumber.ResNumber.out"
```
Convert TAPO format to multiple structure alignment.

Available options as follow:

* < fileInTapoFormat.in>: input file
* < fileOut.msa>: output file.

```
$ cd $YOUR_INSTALL_DIR
$ cd BioApps/tapo-v1.1.2
$ sh bin/run  utils.ExportMsa "fileInTapoFormat.in fileOut.msa"
```

## FAQ?

Can we run TAPO on Window OS?
Maybe. You need to install Cygwin (https://www.cygwin.com/) and rebuild all dependencies.  
