# MSG

This code is the implemention of MSG methods published as following:

Y. Ji, Q. Wei, R. Chen, Q. Wang, R. Tao, B. Li, Integration of multidimensional splicing data and GWAS summary statistics for risk gene discovery. PLoS Genet. 18, e1009814 (2022).


# Hassle-free way to run the code: singularity

singularity > 3.8

```{bash}
conda create -p /path_to_env/msg

conda activate /path_to_env/msg

conda install -c conda-forge singularity

```


## step to run it

1. Download the data.tar from the [link](https://www.dropbox.com/scl/fo/vly3z0mxawa5x0v9k4elh/ADUtksnvW0Iyo71twsQDNGc?rlkey=vi69osfby7ipnjru9wva1tte3&dl=0):

	- gene.500k.id: The file include the coordinates of upstreaming and downstreaming 500k of a gene

	- hg38.vcf: The rsid for each snp in 1k genome.

	- GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz: The genotype data for all the individuals in GTEx. Used for generating x matrix.

	- samples.filtered.used: The samples used for running the MSG methods. Need to be generated for each tissue by users. Users can use all samples in the genotype file. It will be automatically filtered during the excution.

	- LDREF: The linkage disequilibrim file for each chromosome.

	- all_splicing_filled_expand: The splicing expression in GTEx. Used for generating y matrix.
	
	- sumstats: The summary statistics of a specific disease. Users need to format their own file by this way. 

3. Git the code to local direcoty:

	- scripts: include all the scripts for the MTSG method.
	
	- runRmd.R: A script to run rmd file through command line with options can be specified.

	- perGene.process.x.rmd: generate x matrix for each gene
	
	- perGene.process.y.rmd: generate y matrix for each gene

	- generate_db_and_cov.R: The R script to generate cov file for each gene.

	- gtex_comp_MSG_ACAT_GBJ_120522.R: The R script run the MSG method.
	
	- impute_missing_0120.R: Script to impute splcing expressions
	
	- sda_static_linux: Binary to excute tensor decomposition
	
	- others: required scripts to run MTSG
	

4. Put the data and the code in the same dir. If not, please specify absolute path. Run for gene ENSG00000000457:


### step-by-step instruction 
### using conda to install enviroment


```{bash, label = "", linewidth = 85, eval=opt$eval}
## 1. install conda
## install mambaforge - a variant of conda
ENV=/Path_to_Conda/miniconda3

mkdir download/
cd download/
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
mkdir -p $ENV
bash Mambaforge-Linux-x86_64.sh -b -p $ENV

export PATH=$PATH:$ENV/bin

## 2. install R, perl, cpanm, parallel and other r packages.
### activate the env
source $ENV/bin/activate
### install
mamba install -y -c conda-forge r-base perl perl-app-cpanminus parallel r-data.table r-optparse r-magrittr r-reshape2 r-pma r-mvtnorm  r-rcpp r-devtools r-stringi r-rcppeigen r-gbj r-glmnet r-this.path r-dplyr r-devtools

## 3. link compilers
parallel echo ln -s $ENV/bin/x86_64-conda_cos6-linux-gnu-cc $ENV/bin/{=s/.*-//=} ::: $ENV/bin/x86_64-conda_cos6-linux-gnu-cc $ENV/bin/x86_64-conda_cos6-linux-gnu-gcc $ENV/bin/x86_64-conda_cos6-linux-gnu-g++ | bash


## 4. install perl modules
mamba install -y -c bioconda perl-dbi perl-ipc-run perl-mce perl-module-build perl-string-random perl-mce-shared
cpanm Switch
cpanm List::Uniq
cpanm Array/Utils.pm
cpanm MCE/Hobo.pm
cpanm --force Env/Modify.pm
wget https://github.com/crotoc/Bundle-Wrapper/raw/master/Bundle-Wrapper-0.03.tar.gz
cpanm Bundle-Wrapper-0.03.tar.gz

## 5. install tabix
mamba install -y -c bioconda tabix

## 6. install r packages need by MSG
R -e 'library(devtools);install_github("gabraham/plink2R", subdir="plink2R")'

mamba install -y r-bigmemory r-bigmemory.sri  r-bh r-uuid
R -e 'library(devtools);install_github("XingjieShi/TisCoMM")'

mamba install -y r-base r-argparse r-optparse r-ggplot2 r-plyr r-dplyr r-magrittr r-reshape2 r-data.table r-gtools r-rcolorbrewer r-extrafont r-gridextra r-rmarkdown r-knitr r-formatr r-devtools r-readr 
R -e 'library(devtools);install_github("crotoc/myutils",dep=FALSE)'

## activate the ENV
conda activate $ENV
```


### 

```{bash}
## clone scripts
git clone https://github.com/sallymeeyan/MTSG.git

cd MTSG

## !!!! Download the MSG.sif and data dir to MSG_code dir
wget https://www.dropbox.com/scl/fo/vly3z0mxawa5x0v9k4elh/ADUtksnvW0Iyo71twsQDNGc?rlkey=vi69osfby7ipnjru9wva1tte3 -O MSG.zip
unzip MSG.zip

## activate the singularity env
conda activate /path_to_env/msg

## run one gene to test
## generate x matrix for ENSG00000000457
## output: test/xmatrix/ENSG00000000457.select.rsid.vcf.ac.x_all
Rscript scripts/runRmd.R --rmd scripts/perGene.process.x.rmd  --gene ENSG00000000457 --dir_data data/ --dir_out test/xmatrix/ --output ENSG00000000457  --eval TRUE --ext html

## generate cov file for ENSG00000000457
## output: test/cov/ENSG00000000457.select.rsid.vcf.ac.cov.RData
Rscript scripts/generate_db_and_cov.R --input test/xmatrix/ENSG00000000457.select.rsid.vcf.ac.x_all --dir_cov test/cov/ --ref_ld_chr data/LDREF/1000G.EUR.


## generate y matrix for ENSG00000000457
## output: test/ymatrix/ENSG00000000457.all.final.matrix.decomp
Rscript scripts/runRmd.R --rmd scripts/perGene.process.y.rmd  --gene ENSG00000000457 --dir_data data/ --dir_out test/ymatrix --output ENSG00000000457 --eval TRUE --ext html


## do MTSG
Rscript scripts/gtex_comp_MSG_ACAT_GBJ_120522.R --x test/xmatrix/ENSG00000000457.select.rsid.vcf.ac.x_all --y test/ymatrix/ENSG00000000457.all.final.matrix.decomp  --model_training --save_model --cov test/cov/ENSG00000000457.select.rsid.vcf.ac.cov.RData --sumstats data/sumstats/clozukscz.sumstats --dir_out test/MTSG/ --verbose TRUE


```







