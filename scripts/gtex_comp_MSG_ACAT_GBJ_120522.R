

# script:
# real data: /scratch/cgg/jiy1/CCA/MSG_ACAT_GBJ_111321
# simulation:

# goal: take x_all, y_all from gtex files, compute MSG results, combine using haky's, MSG, ACAT
# DATE: 11/13/21

# splicing only (qiang's output yall include expression + splicing  )
# permute y rows (so person id mixed up)

# add Bonf version
# add polled test of all isoforms (chisq)

# need a new list with more than 2 splicing events
# awk '$11 >=2' /fs0/jiy1/CCA/gtex/try/geneid.nomiss.output > geneid.splt2.output
# /fs0/jiy1/CCA/gtex/try/geneid.splt2.output



# gen cov
# awk '{print "Rscript /home/jiy1/CCA_scripts/TisCoMM/gtex/gtex_comp_120120.R --x" , "/scratch/cgg/jiy1/CCA/ba9_data/ba9/"$1"."$8".x_all",  "--y","/scratch/cgg/jiy1/CCA/ba9_data/ba9/"$1"."$8".y_all --ref_ld_chr  /scratch/cgg/jiy1/CCA/BioVU_5000REF/r5000/chr  --sumstats /scratch/cgg/jiy1/ldsc/summstats/MDD/MDD2018_2.sumstats.gz --generate_db_and_cov"}'  /scratch/cgg/jiy1/CCA/geneid.nomiss.output  > myjob.120120


# model training


# awk '{print "Rscript /home/jiy1/CCA_scripts/TisCoMM/gtex/gtex_comp_MSG_ACAT_GBJ_111321.R --x" , "/scratch/cgg/jiy1/CCA/ba9_data/ba9/"$1"."$8".x_all",  "--y","/scratch/cgg/jiy1/CCA/ba9_data/ba9/"$1"."$8".y_all   --sumstats /scratch/cgg/jiy1/ldsc/summstats/MDD/MDD2018_2.sumstats.gz --model_training"}'  /scratch/cgg/jiy1/CCA/geneid.nomiss.output  > myjob.train.111321


# put every 100 commands in a job
# cov jobs
# for c in {0..11556..100};do echo "less /scratch/cgg/jiy1/CCA/mdd_biovuref_120220/myjob.120120 |head -$c |tail -100 | /home/weiq1/script/bingshan/local-run -n 1 -sleep 0.1";done > cov.tmpjob
# sed -i '1d' cov.tmpjob

# echo 'less myjob.120120 | tail -56 | /home/weiq1/script/bingshan/local-run -n 1 -sleep 0.1' >> cov.tmpjob

# less cov.tmpjob| perl /home/weiq1/script/bingshan/ssub_array2 -pd `pwd`  -pm 16000 -ph 8 -email ying.ji@vanderbilt.edu -e

# model jobs
# for c in {0..11556..100};do echo "less /scratch/cgg/jiy1/CCA/mdd_biovuref_120220/myjob.train.120120 |head -$c |tail -100 | /home/weiq1/script/bingshan/local-run -n 1 -sleep 0.1";done > train.tmpjob
# sed -i '1d' train.tmpjob

# echo 'less myjob.train.120120 | tail -56 | /home/weiq1/script/bingshan/local-run -n 1 -sleep 0.1' >> train.tmpjob

# less train.tmpjob| perl /home/weiq1/script/bingshan/ssub_array2 -pd `pwd`  -pm 8000 -ph 40 -email ying.ji@vanderbilt.edu -e



###################
# example
# dir: /scratch/cgg/jiy1/CCA/sCCA_GBJ_110621
####################

# Rscript /home/jiy1/CCA_scripts/TisCoMM/gtex/gtex_comp_MSG_ACAT_GBJ_111321.R --x /scratch/cgg/jiy1/CCA/ba9_data/ba9/ENSG00000137601.NEK1.x_all --y /scratch/cgg/jiy1/CCA/ba9_data/ba9/ENSG00000137601.NEK1.y_all    --model_training

# Rscript /home/jiy1/CCA_scripts/TisCoMM/gtex/gtex_comp_MSG_ACAT_GBJ_111321.R --x /scratch/cgg/jiy1/CCA/ba9_data/ba9/ENSG00000164967.RPP25L.x_all --y /scratch/cgg/jiy1/CCA/ba9_data/ba9/ENSG00000164967.RPP25L.y_all    --model_training

# Rscript /home/jiy1/CCA_scripts/TisCoMM/gtex/gtex_comp_MSG_ACAT_GBJ_111321.R --x /scratch/cgg/jiy1/CCA/ba9_data/ba9/ENSG00000132323.ILKAP.x_all --y /scratch/cgg/jiy1/CCA/ba9_data/ba9/ENSG00000132323.ILKAP.y_all    --model_training




start_time <- Sys.time()

# use R 3.6
#.libPaths("/home/jiy1/R/rlib-3.6.0")
library(mvtnorm)
library(GBJ)
library(glmnet)
library(Rcpp)
library(TisCoMM)
library(optparse)
library(plink2R)
library(magrittr)
library(this.path)
library(data.table)
suppressMessages(require(dplyr))
suppressMessages(require(magrittr))


## mydir <- "/fs0/chenr6/git_projects/MSG_code/scripts"
##mydir  <- "/nobackup/cgg/yany14/git_project/MTSG/scripts"
mydir <- this.dir()
##print(mydir)

source(paste(mydir,"functions.r",sep="/"))
#source("/fs0/jiy1/CCA/TisCoMM/simulation/functions.r")  
# need to compile c code first: R CMD SHLIB optim.c
# and add absolute path to c code to load
source(paste(mydir,"multiEstB.R",sep="/"))
source(paste(mydir,"allele_qc.r",sep = "/"))

sourceCpp(paste(mydir,"PXem_ss.cpp",sep = "/"))

sourceCpp(paste(mydir,"PXem.cpp",sep = "/"))

# a file with cca script
source(paste(mydir,"cca_092930.r",sep = "/"))



option_list <- list(
    ## Type : logical, integer, souble complex, character
	make_option(c("--x"), action="store", #dimension
                type='character', help="x matrix genotype"),
     make_option(c("--y"),action="store", #z share with b2
                type='character', help="y phenotype"),
     make_option(c("--sumstats"), action="store", default="/scratch/cgg/jiy1/ldsc_2020/ldsc/sumstats/clozukscz.sumstats.gz", type='character',help="Path to summary statistics (must have SNP and Z column headers) [required]"),
    make_option(c("--ref_ld_chr"), action="store", 
    #default="/data1/jiy1/LDREF/1000G.EUR.",
    default=NA,
    type='character',help="Prefix to reference LD files in binary PLINK format by chromosome [required]"),
    make_option(c("--cov_dir"), action="store", 
                default="/scratch/cgg/jiy1/CCA/BioVU_5000REF/cov",
                type='character',help="Prefix to reference LD files in RData format by gene id [required]"),
    make_option(c("--cov"), action="store",
                default="x",
                type="character",
                help="cov file [default: %default]"),
    make_option(c("--dir_out"), action="store",
                default="1",
                type="character",
                help=" [default: %default]"),
    make_option("--model_training", action="store_true", default=FALSE),
    make_option("--generate_db_and_cov", action="store_true", default=FALSE),
    make_option("--save_model",action="store_true",default=FALSE),
    make_option("--bonf",action="store_true",default=FALSE),
    make_option(c("--test"), action="store",
                default=FALSE,
                type="logical",
                help=" [default: %default]"),
    make_option(c("--verbose"), action="store",
                default=FALSE,
                type="logical",
                help=" [default: %default]")
)



######################################################
#######        parse arguments ########################
#######################################################


opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

if(opt$test){
    save.image("test.rda")
    q("no")
}


if(opt$model_training & opt$generate_db_and_cov){
  stop('You can only do model training or generating db and cov at a time')
}
if(!(opt$model_training | opt$generate_db_and_cov | opt$bonf )){
  stop('Please use "--model_training" for model training or use "--generate_db_and_cov" for db and cov generating or use "--bonf" for bonferroni correction' )
}



# x.file="/data1/jiy1/gtex/frontal_cortex/ba9/ENSG00000100197.CYP2D6.x_all"
# y.file="/data1/jiy1/gtex/frontal_cortex/ba9/ENSG00000100197.CYP2D6.y_all"


# x.file = "/data1/jiy1/gtex/frontal_cortex/ba9/ENSG00000164967.RPP25L.x_all"

# y.file="/data1/jiy1/gtex/frontal_cortex/ba9/ENSG00000164967.RPP25L.y_all"

## opt$x <- "/scratch/cgg/yany14/MSG/test_out1021/ENSG00000204396.select.rsid.vcf.ac.x_all"
## opt$y <- "/scratch/cgg/yany14/MSG/test_out1021/all.final.matrix.decomp"     
## opt$ref_ld_chr <- "/scratch/cgg/yany14/MSG/data/LDREF/1000G.EUR."
## opt$cov_dir <- "/scratch/cgg/yany14/MSG/test_out1021/"
if(opt$verbose)
    cat("## Reading X")

if(!file.exists(opt$x)){
    cat("## X file not exist",opt$x,"\n")
    q(save = "no")
}

x.file = opt$x
print(x.file)

######################
# generate output loc
######################
if(opt$verbose){
    cat("## generate output loc\n")
}
## Provide the dir name(i.e sub dir) that you want to create under main dir:

output_dir <- paste0(opt$dir_out,"/results/")

if (!dir.exists(opt$dir_out)){
    dir.create(opt$dir_out,recursive = TRUE,showWarnings = F)
}

if (!dir.exists(output_dir)){
    dir.create(output_dir,recursive = TRUE)
} else {
    print("results Dir already exists!")
}

dir.create(paste0(output_dir,"/all"),recursive = TRUE)

if (opt$save_model){
    dir.create(paste0(output_dir,"/models"),recursive = TRUE)
}


######################################################################################
##loading input and filter
######################################################################################
####x_all
if(opt$verbose){
    cat("## loading x input\n")
}

if (is.null(x.file) || x.file == "None" || x.file == "NULL") {
    stop("x file is required");
} else {
    x <- fread(x.file);
}

####################################################################
#######		get some variables
#####################################################################
gid <- gsub("\\.x_all.*","",basename(x.file),perl = T)
chr <- unique(x[[1]])[1]

####################################################################
#######		 make cov files
#####################################################################
if(opt$verbose){
    cat("## make cov files\n")
}

if(opt$generate_db_and_cov){
    print(paste0('INFO generating db and cov files: ',gid))
    ## read in ref panel + sum stats
    ## note: don't type .bim
    ## genos = read_plink("/data1/jiy1/LDREF/1000G.EUR.22",impute="avg")
    ## Load in reference data: 1kg
    ##  print(paste(opt$ref_ld_chr,chr,sep=''))
    genos = read_plink(paste(opt$ref_ld_chr,chr,sep=''),impute="avg") 
    ## subset the common snps of x matrix and genos 
    comm_snps <- intersect(x[[2]],genos$bim[,2]) 
    cur.genos = scale(genos$bed[,match(comm_snps,colnames(genos$bed))])
    LD = var(cur.genos)
    if(!dir.exists(paste0(opt$cov_dir,'/'))){
        dir.create(paste0(opt$cov_dir,'/'))
    }
    outname1 <- paste0(opt$cov_dir,'/',paste(gid,"cov.RData",sep="."))
    print(outname1)
    print(str(genos))
    bim = genos$bim[match(comm_snps,genos$bim[,2]),]
    save(bim,LD,file=outname1)
}

if (!opt$bonf && !opt$model_training ){
    cat("## Don't do any model training or testing")
    q(save = "no")
}	


############
##y_all
if(opt$verbose)
    cat("## Reading Y")
if(!file.exists(opt$y)){
    cat("## Y file not exist",opt$y,"\n")
    q(save = "no")
}
y.file = opt$y
print(y.file)
#output = opt$output

if(opt$verbose){
    cat("## loading y input\n")
}
if (is.null(y.file) || y.file == "None" || y.file == "NULL") {
	stop("y file is required");
} else {
    ##y <- read.table(y.file)
    y <- fread(y.file)
    ## delete expression 
    y <- y[get(names(y)[1]) != "Expression",]
    ## shuffle splicing
    set.seed(42)
    y <- y[sample(nrow(y)),]

    ## identify columns with constant value
    y <- y[apply(y[,-1, with =F],1,sd)!=0,]

    if( length(names(x)[-1:-7]) != length(names(y)[-1]) || !all(names(x)[-1:-7] == names(y)[-1])){
        cat("## Y matrix and X matrix must contain the same samples with the same order\n")
        cat("## filtering...\n")
        overlap_samples <- intersect(names(x)[-1:-7],names(y)[-1])
        x <- x[,c(names(x)[1:7],overlap_samples),with = F]
        y <- y[,c(names(y)[1],overlap_samples),with = F]
    }

    if(all(names(x)[-1:-7] == names(y)[-1])){
        y <- data.table(y[[1]],y[,names(x)[-1:-7],with = F])
    }else{
        cat("## Y matrix and X matrix still don't have same samples with the same order..Quit\n")
        q(save = "no")
    }
    if(nrow(y)<2){
        cat("## Y matrix needs at least two features\n")
        q(save = "no")
    }
    ## transpose y to individual by splicing
    yname <- y[[1]]
    y <- data.table(t(y[,-1,with=F]))
    names(y) <- yname

    ## row scale: so each omic become mean 0 sd 1
    y_all <- round(scale(y),digits=8)
}

#n= length(y_all)

if(opt$verbose){
    cat("## checking x_all row number\n")
}

if(nrow(x)<1){
    cat("## X matrix needs at least one line\n")
    q(save = "no")
}



####################################################################
#######		load saved cov rdata
#####################################################################
if(!opt$model_training){
    cat("## model training is FALSE\n")
    q(save = "no")
}

## load in LD file
if(file.exists(opt$cov)){
    load(opt$cov)
}else{
    ## get all files in the dir
    cov_files <- list.files(opt$cov_dir)

    ## get ENSG id
    gid_ENSG <- gsub("(ENSG\\d+).*","\\1",gid,perl =T)

    ## Sys.glob use the widecard * in the file name
    if(!file.exists(Sys.glob(paste0(opt$cov_dir,'/',paste(gid_ENSG,"*cov.RData",sep=".")), dirmark = FALSE))){
        cat("## cov file not exist", Sys.glob(paste0(opt$cov_dir,'/',paste(gid_ENSG,"*cov.RData",sep=".")), dirmark = FALSE),"\n")
        q(save = "no")
    }

    load(Sys.glob(paste0(opt$cov_dir,'/',paste(gid_ENSG,"*cov.RData",sep=".")), dirmark = FALSE))
}
## Clean LD rows and cols with all NA
vectorIsAllna <- function(x){
    if(length(unique(x))==1 && is.na(unique(x))){
        return(TRUE)
    }else{
        return(FALSE)
    }
}

LD <- LD[!apply(LD,1,function(x){vectorIsAllna(x)}),!apply(LD,2,function(x){vectorIsAllna(x)})]

## replace NA cells of LD to 0
LD[is.na(LD)] <- 0

####################################################################
#######		sum stats
####################################################################
if(opt$verbose){
    cat("## read sum stats\n")
}
## because in this section, there are a lot of filtering steps. To make sure the final matrix is conformable, we need to update all the input matrices after each filter step.

s.out2<-data.frame(GID=gid,CHR=chr,stringsAsFactors=FALSE)

## Load in summary stats
sumstat = fread(opt$sumstats)

## get comm_snp of x_all, LD, sumstat
comm_snps <- intersect(intersect(intersect(x[[2]],bim[,2]),colnames(LD)),sumstat$SNP)

## update sumstat and bim
sumstat <- sumstat[match(comm_snps,SNP),]
bim <- bim[match(comm_snps,bim[,2]),]

## QC / allele-flip the input and output
qc = allele.qc( sumstat$A1 , sumstat$A2 , bim[,5] , bim[,6] )

## Flip Z-scores for mismatching alleles in sumstats
sumstat$Z[ qc$flip ] = -1 * sumstat$Z[ qc$flip ]
sumstat$A1[ qc$flip ] = bim[qc$flip,5]
sumstat$A2[ qc$flip ] = bim[qc$flip,6]

## FILTER: Remove strand ambiguous SNPs (if any)
if ( sum(!qc$keep) > 0 ) {
    bim = bim[qc$keep,]
    sumstat = sumstat[qc$keep,]
}

## UPDATE comm_snps after filtering bim and sumstat
comm_snps <- intersect(intersect(intersect(x[[2]],bim[,2]),colnames(LD)),sumstat$SNP)
x <- x[match(comm_snps,get(names(x)[2])),]
LD <- LD[match(comm_snps,row.names(LD)),match(comm_snps,colnames(LD))]
sumstat <- sumstat[match(comm_snps,SNP),]
bim <- bim[match(comm_snps,bim[,2]),]


####################################################################
#######		filtering X with columns of zero variation
#####################################################################
x <- x[apply(x[,-1:-7,with = F], 1, sd)!=0,]

if(nrow(x)<2 ){
    cat("## X matrix needs at least two columns\n")
    q(save = "no")
}

## UPDATE comm_snps after filtering bim and sumstat
comm_snps <- intersect(intersect(intersect(x[[2]],bim[,2]),colnames(LD)),sumstat$SNP)
x <- x[match(comm_snps,get(names(x)[2])),]
LD <- LD[match(comm_snps,row.names(LD)),match(comm_snps,colnames(LD))]
sumstat <- sumstat[match(comm_snps,SNP),]
bim <- bim[match(comm_snps,bim[,2]),]

comm_snps %>% length
x %>% dim
LD %>% dim
bim %>% dim

## finalize all variable used bellow
## Match up the LDREF SNPs and the summary stats
x_all <- t(x[,-1:-7, with = F])
colnames(x_all) <- x[[2]]
cur.Z <- sumstat$Z
cur.miss <- is.na(cur.Z)
cur.LD <- LD
snps <- x[,1:7,with = F]

####################################################################
#######		single splice enet and bonferroni
#####################################################################
## do Bonferroni
## mcca single splice level
if (opt$bonf){	
    ms.out<-data.frame(GID=character(0),CHR=character(0),isoform=character(0),NWGT=numeric(0),TWAS.Z=numeric(0),TWAS.P=numeric(0),stringsAsFactors=FALSE)
    ##s.out <- matrix(0,ncol=4,nrow=ncol(s.wgt.matrix)) 
    for (w in 1:ncol(wgt.matrix)){ # only those not dropped tissues
        ms.out[w,"model"] = "enet"
        ms.out[w,"GID"] = gid
        ms.out[w,"CHR"] = chr
        ms.out[w,"isoform"] = w
        ms.out[w,"NWGT"] = sum(wgt.matrix[,w]!=0)
        cur.twasz = wgt.matrix[!cur.miss,w] %*% cur.Z[!cur.miss]
        cur.twasr2pred = wgt.matrix[!cur.miss,w] %*% cur.LD %*% wgt.matrix[!cur.miss,w]
        cur.twas = cur.twasz / sqrt(cur.twasr2pred)
        TWAS.Z = cur.twas
        TWAS.P = 2*(pnorm( abs(TWAS.Z) , lower.tail=F))
        print(TWAS.Z)
        print(TWAS.P)
        ms.out[w,"TWAS.Z"] = TWAS.Z
        ms.out[w,"TWAS.P"] = TWAS.P
    }
    print(ms.out)
    outname1 <- paste0("all/",paste(gid,"bysplice","enet",sep="."))
    print(outname1)
    write.table( ms.out, quote=F , row.names=F , sep='\t' , file=outname1 )
}



##########################
## cca multiple canonical vectors
############################
if(opt$verbose){
    cat(paste0('INFO model training: ',gid))
}

if(opt$model_training){
    ## if ncol(y_all) > ncol(x_all), CCA will throw error because the V dimension is the minimium of these two numbers
    K = min(c(ncol(y_all),ncol(x_all))) 
    mod = CCA(x_all,y_all,K=K,typex="standard",typez="standard")
    B_hat1 = mod$u
    print(dim(B_hat1))
    print(dim(snps))
                                        #rownames(B_hat1) = snps[,2]
    print(head(B_hat1))
    sel = which(apply(B_hat1, 2, sd)!= 0)

    B_hat1 <- as.matrix(B_hat1[!cur.miss, sel])
                                        # 333 333					
                                        # this should be test_stats1	
    B_hats = B_hat1
    eta2 = diag(sqrt( t(B_hats) %*% cur.LD %*% B_hats ))
    sig_x2 <- sqrt(diag(cur.LD))
                                        # str(sig_x2)
                                        #  Named num [1:333] 1 1 1 1 1 ...
                                        #  - attr(*, "names")= chr [1:333] "rs132788" "rs132793" "rs8779" "rs132806" ...
    Ztilde <- cur.Z[!cur.miss] 
                                        #Ztilde <- hatmu/hats
    Lambda2 <-  diag(sig_x2) %*% sweep(B_hats, 2, eta2, "/")
    Lambda_sub2 <- Lambda2[, which(eta2 != 0)]
                                        # 333 6
    test_stats2 <- t(Lambda_sub2) %*% Ztilde
    cur.Zs = cur.Z[!cur.miss]	
                                        # now get rhoGE2
                                        # Ge_impute1 <- cur.genos[,!cur.miss]%*%B_hat1
                                        # rhoGE1 <- cor(Ge_impute1)
                                        # rhoGE_svd1 <- svd(rhoGE1)

    rhoGE2 = matrix(0,nrow=ncol(B_hats),ncol=ncol(B_hats))
    for (i in 1:ncol(B_hats)){
        for (j in 1:ncol(B_hats)){
            nom = t(B_hats[,i])%*%cur.LD%*%B_hats[,j]
            ##print(nom)
            denom = sqrt((t(B_hats[,i])%*%cur.LD%*%B_hats[,i])*(t(B_hats[,j])%*%cur.LD%*%B_hats[,j]))
            rhoGE2[i,j] = nom/denom
        }
    }


    rhoGE_svd2 = svd(rhoGE2)
    cond = 30
    ind_top <- which(rhoGE_svd2$d[1]/rhoGE_svd2$d < cond)
    if(length(ind_top) == 0) ind_top <- 1
    u <- rhoGE_svd2$u
    v <- rhoGE_svd2$v
    us <- as.matrix(u[, ind_top])
    vs <- as.matrix(v[, ind_top])
    d <- rhoGE_svd2$d[ind_top]
    if(length(d) > 1) ds <- diag(1/d) else ds = matrix(1/d)
    rhoGE_ginv <- vs %*% ds %*% t(us)
    mcca.chisq2 <- c(t(test_stats2) %*% rhoGE_ginv %*% test_stats2)
    mcca.dof = length(ind_top)
    mcca.p <- pchisq(mcca.chisq2,length(ind_top),lower.tail=F) 

    s.out2[,"mcca.p"] = rep(mcca.p,nrow(s.out2))
    s.out2[,"mcca.chisq"] = rep(mcca.chisq2,nrow(s.out2))
    s.out2[,"mcca.dof"] = rep(mcca.dof,nrow(s.out2))
    print("MSG results")
    print(s.out2)
    
########### GBJ ######################
    
    Lambda2 <-  diag(sig_x2) %*% sweep(B_hats, 2, eta2, "/")
    Lambda_sub <- Lambda2[, which(eta2 != 0)]
    ## 333 6
    test_stats <- t(Lambda_sub) %*% Ztilde
    
    
    ## R: LD matrix
    lam = 0.95
    RR = cur.LD
    p = sum(!cur.miss)
    
    R = RR*lam + (1 - lam)*diag(p)

    
    cov_Z <- t(Lambda_sub) %*% R %*% Lambda_sub
    cov_Z <- round(cov_Z,digits = 10) ## To make sure it's symmetric
    tc.gbj <- GBJ(test_stats=test_stats, cor_mat=cov_Z)
    s.out2[,"MSG.gbj.p"] = rep(tc.gbj$GBJ_pvalue,nrow(s.out2))
    s.out2[,"MSG.gbj"] = rep(tc.gbj$GBJ,nrow(s.out2))
    s.out2[,"MSG.gbj.dof"] = rep(length(test_stats),nrow(s.out2))
    
    print("MSG - GBJ")
    print(s.out2)
    
######### acat ###########
    Pvals = rep(0,ncol(B_hats))
    
    print(paste("cur.z:",length(cur.Z))) # 420
    print(paste("cur.LD:",dim(cur.LD))) # 420
    print(paste("B_hats:",dim(B_hats)))
    
    
    
    for (j in 1:ncol(B_hats)){
        
        cur.twasz = B_hats[,j] %*% cur.Zs
        cur.twasr2pred = B_hats[,j] %*% cur.LD %*% B_hats[,j]
        if ( cur.twasr2pred > 0 ) {
            cur.twas = cur.twasz / sqrt(cur.twasr2pred)
            cur.p = 2*(pnorm( abs(cur.twas) , lower.tail=F))
            Pvals[j] <- cur.p
        } else{
            Pvals[j] <- NA
        }
    }
    
                                        # acat combine everything
    Pvalsnew <- Pvals[!is.na(Pvals)]
    Pvals <- Pvalsnew
    
    Weights<-rep(1/length(Pvals),length(Pvals))    
    
    is.small<-(Pvals<1e-16)
    if (sum(is.small)==0){
        cct.stat<-sum(Weights*tan((0.5-Pvals)*pi))
    }else{
        cct.stat<-sum((Weights[is.small]/Pvals[is.small])/pi)
        cct.stat<-cct.stat+sum(Weights[!is.small]*tan((0.5-Pvals[!is.small])*pi))
    }
#### check if the test statistic is very large.
    if (cct.stat>1e+15){
        pval<-(1/cct.stat)/pi
    }else{
        pval<-1-pcauchy(cct.stat)
    }
                                        #pvalue[i,9] = pval 
    
    


    s.out2[,"MSG.acat.p"] <- pval 

    s.out2[,"MSG.acat.dof"] <- length(Pvals)
    s.out2[,"n_spl"] <- ncol(y_all)
    
    outname1 <- paste0(output_dir,"all/",paste(gid,"results","MSG_GBJ_ACAT.txt",sep="."))
    print(outname1)
    write.table( s.out2 , quote=F , row.names=F , sep='\t' , file=outname1 )

                                        # ind_top <- which(mod$cors > cond)
                                        # 
                                        # 
                                        # if(length(ind_top) == 0) ind_top <- 1
                                        # 
                                        # #u <- rhoGE_svd$u
                                        # #v <- rhoGE_svd$v
                                        # #d <- rhoGE_svd$d[ind_top]
                                        # u <- mod$u
                                        # v <- mod$v
                                        # d <- mod$d
                                        # us <- as.matrix(u[, ind_top])
                                        # vs <- as.matrix(v[, ind_top])
                                        # 
                                        # #if(length(d) > 1) ds <- diag(1/d) else ds = matrix(1/d)
                                        # rhoGE_ginv <- vs %*% ds %*% t(us)
                                        # tc.multixcan.chisq <- c(t(test_stats) %*% rhoGE_ginv %*% test_stats)
                                        # tc.multixcan.dof = length(ind_top)
                                        # tc.multixcan.p <- pchisq(tc.multixcan.chisq,length(ind_top),lower.tail=F) 



##########################
                                        # cca
############################

                                        # simulation code
                                        # cca: from file /fs0/jiy1/CCA/TisCoMM/simulation/cca_092930.R
                                        # pvalue[i,4] <- cca(x1,y,x2,z)$pval



                                        # ccao = ccapa(x_all,y_all)
                                        # 
                                        # 	cca.wgt.matrix <- ccao[[1]]
                                        # 	vout <- ccao[[2]]
                                        # 	penaltyx <- ccao[[3]]
                                        # 	penaltyz <- ccao[[4]]
                                        # 	ccacor <- ccao[[5]]
                                        # 
                                        # 
                                        # 	rownames(cca.wgt.matrix) = snps[,2]
                                        # 	# cur.LD = 
                                        # 	#cur.LD = t(cur.genos[,!cur.miss]) %*% cur.genos[,!cur.miss] / (nrow(cur.genos[,!cur.miss])-1)
                                        # 	cur.twasz = cca.wgt.matrix[!cur.miss,1] %*% cur.Z[!cur.miss]
                                        # 	cur.twasr2pred = cca.wgt.matrix[!cur.miss,1] %*% cur.LD %*% cca.wgt.matrix[!cur.miss,1]
                                        # 
                                        # 	cca.NWGT = sum( cca.wgt.matrix[!cur.miss,1] != 0 )
                                        # 	cur.twas = cur.twasz / sqrt(cur.twasr2pred)
                                        # 	cca.TWAS.Z = cur.twas
                                        # 	cca.TWAS.P = 2*(pnorm( abs(cca.TWAS.Z) , lower.tail=F))
                                        # 
                                        # 
                                        # 
                                        # 	s.out2[,"cca.z"] = rep(cca.TWAS.Z,nrow(s.out2))
                                        # 	s.out2[,"cca.p"] = rep(cca.TWAS.P,nrow(s.out2))
                                        # 	s.out2[,"cca.v"] = paste(round(vout[,1],4),collapse=",")
                                        # 	s.out2[,"cca.ydim"] = sum(abs(vout[,1]) > 0.05)
                                        # 	s.out2[,"cca.cor"] = ccacor
                                        # 	s.out2[,"cca.penalty.x"] = penaltyx
                                        # 	s.out2[,"cca.penalty.y"] = penaltyz




    

    if (opt$save_model) {
        outname2 <- paste0(output_dir,"models/",paste(gid,"models","RData",sep="."))
        print(outname2)
        if(exists("wgt.matrix")){
            save(wgt.matrix, B_hats, rhoGE_svd2, test_stats2, mcca.dof, rhoGE_ginv, ind_top, mcca.chisq2, snps, file=outname2)
        } else{
            save(B_hats,rhoGE_svd2, test_stats2, mcca.dof, rhoGE_ginv, ind_top, mcca.chisq2, snps, file=outname2)
        }
    }



}



end_time <- Sys.time()
end_time - start_time


# test one: scz summstats
# Rscript /home/jiy1/CCA_scripts/TisCoMM/gtex/gtex_comp_savemod_020221.R --x /data1/jiy1/gtex/frontal_cortex/ba9/ENSG00000137601.NEK1.x_all --y /data1/jiy1/gtex/frontal_cortex/ba9/ENSG00000137601.NEK1.y_all    --model_training --save_model
