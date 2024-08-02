start_time <- Sys.time()

library(mvtnorm)
library(PMA)
library(data.table)
library(optparse)
library(this.path)
suppressMessages(require(dplyr))
suppressMessages(require(magrittr))


option_list <- list(
    ## Type : logical, integer, souble complex, character
	make_option(c("--x"), action="store", #dimension
                type='character', help="x matrix genotype"),
     make_option(c("--y"),action="store", #z share with b2
                type='character', help="y phenotype"),
    make_option(c("--cov_dir"), action="store", 
                default="/nobackup/cgg/yany14/MSG/dir_cov_1219",
                type='character',help="Prefix to reference LD files in RData format by gene id [required]"),
    make_option(c("--cov"), action="store",
                default="x",
                type="character",
                help="cov file [default: %default]"),
    make_option(c("--nsims", action="store",
                type='integer',help="number of simulations")),
    make_option(c("--dir_out"), action="store",
                default="1",
                type="character",
                help=" [default: %default]"),
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

nsims <- opt$nsims

######################################################################################
##loading input and filter
######################################################################################
####x_all
if(opt$verbose)
    cat("## Reading X")

if(!file.exists(opt$x)){
    cat("## X file not exist",opt$x,"\n")
    q(save = "no")
}

x.file = opt$x
print(x.file)


####### get some variables

gid <- gsub("\\.x_all.*","",basename(x.file),perl = T)


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


############
##cov file
if(opt$verbose)
    cat("## Reading cov")
if(!file.exists(opt$cov)){
    cat("## cov file not exist",opt$y,"\n")
    q(save = "no")
}
cov.file <- opt$cov
print(cov.file)


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



# x: snp*subj
#x.file <- fread("~/Box Sync/Vandy_forth_year/LiLab/MSG/MSG_public-main/example_data/ENSG00000100197.CYP2D6.toy.x_all")
#y.file <- fread("~/Box Sync/Vandy_forth_year/LiLab/MSG/MSG_public-main/example_data/ENSG00000100197.CYP2D6.y_all")

#x.file <- fread("/Users/yanyan/Box Sync/Vandy_forth_year/LiLab/Splicing/splicing/x_y_matrices/ENSG00000007255.select.rsid.vcf.ac.x_all")

#y.file <- fread("/Users/yanyan/Box Sync/Vandy_forth_year/LiLab/Splicing/splicing/x_y_matrices/ENSG00000007255.all.final.matrix.decomp")

# r object, load in 'bim' and 'LD'
#cov <- load("/Users/yanyan/Box Sync/Vandy_forth_year/LiLab/Splicing/splicing/x_y_matrices/ENSG00000007255.select.rsid.vcf.ac.cov.RData")

#if (is.null(x.file) || x.file == "None" || x.file == "NULL") {
if (is.null(x.file)) {
  stop("x file is required");
} else {
  #x <- read.table(x.file,as.is=T, head = F);
  x <- fread(x.file)
  x1<-x[,c(1,2,3,4,5,6)] #basic info: CHROM,rsid,QUAL,POS,REF,ALT 
  # remove CHROM,rsid,QUAL,POS,REF,ALT,INFO
  x2<-as.matrix(x[,c(-1,-2,-3,-4,-5,-6, -7)]) # snp*subjects
  x2 <- t(x2)	# subjects * snps
  rm(x)
  
  p <- ncol(x2) # number of SNPs
  
  x.na <- sapply(1:p, function(j) {any(is.na(x2[, j]))}) # if any SNPs is NA
  
  if (any(x.na)) {  # remove the SNP if any NA
    x2 <- x2[, !x.na]
    p <- ncol(x2)	  
    x1 <- x1[!x.na,]    # used to extract the LD, snp on the row
  }
}

x_all <- x2

#y_all
#if (is.null(y.file) || y.file == "None" || y.file == "NULL") {
if (is.null(y.file)) {
  stop("y file is required");
} else {
  y <- fread(y.file)
  #y <- read.table(y.file,as.is=T,head=F)
  #y3 <- y[!(y[,1] == "Expression"),] # delete expression 
  y2 <- copy(y)
  
  # remove the splicing index
  y2$Expression <- NULL 
  

  # permute rows of this 
  set.seed(42)
  rows <- sample(nrow(y2))
  y2 <- y2[rows,] # splicing * subjects
  rm(y)
  
  # transpose y2, 
  y2 <- as.matrix(t(y2)) # subjects * splicing
  
  # row scale: so each omic become mean 0 sd 1
  #	y_all <- lapply(1:nrow(y2), function(r) { y2 <- scale(as.matrix(y2[r,])) })
  y_all <- scale(y2)
  
  rm(y2)
  
}



x1 <- x_all # 581*372, genotype
y <- y_all # 581*7, splicing

#load cov file
load(cov.file)
cur.LD <- LD # 372*372

#library(PMA)
mod = CCA(x1,y,K=ncol(y), typex = "standard", typez = "standard")
B_hat1 = mod$u # weights for x1


#library(mvtnorm)
######## mcca function################
  mccasum <- function(B_hat1,cur.LD, n_p = 5000){
  # genotypes in GWAS
  x3 <- rmvnorm(n_p, mean=rep(0, p), sigma = cur.LD)
  x3mean <- matrix(rep(colMeans(x3), n_p), ncol=p, byrow = T)
  #x3sd <- matrix(rep(apply(x3, 2, sd), n_p), ncol=p, byrow = T)
  x3 <- (x3 - x3mean)#/x3sd/sqrt(n_p)
  
  z <- rnorm(n_p, 0, sqrt(3))
  
  # summary statisitcs for GWAS
  ## hatmu = matrix(0, p, 1)
  ## hats  = matrix(0, p, 1)
  ## hatp = matrix(0,p,1)
  hatz = matrix(0,p,1) # add z score 
  for (j in 1:p){
    fm <- lm(z ~ 1 + x3[, j]);
    ## hatmu[j] <- summary(fm)$coefficients[2,1]
    ## hats[j]  <- summary(fm)$coefficients[2,2]
    ## hatp[j] <- summary(fm)$coefficients[2,4]
    hatz[j] <- summary(fm)$coefficients[2,3]
  }
  
  if (sum(abs(B_hat1)) > 1e-6) {
    B_hats <- as.matrix(B_hat1[, apply(B_hat1, 2, sd)!= 0]) 
    #Ge_impute1 <- x3%*%B_hat1 # weight from cca * SNPs
    #eta <- apply(Ge_impute1, 2, sd) # column wise sd
    eta <- diag(sqrt( t(B_hats) %*% cur.LD %*% B_hats )) 
   
    
    #sig_x <- sqrt(diag(var(x3))) # sigma of the ref panel
    sig_x <- sqrt(diag(cur.LD))
    
    #Ztilde <- hatmu/hats # Z-stat from summary
    #Ztilde <- hatmu/hats
    
    #simulate the summary statistics instead
    Ztilde <- hatz
    Lambda1 <-  diag(sig_x) %*% sweep(B_hats, 2, eta, "/")
    Lambda_sub1 <- Lambda1[, which(eta != 0)]
    # 
    test_stats1 <- t(Lambda_sub1) %*% Ztilde
    #rhoGE1 <- cor(Ge_impute1)
    rhoGE1 <-  matrix(0,nrow=ncol(B_hats),ncol=ncol(B_hats))
    
    for (i in 1:ncol(B_hats)) {
      for (j in 1:ncol(B_hats)) {
        nom <-  t(B_hats[,i])%*%cur.LD%*%B_hats[,j]
        denom <-  sqrt((t(B_hats[,i])%*%cur.LD%*%B_hats[,i])*(t(B_hats[,j])%*%cur.LD%*%B_hats[,j]))
        rhoGE1[i,j]  <-  nom/denom
      }
    }
    
    rhoGE_svd1 <- svd(rhoGE1)
    cond <- 30
    ind_top <- which(rhoGE_svd1$d[1]/rhoGE_svd1$d < cond)
    if(length(ind_top) == 0) ind_top <- 1
    u <- rhoGE_svd1$u
    v <- rhoGE_svd1$v
    us <- as.matrix(u[, ind_top])
    vs <- as.matrix(v[, ind_top])
    d <- rhoGE_svd1$d[ind_top]
    if(length(d) > 1) ds <- diag(1/d) else ds = matrix(1/d)
    rhoGE_ginv <- vs %*% ds %*% t(us)
    mcca.chisq <- c(t(test_stats1) %*% rhoGE_ginv %*% test_stats1)
    mcca.dof  <-  length(ind_top)
    mcca.p <- pchisq(mcca.chisq,length(ind_top),lower.tail=F)
    
    return(mcca.p)
  } else{
    return(NA)
  }
}

#mccasum(B_hat1 =B_hat1, cur.LD = cur.LD)

rlt <- vector()
for (i in 1:nsims) {
  rlt[i] <- suppressWarnings(mccasum(B_hat1 =B_hat1, cur.LD = cur.LD))
}

out <- data.table(gene = gid, pval = rlt)

outname1 <- paste0(output_dir,"all/",paste(gid,"results","MSG_GBJ_ACAT.txt",sep="."))
print(outname1)
write.table(out, quote=F , row.names=F , sep='\t' , file=outname1 )

end_time <- Sys.time()
end_time - start_time





