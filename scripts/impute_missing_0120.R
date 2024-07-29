#!/usr/bin/env Rscript
suppressMessages(require(optparse))
if(!suppressWarnings(suppressMessages(require("myutils",character.only=TRUE)))==FALSE){
    suppressWarnings(suppressMessages(require(myutils)))
}else{
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/addGeom.R" ,sep = ""))
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/formatPvalue.R" ,sep = ""))
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/opendev.R" ,sep = ""))
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/options.R" ,sep = ""))
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/plotPLabelFormat.R" ,sep = ""))
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/prettyPaginate.R" ,sep = ""))
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/using.R" ,sep = ""))
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/VAlignPlots.R" ,sep = ""))
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/view.R" ,sep = ""))
}



option_list <- list(
    make_option(c("--gid"), action="store",
                default="ENSG00000000457",
                type="character",
                help=" [default: %default]"),
    make_option(c("--impute_mode"), action="store",
                default=TRUE,
                type="logical",
                help=" [default: %default]")
)

option <- myoption()

option_list <- c(option_list,
                 option$table,
                 option$plot,
                 option$heatmap,
                 ## option$gwas,
                 ## option$summary,
                 ## option$highc,
                 ## option$enhancer,
                 ## option$gibbs,
                 option$general
                 )

opt <- parse_args(OptionParser(option_list=option_list))

if(opt$test){
    save.image(file = paste(opt$output,".rda",sep = ""));
    stop("test rda saved...")
}


suppressMessages(require(data.table))
suppressMessages(require(readr))
suppressMessages(require(dplyr))
suppressMessages(require(magrittr))
##suppressMessages(require(mice))
suppressMessages(require(Amelia))


lf <- list.files(opt$input)
idx <- intersect(grep("^splice_",lf,perl = T),grep(opt$gid,lf,perl = T))
print(idx)
if(length(idx)==0){
    stop("no splice files")
}
lf <- lf[idx]
print(lf)

alldat <- lapply(lf,
                 FUN=function(x){
                     fread(cmd=paste("cat",paste(opt$input, x, sep = "/"),sep=" "))}
                     )
                 



## Amelia

cat("impute...\n")
set.seed(100)
all.splicings3 <- lapply(1:length(alldat), function(i){
    cat(lf[i])
    ##splicing matrix of each tissue
    dat <- alldat[[i]]
    sp <- dat[[5]][1]
    cnames <- colnames(dat)[6:ncol(dat)]
    dat2 <- dat[,-(2:5)] %>% t() 
    tissues <- as.character(dat2[1,])
    dat2 <- data.table(dat2[-1,],keep.rownames = T)
    names(dat2)[-1] <- tissues
    dat2 <- dat2[,(names(dat2)[-1]):=lapply(.SD,as.numeric),.SDcols = names(dat2)[-1]]

    ## identify the un-vary columns, these columns should not be imputed
    missids <- apply(dat2[,-1,with = F],2,function(x){
        summ <- summary(x) %>% unname()
        m <- ifelse(identical(summ[1],summ[3],summ[5]),1,0)
        return(m)
    })

    ## row name col should not be imputed as well
    missids <- c(1,missids)

    missidx <-  which(missids==1)
    misscol <-  length(missidx)-1 ## tissue col is always 1 so should be minused
    print(misscol)

    ## Data do not needs imputation
    dat3 <- dat2[,missidx,with=F]

    ## Data needs imputation
    dat4 <- dat2[,-missidx,with=F]
    if(misscol <= 20){
        if(opt$impute_mode){
            sink("/dev/null") ## suppress amelia's message
            dat.imp <- amelia(dat4,m=2, parallel = "multicore", ncpus = opt$threads)
            sink(NULL) ## ending sending the message to /dev/null
            imp <- dat.imp$imputations[[2]]
            ##combine imputed and unimputed cols
            dat5 <- data.table(dat3,imp)[,names(dat2),with=F]

            ## fill in NA with minimium after imputation
            if(misscol >0){
                mm <- min(dat2[,-1,with=F],na.rm = T)
                ##fill in the unvarying columns 
                dat5[is.na(dat5)] <- mm
            }
            return(dat5)
        }else{
            return(dat2)
        }
    }
})

names(all.splicings3) <- lf
all.splicings4 <- lapply(1:length(all.splicings3),
                         function(i){
                             if(is.null(all.splicings3[[i]])){NULL}
                             else{
                                 data.table(splice=rep(names(all.splicings3)[i],nrow(all.splicings3[[i]])),all.splicings3[[i]])}})

cat("combine results...\n")
## dim: number of splicing * number of person = rownumber
## dim: colnumber = number of tissues 
## imps <- do.call(rbind,all.splicings3)

## available splicing events index
spid <- lapply(all.splicings3,is.null) %>% do.call(c,.) %>% as.numeric %>% ifelse(.,0,1)
print(spid)
if(sum(spid)==0){
    write.table(spid,paste0(opt$dir_out,"/",opt$gid,".all.final.matrix"),row.names = T, col.names = T)
    q(save = "no")
}

## use one of the matrix for individual IDs and number of individuals
## spid1 <- which(spid==1)[1]
##rownames - individual IDs
## snglrownames <- all.splicings3[[spid1]][[1]]
## num_sample <- dim(all.splicings3[[spid1]])[1]

##re-arrange the matrix- each tissue format a block with dims: number of individual* number of sp events
##because the input for SDA required to be stacking of each tissue 
## imps2 <- lapply(2:ncol(imps),function(i){
##     matrix(imps[[i]],ncol = dim(imps)[1]/num_sample, nrow = num_sample)
## }) %>% do.call(rbind,.)

## assign rownames to the reshaped data
## cat("write files...\n")
## imps2 <- data.table(sample = rep(snglrownames,ncol(imps)-1),imps2)


## dim: number of splicing * number of person = rownumber
## dim: colnumber = number of tissues 
comb_data <- do.call(rbind,all.splicings4)

comb_data_melt <- melt(comb_data,id.vars = names(comb_data)[1:2])

f <- paste(names(comb_data_melt)[3],"+",names(comb_data_melt)[2],"~", names(comb_data_melt)[1])
comb_data2 <- dcast(comb_data_melt,as.formula(f))

comb_data2 <- comb_data2[order(get(names(comb_data2)[1]),get(names(comb_data2)[2])),]


## write.table(imps2,paste0(opt$dir_out,"/",opt$gid,".all.final.matrix.wt.names"),row.names = T, col.names = T)
## write.table(imps2,paste0(opt$dir_out,"/",opt$gid,".all.final.matrix"),row.names = T, col.names = T, quote=F)
write.table(comb_data2,paste0(opt$dir_out,"/",opt$gid,".all.final.matrix"),row.names = F, col.names = T, quote=F)
