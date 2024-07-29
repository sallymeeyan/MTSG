#!/usr/bin/env Rscript
suppressMessages(require(optparse))

option_list <- list(
    make_option(c("--input"), action="store",
                default=NULL,
                type="character",
                help=" [default: %default]"),
    make_option(c("--ref_ld_chr"), action="store", 
                default="/data1/jiy1/LDREF/1000G.EUR.",
                #default=NA,
                type='character',help="Prefix to reference LD files in binary PLINK format by chromosome [required]"),
    make_option(c("--dir_cov"), action="store", 
                default="/scratch/cgg/jiy1/CCA/BioVU_5000REF/cov",
                type='character',help="Prefix to reference LD files in RData format by gene id [required]"),
    make_option(c("--test"), action="store",
                default=FALSE,
                type="logical",
                help=" [default: %default]")

)

option_list <- c(option_list)

opt <- parse_args(OptionParser(option_list=option_list))

if(opt$test){
    save.image("test.rda")
    q("no")
}

suppressMessages(require("magrittr"))
suppressMessages(require("reshape2"))
suppressMessages(require("data.table"))
suppressMessages(require(plink2R))

x <- fread(cmd=paste("cat",opt$input,sep=" "))

gid <- gsub("\\.x_all.*","",basename(opt$input),perl = T)
chr <- gsub("chr","",unique(x[[1]])[1],perl = T)

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
if(!dir.exists(paste0(opt$dir_cov,'/'))){
    dir.create(paste0(opt$dir_cov,'/'))
}
outname1 <- paste0(opt$dir_cov,'/',paste(gid,"cov.RData",sep="."))
print(outname1)
print(str(genos))
bim = genos$bim[match(comm_snps,genos$bim[,2]),]
save(bim,LD,file=outname1)


## if(opt$progress_bar){
##     pb <- txtProgressBar(min = 0,max = opt$number, style = 3)
## }
## if(opt$progress_bar){
##     setTxtProgressBar(pb, i)
## }


#####################################################
############# title of figure
## opt$ggtitle <- "1"
## x_plot <- x_plot  + ggtitle(opt$ggtitle)




#####################################################
############# label using for gridExtra

## fs <- 15
## A.grob <- textGrob(
##     label = "a",
##     x = unit(0, "lines"),
##     y = unit(0, "lines"),
##     hjust = 0, vjust = 0,
##     gp = gpar(fontsize = fs))
## B.grob <- textGrob(
##     label = "b",
##     x = unit(0, "lines"),
##     y = unit(0, "lines"),
##     hjust = 0, vjust = 0,
##     gp = gpar(fontsize = 14))
## C.grob <- textGrob(
##     label = "c",
##     x = unit(0, "lines"),
##     y = unit(0, "lines"),
##     hjust = 0, vjust = 0,
##     gp = gpar(fontsize = 14))



#####################################################
###########grid.arrange

## lay <- rbind(c(1,1,2,2,2),
##              c(1,1,2,2,2),
##              c(3,3,2,2,2),
##              c(3,3,2,2,2),
##              c(3,3,2,2,2)
##              )

## opt$nrow <- 1
## opt$width <- 7.2
## opt$height <- 8
## dev.new(width=opt$width,height=opt$height)
## grid.arrange(arrangeGrob(top=A.grob,outplot[[1]]),arrangeGrob(top=B.grob,outplot[[2]]),arrangeGrob(top=C.grob,outplot[[3]]),layout_matrix=lay)


#####################################################
###########show/save figure

### if(is.na(opt$output)){
###	x11() 
### 	cat("Plotting QQ plot ...\n") else cat("Saving QQ plot to ", opt$output, " ...\n");	
### 	p
### 	message("Click the graph to exit..."); loc=locator(1); 
### 	loc=locator(1);
### 	res=dev.off();
### } else {
###   	 #Save with ggsave
###   	 ggsave(x_plot,file=paste(opt$output,opt$ext,sep=""),width=opt$width,height=opt$height)
###   	 #Save with my function
### 	 opendev(opt)
###	 dev.off()
### }



