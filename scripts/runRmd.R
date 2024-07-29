#!/usr/bin/env Rscript
suppressMessages(require(optparse))
if(!suppressWarnings(require("myutils",character.only=TRUE)==FALSE)){
    suppressMessages(require(myutils))
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

suppressMessages(using("myutils"))

option_list <- list(
    make_option(c("--gene"), action="store",
                default=,
                type="character",
                help=" [default: %default]"),
    make_option(c("--dir_data"), action="store",
                default="/panfs/accrepfs.vampire/nobackup/cgg/yany14/git_project/MTSG/data",
                type="character",
                help=" [default: %default]"),
    make_option(c("--dir_script"), action="store",
                default="/panfs/accrepfs.vampire/nobackup/cgg/yany14/git_project/MTSG/scripts",
                type="character",
                help=" [default: %default]"),
    make_option(c("--dir_intermediates"), action="store",
                default="",
                type="character",
                help=" [default: %default]"),
    make_option(c("--dir_knit_root"), action="store",
                default="",
                type="character",
                help=" [default: %default]"),
    make_option(c("--sumstats"), action="store",
                default="1",
                type="character",
                help=" [default: %default]"),
    make_option(c("--feature_col"), action="store",
                default=2,
                type="integer",
                help=" [default: %default]"),
    make_option(c("--double"), action="store",
                default=0.01,
                type="double",
                help=" [default: %default]"),
    make_option(c("--optFile"), action="store",
                default="/home/yany14/Rscripts/rmd/cca_multiple/cca.rmdOpt.R",
                type="character",
                help=" [default: %default]"),
    make_option(c("--run_pandoc"), action="store",
                default=TRUE,
                type="logical",
                help=" [default: %default]")
)

option <- myoption()
option_list <- c(option_list,
                 option$table,
                 option$plot,
                 option$singlecell,
                 ## option$gwas,
                 ## option$summary,
                 ## option$highc,
                 ## option$enhancer,
                 ## option$gibbs,
                 option$general,
                 option$heatmap
                 )

opt <- parse_args(OptionParser(option_list=option_list))


suppressMessages(using("ggplot2"))
suppressMessages(library(plyr))
suppressMessages(using("dplyr"))
suppressMessages(using("magrittr"))
suppressMessages(using("reshape2"))
suppressMessages(using("data.table"))
suppressMessages(using("gtools"))
suppressMessages(using("RColorBrewer"))
suppressMessages(using("extrafont"))
suppressMessages(using("gridExtra"))
suppressMessages(using("grid"))
suppressMessages(using("rmarkdown"))

##preload functions

loadfonts()
theme_figure <- theme_classic() + theme(panel.background = element_rect(fill = 'transparent', colour = NA),legend.title=element_blank(),text=element_text(family="Arial"),axis.title.x = element_blank(),plot.title = element_text(hjust=0.5, face="bold",size=12),axis.line = element_blank())


## if(opt$head)
## {
##     x <- fread(cmd=paste("cat",opt$input,sep=" "))
## }else{
##     x <- fread(cmd=paste("cat",opt$input,sep=" "),head=F,fill=TRUE)
## }

if(grepl(",",opt$input,perl = TRUE)){
    input <- strsplit(opt$input,split=",",perl=T) %>% unlist
    input <- sapply(input,normalizePath)
    opt$input <- paste(input,collapse = ",")
}else{
    opt$input <- normalizePath(opt$input)
}

print(opt$input)

if(length(opt$input2)){
    if(grepl(",",opt$input2,perl = TRUE)){
        input <- strsplit(opt$input2,split=",",perl=T) %>% unlist
        input <- sapply(input,normalizePath)
        opt$input2 <- paste(input,collapse = ",")
    }else{
        opt$input2 <- normalizePath(opt$input2)
    }
    print(opt$input2)
}




if(opt$rds!="1"){
    opt$rds <- normalizePath(opt$rds)
}
opt$rmd <- normalizePath(opt$rmd)
dir.create(opt$dir_out,recursive = T,showWarnings = F)
opt$dir_out <- normalizePath(opt$dir_out)

opt$dir_working <- getwd()
print(opt$dir_working)

print(opt$rds)
print(opt$rmd)
print(opt$dir_out)
print(opt$output)
if(opt$test){
    save.image(file="test.rda")
    quit("no")
}


dir.create(paste0(opt$dir_out,"/knit_root_",opt$output),recursive =T, showWarnings =F)
dir.create(paste0(opt$dir_out,"/intermediates_",opt$output),recursive =T, showWarnings =F)

rmarkdown::render(input = opt$rmd, output_file = paste(opt$output,opt$ext,sep = ".") , output_dir = opt$dir_out,knit_root_dir=paste0(opt$dir_out,"/knit_root_",opt$output), intermediates_dir = paste(opt$dir_out,"/intermediates_",opt$output,sep = ""), clean = opt$clean,run_pandoc=opt$run_pandoc)

## remove intermediate dir
unlink(paste0(opt$dir_out,"/knit_root_",opt$output),recursive = T)
unlink(paste(opt$dir_out,"/intermediates_",opt$output,sep = ""),recursive = T)

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



