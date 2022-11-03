#how drastically does sample number (n) affect DEG results?
#perform downsampling to find out
#input is a sample info file and a count matrix - expects raw (not normalized) read counts - further reading in deseq2 vignette

rm(list = ls())
library(dplyr)
library(DESeq2)

#     !user defines!      #######################################################################################################################################################################
setwd("/mnt/hdd/atrostle/example_downsample/") # your dir
reps <- 15  #number of repetitions - each call to deseq takes ~1 minute - depends on your system, multiplied by cut. Each result as text file is ~10MB. Depends on size of imput data
cut <- 4 #number of samples to cut - test from your full n down to (n-cut) if your data is asymetrical on condition, goes from smaller group
num <- "null" #the conditon in numerator (Experimental condition)
denom <- "WT" #the condition in the denominator (Control condition)
Counts <- read.delim("/mnt/hdd/atrostle/example_downsample/countmatrix.txt", stringsAsFactors=F, sep="") # count matrix - genes (or features) are rows, sample names are column names
id.info <- read.delim("/mnt/hdd/atrostle/example_downsample/SampleInfo.txt", stringsAsFactors=F, sep="") #sampleinfo - requires names of relevant fields as "Sample" (sample name/ID) and "Genotype" (condition/phenotype)
#################################################################################################################################################################################################

#make directory to store the results. expects this name but feel free to tweak here and lower if you must
dir.create("results")
sample.names <- (id.info[,"Sample"])
Genotype <- factor(id.info[,"Genotype"])

#get smaller n number to use for downsample
effectiveN <- min(length(which(id.info$Genotype == num)),length(which(id.info$Genotype == denom)))

#make empty results df to fill #we will also write individual deseq results out to files for reference and flexibility of use
downsample_df <- data.frame(gene=row.names(Counts))

#loop for number of samples to cut (cut) and then number of reps to perform (rep)
for(i in 1:cut){
  for(j in 1:reps){
    #format sampleinfo
    colData.counts <- data.frame(genotype=id.info[,"Genotype"])
    rownames(colData.counts) <- id.info$Sample
    colData.counts$index <- row_number(colData.counts)
    
    #get how many n to use rather than how many to cut
    keep = (effectiveN - i)
    
    #get relevant index and keep those
    numindex <- which(colData.counts$genotype %in% as.factor(num))
    denomindex <- which(colData.counts$genotype %in% as.factor(denom))
    #pick which to keep
    wtc <- sample(numindex, keep, replace=FALSE)
    koc <- sample(denomindex, keep, replace=FALSE)
    keepc <- c(wtc,koc)
    cut_colData.counts <- colData.counts[colData.counts$index %in% keepc, ]
    
    cutCounts <- Counts
    cutCounts = cutCounts[,(names(cutCounts)) %in% rownames(cut_colData.counts)]
    
    ddscounts <- DESeqDataSetFromMatrix(countData = cutCounts,
                                        colData = cut_colData.counts,
                                        design = ~ genotype)
    #dim(ddscounts)
    
    # filter out very low read counts
    #10 in half or more
    ddscounts.good = ddscounts[rowSums(counts(ddscounts) >= 10) >= (0.5 * ncol(ddscounts)), ]
    dim(ddscounts.good)
    dds.good <-DESeq(ddscounts.good)
    
    res.nullVSWT_v1 <- results(dds.good, contrast=c("genotype",num,denom))
    name <- paste0(num,"VS",denom)
    colnames(res.nullVSWT_v1) <- paste(name, colnames(res.nullVSWT_v1), sep = "_")
    grp.mean <- sapply(levels(dds.good$genotype), function(lvl) rowMeans(counts(dds.good, normalized=TRUE)[,dds.good$genotype == lvl] ) )
    colnames(grp.mean)
    
    norm.counts <- counts(dds.good, normalized=TRUE)
    res.nullVSWT <- data.frame(res.nullVSWT_v1)
    all <- data.frame(res.nullVSWT, grp.mean, norm.counts)
    all <- data.frame(gene=rownames(all), all)
    
    #make sure to rename columns to work smoothly with tool
    write.table(all, file=paste0("./results/DEG_n_of_",keep,"_rep_",j,".txt"), row.names = F, sep="\t",quote=F)
    
    #only keep a few relevant columns
    mervec <- all[c(1,3,7)]
    names(mervec) <- gsub(name,paste0("n",keep,"_rep",j),names(mervec))
    downsample_df <- merge(mervec,downsample_df,by="gene",all=T)
  }   
}

# hang onto image for convenience
save.image(file = "downsample_results.Rdata")



#bonus code to make a quick deg number plot - only padj filtering, you can also change minimally to filter on effect size
#padjusted = FDR
FDR <- 0.01 #  <-------------- change if you want!
#get only padj columns
downsample_df_padj <- downsample_df[grep("padj", names(downsample_df), value = TRUE)]
#count sig per column, then format
fun <- function(x){return (length(which(x < FDR)))}
plotRes <- data.frame(lapply(downsample_df_padj,fun))
plotRes <- melt(plotRes)
plotRes$n <- as.factor(gsub("_rep.*","",as.character(plotRes$variable))) #could be more elegant

#plot!
f=14
ggplot(plotRes,aes(x=n,y=value)) + 
  theme_bw() +
  geom_jitter(shape=16,size=0.8, alpha=0.8, color="gray50", position=position_jitter(0.2)) +
  geom_boxplot(outlier.alpha = 0, width = .6,alpha=.3,aes(fill=n)) + 
  labs(title=paste0("Downsample - FDR < ",FDR), x ="# of samples (Each Condition)", y = "# of DEG After Removal")+
  theme(axis.text.x = element_text(size=f, face="bold"),axis.title.x = element_text(size=f, face="bold"),
        axis.title.y = element_text(size=f, face="bold"),axis.text.y = element_text(size=f, face="bold"),
        plot.title = element_text(size=f+2, face="bold"), legend.title = element_text(size=f-2, face="bold"), 
        legend.text= element_text(size=f-2, face="bold"))






