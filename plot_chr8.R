library(getopt)
suppressMessages(library(dplyr))
chr <- c(1:22, "X", "Y")
chrs = paste("chr", chr, sep = "")
opt_spec <- matrix(c('help', 'h', 0, 'logical', 'help manual',
                     'input', 'i', 1, 'character', 'filepath for *.BAF.txt',
                     'ratioGraph', 'r', 0, 'logical', 'if RatioGraph should be drawn',
                     'BAFGraph', 'f', 0, 'logical', 'if BAFGraph should be drawn',
                     'output', 'o', 1, 'character', 'filepath for output figures',
                     'gene', 'g', 1, 'character', 'gene/transcript symbol for output figures',
                     'chr', 'c', 1, 'character', 'chrmosome for output figures',
                     'cytoBand', 'y', 1, 'character', 'UCSC cytoBand for output figures',
                     'bed', 'b', 1, 'character', 'panel bed file'),
                   byrow=TRUE, ncol=5)
opt = getopt(opt_spec, commandArgs(TRUE))

# test params
# opt <- list()
# opt$input <- "freec3/20110533"
# opt$output <- "chr8_plots"
# opt$ratioGraph <- T
# opt$BAFGraph <- T
# opt$bed <- "~/pt011/data/bed/SureSelect_Human_All_Exon_V6_r2_hs_hg38/S07604514_Covered.bed"
# opt$chr <- "chr8"
# opt$cytoBand <- "~/pt011/data/ucsc/human/cytoBand.txt.gz"
# opt$gene <- "FGFR1"

if (!is.null(opt$help)){
  cat(getopt(opt_spec, usage=TRUE))
  q(save='no', status=1)
}
stopifnot(file.exists(opt$input, opt$bed))

header <- readLines(opt$bed, n = 10)
nskip <- sum(grepl("^track|^browser|^#", header, ignore.case = T))
bed <- read.delim(opt$bed,  header = F, stringsAsFactors = F, skip = nskip)
colnames(bed) <- c("chrom", "start", "end", "gene")

cyto <- NULL

# restrict bed to target region
if (!is.null(opt$chr) && opt$chr %in% chrs) {
  bed <- bed %>% filter(chrom == opt$chr)
  if (!is.null(opt$cytoBand)) {
    cyto <- read.table(opt$cytoBand, sep = "\t")
    colnames(cyto) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
    cyto <- cyto %>% mutate(arm = substr(name, 1, 1)) %>%
      filter(chrom == opt$chr)
  }
  xlabel <- opt$chr
} else if (!is.null(opt$gene)) {
  bed <- bed %>% filter(grepl(paste0("ref\\|", opt$gene, ","), gene))
  xlabel <- opt$gene
  # stop if one gene found in multiple chroms
  stopifnot(length(unique(bed$chrom)) == 1)
} else {
  stop("gene & chr params are both invalid!")
}
lim <- sum(bed$end - bed$start)


dir.create(opt$output, showWarnings = F, recursive = T)
if(is.null(opt$BAFGraph) == F){
  baflist <- list.files(path = opt$input, pattern = "BAF.txt$")
  for(v in 1:length(baflist)){
    cat(paste(baflist[v], "\n"))
    patient <- unlist(strsplit(baflist[v], split = "_BAF"))[1]
    baf <- read.delim(paste(opt$input, baflist[v], sep = "/"), header = T, stringsAsFactors = F)
    
    bed_new <-  NULL
    baf_new <-  NULL

    # filter baf to bed region
    baf_chr <- baf %>% mutate(Chromosome = paste0("chr", Chromosome)) %>%
      filter(Chromosome == bed$chrom[1],
             Position > bed$start[1],
             Position <= max(bed$end))
    bed_chr <- bed

    # re-cal start_norm end_norm in x axis
    for(j in 1:nrow(bed_chr)){
      bed_chr$end_norm[j] <- sum(bed_chr$end[1:j] - bed_chr$start[1:j])
      bed_chr$start_norm[j] <- bed_chr$end_norm[j] - (bed_chr$end[j] - bed_chr$start[j])
    }

    for(j in 1:nrow(baf_chr)){
      n <- which(bed_chr$start < baf_chr$Position[j] & bed_chr$end >= baf_chr$Position[j])
      baf_chr$pos_norm[j] <- baf_chr$Position[j] - bed_chr$start[n] + bed_chr$start_norm[n]
    }
    
    bed_new <- bed_chr

    baf_new <- baf_chr

    # tiff(filename = paste(opt$output, "/", gsub("txt", "tiff", baflist[v]), sep = ""),
    #      width = 1500, height = 400,
    #      units = "px", pointsize = 14, bg = "white", res = NA)
    title <- paste(patient, xlabel, sep = " - ")
    tiff(filename = paste0(opt$output, "/", gsub(" - ", "_", title, xlabel), "_BAF.tiff"),
         width = 3000, height = 800,
         units = "px", pointsize = 10, bg = "white", res = 300)
    # par(mar = c(5,5.5,4,2))
    par(mar = c(2,5,2,1))
    plot(baf_new$pos_norm,baf_new$BAF, xlim = c(0,lim), ylim = c(-0.1,1.1),
         xlab = xlabel, cex.main = 1.1, cex.axis = 1, cex.lab = 1.2, las = 1,
         ylab = "BAF", main = paste(patient, xlabel, sep = " - "),
         pch = ".", col = colors()[1], xaxs="i", xaxt = "n")
    # 将BAF 0.5的位点用橙色标记
    tt <- which(baf_new$A == 0.5)		
    points(baf_new$pos_norm[tt], baf_new$BAF[tt], pch = 20, cex = 0.6, col = "darkorange2")
    tt <- which(baf_new$A != 0.5 & baf_new$A >= 0)
    points(baf_new$pos_norm[tt], baf_new$BAF[tt], pch = 20, cex = 0.6, col = "royalblue2")
    tt <- 1
    pres <- 1
    
    if (length(baf_new$A) > 4) {
      for (j in c(2:(length(baf_new$A) - pres - 1))) {
        if (baf_new$A[j] == baf_new$A[j+pres]) {	
          tt[length(tt)+1] <- j 
        }
      }
      tt <- intersect(tt, which((baf_new$A == 0.5 | baf_new$A == 1 | baf_new$A == 0) & baf_new$A >= 0))
      points(baf_new$pos_norm[tt], baf_new$A[tt], pch = ".", col = "black", cex=3)
      points(baf_new$pos_norm[tt], baf_new$B[tt], pch = ".", col = "black", cex=3)	
    }
    if (length(baf_new$A) > 4) {
      for (j in c(2:(length(baf_new$A) - pres - 1))) {
        if (baf_new$A[j] == baf_new$A[j+pres]) {	
          tt[length(tt)+1] <- j 
        }
      }
      tt <- intersect(tt, which(baf_new$A != 0.5 & baf_new$A != 0 & baf_new$A != 1 & baf_new$A >= 0))
      points(baf_new$pos_norm[tt], baf_new$A[tt], pch = 20, col = "magenta3", cex=0.4)
      points(baf_new$pos_norm[tt], baf_new$B[tt], pch = 20, col = "magenta3", cex=0.4)	
    }
    
    if (!is.null(cyto)) {
      p_pos <- max(cyto$chromEnd[cyto$arm == 'p'])
      centromere <- bed_new %>% filter(start >= p_pos, lag(end) < p_pos) %>% pull(start_norm)
      abline(v = centromere, col = "darkorchid", lwd = 2, lty = 2)
    }
    
    dev.off()
  }
}


if(is.null(opt$ratioGraph) == F){
  ratiolist <- list.files(path = opt$input, pattern = "ratio.txt$")
  cnvlist <- list.files(path = opt$input, pattern = "CNVs$")
  for(v in 1:length(ratiolist)){
    cat(paste(ratiolist[v], "\n"))
    patient <- unlist(strsplit(ratiolist[v], split = "_"))[1]
    ratio <- read.delim(paste(opt$input, ratiolist[v], sep = "/"), header = T, stringsAsFactors = F)
    ratio$Ratio[which(ratio$Ratio <= 0.01)] <- 0.01
    ratio$End <- as.integer(unlist(strsplit(ratio$Gene, split = "-"))[seq(2, nrow(ratio)*2, 2)])
    ratio$Start <- ratio$Start - 1
    cnv <- read.delim(paste(opt$input, cnvlist[v], sep = "/"), header = F, stringsAsFactors = F)
    colnames(cnv)[1:4] <- c("Chromosome", "Start", "End", "CopyNumber")
    cnv$CopyNumber[which(cnv$CopyNumber <= 0.01)] <- 0.01
    
    bed_new <-  NULL
    ratio_new <-  NULL
    cnv_new <-  NULL
    
    bed_chr <- bed
    ratio_chr <- ratio %>% mutate(Chromosome = paste0("chr", Chromosome)) %>%
      filter(Chromosome == bed$chrom[1],
             Start > bed$start[1],
             Start <= max(bed$end))
    cnv_chr <- cnv %>% mutate(Chromosome = paste0("chr", Chromosome)) %>%
      filter(Chromosome == bed$chrom[1],
             Start >= bed$start[1],
             End <= max(bed$end))

    for(j in 1:nrow(bed_chr)){
      bed_chr$end_norm[j] <- sum(bed_chr$end[1:j] - bed_chr$start[1:j])
      bed_chr$start_norm[j] <- bed_chr$end_norm[j] - (bed_chr$end[j] - bed_chr$start[j])
    }

    for(j in 1:nrow(ratio_chr)){
      n <- which(bed_chr$start <= ratio_chr$Start[j] & bed_chr$end >= ratio_chr$Start[j])
      ratio_chr$start_norm[j] <- ratio_chr$Start[j] - bed_chr$start[n] + bed_chr$start_norm[n]
      ratio_chr$end_norm[j] <- ratio_chr$End[j] - bed_chr$start[n] + bed_chr$start_norm[n]
    }
    for(j in 1:nrow(cnv_chr)){
      n <- which(bed_chr$start <= cnv_chr$Start[j] & bed_chr$end >= cnv_chr$Start[j])
      cnv_chr$start_norm[j] <- cnv_chr$Start[j] - bed_chr$start[n] + bed_chr$start_norm[n]
      n <- which(bed_chr$start <= cnv_chr$End[j] & bed_chr$end >= cnv_chr$End[j])
      cnv_chr$end_norm[j] <- cnv_chr$End[j] - bed_chr$start[n] + bed_chr$start_norm[n]
    }
    bed_new <- bed_chr
    ratio_new <- ratio_chr
    cnv_new <- cnv_chr

    title <- paste(patient, xlabel, sep = " - ")
    tiff(filename = paste0(opt$output, "/", gsub(" - ", "_", title, xlabel), "_ratio.tiff"),
         width = 3000, height = 800,
         units = "px", pointsize = 10, bg = "white", res = 300)
    par(mar = c(2,5,2,1))
    plot(ratio_new$start_norm,ratio_new$Ratio,
         ylim = c(min(log2(min(ratio_new$Ratio)),-2), max(log2(max(ratio_new$Ratio)),2)),
         xlim = c(0,lim),
         ylab = "Log2 Copy Number",
         xlab = xlabel, cex.main = 1.5, xaxs="i", 
         cex.axis = 1, cex.lab = 1, las = 2,
         main = title,
         pch = ".", col = colors()[1], xaxt = "n")
    abline(h = 0, col = "grey", lwd = 1, lty = 2)
    abline(h = 1, col = "pink", lwd = 1, lty = 2)
    abline(h = -1, col = "lightblue", lwd = 1, lty = 2)
    # mark CopyNumber 2 Ratios as green
    tt <- which(ratio_new$CopyNumber == 2 )
    if(length(tt) != 0){
      for(j in 1:length(tt)){
        aa <- ratio_new$start_norm[tt[j]]:ratio_new$end_norm[tt[j]]
        points(aa,rep(log2(ratio_new$Ratio[tt[j]]), length(aa)),
               pch = 20, cex = 0.4,col = "darkolivegreen3")
      }
    }
    # mark CopyNumber > 2 Ratios as firebrick3 color
    tt <- which(ratio_new$CopyNumber > 2 )
    if(length(tt) != 0){
      for(j in 1:length(tt)){
        aa <- ratio_new$start_norm[tt[j]]:ratio_new$end_norm[tt[j]]
        points(aa,rep(log2(ratio_new$Ratio[tt[j]]), length(aa)),pch = 20,
               cex = 0.4,col = colors()[136])
      }
    }
    # mark CopyNumber < 2 Ratios in mediumblue color
    tt <- which(ratio_new$CopyNumber < 2 & ratio_new$CopyNumber!= -1)
    if(length(tt) != 0){
      for(j in 1:length(tt)){
        aa <- ratio_new$start_norm[tt[j]]:ratio_new$end_norm[tt[j]]
        points(aa,rep(log2(ratio_new$Ratio[tt[j]]), length(aa)),
               pch = 20, cex = 0.4,col = "mediumblue")
      }
    }
    # mark CopyNumber = 2 in black color
    tt <- which(cnv_new$CopyNumber == 2 )
    if(length(tt) != 0){
      for(j in 1:length(tt)){
        aa <- cnv_new$start_norm[tt[j]]:cnv_new$end_norm[tt[j]]
        points(aa,rep(log2(cnv_new$CopyNumber[tt[j]]/2), length(aa)),
               col = "black", pch = ".", cex = 3)
      }
    }
    # mark CopyNumber > 2 in violet color
    tt <- which(cnv_new$CopyNumber > 2 )
    if(length(tt) != 0){
      for(j in 1:length(tt)){
        aa <- cnv_new$start_norm[tt[j]]:cnv_new$end_norm[tt[j]]
        points(aa,rep(log2(cnv_new$CopyNumber[tt[j]]/2), length(aa)),
               col = "violet", pch = ".", cex = 3)
      }
    }
    # mark CopyNumber < 2 in skyblue color
    tt <- which(cnv_new$CopyNumber < 2 )
    if(length(tt) != 0){
      for(j in 1:length(tt)){
        aa <- cnv_new$start_norm[tt[j]]:cnv_new$end_norm[tt[j]]
        points(aa,rep(log2(cnv_new$CopyNumber[tt[j]]/2), length(aa)),
               col = "skyblue", pch = ".", cex = 3)
      }
    }
    
    if (!is.null(cyto)) {
      p_pos <- max(cyto$chromEnd[cyto$arm == 'p'])
      centromere <- bed_new %>% filter(start >= p_pos, lag(end) < p_pos) %>% pull(start_norm)
      abline(v = centromere, col = "darkorchid", lwd = 2, lty = 2)
    }
    
    dev.off()
  }
}