library(BioCircos)
library(pafr)
library(ggplot2)
library(gridExtra)

library(Repitools)
library(StructuralVariantAnnotation)
library(gTrack)
library(gGnome)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(dplyr)
library(parallel)
library(plyranges)
library(diffloop)

make_bigwig_gTrack = function(fname, name = "NoName", y0=0, y1=1000){
  this_track = read_bigwig(fname)
  this_track = addchr(this_track)
  seqlevels(this_track , pruning.mode="coarse") = standardchrs
  this_track = binnedAverage(window_bins, coverage(this_track, weight=this_track$score), "score")
  this_gt = gTrack(this_track, y.field="score", bar=TRUE, name = name, y0=y0, y1=y1) #, 
  return(this_gt)
}

draw_amplicon = function(file, covfile, window, cov, gencode){

  sv_curated = jJ(file, geno = FALSE, chr.convert = FALSE)
  gg_curated = gG(juncs = sv_curated)
  gt.cov = make_bigwig_gTrack(covfile, name = "coverage", y1=cov)
  gg_curated$set(name = 'reconstruction')
  
  print(gg_curated$nodes)
  
  p<-plot(c(gencode, gt.cov, gg_curated$gt), window, links = gg_curated$grl)
  return(p)
}

draw_amplicon_manual = function(nodesfile, edgesfile, covfile, window, cov, gencode){
  
  # sv_curated = jJ(file, geno = FALSE, chr.convert = FALSE)
  # gg_curated = gG(juncs = sv_curated)
  
  reconstruct_coordinates <- read.table(nodesfile, header=FALSE, sep="\t")
  colnames(reconstruct_coordinates) <- c("chr", "start","end", "strand",	"target",	"coverage",	"structure","fragment")
  nodes<-makeGRangesFromDataFrame(reconstruct_coordinates,
                                  keep.extra.columns=FALSE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field=c("seqnames", "seqname",
                                                   "chromosome", "chrom",
                                                   "chr", "chromosome_name",
                                                   "seqid"),
                                  start.field="start",
                                  end.field=c("end", "stop"),
                                  strand.field="strand",
                                  starts.in.df.are.0based=FALSE)
  
  edges<-read.table(edgesfile, header=FALSE, sep="\t")
  colnames(edges)<- c("n1","n2","n1.side","n2.side")
  gg_curated <- gG(nodes = nodes, edges = edges)
  
  gt.cov = make_bigwig_gTrack(covfile, name = "coverage", y1=cov)
  gg_curated$set(name = 'reconstruction')
  
  p<-plot(c(gencode, gt.cov, gg_curated$gt), window, links = gg_curated$grl, '3')
  return(p)
}

draw_circle = function(file, sample){
  
  df <- read.csv2(file, header=FALSE, sep="\t", row.names = NULL, skip = 1)
  colnames(df) <- c("chr", "start","stop","strand","target","coverage","structure", "fragment")
  df$size <- df$stop - df$start
  
  genomeChr = df$fragment
  lengthChr = df$size
  names(lengthChr) <- genomeChr
  chrLabels <- df$chr
  
  # arcs_chromosomes = c(df$fragment[1,1]) # Chromosomes on which the arcs should be displayed
  # arcs_begin = c(1)
  # arcs_end = c(sum(lengthChr))
  
  tracklist = BioCircosTextTrack('structure',sample, size = "2em", opacity = 0.5,  x = -0.3, y = 0)
  # tracklist = tracklist + BioCircosArcTrack('myArcTrack', arcs_chromosomes, arcs_begin, arcs_end,
  #                     minRadius = 0.9, maxRadius = 0.95, opacities = c(0.4, 0.4, 1, 0.8), color="gray")
  # tracklist = tracklist + BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0.9, maxRadius = 0.95,
  #                          borderColors = "blue", borderSize = 0.6, fillColors = "blue", labels=chrLabels)  
  
  p<-BioCircos(tracklist, genome = lengthChr,
            genomeFillColor = c("tomato2", "darkblue"),
            genomeTicksDisplay = FALSE,
            genomeLabelTextSize = "8pt",
            chrPad = 0.05)
  return(p)
}

draw_dotplot = function(file,sample){
  
  ali <- read_paf(file)
  # x=query, y=target
  p<-dotplot(ali, label_seqs=TRUE, alignment_colour="blue", xlab=sample, ylab="reference") + theme_bw()
  return(p)
}

root <- "/home/mada/Projects/henssen_lab/Repositories/reconstruct/simulation/data/process"
sample <- "CB2120_initial"

bedfile <- paste(root, sample,"reads/circle-template.bed",sep="/")
paffile1 <-paste(root, sample,"reconstruct/ref-to-reconstruct.paf",sep="/")
paffile2 <-paste(root, sample,"reconstruct/ref-to-sim.paf",sep="/")
fas <- c(Decoil=paste(root, sample,"reconstruct/reconstruct.fasta",sep="/"),
         Sim_Truth=paste(root, sample,"reads/circle-template.fa",sep="/"),
         Shasta=paste(root, sample,"assembly/Assembly.fasta",sep="/"))
fsnv <- paste(root, sample, "sv/sv.sniffles.bedpe", sep="/")
fcov <- paste(root, sample, "mapping/coverage.bw", sep="/")

# --------------------- print simulated circle --------------------------------- #

par(mfrow = c(1, 1))
draw_circle(bedfile,sample)
draw_dotplot(paffile1, paste(sample,"simulated",sep="_")) # align reference fragments to simulation (truth)
draw_dotplot(paffile2,  paste(sample,"reconstruct",sep="_")) # align reference fragments to reconstruction (reconstruct)


window <- GRanges("chr2", IRanges(15631400,16578682))
root <- "/mnt/henssen_lab/data/Nanopore/Evolution/Process"
sample <- "CB2120-initial_26042022"
fsnv <- paste(root, sample, "hg19/ngmlr_hg19.sniffles.filtered.high_curated.bedpe", sep="/")
fcov <- paste(root, sample, "hg19/coverage_hg19.bw", sep="/")
initial<-draw_amplicon(fsnv, fcov, window, 1800, gencode)

reconstruct_coordinates <- read.table("/home/mada/Projects/henssen_lab/Repositories/reconstruct/simulation/data/raw/CB2120_initial_nodes.bed", header=FALSE, sep="\t")
colnames(reconstruct_coordinates) <- c("chr", "start","end", "strand",	"target",	"coverage",	"structure","fragment")
nodes<-makeGRangesFromDataFrame(reconstruct_coordinates,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
edges<-read.table("/home/mada/Projects/henssen_lab/Repositories/reconstruct/simulation/data/raw/CB2120_initial_connections.bed", header=FALSE, sep="\t")
colnames(edges)<- c("n1","n2","n1.side","n2.side")

file_nodes <- "/home/mada/Projects/henssen_lab/Repositories/reconstruct/simulation/data/raw/CB2120_initial_nodes.bed"
file_egdes <- "/home/mada/Projects/henssen_lab/Repositories/reconstruct/simulation/data/raw/CB2120_initial_connections.bed"
inital_manual<-draw_amplicon_manual(file_nodes, file_egdes, fcov, window, 1800, gencode)



# library(DECIPHER)

# specify the path to each FASTA file (in quotes)
# each genome must be given a unique identifier here
# for example: Genome1, Genome2, etc.
fas <- c(Decoil="/home/mada/Projects/henssen_lab/Repositories/reconstruct/simulation/data/process/AnBC/reconstruct/reconstruct.fasta",
         Sim_Truth="/home/mada/Projects/henssen_lab/Repositories/reconstruct/simulation/data/process/AnBC/reads/circle-template.fa",
         Shasta="/home/mada/Projects/henssen_lab/Repositories/reconstruct/simulation/data/process/AnBC/assembly/Assembly.fasta")

# specify where to create the new sequence database
db <- "db4"

# load the sequences from the file in a loop
for (i in seq_along(fas)) {
  Seqs2DB(fas[i], "FASTA", db, names(fas[i]))
}

# map the syntenic regions between each genome pair
synteny <- FindSynteny(db)

# print a summary of the syntenic map (optional)
synteny

# view the syntenic regions (optional)
pairs(synteny) # displays a dot plot of all pairs
plot(synteny) # displays a bar plot of adjacent pairs

# perform alignments of all pairs of genomes
DNA <- AlignSynteny(synteny, db)

# view the aligned syntenic blocks for each pair
unlist(DNA[[1]]) # Genomes 1 and 2
unlist(DNA[[2]]) # Genomes 1 and 3
unlist(DNA[[3]]) # Genomes 2 and 3


