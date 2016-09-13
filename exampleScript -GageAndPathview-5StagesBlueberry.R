#
#  Goal to run pathview on some mapped contigs of blueberry to KEGG Ids. 
#   KEGG IDs were used.
#

rm(list = ls())
setwd("/Users/rreid2/Documents/blueberry/bbgage/R/withNULL")

source("http://bioconductor.org/biocLite.R")
biocLite("pathview")
biocLite("gage")

library(pathview)
library(gage)
#data(gse16873.d)

### To see list of ID types suitable for pathview:
#gene.idtype.list
data(rn.list)
names(rn.list)

#### Get the data read in.  hhas header, 1st column is KEGG ID.
bbfpkm <- read.delim("/Users/rreid2/Documents/R/heatmaps/MASTERTABLE-kaas-fpkm-just15-withNULLS.txt")
head(bbfpkm[1,])

### Log Transform
bblog2 <- log2(bbfpkm[,-1]+1)

#### Get the data frame into a datamatrix !!! (note the [,-1] will ignore the first column !!!)
dmat2=as.matrix(bblog2)

## Check dmatrix is values
str(dmat2)
#### Make rownames for dmatrix (1st column)
rownames(dmat2)=bbfpkm[,1]

### getting difference between Condition1 ([,1:3])   and condition5  ([,13:15])
dmatleaf <- dmat2[,1:3]-dmat2[,4:6]
avgpad <- rowMeans(dmat2[,1:3])
avgcup <- rowMeans(dmat2[,4:6])
avggreen <- rowMeans(dmat2[,7:9])
avgpink <- rowMeans(dmat2[,10:12])
avgripe <- rowMeans(dmat2[,13:15])
avgNULL <- rowMeans(dmat2[,16:18])
avgripeminuspad <- avgripe-avgpad
padnull <- avgpad - avgNULL
cupnull <- avgcup - avgNULL
greennull <- avggreen - avgNULL
pinknull <- avgpink - avgNULL
ripenull <- avgripe - avgNULL

### Give names to columns for for categories
padnames <- grep('pad',colnames(bblog2[1,]))
cupnames <- grep('cup',colnames(bblog2[1,]))
greennames <- grep('green',colnames(bblog2[1,]))
pinknames <- grep('pink',colnames(bblog2[1,]))
ripenames <- grep('ripe',colnames(bblog2[1,]))
nullnames <- grep('null',colnames(bblog2[1,]))

### Need to define KeggSET for gage (in this case gneric Kegg id is "KO" entered as a species.)
kg.ko=kegg.gsets(species="ko")
kg.sigmet=kg.ko$kg.sets[kg.ko$sigmet.idx]

##### Running GAGE
test.g <- gage(dmat2, gsets =kg.sigmet, ref = padnames, samp = ripenames, compare = "unpaired", same.dir = F)
head(test.g)
write.table(test.g$greater,file="test.greater.txt",row.names=T,sep='\t',quote=F) #write the $ greater table
write.table(test.g$stats,file="test.stats.txt",row.names=T,sep='\t',quote=F)

## Testing each stage versus pads to see what pathways change.
gagenull2Pad <- gage(dmat2, gsets =kg.sigmet, ref = nullnames, samp = padnames, compare = "unpaired", same.dir = T)
gagenull2Cup <- gage(dmat2, gsets =kg.sigmet, ref = nullnames, samp = cupnames, compare = "unpaired", same.dir = T)
gagenull2Green <- gage(dmat2, gsets =kg.sigmet, ref = nullnames, samp = greennames, compare = "unpaired", same.dir = T)
gagenull2Pink <- gage(dmat2, gsets =kg.sigmet, ref = nullnames, samp = pinknames, compare = "unpaired", same.dir = T)
gagenull2Ripe <- gage(dmat2, gsets =kg.sigmet, ref = nullnames, samp = ripenames, compare = "unpaired", same.dir = T)

write.table(gagenull2Pad$greater,file="null2pad.greater.txt",row.names=T,sep='\t',quote=F) #write the $ greater table
write.table(gagenull2Cup$greater,file="null2cup.greater.txt",row.names=T,sep='\t',quote=F) #write the $ greater table
write.table(gagenull2Green$greater,file="null2green.greater.txt",row.names=T,sep='\t',quote=F) #write the $ greater table
write.table(gagenull2Pink$greater,file="null2pink.greater.txt",row.names=T,sep='\t',quote=F) #write the $ greater table
write.table(gagenull2Ripe$greater,file="null2ripe.greater.txt",row.names=T,sep='\t',quote=F) #write the $ greater table

##Testing pink and ripe against green
gageGreen2Pink <- gage(dmat2, gsets =kg.sigmet, ref = greennames, samp = pinknames, compare = "unpaired", same.dir = T)
gageGreen2Ripe <- gage(dmat2, gsets =kg.sigmet, ref = greennames, samp = ripenames, compare = "unpaired", same.dir = T)

write.table(gageGreen2Pink$greater,file="green2pink.greater.txt",row.names=T,sep='\t',quote=F) #write the $ greater table
write.table(gageGreen2Ripe$greater,file="green2ripe.greater.txt",row.names=T,sep='\t',quote=F) #write the $ greater table

##Testing pad versus OTHER 4 stages
gagepad2Cup <- gage(dmat2, gsets =kg.sigmet, ref = padnames, samp = cupnames, compare = "unpaired", same.dir = T)
gagepad2Green <- gage(dmat2, gsets =kg.sigmet, ref = padnames, samp = greennames, compare = "unpaired", same.dir = T)
gagepad2Pink <- gage(dmat2, gsets =kg.sigmet, ref = padnames, samp = pinknames, compare = "unpaired", same.dir = T)
gagepad2Ripe <- gage(dmat2, gsets =kg.sigmet, ref = padnames, samp = ripenames, compare = "unpaired", same.dir = T)

write.table(gagepad2Cup$greater,file="pad2cup.greater.txt",row.names=T,sep='\t',quote=F) #write the $ greater table
write.table(gagepad2Green$greater,file="pad2green.greater.txt",row.names=T,sep='\t',quote=F) #write the $ greater table
write.table(gagepad2Pink$greater,file="pad2pink.greater.txt",row.names=T,sep='\t',quote=F) #write the $ greater table
write.table(gagepad2Ripe$greater,file="pad2ripe.greater.txt",row.names=T,sep='\t',quote=F) #write the $ greater table
###########  PATHVIEW     ############

### Subtract away the mean of the other condition....
ripeMinusPad=dmat2[,ripenames]-rowMeans(dmat2[,padnames])
## pathview the sulfur
pv.out <- pathview(gene.data = ripeMinusPad, pathway.id = "00941", out.suffix = "flavonoid-ripe-pad-averaged", kegg.native = T, species = "ko")
## pathview - photosynthesis antenna proteins
pv.out <- pathview(gene.data = ripeMinusPad, pathway.id = "00196", out.suffix = "photosynth_antenna-ripe-avgPad", kegg.native = T, species = "ko")
## porphyrin & chlorophyl
pv.out <- pathview(gene.data = ripeMinusPad, pathway.id = "00860", out.suffix = "porphorynChlorophyl-ripe-avgPad", kegg.native = T, species = "ko")
## Starch and sucrose  00500
pv.out <- pathview(gene.data = ripeMinusPad, pathway.id = "00500", out.suffix = "starchSucrose-ripe-avgPad", kegg.native = T, species = "ko")

##############################################
###### HEATMAP of GAGE results   #############
##############################################
## gage only forward dir
gagepad2Ripe <- gage(dmat2, gsets =kg.sigmet, ref = padnames, samp = ripenames, compare = "unpaired", same.dir = T)
## gage BOTH Dir
gagepad2Ripe.bothdir <- gage(dmat2, gsets =kg.sigmet, ref = padnames, samp = ripenames, compare = "unpaired", same.dir = F)
## We make a UP
gagepad2Ripe.kegg.esg.up <-  esset.grp(gagepad2Ripe$greater, dmat2, gsets = kg.sigmet, ref = padnames, samp = ripenames,  test4up = T, output = T, outname = "gagepad2Ripe.kegg.up", make.plot = T)
## We make a DOWN
gagepad2Ripe.kegg.esg.down <-  esset.grp(gagepad2Ripe$less, dmat2, gsets = kg.sigmet, ref = padnames, samp = ripenames,  test4up = F, output = T, outname = "gagepad2Ripe.kegg.down", make.plot = F)

## Make heatmap
gs=unique(unlist(kg.sigmet[rownames(gagepad2Ripe$greater)[1:3]]))
essData=essGene(gs, dmat2, ref =padnames, samp =ripenames)
head(essData)
ref1=1:3
samp1=4:6

## Individual Gene Heat Map  (doing top 3 in this case....)
my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
for (gs in rownames(gagepad2Ripe$greater)[1:3]) {
outname = substr(gs, 0, 20)
geneData(genes = kg.sigmet[[gs]], exprs = essData, ref = ref1,outname=outname,cols = my_palette,
samp = samp1, txt = TRUE, heatmap = TRUE,Colv = FALSE, Rowv = FALSE,scale = "row", dendrogram = "none", limit = 3, scatterplot = TRUE)
}

## Whole Set heatmap
pad2Ripe.both.2d.sig<-sigGeneSet(gagepad2Ripe.bothdir, outname="pad2Ripe-both-heatmap.kegg",heatmap=TRUE,pdf.size = c(11,8))


########   OLDER STUFF    ##################

#### Read in the IDmap for soltub and refseq
soltubID <- read.delim("
                       ",header=-FALSE,sep="\t")
soltubID=as.matrix(soltubID)
write.table(soltubID,file="soltub.txt",row.names=T,sep='\t',quote=F) #write the table
## Check that soltub is matrix (will be in " "  )
str(soltubID)









#### Now can run mol.sum to map the refseq ids to kegg IDs
dmat2=mol.sum(mol.data=dmatrix, id.map=soltubID[,2:1])
#### Running pathview on overall data and sulfur as the pathway
#pv.out <- pathview(gene.data = dmat2, pathway.id = "00920",species = "sot", out.suffix = "Test1-sulfur", kegg.native = T)

### getting difference between Treatment ([,1:3])   and control  ([,4:6])
dmatleaf <- dmat2[,1:3]-dmat2[,4:6]
avgleafC <- rowMeans(dmat2[,4:6])
avgleafT <- rowMeans(dmat2[,1:3])
avgleafdiff <- avgleafT-avgleafC

## pathview the sulfur
pv.out <- pathview(gene.data = avgleafdiff, pathway.id = "00920",species = "sot", out.suffix = "S10-leaf-sulfur-averaged", kegg.native = T)
pv.out <- pathview(gene.data = dmatleaf, pathway.id = "00920",species = "sot", out.suffix = "S10-leaf-sulfur", kegg.native = T)
## Do for root
dmatroot <- dmat2[,7:9]-dmat2[,10:12]
avgrootC <- rowMeans(dmat2[,10:12])
avgrootT <- rowMeans(dmat2[,7:9])
avgrootdiff <- avgrootT-avgrootC
pv.out <- pathview(gene.data = avgrootdiff, pathway.id = "00920",species = "sot", out.suffix = "S10-root-sulfur-averaged", kegg.native = T)
pv.out <- pathview(gene.data = dmatroot, pathway.id = "00920",species = "sot", out.suffix = "S10-root-sulfur", kegg.native = T)

### Pathview for Photosynthesis
pv.out <- pathview(gene.data = avgrootdiff, pathway.id = "00195",species = "sot", out.suffix = "S10-root-photosynthsis-averaged", kegg.native = T)
pv.out <- pathview(gene.data = avgleafdiff, pathway.id = "00195",species = "sot", out.suffix = "S10-leaf-photosynthsis-averaged", kegg.native = T)


### Gage
### Give names to columns for for categories
leafTnames <- grep('TL',colnames(s10log2[1,]))
leafCnames <- grep('CL',colnames(s10log2[1,]))
rootTnames <- grep('TR',colnames(s10log2[1,]))
rootCnames <- grep('CR',colnames(s10log2[1,]))

### Load in SOT pathway table
soltubpathways <- read.delim("/Users/rreid2/Documents/R/pathview/kegg-apis/kegg-pathway-sot.txt",header=-FALSE,sep="\t")
head=as.matrix(soltubpathways)
spaths=soltubpathways
spaths[,1]=gsub("sot:", "", spaths[,1])
spaths[,2]=gsub("path:", "", spaths[,2])
sot.gs=split(spaths[,1], spaths[,2])
lapply(head(sot.gs), head)

### Run Gage with local SOT pathway pulled down above
root.g <- gage(dmat2, gsets =sot.gs, ref = rootCnames, samp = rootTnames, compare = "unpaired")
write.table(root.g,file="root.g.txt",row.names=T,sep='\t',quote=F) #write the table

### Run Gage where SOT pathway is pulled from web.
#root.g <- gage(dmat2, gsets = kegg.gsets(species = "sot"), ref = rootCnames, samp = rootTnames, compare = "unpaired")
head(root.g$greater)

## Greater
write.table(root.g$greater, file="root-gage-greater.txt",row.names=T,sep='\t',quote=F) #write the table
## Just significant ones.
grapesig <-sigGeneSet(root.g, outname ="sig-root",cutoff = 0.8,heatmap=TRUE)







### Pathview on the Gage results of highest abundance
## 04712
pv.out <- pathview(gene.data = avgrootdiff, pathway.id = "04712",species = "sot", out.suffix = "S10-root-circRhythm-averaged", kegg.native = T)
pv.out <- pathview(gene.data = avgleafdiff, pathway.id = "04712",species = "sot", out.suffix = "S10-leaf-circRhythm-averaged", kegg.native = T)

### DNA replication
###03030
pv.out <- pathview(gene.data = avgrootdiff, pathway.id = "03030",species = "sot", out.suffix = "S10-root-DNArep-averaged", kegg.native = T)
pv.out <- pathview(gene.data = avgleafdiff, pathway.id = "03030",species = "sot", out.suffix = "S10-leaf-DNARep-averaged", kegg.native = T)

###00650
### Butanoate metabolism
pv.out <- pathview(gene.data = avgrootdiff, pathway.id = "00650",species = "sot", out.suffix = "S10-root-butonoate-averaged", kegg.native = T)
pv.out <- pathview(gene.data = avgleafdiff, pathway.id = "00650",species = "sot", out.suffix = "S10-leaf-butonoate-averaged", kegg.native = T)

###03440
pv.out <- pathview(gene.data = avgrootdiff, pathway.id = "03440",species = "sot", out.suffix = "S10-root-homrecomb-averaged", kegg.native = T)
pv.out <- pathview(gene.data = avgleafdiff, pathway.id = "03440",species = "sot", out.suffix = "S10-leaf-homrecomb-averaged", kegg.native = T)

#### To get gene names out, to see gene names mapped.
head(pv.out$plot.data.gene)
### Writing these to file
write.table(pv.out$plot.data.gene,file="sulfur-kegg-pathwaview.txt",row.names=T,sep='\t',quote=F) #write the table

###Getting Genes for oxidative phosphorylation
pv.oxidpho.root.out <- pathview(gene.data = avgrootdiff, pathway.id = "00190",species = "sot", out.suffix = "S10-root-oxidpho-averaged", kegg.native = T)
pv.oxidpho.leaf.out <- pathview(gene.data = avgleafdiff, pathway.id = "00190",species = "sot", out.suffix = "S10-leaf-oxidpho-averaged", kegg.native = T)
write.table(pv.oxidpho.root.out$plot.data.gene,file="oxidpho-kegg-genes-root.txt",row.names=T,sep='\t',quote=F) #write the table
write.table(pv.oxidpho.leaf.out$plot.data.gene,file="oxidpho-kegg-genes-leaf.txt",row.names=T,sep='\t',quote=F) #write the table

###Getting Genes for terpenoid backbone biosynthesis
pv.terp.root.out <- pathview(gene.data = avgrootdiff, pathway.id = "00900",species = "sot", out.suffix = "S10-root-terp-averaged", kegg.native = T)
pv.terp.leaf.out <- pathview(gene.data = avgleafdiff, pathway.id = "00900",species = "sot", out.suffix = "S10-leaf-terp-averaged", kegg.native = T)
write.table(pv.terp.root.out$plot.data.gene,file="terp-kegg-genes-root.txt",row.names=T,sep='\t',quote=F) #write the table
write.table(pv.terp.leaf.out$plot.data.gene,file="terp-kegg-genes-leaf.txt",row.names=T,sep='\t',quote=F) #write the table

