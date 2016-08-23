#
#  Goal to run UpSetR to create a cool looking Venn diagram for the master blueberry pathway list made
#   
#

rm(list = ls())
setwd("/Users/rreid2/Documents/blueberry/bbgage/R/upsetR")

install.packages("UpSetR")
library(UpSetR)


qvalTable <- read.delim("/Users/rreid2/Documents/blueberry/bbgage/R/upsetR/masterTable-keggid.new",header=TRUE,sep="\t")
head(qvalTable)

listinput <- list(pad = qvalTable[,3], cup = qvalTable[,4], green = qvalTable[,5], pink = qvalTable[,6], ripe = qvalTable[,7])
upset(fromList(listinput))
head(qvalTable[])
