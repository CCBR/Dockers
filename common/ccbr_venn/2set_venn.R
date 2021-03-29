#!/usr/bin/Rscript
rm(list=ls())
library(ggplot2)
library(argparse)
library(eulerr)

parser <- ArgumentParser()

parser$add_argument("-l", "--list1", type="character", required=TRUE,
                    help="1st list")

parser$add_argument("-r", "--list2", type="character", required=TRUE,
                    help="2nd list")

parser$add_argument("-p", "--plot", type="character", required=TRUE,
                    help="output Venn PNG file")

parser$add_argument("-m", "--list1only", type="character", required=TRUE,
                    help="items unique to 1st list")

parser$add_argument("-s", "--list2only", type="character", required=TRUE,
                    help="items unique to 2nd list")

parser$add_argument("-c1", "--category1", type="character", required=TRUE,
                    help="name of list1")

parser$add_argument("-c2", "--category2", type="character", required=TRUE,
                    help="name of list2")

parser$add_argument("-c", "--common", type="character", required=TRUE,
                    help="common list of circRNAs")

parser$add_argument("-t", "--plottitle", type="character", required=TRUE,
                    help="title for the plot")

args <- parser$parse_args()


lefttable=read.csv(args$list1,header = FALSE,sep = "\t")
righttable=read.csv(args$list2,header = FALSE, sep = "\t")

# u=union(lefttable$V1,righttable$V1)
i=intersect(lefttable$V1,righttable$V1)
l=setdiff(lefttable$V1,righttable$V1)
r=setdiff(righttable$V1,lefttable$V1)

xx=list(left=lefttable$V1,right=righttable$V1)

png(args$plot)
# ggVennDiagram(xx,category.names = c(args$category1,args$category2))+ scale_fill_gradient(low="blue",high = "red") + ggtitle(args$title)
x=euler(c('A'=length(l),'B'=length(r),'A&B'=length(i)),shape="ellipse")
plot(x,quantities = list(type=c("counts","percent")),fills = list(fill = c("red", "steelblue4"), alpha = 0.6), labels=c(args$category1,args$category2),main = args$plottitle)
dev.off()

write(l,args$list1only)
write(r,args$list2only)
write(i,args$common)

