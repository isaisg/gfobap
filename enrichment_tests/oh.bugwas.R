# (C) Copyright 2017 Isai Salas Gonzalez
#
#    This file is part of gfobap.
#
#    gfobap is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    gfobap is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with gfobap.  If not, see <http://www.gnu.org/licenses/>.

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("Two argument must be supplied a matrix ", call.=FALSE)
}
library(ape)
library(phytools)
mat<-args[1]

#####PA versus NPA#####
#Read the Map #Change the path if needed
#Map<-read.table("/pine/scr/i/s/isai/3837_march_2017/metadata_3837_genomes_june2016_complete.tsv",header=T,sep="\t")
Map<-read.table("metadata_3837_genomes_june2016_complete.tsv",header=T,sep="\t")
#Read the matrix file passed
Tab<-read.table(mat,header=T,sep="\t",row.names=1,check.names=F)
Tab<-t(Tab)
#The logic is to subset first using the columns of the matrix
Map<-Map[match(rownames(Tab),Map$taxon_oid),]
#Here we are  removing soil
Map<-Map[which(Map$Classification!="soil"),]
Map<-droplevels(Map)
Tab<-Tab[match(Map$taxon_oid,rownames(Tab)),]
#Remove zero that sum to zero
Tab<-Tab[,which(colSums(Tab)!=0)]
#Binarize
Tab[Tab>1]<-1
#Remove singletons and core genes
Tab<-Tab[,which(colSums(Tab)!=1)]
Tab<-Tab[,which(colSums(Tab)!=nrow(Tab))]

Tab[Tab==1]<-"A"
Tab[Tab==0]<-"G"


Tab<-t(Tab)

trait<-as.numeric(Map$Classification)-1
df_traits<-data.frame(ID=Map$taxon_oid,pheno=trait)

#Create the dataframe with the 
tab_df<-data.frame(ps=rownames(Tab))
rownames(Tab)<-NULL
tab_df<-cbind(tab_df,Tab)


#Read the phylogenetic tree and root if 
#tree<-read.tree("/pine/scr/i/s/isai/3837_march_2017/3837_genomes_31scg_june2016.newick")
tree<-read.tree("3837_genomes_31scg_june2016.newick")
#Root the tree first using MRCA
#Use the firmicutes
#Paenibacillus Isolate 2517572151
#Bacillus Isolate 2623620997
root_node<-phytools::findMRCA(tree=tree,tips=c("2517572151","2623620997"))
tree<-reroot(tree,node.number=root_node)

#Subset the phylogenetic tree
subtree<-drop.tip(phy=tree,tip=which(!(tree$tip.label%in%Map$taxon_oid)))


#Create directory to contain the split files
outdir="bugwas_pa_npa"
dir.create(outdir)

setwd('bugwas_pa_npa')
write.tree(file="tree.txt",phy=subtree)
write.table(x=tab_df,file="geno.txt",sep="\t",row.names=F,col.names=T,quote=F,append=F)
write.table(x=df_traits,file="pheno.txt",sep="\t",row.names=F,col.names=T,quote=F,append=F)



##Run Bugwas
library(bugwas)
gem_path="~/bin/bugwas/gemma/gemma.0.93b"
prefix="pa_npa"
data <- lin_loc(gen = "geno.txt", pheno = "pheno.txt", phylo = "tree.txt", prefix = prefix, gem.path = gem_path)
saveRDS(data,file="lin_loc.rds")

######RA versus roil######
setwd('../')
rm(data,gem_path,data,subtree,outdir,tree,root_node,tab_df,df_traits,trait,Tab,Map)
#Read the Map #Change the path if needed
#Map<-read.table("/pine/scr/i/s/isai/3837_march_2017/metadata_3837_genomes_june2016_complete.tsv",header=T,sep="\t")
Map<-read.table("metadata_3837_genomes_june2016_complete.tsv",header=T,sep="\t")
#Read the matrix file passed
Tab<-read.table(mat,header=T,sep="\t",row.names=1,check.names=F)
Tab<-t(Tab)
#The logic is to subset first using the columns of the matrix
Map<-Map[match(rownames(Tab),Map$taxon_oid),]
#Here we are  removing soil
Map_soil<-as.character(Map$taxon_oid[which(Map$Classification=="soil")])
Map_ra<-as.character(Map$taxon_oid[which(Map$Root.NotRoot=="RA")])
myids<-c(Map_soil,Map_ra)
Map<-Map[match(myids,Map$taxon_oid),]
Map<-droplevels(Map)
Tab<-Tab[match(Map$taxon_oid,rownames(Tab)),]
#Remove cols that sum to zero
Tab<-Tab[,which(colSums(Tab)!=0)]
#Binarize
Tab[Tab>1]<-1
#Remove singletons and core genes
Tab<-Tab[,which(colSums(Tab)!=1)]
Tab<-Tab[,which(colSums(Tab)!=nrow(Tab))]

Tab[Tab==1]<-"A"
Tab[Tab==0]<-"G"


Tab<-t(Tab)

trait<-c(rep(0,length(Map_soil)),rep(1,length(Map_ra)))
df_traits<-data.frame(ID=Map$taxon_oid,pheno=trait)

#Create the dataframe with the 
tab_df<-data.frame(ps=rownames(Tab))
rownames(Tab)<-NULL
tab_df<-cbind(tab_df,Tab)


#Read the phylogenetic tree and root if 
#tree<-read.tree("/pine/scr/i/s/isai/3837_march_2017/3837_genomes_31scg_june2016.newick")
tree<-read.tree("3837_genomes_31scg_june2016.newick")

#Root the tree first using MRCA
#Use the firmicutes
#Paenibacillus Isolate 2517572151
#Bacillus Isolate 2623620997
root_node<-phytools::findMRCA(tree=tree,tips=c("2517572151","2623620997"))
tree<-reroot(tree,node.number=root_node)

#Subset the phylogenetic tree
subtree<-drop.tip(phy=tree,tip=which(!(tree$tip.label%in%Map$taxon_oid)))




outdir="bugwas_ra_soil"
dir.create(outdir)

setwd('bugwas_ra_soil')
write.tree(file="tree.txt",phy=subtree)
write.table(x=tab_df,file="geno.txt",sep="\t",row.names=F,col.names=T,quote=F,append=F)
write.table(x=df_traits,file="pheno.txt",sep="\t",row.names=F,col.names=T,quote=F,append=F)

##Run Bugwas
library(bugwas)
gem_path="~/bin/bugwas/gemma/gemma.0.93b"
prefix="ra_soil"
data <- lin_loc(gen = "geno.txt", pheno = "pheno.txt", phylo = "tree.txt", prefix = prefix, gem.path = gem_path)
saveRDS(data,file="lin_loc.rds")
