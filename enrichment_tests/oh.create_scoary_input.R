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

args=commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("One argument must be supplied a matrix file", call.=FALSE)
}
mat<-args[1]

library(ape)
library(phytools)

######PA versus NPA######
#Map<-read.table("/pine/scr/i/s/isai/3837_march_2017/metadata_3837_genomes_june2016_complete.tsv",header=T,sep="\t")
Map<-read.table("metadata_3837_genomes_june2016_complete.tsv",header=T,sep="\t")
Tab<-read.table(mat,header=T,sep="\t",row.names=1,check.names=F)
Tab<-t(Tab)
#The logic is to subset first using the columns of the matrix
Map<-Map[match(rownames(Tab),Map$taxon_oid),]
Map<-Map[which(Map$Classification!="soil"),]
Map<-droplevels(Map)
Tab<-Tab[match(Map$taxon_oid,rownames(Tab)),]
#Remove cols that sum to zero
Tab<-Tab[,which(colSums(Tab)!=0)]
Tab<-t(Tab)

trait<-as.numeric(Map$Classification)-1
df_traits<-data.frame(TraitY=trait)
rownames(df_traits)<-Map$taxon_oid

df_Tab<-data.frame(genes=rownames(Tab))
df_Tab_vals<-as.data.frame(Tab)
df_Tab<-cbind(df_Tab,df_Tab_vals)

#Read the general phylogenetic tree
#tree<-read.tree("/pine/scr/i/s/isai/3837_march_2017/3837_genomes_31scg_june2016.newick")
tree<-read.tree("3837_genomes_31scg_june2016.newick")
#Root the tree first using MRCA
#Use the firmicutes
#Paenibacillus Isolate 2517572151
#Bacillus Isolate 2623620997
root_node<-phytools::findMRCA(tree=tree,tips=c("2517572151","2623620997"))
tree<-reroot(tree,node.number=root_node)

subtree<-drop.tip(phy=tree,tip=which(!(tree$tip.label%in%Map$taxon_oid)))


write.table(df_traits,"traits_scoary_pa_npa.tsv",row.names=T,col.names=NA,sep=",",append=F,quote=F)
write.table(df_Tab,"matrix_scoary_pa_npa.tsv",row.names=F,col.names=T,sep=",",append=F,quote=F)
write.tree(subtree,"tree_scoary_pa_npa.newick")

######Soil versus Root######
rm(Map,Tab,trait,df_traits,df_Tab,df_Tab_vals,tree,root_node,subtree)
#Map<-read.table("/pine/scr/i/s/isai/3837_march_2017/metadata_3837_genomes_june2016_complete.tsv",header=T,sep="\t")
Map<-read.table("metadata_3837_genomes_june2016_complete.tsv",header=T,sep="\t")
Tab<-read.table(mat,header=T,sep="\t",row.names=1,check.names=F)
Tab<-t(Tab)
#The logic is to subset first using the columns of the matrix
Map<-Map[match(rownames(Tab),Map$taxon_oid),]
#Subset soil and root at the same time
Map_soil<-as.character(Map$taxon_oid[which(Map$Classification=="soil")])
Map_ra<-as.character(Map$taxon_oid[which(Map$Root.NotRoot=="RA")])
myids<-c(Map_soil,Map_ra)
Map<-Map[match(myids,Map$taxon_oid),]
Map<-droplevels(Map)
Tab<-Tab[match(Map$taxon_oid,rownames(Tab)),]
#Remove cols that sum to zero
Tab<-Tab[,which(colSums(Tab)!=0)]
Tab<-t(Tab)

#Detine the trait RA in this case should be 1 and soil 0
trait<-c(rep(0,length(Map_soil)),rep(1,length(Map_ra)))
df_traits<-data.frame(TraitY=trait)
rownames(df_traits)<-Map$taxon_oid

df_Tab<-data.frame(genes=rownames(Tab))
df_Tab_vals<-as.data.frame(Tab)
df_Tab<-cbind(df_Tab,df_Tab_vals)

#Read the general phylogenetic tree
#tree<-read.tree("/pine/scr/i/s/isai/3837_march_2017/3837_genomes_31scg_june2016.newick")
tree<-read.tree("3837_genomes_31scg_june2016.newick")
#Root the tree first using MRCA
#Use the firmicutes
#Paenibacillus Isolate 2517572151
#Bacillus Isolate 2623620997
root_node<-phytools::findMRCA(tree=tree,tips=c("2517572151","2623620997"))
tree<-reroot(tree,node.number=root_node)

subtree<-drop.tip(phy=tree,tip=which(!(tree$tip.label%in%Map$taxon_oid)))


write.table(df_traits,"traits_scoary_ra_soil.tsv",row.names=T,col.names=NA,sep=",",append=F,quote=F)
write.table(df_Tab,"matrix_scoary_ra_soil.tsv",row.names=F,col.names=T,sep=",",append=F,quote=F)
write.tree(subtree,"tree_scoary_ra_soil.newick")

