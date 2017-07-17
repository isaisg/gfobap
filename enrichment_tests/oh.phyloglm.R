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

#Print on the matrix we are working on 
print(paste("Working on ",mat,sep=""))

suppressMessages(library(phytools))
suppressMessages(library(phylolm))
Tab<-read.table(mat,header=T,sep="\t",row.names=1,check.names=F)

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
subtree<-drop.tip(phy=tree,tip=which(!(tree$tip.label%in%rownames(Tab))))


#Loop to perform phylogLM
start=2
y<-Tab[,1]
names(y)<-rownames(Tab)
Res <- NULL

for(column in start:ncol(Tab)){
	print(column)
	x<-Tab[,column]
	orthogroup<-colnames(Tab)[column]
	names(x)<-rownames(Tab)
	dat<-as.data.frame(cbind(y,x))
	#Do a first step running phyloglm without bootstrapping
	#Probably it would be worth to rerun phyloglm with bootstraping afterwards 
	#Add part to return res object
	m1 <- tryCatch(phyloglm(formula=y~x, data = dat,phy = subtree,method = "logistic_IG10"), error = function(e) list(coefficients = NA))
    if(is.na(coef(m1)[1])){
      res <- data.frame(orthogroup.id = orthogroup,
                        Estimate = NA,
                        SE = NA,
                        z.value = NA,
                        p.value = NA)
    }else{
      m1.sum <- summary(m1)
      res <- data.frame(orthogroup.id = orthogroup,
                        Estimate = m1.sum$coefficients["x",1],
                        SE = m1.sum$coefficients["x",2],
                        z.value = m1.sum$coefficients["x",3],
                        p.value = m1.sum$coefficients["x",4])
    }
    rm(m1,m1.sum)
    Res <- rbind(Res,res)
}

#Create outfile
outfile<-paste(mat,".phyloglm.tab",sep="")
write.table(x=Res,file=outfile,col.names=T,row.names=F,quote=F,append=F,sep="\t")








