# (C) Copyright 2017 Isai Salas Gonzalez
#
#    This file is part of gdobap.
#
#    gdobap is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    gdobap is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with gdobap.  If not, see <http://www.gnu.org/licenses/>.

args=commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("One argument must be supplied a matrix file", call.=FALSE)
}
prefix<-args[1]
name_taxon<-args[2]

#Print on the m
res<-list.files(pattern=".tab$")
master_res<-NULL
for(file in res){
  Tab<-read.table(file,sep="\t",header=T,check.names=F)
  master_res<-rbind(master_res,Tab)
}

colnames(master_res)<-paste(prefix,colnames(master_res),sep="_")

outfile<-paste("phyloglm_",name_taxon,"_results.tsv",sep="")
write.table(x=master_res,file=outfile,append=F,sep="\t",quote=F,row.names=F,col.names=T)
