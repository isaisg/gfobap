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


library(RColorBrewer)
#library(ggtree)
library(ape) 
library(phytools)
library(pROC)
library(ggplot2)
library(labdsv)
library(vegan)


##Define colors for the analysis
#RGB colors 
#PA 22,118,15 #16760f
#npa 255,177,0 #ffb100
#soil: 102,51,0 #663300
#RA: 65,105,225 #4169e1
pal_panpa<-c("#ffb100","#16760f")
pal_rasoil<-c("#4169e1","#663300")
pal_def<-brewer.pal(name = "Set1",n = 9)
large_colors<-colorRampPalette(colors = brewer.pal(n = 9,name = "Set1"))(17)

#Function to save
oh.ggsave.svg<-function(ggobject=p,outname="",outdir="figures/",width=20,height=15,device="svg",dpi=600){
  dir.create(outdir, showWarnings = FALSE)
  myfilename<-paste(outdir,outname,sep="")
  ggsave(filename = myfilename,plot = ggobject,width =width,height = height,device = device,dpi = dpi)
}

##Functions taken from ohchibi own package to perform Pcoa and plot Pcoa objects
#Function to perform pco and plot the pco object
oh.pco<-function(Tab=Tab,Map=Map,ndim=3,eig=T,distfun=distfun,id_var="Id"){
  Tab_dist <- distfun(Tab)
  res <- cmdscale(Tab_dist, k = ndim, eig = eig)
  perc_var <- (100*res$eig / sum(res$eig[res$eig > 0 ]))
  perc_var<-perc_var[1:ndim]
  points<-res$points
  #Define name for the columns in the points object
  newcols<-NULL
  for (i in 1:ndim){
   temp<-paste("PCo",i,sep="")
   newcols<-c(newcols,temp)
  }
  colnames(points)<-newcols
  names(perc_var)<-newcols
  perc_var<-round(perc_var,digits = 2)
  points<-points[match(Map[,which(colnames(Map)==id_var)],rownames(points)),]
  Map_pco <-cbind(Map,points)
  toret=list(Map_pco = Map_pco, variance_explained=perc_var)
  return(toret)
}

oh.plot.pco<-function(list_ohpco=NULL,col_val=NULL,comp_a="PCo1",comp_b="PCo2",mypch=22,size=25,alpha=0.7,stroke=1.5){
  Map_pco<-list_ohpco$Map_pco
  myvar=list_ohpco$variance_explained
  p <- ggplot(data = Map_pco,aes_string(x = comp_a,y = comp_b))+
    geom_point(size = size,alpha=alpha,pch=mypch,colour="#414141",stroke = stroke,aes(fill = get(col_val)))+
    #scale_fill_manual(values = mix.colors)+
    xlab(label = paste(comp_a,"(",myvar[which(names(myvar)==comp_a)],"%)",sep="")) + 
    ylab(label = paste(comp_b,"(",myvar[which(names(myvar)==comp_b)],"%)",sep="")) +
    guides(fill = guide_legend(override.aes = list(size=12))) +
    theme(axis.line = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.grid.major = element_line(colour =   "#D9D9D9"),
          panel.grid.minor = element_line(colour = "#D9D9D9"),
          panel.border = element_rect(fill=NA,color =  "#414141",size = 1),
          axis.ticks = element_line(colour = "black",size = 2.5),
          axis.text.x = element_text(family = "AvantGarde",face = "plain",size =20,colour="#414141"),
          axis.text.y = element_text(family = "AvantGarde",face="plain",size=20,colour="#414141"),
          axis.title.x = element_text(family = "AvantGarde",face="plain",size = 30,colour = "#414141"),
          axis.title.y = element_text(family = "AvantGarde",face="plain",size=30,colour="#414141"),
          legend.background = element_blank(),legend.key.size = unit(2,"point"),
          legend.title=element_blank(),legend.key = element_blank(), 
          legend.text = element_text(size=20,
          family = "AvantGarde",face = "plain",colour = "#414141"),
          legend.position ="right") 
  p<-p + guides(fill=guide_legend(keywidth=0.5,keyheight=0.5,default.unit="inch",override.aes = list(size=10,stroke=stroke,shape=mypch,alpha=alpha)))
  return(p)
}



#Function to perform PCOA using and enriched list
wrapper.pcoa.glm.enrichment.pa.npa<-function(file=NULL,orthogroups=NULL){
  different_names<-orthogroups
  Tab<-read.table(file,header=T,row.names = 1,sep="\t",check.names = F)
  Map<-read.table("../metadata_3837_genomes_june2016_complete.tsv",header = T,sep="\t")
  Map<-Map[match(colnames(Tab),Map$taxon_oid),]
  Map<-Map[which(Map$Classification!="soil"),]
  Map<-droplevels(Map)
  Tab<-Tab[,match(Map$taxon_oid,colnames(Tab))]
  Tab<-Tab[which(rowSums(Tab)!=0),]
  Tab<-Tab[which(rownames(Tab)%in%different_names),]
  original_num<-ncol(Tab)
  #Removegenomes that go to zero
  Tab<-Tab[,which(colSums(Tab)!=0)]
  if(ncol(Tab)!=original_num){
    mypco<-NULL
    mypco$Map_pco<-NA
  }else{
    Map<-Map[match(colnames(Tab),Map$taxon_oid),]
    #Use canberra function
    distfun<-function(x)vegdist(x = x,method = "canberra")
    mypco<-tryCatch(oh.pco(Tab = t(Tab),Map = Map,distfun = distfun,id_var = "taxon_oid"),error = function(e) list(Map_pco = NA))
  }
  if(is.na(mypco$Map_pco)){
    lista=list(pco_res=NA,glm_res=NA,pco_plot=NA,roccurve=NA,summary=NA)
    return(lista)
    }else{
    #Perform glm using the two first components
    m1<-glm(formula = Classification~PCo1 + PCo2,family = binomial(link = 'logit'),data = mypco$Map_pco, maxit = 1000)
    slope <- coef(m1)["PCo1"] / (-coef(m1)["PCo2"])
    intercept <- coef(m1)["(Intercept)"] / (-coef(m1)["PCo2"])
    #Create the classification plot
    p<-oh.plot.pco(list_ohpco = mypco,col_val = "Classification",
                   comp_a = "PCo1",comp_b = "PCo2",mypch = 21,size = 10) + scale_fill_manual(values = pal_panpa)
    p<-p + geom_abline(intercept = intercept,slope = slope,color="#414141",size=1)
    #RocCurve
    predpr<-predict(m1,type = c("response"))
    roccurve<-roc(mypco$Map_pco$Classification ~ predpr,percent=F,ci=F,plot=F)
    sum<-data.frame(num_genomes=nrow(mypco$Map_pco),aic=m1$aic,auc=roccurve$auc,num_orthogroups=nrow(Tab))
    lista=list(pco_res=mypco,glm_res=m1,pco_plot=p,roccurve=roccurve,summary=sum)
    return(lista)
  }
}

#Perfrom PCOA enrichment random
wrapper.pcoa.glm.enrichment.pa.npa.control<-function(file=NULL,orthogroups=NULL){
  different_names<-orthogroups
  Tab<-read.table(file,header=T,row.names = 1,sep="\t",check.names = F)
  Map<-read.table("../metadata_3837_genomes_june2016_complete.tsv",header = T,sep="\t")
  Map<-Map[match(colnames(Tab),Map$taxon_oid),]
  Map<-Map[which(Map$Classification!="soil"),]
  Map<-droplevels(Map)
  Tab<-Tab[,match(Map$taxon_oid,colnames(Tab))]
  Tab<-Tab[which(rowSums(Tab)!=0),]
  #Select random genes
  taxon<-unique(gsub(pattern = "_.*",replacement = "",x = rownames(Tab),perl = T))
  different_names<-grep(pattern = taxon,x = different_names,value = T)
  chosen_rows<-sample(x = rownames(Tab),size = length(different_names),replace = T)
  Tab<-Tab[match(chosen_rows,rownames(Tab)),]
  Tab<-Tab[,which(colSums(Tab)!=0)]
  Map<-Map[match(colnames(Tab),Map$taxon_oid),]
  distfun<-function(x)vegdist(x = x,method = "canberra")
  mypco<-oh.pco(Tab = t(Tab),Map = Map,distfun = distfun,id_var = "taxon_oid")
  #Perform glm using the two first components
  m1<-glm(formula = Classification~PCo1 + PCo2,family = binomial(link = 'logit'),data = mypco$Map_pco)
  slope <- coef(m1)["PCo1"] / (-coef(m1)["PCo2"])
  intercept <- coef(m1)["(Intercept)"] / (-coef(m1)["PCo2"])
  #Create the classification plot
  p<-oh.plot.pco(list_ohpco = mypco,col_val = "Classification",
  comp_a = "PCo1",comp_b = "PCo2",mypch = 21,size = 10) + scale_fill_manual(values = pal_panpa)
  p<-p + geom_abline(intercept = intercept,slope = slope,color="#414141",size=1)
  predpr<-predict(m1,type = c("response"))
  roccurve<-roc(mypco$Map_pco$Classification ~ predpr,percent=T,ci=T,plot=F)
  sum<-data.frame(num_genomes=nrow(mypco$Map_pco),aic=m1$aic,auc=roccurve$auc,num_orthogroups=nrow(Tab))
  lista=list(pco_res=mypco,glm_res=m1,pco_plot=p,roccurve=roccurve,summary=sum)
  return(lista)
}

##Subsampling 20 % of the orthogroups at random in the selected orthogroups
wrapper.pcoa.glm.enrichment.pa.npa.control.subsample<-function(file=NULL,orthogroups=NULL){
  different_names<-orthogroups
  Tab<-read.table(file,header=T,row.names = 1,sep="\t",check.names = F)
  Map<-read.table("../metadata_3837_genomes_june2016_complete.tsv",header = T,sep="\t")
  Map<-Map[match(colnames(Tab),Map$taxon_oid),]
  Map<-Map[which(Map$Classification!="soil"),]
  Map<-droplevels(Map)
  Tab<-Tab[,match(Map$taxon_oid,colnames(Tab))]
  Tab<-Tab[which(rowSums(Tab)!=0),]
  Tab<-Tab[rownames(Tab)%in%different_names,]
  #Select random genes
  chosen_rows<-sample(x =rownames(Tab) ,size = ceiling(nrow(Tab)*0.4),replace = T)
  Tab<-Tab[match(chosen_rows,rownames(Tab)),]
  Tab<-Tab[,which(colSums(Tab)!=0)]
  Map<-Map[match(colnames(Tab),Map$taxon_oid),]
  distfun<-function(x)vegdist(x = x,method = "canberra")
  mypco<-oh.pco(Tab = t(Tab),Map = Map,distfun = distfun,id_var = "taxon_oid")
  #Perform glm using the two first components
  m1<-glm(formula = Classification~PCo1 + PCo2,family = binomial(link = 'logit'),data = mypco$Map_pco,maxit=1000)
  slope <- coef(m1)["PCo1"] / (-coef(m1)["PCo2"])
  intercept <- coef(m1)["(Intercept)"] / (-coef(m1)["PCo2"])
  #Create the classification plot
  p<-oh.plot.pco(list_ohpco = mypco,col_val = "Classification",
                 comp_a = "PCo1",comp_b = "PCo2",mypch = 21,size = 10) + scale_fill_manual(values = pal_panpa)
  p<-p + geom_abline(intercept = intercept,slope = slope,color="#414141",size=1)
  predpr<-predict(m1,type = c("response"))
  roccurve<-roc(mypco$Map_pco$Classification ~ predpr,percent=F,ci=F,plot=F)
  sum<-data.frame(num_genomes=nrow(mypco$Map_pco),aic=m1$aic,auc=roccurve$auc,num_orthogroups=nrow(Tab))
  lista=list(pco_res=mypco,glm_res=m1,pco_plot=p,roccurve=roccurve,summary=sum)
  return(lista)
}

wrapper.pcoa.glm.enrichment.ra.soil<-function(file=NULL,orthogroups=NULL){
  different_names<-orthogroups
  Tab<-read.table(file,header=T,row.names = 1,sep="\t",check.names = F)
  Map<-read.table("../metadata_3837_genomes_june2016_complete.tsv",header = T,sep="\t")
  Map<-Map[match(colnames(Tab),Map$taxon_oid),]
  Map_soil<-as.character(Map$taxon_oid[which(Map$Classification=="soil")])
  Map_ra<-as.character(Map$taxon_oid[which(Map$Root.NotRoot=="RA")])
  myids<-c(Map_soil,Map_ra)
  Map<-Map[match(myids,Map$taxon_oid),]
  Map<-droplevels(Map)
  Tab<-Tab[,match(Map$taxon_oid,colnames(Tab))]
  Tab<-Tab[which(rowSums(Tab)!=0),]
  Tab<-Tab[which(rownames(Tab)%in%different_names),]
  #Removegenomes that go to zero
  Tab<-Tab[,which(colSums(Tab)!=0)]
  Map<-Map[match(colnames(Tab),Map$taxon_oid),]
  #Use canberra function
  distfun<-function(x)vegdist(x = x,method = "canberra")
  mypco<-tryCatch(oh.pco(Tab = t(Tab),Map = Map,distfun = distfun,id_var = "taxon_oid"),error = function(e) list(Map_pco = NA))
  if(is.na(mypco$Map_pco)){
    lista=list(pco_res=NA,glm_res=NA,pco_plot=NA,roccurve=NA,summary=NA)
    return(lista)
  }else{
    #Perform glm using the two first components
    m1<-glm(formula = Classification~PCo1 + PCo2,family = binomial(link = 'logit'),data = mypco$Map_pco, maxit = 1000)
    slope <- coef(m1)["PCo1"] / (-coef(m1)["PCo2"])
    intercept <- coef(m1)["(Intercept)"] / (-coef(m1)["PCo2"])
    #Create the classification plot
    p<-oh.plot.pco(list_ohpco = mypco,col_val = "Root.NotRoot",
                   comp_a = "PCo1",comp_b = "PCo2",mypch = 21,size = 10) + scale_fill_manual(values = pal_rasoil)
    p<-p + geom_abline(intercept = intercept,slope = slope,color="#414141",size=1)
    #RocCurve
    predpr<-predict(m1,type = c("response"))
    roccurve<-roc(mypco$Map_pco$Classification ~ predpr,percent=F,ci=F,plot=F)
    sum<-data.frame(num_genomes=nrow(mypco$Map_pco),aic=m1$aic,auc=roccurve$auc,num_orthogroups=nrow(Tab))
    lista=list(pco_res=mypco,glm_res=m1,pco_plot=p,roccurve=roccurve,summary=sum)
    return(lista)
  }
}

wrapper.pcoa.glm.enrichment.ra.soil.control<-function(file=NULL,orthogroups=NULL){
  different_names<-orthogroups
  Tab<-read.table(file,header=T,row.names = 1,sep="\t",check.names = F)
  Map<-read.table("../metadata_3837_genomes_june2016_complete.tsv",header = T,sep="\t")
  Map<-Map[match(colnames(Tab),Map$taxon_oid),]
  Map_soil<-as.character(Map$taxon_oid[which(Map$Classification=="soil")])
  Map_ra<-as.character(Map$taxon_oid[which(Map$Root.NotRoot=="RA")])
  myids<-c(Map_soil,Map_ra)
  Map<-Map[match(myids,Map$taxon_oid),]
  Map<-droplevels(Map)
  Tab<-Tab[,match(Map$taxon_oid,colnames(Tab))]
  Tab<-Tab[which(rowSums(Tab)!=0),]
  #Select random genes
  taxon<-unique(gsub(pattern = "_.*",replacement = "",x = rownames(Tab),perl = T))
  different_names<-grep(pattern = taxon,x = different_names,value = T)
  chosen_rows<-sample(x = rownames(Tab),size = length(different_names),replace = T)
  Tab<-Tab[match(chosen_rows,rownames(Tab)),]
  Tab<-Tab[,which(colSums(Tab)!=0)]
  Map<-Map[match(colnames(Tab),Map$taxon_oid),]
  distfun<-function(x)vegdist(x = x,method = "canberra")
  mypco<-oh.pco(Tab = t(Tab),Map = Map,distfun = distfun,id_var = "taxon_oid")
  #Perform glm using the two first components
  m1<-glm(formula = Classification~PCo1 + PCo2,family = binomial(link = 'logit'),data = mypco$Map_pco)
  slope <- coef(m1)["PCo1"] / (-coef(m1)["PCo2"])
  intercept <- coef(m1)["(Intercept)"] / (-coef(m1)["PCo2"])
  #Create the classification plot
  p<-oh.plot.pco(list_ohpco = mypco,col_val = "Root.NotRoot",
                 comp_a = "PCo1",comp_b = "PCo2",mypch = 21,size = 10) + scale_fill_manual(values = pal_rasoil)
  p<-p + geom_abline(intercept = intercept,slope = slope,color="#414141",size=1)
  predpr<-predict(m1,type = c("response"))
  roccurve<-roc(mypco$Map_pco$Classification ~ predpr,percent=T,ci=T,plot=F)
  sum<-data.frame(num_genomes=nrow(mypco$Map_pco),aic=m1$aic,auc=roccurve$auc,num_orthogroups=nrow(Tab))
  lista=list(pco_res=mypco,glm_res=m1,pco_plot=p,roccurve=roccurve,summary=sum)
  return(lista)
}
 
#Resample the function x given times
wrapper.boot.pcoa.glm.control<-function(file=NULL,orthogroups=NULL,n=1000,mfunction){
  res<-NULL
  for(i in 1:n){
    print(i)
    res_lista<-mfunction(file=file,orthogroups=orthogroups)
    isum<-res_lista$summary
    res<-rbind(res,isum)
  }
  return(res)
}

 heatmap.ggtree.binary.pa.npa<-function(file=NULL,orthogroups=NULL,tree=tree){
   different_names<-orthogroups
   Tab<-read.table(file,header=T,row.names = 1,sep="\t",check.names = F)
   Map<-read.table("../metadata_3837_genomes_june2016_complete.tsv",header = T,sep="\t")
   Map<-Map[match(colnames(Tab),Map$taxon_oid),]
   Map<-Map[which(Map$Classification!="soil"),]
   Map<-droplevels(Map)
   Tab<-Tab[,match(Map$taxon_oid,colnames(Tab))]
   Tab<-Tab[which(rowSums(Tab)!=0),]
   Tab<-Tab[which(rownames(Tab)%in%different_names),]
   #Remove 6 genomes that go to zero
   Tab<-Tab[,which(colSums(Tab)!=0)]
   Map<-Map[match(colnames(Tab),Map$taxon_oid),]
   #Subset the tree
   subtree<-drop.tip(phy = tree,tip = tree$tip.label[which(!(tree$tip.label%in%Map$taxon_oid))])
   binary<-Tab
   binary[binary>1]<-1
   p <- ggtree(subtree)
   p<-p %<+% Map + geom_tippoint(aes(size=3,shape=Classification, color=Classification), alpha=0.6) + scale_color_manual(values = pal_panpa)
   p<-gheatmap(p = p,data = t(binary),  colnames=F,low = "lightgrey",high = "black")
   return(p)
}

##Function to wrap and compare individual methods plus combinations of two 
wrapper.pcoa.glm.enrichment.methods<-function(myfun=NULL,list_comparisons=NULL,filename=NULL){
  #Create a df to contain the results
  res<-NULL
  #Test singleton tests
  for(method in list_comparisons){
    print(method)
    myres<-myfun(file = filename,orthogroups = get(method))
    df<-myres$summary
    print(df)
    #Append the summary to the dataframe object
    res<-rbind(res,df)
  }
  #Add the name to the rows
  rownames(res)<-list_comparisons
  selected_tests<-list_comparisons
  for(j in 2:2){
    cblock<-combn(selected_tests,j)
    #Loop over the combination block
    for(k in 1:ncol(cblock)){
      temp_c<-cblock[,k]
      names_and_comb<-paste(temp_c[1],sep = "")
      names_or_comb<-paste(temp_c[1],sep = "")
      vec_and<-get(temp_c[1])
      vec_or<-get(temp_c[1])
      #In this part i should create the combination structure
      for(m in 2:length(temp_c)){
        names_and_comb<-paste(names_and_comb,"_and_",temp_c[m],sep = "")
        names_or_comb<-paste(names_or_comb,"_or_",temp_c[m],sep = "")
        vec_and<-intersect(vec_and,get(temp_c[m]))
        vec_or<-union(vec_or,get(temp_c[m]))

      }
      #Perfrom the wrapper call
      list_comparisons<-c(list_comparisons,names_and_comb,names_or_comb)
      myres<-wrapper.pcoa.glm.enrichment.pa.npa(file = filename,orthogroups = vec_and)
      df<-myres$summary
      print(df)
      res<-rbind(res,df)
      myres<-wrapper.pcoa.glm.enrichment.pa.npa(file = filename,orthogroups = vec_or)
      df<-myres$summary
      print(df)
      res<-rbind(res,df)
    }
  }
  rownames(res)<-list_comparisons
  return(res)
}
