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


setwd('./')
source('functions_analysis_pcoa.R')

#Read the significant orthogroups
ogs_sig<-read.table("sigOGs.txt")
selected_ogs<-as.character(ogs_sig$V4)

#Change some labels to match the matrices labels
selected_ogs<-gsub(pattern = "Actinobacteria1",replacement = "Actinobacteria_one",x = selected_ogs)
selected_ogs<-gsub(pattern = "Actinobacteria2",replacement = "Actinobacteria_two",x = selected_ogs)
ogs_sig$V4<-factor(selected_ogs)

#Define columns 
colnames(ogs_sig)<-c("Taxon","EnrichmentDirection","Method","Orthogroup")


#We have five different methods benchmarked
#hypergbin hypergcn phyloglmbin phyloglmcn scoary

#We need to PCO-GLM analysis for the 9 taxon analyzed
matrices_files<-list.files(pattern = "matrix_")

r<-NULL
r_inter<-NULL
r_var<-NULL
for(file in matrices_files){
  suffix<-file
  suffix<-gsub(pattern = "\\.tsv",replacement = "",x = gsub(pattern = "matrix_raw_orthogroups_",replacement = "",x = suffix))
  for(met in levels(ogs_sig$Method)){
    myres<-wrapper.pcoa.glm.enrichment.pa.npa(file = file,
                                              orthogroups = subset(ogs_sig,Method==met)[,4])
    if(is.na(myres$pco_res)){
      r1<-data.frame(Classification=NA,Taxonomic_Group=suffix,PCo1=NA,PCo2=NA,PCo3=NA,Method=met)
      r<-rbind(r,r1)
    }else{
      temp_Map<-myres$pco_res$Map_pco
      temp_Map$Method<-rep(met,nrow(temp_Map))
      r1<-temp_Map[,c(3,5,12:15)]
      r<-rbind(r,r1)
      #Deal with the variance
      r1_var<-data.frame(PCo1=as.vector(myres$pco_res$variance_explained[1]),PCo2=as.vector(myres$pco_res$variance_explained[2]),
                         PCo3=as.vector(myres$pco_res$variance_explained[3]))
      r1_var$Taxonomic_Group<-as.vector(unique(r1[,2]))
      r1_var$Method<-met
      r_var<-rbind(r_var,r1_var)
      #Deal with the glm result
      data.frame(r1_var)
      m1<-myres$glm_res
      myaic<-AIC(m1)
      slope <- coef(m1)["PCo1"] / (-coef(m1)["PCo2"])
      intercept <- coef(m1)["(Intercept)"] / (-coef(m1)["PCo2"])
      r1_inter<-data.frame(slope=slope,intercept=intercept,Method=met,Taxonomic_Group=unique(r1[,2]),AIC=myaic)
      rownames(r1_inter)<-NULL
      r_inter<-rbind(r_inter,r1_inter)
    }
  }
}

##As graphical comparison I will run the same analysis over the whole matrix and N number of random genomes
for(file in matrices_files){
  met<-"Full_Matrix"
  suffix<-file
  suffix<-gsub(pattern = "\\.tsv",replacement = "",x = gsub(pattern = "matrix_raw_orthogroups_",replacement = "",x = suffix))
  temp_Tab<-read.table(file = file,header = T,sep = "\t",row.names = 1,check.names = F)
  ogs<-rownames(temp_Tab)
  myres<-wrapper.pcoa.glm.enrichment.pa.npa(file = file,
                                              orthogroups = ogs)

  temp_Map<-myres$pco_res$Map_pco
  temp_Map$Method<-rep(met,nrow(temp_Map))
  r1<-temp_Map[,c(3,5,12:15)]
  r<-rbind(r,r1)
  #Deal with the variance
  r1_var<-data.frame(PCo1=as.vector(myres$pco_res$variance_explained[1]),PCo2=as.vector(myres$pco_res$variance_explained[2]),
                     PCo3=as.vector(myres$pco_res$variance_explained[3]))
  r1_var$Taxonomic_Group<-as.vector(unique(r1[,2]))
  r1_var$Method<-met
  r_var<-rbind(r_var,r1_var)
  #Deal with the glm result
  data.frame(r1_var)
  m1<-myres$glm_res
  myaic<-AIC(m1)
  slope <- coef(m1)["PCo1"] / (-coef(m1)["PCo2"])
  intercept <- coef(m1)["(Intercept)"] / (-coef(m1)["PCo2"])
  r1_inter<-data.frame(slope=slope,intercept=intercept,Method=met,Taxonomic_Group=unique(r1[,2]),AIC=myaic)
  rownames(r1_inter)<-NULL
  r_inter<-rbind(r_inter,r1_inter)
}

#Subset random genes and compute the same pcoa
# df_average_ogs<-data.frame(Taxon=levels(ogs_sig$Taxon),aver=round(colMeans(table(ogs_sig$Method,ogs_sig$Taxon))))
# rownames(df_average_ogs)<-NULL
# df_average_ogs$Taxon<-factor(c("Acinetobacter","Actinobacteria_one","Actinobacteria_two","Alphaproteobacteria","Bacillales","Bacteroidetes","Burkholderiales","Pseudomonas","Xanthomonadaceae"))
# 
# for(file in matrices_files){
#   met<-"Random"
#   suffix<-file
#   suffix<-gsub(pattern = "\\.tsv",replacement = "",x = gsub(pattern = "matrix_raw_orthogroups_",replacement = "",x = suffix))
#   temp_Tab<-read.table(file = file,header = T,sep = "\t",row.names = 1,check.names = F)
#   ogs<-rownames(temp_Tab)
#   num<-subset(df_average_ogs,Taxon==suffix)[,2]
#   ogs<-sample(ogs,num)
#   myres<-wrapper.pcoa.glm.enrichment.pa.npa(file = file,
#                                             orthogroups = ogs)
#   
#   temp_Map<-myres$pco_res$Map_pco
#   temp_Map$Method<-rep(met,nrow(temp_Map))
#   r1<-temp_Map[,c(3,5,12:15)]
#   r<-rbind(r,r1)
#   #Deal with the variance
#   r1_var<-data.frame(PCo1=as.vector(myres$pco_res$variance_explained[1]),PCo2=as.vector(myres$pco_res$variance_explained[2]),
#                      PCo3=as.vector(myres$pco_res$variance_explained[3]))
#   r1_var$Taxonomic_Group<-as.vector(unique(r1[,2]))
#   r1_var$Method<-met
#   r_var<-rbind(r_var,r1_var)
#   #Deal with the glm result
#   data.frame(r1_var)
#   m1<-myres$glm_res
#   slope <- coef(m1)["PCo1"] / (-coef(m1)["PCo2"])
#   intercept <- coef(m1)["(Intercept)"] / (-coef(m1)["PCo2"])
#   r1_inter<-data.frame(slope=slope,intercept=intercept,Method=met,Taxonomic_Group=unique(r1[,2]))
#   rownames(r1_inter)<-NULL
#   r_inter<-rbind(r_inter,r1_inter)
# }

#Readjust the levels for method
r$Method<-factor(r$Method,levels = c("Full_Matrix","hypergbin","hypergcn","phyloglmbin","phyloglmcn","scoary"))
r$Taxonomic_Group<-factor(gsub(pattern = "Actinobacteria_two",replacement = "Actinobacteria 2",x = gsub(pattern = "Actinobacteria_one",replacement = "Actinobacteria 1",x = r$Taxonomic_Group)))
r_inter$Taxonomic_Group<-factor(gsub(pattern = "Actinobacteria_two",replacement = "Actinobacteria 2",x = gsub(pattern = "Actinobacteria_one",replacement = "Actinobacteria 1",x = r_inter$Taxonomic_Group)))

#Create the figure ussing ggplot
r_temp<-droplevels(subset(r,Taxonomic_Group=="Acinetobacter" | Taxonomic_Group =="Actinobacteria 1" | Taxonomic_Group=="Actinobacteria 2"))
r_var_temp<-droplevels(subset(r_var,Taxonomic_Group=="Acinetobacter" | Taxonomic_Group =="Actinobacteria 1" | Taxonomic_Group=="Actinobacteria 2" ))
r_inter_temp<-droplevels(subset(r_inter,Taxonomic_Group=="Acinetobacter" | Taxonomic_Group =="Actinobacteria 1" | Taxonomic_Group=="Actinobacteria 2"))

r_inter_temp$x<-rep(-0.55,nrow(r_inter_temp))
r_inter_temp$y<-rep(0.15,nrow(r_inter_temp))
r_inter_temp[2,7]<--0.2
r_inter_temp[11,6]<-0.2
r_inter_temp[6,6]<-0.2
r_inter_temp[6,7]<-0.2
r_inter_temp[7,6]<-0.2
r_inter_temp[7,7]<-0.2
r_inter_temp[8,6]<-0.2
r_inter_temp[8,7]<-0.5



p <- ggplot(data = r_temp,aes_string(x = "PCo1",y = "PCo2"))+ 
  geom_point(size = 4,alpha=1,pch=21,colour="#414141",stroke = 1.5,aes(fill = get("Classification")))+
  facet_wrap(~ Taxonomic_Group + Method, ncol = 6,scales = "free_y",labeller = label_wrap_gen(multi_line=T,width = 100)) +
  xlab(label = paste("PCo1")) + 
  ylab(label = paste("PCo2"))+
  # xlab(label = paste("PCo1","(","variance_explained"[which(names("variance_explained")=="PCo1")],"%)",sep="")) + 
  # ylab(label = paste("PCo2","(","variance_explained"[which(names("variance_explained")=="PCo2")],"%)",sep="")) +
  guides(fill = guide_legend(override.aes = list(size=12))) +
  theme(axis.line = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour =   "#D9D9D9"),
        panel.grid.minor = element_line(colour = "#D9D9D9"),
        panel.border = element_rect(fill=NA,color =  "#414141",size = 1),
        axis.ticks = element_line(colour = "black",size = 2.5),
        axis.text.x = element_text(family = "AvantGarde",face = "plain",size =20,colour="#414141",angle = 90),
        axis.text.y = element_text(family = "AvantGarde",face="plain",size=20,colour="#414141"),
        axis.title.x = element_text(family = "AvantGarde",face="plain",size = 30,colour = "#414141"),
        axis.title.y = element_text(family = "AvantGarde",face="plain",size=30,colour="#414141"),
        legend.background = element_blank(),legend.key.size = unit(2,"point"),
        legend.title=element_blank(),legend.key = element_blank(), 
        legend.text = element_text(size=20,
        family = "AvantGarde",face = "plain",colour = "#414141"),
        legend.position ="right",
        strip.background = element_rect(fill = "#D9D9D9",color = "#414141"),strip.text = element_text(size = 13,)) +
  scale_fill_manual(breaks=c("NPA","PA"),values = pal_panpa)
p<-p + geom_abline(aes(intercept = intercept,slope = slope),data = r_inter_temp)
p<-p + geom_label(data = r_inter_temp,mapping = aes(x=x,y=y,label=round(AIC,1)),size=7,inherit.aes = F)
oh.ggsave.svg(ggobject = p,outname = "pcoa_orthogroups_acitenobacter_actinobacteria1_actinobacteria2.svg")

#Repeat for the second part 
r_temp<-droplevels(subset(r,Taxonomic_Group=="Alphaproteobacteria" | Taxonomic_Group=="Bacillales" | Taxonomic_Group =="Bacteroidetes" ))
r_var_temp<-droplevels(subset(r_var,Taxonomic_Group=="Alphaproteobacteria" | Taxonomic_Group=="Bacillales" | Taxonomic_Group =="Bacteroidetes" ))
r_inter_temp<-droplevels(subset(r_inter,Taxonomic_Group=="Alphaproteobacteria" | Taxonomic_Group=="Bacillales" | Taxonomic_Group =="Bacteroidetes"))

r_inter_temp$x<-rep(0.55,nrow(r_inter_temp))
r_inter_temp$y<-rep(-0.1,nrow(r_inter_temp))

r_inter_temp[9,7]<-0
r_inter_temp[10,7]<--0.3
r_inter_temp[11,7]<-0
r_inter_temp[3,7]<-0
r_inter_temp[5,7]<--0.25
r_inter_temp[15,6]<--0.4
r_inter_temp[15,7]<--0.15
r_inter_temp[6,6]<--0.42
r_inter_temp[6,7]<--0.3
r_inter_temp[8,7]<--0.22
r_inter_temp[9,6]<--0.4
r_inter_temp[9,7]<--0.35
r_inter_temp[11,6]<--0.4
r_inter_temp[11,7]<-0.35




p <- ggplot(data = r_temp,aes_string(x = "PCo1",y = "PCo2"))+
  geom_point(size = 4,alpha=1,pch=21,colour="#414141",stroke = 1.5,aes(fill = get("Classification")))+
  facet_wrap(~ Taxonomic_Group + Method, ncol = 6,scales = "free_y",labeller = label_wrap_gen(multi_line=T,width = 100)) +
  xlab(label = paste("PCo1")) + 
  ylab(label = paste("PCo2"))+
  # xlab(label = paste("PCo1","(","variance_explained"[which(names("variance_explained")=="PCo1")],"%)",sep="")) + 
  # ylab(label = paste("PCo2","(","variance_explained"[which(names("variance_explained")=="PCo2")],"%)",sep="")) +
  guides(fill = guide_legend(override.aes = list(size=12))) +
  theme(axis.line = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour =   "#D9D9D9"),
        panel.grid.minor = element_line(colour = "#D9D9D9"),
        panel.border = element_rect(fill=NA,color =  "#414141",size = 1),
        axis.ticks = element_line(colour = "black",size = 2.5),
        axis.text.x = element_text(family = "AvantGarde",face = "plain",size =20,colour="#414141",angle=90),
        axis.text.y = element_text(family = "AvantGarde",face="plain",size=20,colour="#414141"),
        axis.title.x = element_text(family = "AvantGarde",face="plain",size = 30,colour = "#414141"),
        axis.title.y = element_text(family = "AvantGarde",face="plain",size=30,colour="#414141"),
        legend.background = element_blank(),legend.key.size = unit(2,"point"),
        legend.title=element_blank(),legend.key = element_blank(), 
        legend.text = element_text(size=20,
                                   family = "AvantGarde",face = "plain",colour = "#414141"),
        legend.position ="right",
        strip.background = element_rect(fill = "#D9D9D9",color = "#414141"),strip.text = element_text(size = 13,)) +
  scale_fill_manual(breaks=c("NPA","PA"),values = pal_panpa)
p<-p + geom_abline(aes(intercept = intercept,slope = slope),data = r_inter_temp)
p<-p + geom_label(data = r_inter_temp,mapping = aes(x=x,y=y,label=round(AIC,1)),size=7,inherit.aes = F)

oh.ggsave.svg(ggobject = p,outname = "pcoa_orthogroups_alphaproteobacteria_bacillales_bacteroidetes.svg")

#Repeat for the third part 
r_temp<-droplevels(subset(r,Taxonomic_Group=="Burkholderiales" | Taxonomic_Group=="Pseudomonas" | Taxonomic_Group=="Xanthomonadaceae"))
r_var_temp<-droplevels(subset(r_var,Taxonomic_Group=="Burkholderiales" | Taxonomic_Group=="Pseudomonas" | Taxonomic_Group=="Xanthomonadaceae" ))
r_inter_temp<-droplevels(subset(r_inter, Taxonomic_Group=="Burkholderiales" | Taxonomic_Group=="Pseudomonas" | Taxonomic_Group=="Xanthomonadaceae"))

r_inter_temp$x<-rep(0.5,nrow(r_inter_temp))
r_inter_temp$y<-rep(-0.2,nrow(r_inter_temp))

r_inter_temp[2,7]<-0.1
r_inter_temp[3,7]<--0.25
r_inter_temp[8,6]<-0.1
r_inter_temp[8,7]<--0.3




p <- ggplot(data = r_temp,aes_string(x = "PCo1",y = "PCo2"))+
  geom_point(size = 4,alpha=1,pch=21,colour="#414141",stroke = 1.5,aes(fill = get("Classification")))+
  facet_wrap(~ Taxonomic_Group + Method, ncol = 6,scales = "free_y",labeller = label_wrap_gen(multi_line=T,width = 100)) +
  xlab(label = paste("PCo1")) + 
  ylab(label = paste("PCo2"))+
  # xlab(label = paste("PCo1","(","variance_explained"[which(names("variance_explained")=="PCo1")],"%)",sep="")) + 
  # ylab(label = paste("PCo2","(","variance_explained"[which(names("variance_explained")=="PCo2")],"%)",sep="")) +
  guides(fill = guide_legend(override.aes = list(size=12))) +
  theme(axis.line = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour =   "#D9D9D9"),
        panel.grid.minor = element_line(colour = "#D9D9D9"),
        panel.border = element_rect(fill=NA,color =  "#414141",size = 1),
        axis.ticks = element_line(colour = "black",size = 2.5),
        axis.text.x = element_text(family = "AvantGarde",face = "plain",size =20,colour="#414141",angle=90),
        axis.text.y = element_text(family = "AvantGarde",face="plain",size=20,colour="#414141"),
        axis.title.x = element_text(family = "AvantGarde",face="plain",size = 30,colour = "#414141"),
        axis.title.y = element_text(family = "AvantGarde",face="plain",size=30,colour="#414141"),
        legend.background = element_blank(),legend.key.size = unit(2,"point"),
        legend.title=element_blank(),legend.key = element_blank(), 
        legend.text = element_text(size=20,
                                   family = "AvantGarde",face = "plain",colour = "#414141"),
        legend.position ="right",
        strip.background = element_rect(fill = "#D9D9D9",color = "#414141"),strip.text = element_text(size = 13,)) +
  scale_fill_manual(breaks=c("NPA","PA"),values = pal_panpa)
p<-p + geom_abline(aes(intercept = intercept,slope = slope),data = r_inter_temp)
p<-p + geom_label(data = r_inter_temp,mapping = aes(x=x,y=y,label=round(AIC,1)),size=7,inherit.aes = F)

oh.ggsave.svg(ggobject = p,outname = "pcoa_orthogroups_burkholderiales_pseudomonas_xanthomonadaceae.svg")


