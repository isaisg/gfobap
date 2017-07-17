##Script to benchmark the orthogroups based on phyla apmhorta output
#Load libraries
library(reshape2)
library(ggplot2)
library(gridExtra)

#Color palette to utilize 
paleta<-c("#2988BC","#ED8C72")

setwd('/home/isai/Documents/results/3837_pa_npa_soil_june_2016_resubmit/phyla_amphora_final/')

##Functions###
#Create a function to take a file and compute a melted structure 
create_melted<-function(filename=NULL,unique_name=NULL,general_name=NULL,specific=NULL){
  temp<-read.table(filename,sep="\t")
  colnames(temp)<-c("Marker","Total_SC","Total_Found","Total_Missing","Purity","Fragmentation")
  #melted<-melt(temp,id.vars = c(1,5,6),measure.vars = NULL)
  melted<-melt(temp,id.vars=c(1),measure.vars = c(5,6))
  type<-rep(unique_name,nrow(melted))
  melted$Type<-type
  type_specific<-paste(melted$Marker,melted$variable,general_name,sep="-")
  melted$Marker_specific<-type_specific
  #Subset the markers to cover only specific taxon markers
  melted<-droplevels(melted[grep(pattern = specific,x = melted$Marker),])
  return(melted)
}

#Create a function to save the figures in svg format
oh.ggsave.svg<-function(ggobject=p,outname="",outdir="figures/",width=20,height=15,device="svg",dpi=600){
  dir.create(outdir, showWarnings = FALSE)
  myfilename<-paste(outdir,outname,sep="")
  ggsave(file = myfilename,plot = ggobject,width =width,height = height,device = device,dpi = dpi)
}

#Alphaproteobacteria
uclust_alpha<-create_melted("phylo_amphora_uclust_Alphaproteobacteria.tsv",unique_name = "Uclust_Alphaproteobacteria",general_name = "Alphaproteobacteria",specific="Alpha")
of_alpha<-create_melted(filename = "phylo_amphora_orthofinder_Alphaproteobacteria.tsv",unique_name = "OrthoFinder_Alphaproteobacteria",general_name = "Alphaproteobacteria",specific="Alpha")

#Bacteroidetes
uclust_bacteroidetes<-create_melted("phylo_amphora_uclust_Bacteroidetes.tsv",unique_name = "Uclust_Bacteroidetes",general_name = "Bacteroidetes",specific = "Bacteroidetes")
of_bacteroidetes<-create_melted("phylo_amphora_orthofinder_Bacteroidetes.tsv",unique_name = "OrthoFinder_Bacteroidetes",general_name = "Bacteroidetes",specific="Bacteroidetes")

#Burkholderiales
uclust_burkhol<-create_melted("phylo_amphora_uclust_Burkholderiales.tsv",unique_name = "Uclust_Burkholderiales",general_name = "Burkholderiales",specific="Beta")
of_burkhol<-create_melted("phylo_amphora_orthofinder_Burkholderiales.tsv",unique_name = "OrthoFinder_Burkholderiales",general_name = "Burkholderiales",specific="Beta")

#Xantho
uclust_xantho<-create_melted("phylo_amphora_uclust_Xanthomonadaceae.tsv",unique_name = "Uclust_Xanthomonadaceae",general_name = "Xanthomonadaceae",specific="Gamma")
of_xantho<-create_melted("phylo_amphora_orthofinder_Xanthomonadaceae.tsv",unique_name = "OrthoFinder_Xanthomonadaceae",general_name = "Xanthomonadaceae",specific="Gamma")

#Bacilales
uclust_baci<-create_melted("phylo_amphora_uclust_Bacillales.tsv",unique_name="Uclust_Bacillales",general_name = "Bacillales",specific = "Firmicutes")
of_baci<-create_melted("phylo_amphora_orthofinder_Bacillales.tsv",unique_name="OrthoFinder_Bacillales",general_name = "Bacillales",specific="Firmicutes")

#Pseudomonas
uclust_pseudomonas<-create_melted("phylo_amphora_uclust_Pseudomonas.tsv",unique_name = "Uclust_Pseudomonas",general_name = "Pseudomonas",specific="Gamma")
of_pseudomonas<-create_melted("phylo_amphora_orthofinder_Pseudomonas.tsv",unique_name = "OrthoFinder_Pseudomonas",general_name = "Pseudomonas",specific="Gamma")

#Actino one
uclust_actino_one<-create_melted("phylo_amphora_uclust_Actinobacteria_one.tsv",unique_name = "Uclust_Actinobacteria_one",general_name = "Actinobacteria_one",specific = "Actino")
of_actino_one<-create_melted("phylo_amphora_orthofinder_Actinobacteria_one.tsv",unique_name = "OrthoFinder_Actinobacteria_one",general_name = "Actinobacteria_one",specific="Actino")


#Actino two
uclust_actino_two<-create_melted("phylo_amphora_uclust_Actinobacteria_two.tsv",unique_name = "Uclust_Actinobacteria_two",general_name = "Actinobacteria_two",specific="Actino")
of_actino_two<-create_melted("phylo_amphora_orthofinder_Actinobacteria_two.tsv",unique_name = "OrthoFinder_Actinobacteria_two",general_name = "Actinobacteria_two",specific="Actino")


##Acinetobacter
uclust_acinetobacter<-create_melted("phylo_amphora_uclust_Acinetobacter.tsv",unique_name = "Uclust_Acinetobacter",general_name = "Acinetobacter",specific = "Gamma")
of_acinetobacter<-create_melted("phylo_amphora_orthofinder_Acinetobacter.tsv",unique_name = "OrthoFinder_Acinetobacter",general_name = "Acinetobacter",specific="Gamma")


##Merge all the dataframes
all_uclust<-rbind(uclust_acinetobacter,uclust_actino_one,uclust_actino_two,uclust_alpha,uclust_baci,uclust_bacteroidetes,uclust_burkhol,uclust_pseudomonas,uclust_xantho)
all_of<-rbind(of_acinetobacter,of_actino_one,of_actino_two,of_alpha,of_baci,of_bacteroidetes,of_burkhol,of_pseudomonas,of_xantho)
all<-rbind(all_uclust,all_of)

all$Method<-factor(gsub(pattern = "_.*",replacement = "",x = all$Type))
all$Taxon<-factor(gsub(pattern = "two",replacement = "Actinobacteria_two",x = gsub(pattern ="one",replacement = "Actinobacteria_one",x = gsub(pattern = "[A-z]+?_",replacement = "",x = all$Type,perl=T))))

all_num<-subset(all,variable=="Fragmentation")
all_prop<-subset(all,variable!="Fragmentation")
all_num$value<-1/all_num$value

all<-rbind(all_prop,all_num)

all$Taxon<-factor(gsub("Actinobacteria_two",replacement = "Actinobacteria 2",x = gsub("Actinobacteria_one",replacement = "Actinobacteria 1",x = all$Taxon)))

p<-ggplot(data = all,aes(x=Marker,y=value,col=Method)) + 
  facet_wrap(~Taxon + variable,scales = "free",ncol = 2,labeller = label_wrap_gen(multi_line=FALSE,width = 100)) +
   geom_point()
p<-p + scale_fill_manual(values = paleta) + 
  theme(axis.line = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour =   "#D9D9D9"),
        panel.grid.minor = element_line(colour = "#D9D9D9"),
        panel.border = element_rect(fill=NA,color =  "#414141",size = 1),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour = "black",size = 1.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(family = "AvantGarde",face="plain",size=15,colour="#414141"),
        axis.title.x = element_text(family = "AvantGarde",face="plain",size = 30,colour = "#414141"),
        axis.title.y = element_blank(),
        legend.background = element_blank(),legend.key.size = unit(2,"line"),
        legend.title=element_blank(),legend.key = element_blank(), 
        legend.text = element_text(size=25,family = "AvantGarde",face = "plain",colour = "#414141"),
        legend.position ="right",strip.text = element_text(family = "AvantGarde",colour = "#414141",size = 20),
        strip.background = element_rect(fill = "#D9D9D9",color = "#414141"))



oh.ggsave.svg(ggobject = p,outname = "phyla_amphora_benchmark.svg")

#Create table of markers for legend
all_o<-subset(all,Method=="OrthoFinder")
df<-data.frame(table(all_o$Taxon))
colnames(df)<-c("Taxon","Number of Markers")
grid.table(df,rows=NULL)

##Try another figure
all_sum<-aggregate(value ~ variable + variable + Method+Taxon,all,mean)
colnames(all_sum)[4]<-"Index"
p<-ggplot(data = all_sum,aes(x=variable,y = Index,fill=Method)) + geom_bar(stat = "identity",position = "dodge") +
  facet_wrap(~Taxon) + scale_fill_manual(values = paleta) + 
  theme(axis.line = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour =   "#D9D9D9"),
        panel.grid.minor = element_line(colour = "#D9D9D9"),
        panel.border = element_rect(fill=NA,color =  "#414141",size = 1),
        axis.ticks.x =  element_line(colour = "black",size = 1.5),
        axis.ticks.y = element_line(colour = "black",size = 1.5),
        axis.text.x = element_text(family = "AvantGarde",face="plain",size=20,colour="#414141"),
        axis.text.y = element_text(family = "AvantGarde",face="plain",size=15,colour="#414141"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family = "AvantGarde",face="plain",size=20,colour="#414141"),
        legend.background = element_blank(),legend.key.size = unit(2,"line"),
        legend.title=element_blank(),legend.key = element_blank(), 
        legend.text = element_text(size=25,family = "AvantGarde",face = "plain",colour = "#414141"),
        legend.position ="right",strip.text = element_text(family = "AvantGarde",colour = "#414141",size = 20),
        strip.background = element_rect(fill = "#D9D9D9",color = "#414141"))
oh.ggsave.svg(ggobject = p,outname = "phyla_amphora_benchmark_mean.svg")
