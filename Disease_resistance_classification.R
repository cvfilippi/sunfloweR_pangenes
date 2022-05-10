#Disease-resistance proteins classification based on domain combinations
#This R scrip requires as input a interproscan csv file (Jones et al., 2014) and a TMHMM 2.0 short output file (Krogh et al, 2001, output format: 'One line per protein')
#Both files enerated for the same protein FASTA file
#The script uses R base functions and the package data.table (https://CRAN.R-project.org/package=data.table)
setwd(..)
library(data.table)

#load the interproscan csv file
interpro_file<-read.csv("interpro_file.tsv", header=FALSE, sep="\t", na.strings = " ", fill = TRUE)
colnames(interpro_file)[1]<-"ID"
interpro_file$V9<-as.numeric(as.character(interpro_file$V9)) 
interpro_file_significant<-interpro_file[interpro_file$V9 < 0.05, ] 

#set disease-resistance characteristic domains
LRR<-c("IPR032675","PF13855","SSF52058","SM00369","SM00365","PS51450","PF13306","PF13516","PF12799","PF08263","PF07725","PF07723","PF00560", "PF14580", "IPR013210","IPR011713","IPR003591","IPR001611") ##PF12799 (2 copias) y (PF13306, PF14580, PF13855) tienen varias copias
CC<-c("Coil", "G3DSA:1.20.5.170", "SSF57959")
LysM<-c("PF01476")
NBS<-c("PF00931", "IPR002182", "PF12061")
RPW8<-c("IPR008808", "PS51153", "PF05659")
KIN<-c("SM00219", "G3DSA:3.30.200.20", "IPR000719", "IPR001245", "IPR008271", "IPR011009", "IPR017441", "IPR020635", "IPR021820", "IPR022126", "PF00069", "PF07714", "PF08488", "PS00107", "PS00108", "PS50011", "SM00219", "SM00220", "SSF56112")
TIR<-c("IPR000157", "PF01582", "PF13676","PS50104", "SM00255", "SSF52200")

##Identify proteins containing diseaserresistance characteristic domains 
#LRR
prot_with_LRR<-data.frame(unique(interpro_file_significant[interpro_file_significant$V5 %in% LRR, c("ID", "V4", "V7", "V8")])) 
prot_with_LRR$class<-"LRR" 

#CC
prot_with_CC<-data.frame(unique(interpro_file_significant[interpro_file_significant$V5 %in% CC, c("ID", "V4", "V7", "V8")])) 
prot_with_CC$class<-"CC" ##5820

#LysM
prot_with_LysM<-data.frame(unique(interpro_file_significant[interpro_file_significant$V5 %in% LysM, c("ID", "V4", "V7", "V8")])) 
prot_with_LysM$class<-"LysM" ##26

#NBS
prot_with_NBS<-data.frame(unique(interpro_file_significant[interpro_file_significant$V5 %in% NBS, c("ID", "V4","V7", "V8")]))
prot_with_NBS$class<-"NBS" ##347

#RPW8
prot_with_RPW8<-data.frame(unique(interpro_file_significant[interpro_file_significant$V5 %in% RPW8, c("ID", "V4","V7", "V8")]))
prot_with_RPW8$class<-"RPW8" ##13

#KIN
prot_with_KIN<-data.frame(unique(interpro_file_significant[interpro_file_significant$V5 %in% KIN, c("ID", "V4","V7", "V8")]))
prot_with_KIN$class<-"KIN" ##2703

#TIR
prot_with_TIR<-data.frame(unique(interpro_file_significant[interpro_file_significant$V5 %in% TIR, c("ID", "V4","V7", "V8")]))
prot_with_TIR$class<-"TIR" ##107

#TM
TMHMM_output<-read.table("TMHMM_one_line_output.txt", header=FALSE)[c(1,4,5,6)]
prot_with_TM<-TMHMM_output[TMHMM_output$V5!="PredHel=0",]
colnames(prot_with_TM)<-c("ID", "V4","V7", "V8")
prot_with_TM$class<-"TM"

#merge all proteins containing at least one disease-resistance characteristic domain in a single data frame
motif_complete<-unique(rbind(prot_with_LRR, prot_with_CC, prot_with_LysM, prot_with_NBS, prot_with_RPW8, prot_with_KIN, prot_with_TIR, prot_with_TM)[,c(1,5)])

#identify, for each protein, the combination of disease-resistance characteristic domains 
motif_DT<-setDT(motif_complete)[, lapply(.SD, function(x) toString(na.omit(x))), by = ID]

#Obtain the domain combination observed for the given proteome, as well as the number of proteins with each combination
table(motif_DT$class)
