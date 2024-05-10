library(stringr)
library(dplyr)
library(matlib)


rm(list=ls())
#filtered based on imputation R2
SNPdata<- read.table(file="all_chr_imp_GT0.4",header=FALSE)
colnames(SNPdata)<-c("CHR","SNP","Rsq")
SNPdata$Rsq<-as.numeric(SNPdata$Rsq)
merged<-SNPdata %>% select (SNP,Rsq)


#reading  GWAS files
mydir="/scratch/user/fertility_sequenced_data"
files<-list.files(path=mydir, pattern=".fastGWA" , full.names = TRUE)


##merging and calculating the t-values for SNPs
name<-c()

for(j in 1:length(files)) {
tmpfile<-read.table(files[[j]],header=TRUE, stringsAsFactors = FALSE)

res <- str_match(files[[j]], "Output.\\s*(.*?)\\s*.fastGWA")
name[j]<-res[,2]

filtered<-merge(SNPdata,tmpfile,by="SNP")

filtered$t_value<-filtered$BETA/filtered$SE

tmp<- filtered %>% arrange (SNP) %>% select (SNP,t_value)
colnames(tmp)<-c("SNP", paste(name[j]))
merged<-merge(merged,tmp,by="SNP")
str(merged)
}

##cal correlation
mydata.cor = cor(merged[,3:6],use = "p")

inverted.cor <- inv(mydata.cor)
A<-inverted.cor

##chi square statistics
multiply <- function(x) { x %*% A %*% x}
multi_sig<-apply(as.matrix(merged[,3:6]),1, FUN =multiply)

##Pvalue estimates
pvalue_multi<-pchisq(multi_sig, df=4, lower.tail=FALSE)
merged$P<-pvalue_multi


