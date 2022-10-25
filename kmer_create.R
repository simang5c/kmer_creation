suppressPackageStartupMessages({
library(Biostrings)
library(dplyr)
library(doParallel)
library(data.table)
library(seqRFLP)
library(stringr)
library(factoextra)
library(FactoMineR) })

args <- commandArgs(trailingOnly = TRUE)
kmersize <- as.integer(args[1])

#For example if the fasta file name is Aeembo.fasta the folder name will be Aembo and inside it the subsequences
subsequence<-function(fasta,folder){
  df=readBStringSet(fasta)
  #converting biostring object to dataframe
  dss2df <- function(dff) data.frame(width=width(dff), seq=as.character(dff), names=names(dff))
  #calling the above function
  df1=dss2df(df)
  df11=df1 %>% filter(width>=kmersize) %>% distinct()
  num_seq=dim(df11)[1]
  number_each=as.integer(num_seq/number_of_subfasta) # number of fasta in each sub file=total number of sequences/number of files to be created
  if(number_each<1 || num_seq==1 ){
    number_each=1
    number_of_subfasta=1
  }
  
  #CREATING A FOLDER
  
  dir.create(folder)
  df11=df11 %>% dplyr::select(names,seq) #selects only sequence identifier and sequence
  start1=1
  a=0
  for (i in 1:number_of_subfasta) {
    start1=start1+a
    end1=start1+number_each-1
    if(i<number_of_subfasta){
      test=df11 %>% dplyr::slice(start1:end1)
      filname_out=paste0(i,".fasta")
      filename1=paste0(folder,"/",filname_out)
      df.fasta = dataframe2fas(test, file=filename1)
      
    } else {
      test=df11 %>% dplyr::slice(start1:num_seq)
      filname_out=paste0(i,".fasta")
      filename1=paste0(folder,"/",filname_out)
      df.fasta = dataframe2fas(test, file=filename1)
    }
    start1=end1
    a=1
    
  }
}

#this function creates the subsequences based on the kmer length information
create_subseq <- function(rwnm,seq,kmer) {
  #rng=1:(nchar(seq)-kmer)
  rng=1:nchar(seq)
  tx=seq
  midpoint=as.integer(nchar(tx)/2)
  df_temp=data.frame()
  df_temp2=data.frame()
  subseq1=data.frame()
  c=0
  #df_temp=data.frame()
  #subseq1=data.frame()
  for (i in rng) {
    start=i
    end1=i+kmer #change according to the hexamer
    start2=midpoint + c
    end2=midpoint+kmer +c
    
    if(start < midpoint){
      seq1=substr(tx,start,end1)
      if(nchar(seq1)!=kmersize){ next }
      df_temp=data.frame(rwnm,seq1)
      colnames(df_temp)=c("seqid","kmer")
      subseq1=rbind(subseq1,df_temp)
    } else {
      next
    }
    if(end2<= nchar(tx)){
      seq2=substr(tx,start2,end2)
      if(nchar(seq2)!=kmersize){ next	}
      df_temp2=data.frame(rwnm,seq2)
      colnames(df_temp2)=c("seqid","kmer")
      subseq1=rbind(subseq1,df_temp2)
      c=c+1
    } else {
      next
    }
  }
  return(subseq1)
}



#this function will call the above function after functioning one sequence after another for each organism and sequences separately.
process_sequence<-function(df_in, kmer_len){
  fin=data.frame()
  df_fin=data.frame()
  rwname=df_in$names
  for(i in 1:length(rwname)){
    rw=rwname[i]
    DT <- data.table(df_in)
    sequence=DT[names==rw]$seq
    size=kmer_len-1
    fin=create_subseq(rw,sequence,size) #hexamer number of loops =sequence length-5
    df_fin=rbind(df_fin,fin)
     
  }
  return(df_fin)
}

#kmer generation and writing the kmer count output.
kmergeneration<-function(folder_list, kmersize){
  for (f in folder_list){
    #creating a folder path
    temp=paste0(f,"/")
    fii=list.files(path= temp, pattern="\\.fasta$")
    df_fin=data.frame()
    df_fin1=data.frame()
    df_fin2=data.frame()
    #creating kmers
    df_fin1<-foreach(i=1:length(fii),.combine = rbind) %dopar% {
      tt1=data.frame()
      fpn=paste0(temp,fii[i])
      #cat(sprintf('tasks completed: %s\n', fii[i]))
      df=readAAStringSet(fpn)
      #convertig biostring object to dataframe
      dss2df <- function(dff) data.frame(width=width(dff), seq=as.character(dff), names=names(dff))
      #callign the above function
      df1=dss2df(df)
      df11=df1 %>% distinct()
      tt1=process_sequence(df11, kmersize)
    }
    output1=paste0(temp,f,"_kmer.txt")
    df_fin1=df_fin1 %>% mutate(organism=f)
    df_fin2=df_fin1 %>% group_by(seqid) %>% dplyr::count(kmer)
    df_fin1=df_fin1 %>% group_by(organism) %>% dplyr::count(kmer) #alter this step to cluster sequence or organism
    write.table(df_fin1,output1,row.names = F,sep = "\t")
    output2=paste0(temp,f,"_sequence_kmer.csv")
    write.table(df_fin2,output2,row.names = F,sep = ",")
    
  }
}

#setting number of cores to be used for processing
numCores <- detectCores() - 2  
registerDoParallel(cores=numCores)  
cl <- makeCluster(numCores, type="FORK")

#listing genome/protein fasta files in the directory
fii=list.files(path= ".", pattern="\\.fasta$")
number_of_subfasta=50 
#main parent program
for (i in fii){
  folder=str_replace(i,"\\.fasta","")
  subsequence(i,folder) ###creating folders and splitting the fasta files to smaller chunks
}
dirs=list.dirs(path=".",full.names = TRUE, recursive = F)
folder_list <- gsub("\\./", "", dirs) 


#calling the function for generating kmers
kmergeneration(folder_list, kmersize)
stopCluster(cl)
final_out=data.frame()
for(f in folder_list){
  temp0=paste0(f,"/")
  fls=list.files(path= temp0, pattern="\\.txt$")
  temp1=paste0(temp0,fls)
  #temp1
  df=read.table(temp1,"\t",header = T)
  final_out=rbind(final_out,df)
}
otpt1=paste0("final_",kmersize,"_mer_count.txt")
write.table(final_out,otpt1,sep = "\t",row.names = F)
