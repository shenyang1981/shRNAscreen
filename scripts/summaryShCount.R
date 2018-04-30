#!/usr/bin/env Rscript
library("stringdist")

args = commandArgs(T);

args = args[-1];

infile = args[1]; # collapsed reads

hpfile = args[2]; # annotation of  shRNAs

outfile = args[3]; # results

distcutoff = as.numeric(args[4]); # cutoff for editing distance
distcutoff = ifelse(is.na(distcutoff),2,distcutoff);

cpunum = as.numeric(args[5]); # number of threads

batchsize = as.numeric(args[6]); # number of collapsed reads per batch, e.g. 5000

# read 19mers
myseqs = read.table(infile,header=F,sep="\t",quote='"');
colnames(myseqs) = c("lib","sequence","count");

# read hairpins
# hairpins has three columns: hpid,hpsequence,gene
hp = read.table(hpfile,header=F,sep="\t")
colnames(hp) = c("hpid","hpseq","gene");

i = 1;
allresult = NULL;
while(i < length(myseqs[,2])){
	# prepare collapsed reads
	is = i;
	ie = ifelse(i+batchsize - 1 > length(myseqs[,2]),length(myseqs[,2]),i+batchsize - 1);
	tmpmyseqs = myseqs[is:ie,"sequence"];

	# calculate editing distance for each collapsed read against candidate hairpin shRNA
	mydist = stringdistmatrix(hp[,"hpseq"],tmpmyseqs,nthread = cpunum,method='lv');
	# hit for editing distance < cutoff
	cands = which(mydist<=distcutoff,arr.ind=T);
	# put all hits into data.frame
	result = data.frame("hpid" = hp[cands[,1],"hpid"], "hpseq" = hp[cands[,1],"hpseq"], "seq" = myseqs[cands[,2] + is - 1,"sequence"],"count" = myseqs[cands[,2] + is - 1,"count"],"dist" = mydist[cands],"lib"=myseqs[cands[,2] + is - 1,"lib"],"gene"=hp[cands[,1],"gene"]);
	if(is.null(allresult)){
		allresult = result;
	} else {
		allresult = rbind(allresult,result);
	}
	i = ie+1;
	print(ie);
}

# write results
allresultUniq = allresult[order(allresult[,"seq"],allresult[,"dist"]),];
allresultUniq = allresultUniq[!duplicated(allresultUniq[,c("seq","lib")]),];

write.table(allresultUniq,file=outfile,col.names=T,row.names=F,sep="\t",quote=F);
write.table(allresult,file=paste(outfile,".all",sep=""),col.names=T,row.names=F,sep="\t",quote=F);

usedCols = colnames(allresultUniq);
colnames(myseqs) = c("lib","seq","count");
myseqHpResult = merge(myseqs[,c("seq","count","lib")],allresultUniq[,!grepl("count",usedCols) & !grepl("lib",usedCols)],by="seq",all.x=T)
myseqHpResult = myseqHpResult[order(myseqHpResult[,"count"],decreasing = T),];
write.table(myseqHpResult[,usedCols],file=paste(outfile,".allmature",sep=""),col.names=T,row.names=F,sep="\t",quote=F);
