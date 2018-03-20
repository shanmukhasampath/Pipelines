options(stringsAsFactors=FALSE)
plinkdir <- "~/Software/plink2/plink"
gatk <- '~/Software/GATK-3.8/GenomeAnalysisTK.jar'
reffasta <- '../../human_g1k_v37.fasta'
annovar <- '~/Software/annovar'
shapeit <- '~/Software/shapeit/bin/shapeit'
ceudir <- '../data/CEU'

ifile <- 'variant_calls.vcf.gz'

# Filter the variants by Chromsome as per the Quality Filters
for(i in 1:22){
  # Filter the Vairants
	system(sprintf('java -jar %s -T SelectVariants -R %s -o chrvcf/test%s.vcf.gz --variant %s -L %s -select "AB > 0.75 || QD < 5 || SB > -0.10 || QUAL < 30 || DP < 10" -invertSelect',gatk,reffasta,i,ifile,i))
  
  # Remove Vairants that lie in PseudoGenes and Output to PLINK format
	system(sprintf('%s --vcf chrvcf/test%s.vcf.gz --exclude range %s/PseudoList_hg19 --double-id --allow-no-sex --make-bed --out chrvcf/test_chr%s --biallelic-only strict --snps-only just-acgt',plinkdir,i,psuedolistdir,i))
  
  # Rename the variant IDs
	data <- read.table(sprintf('chrvcf/test_chr%s.bim',i),h=F,stringsAsFactors=F)
	data$V2 <- as.vector(paste(data$V1,data$V4,sep=':'))
	write.table(data,sprintf('chrvcf/test_chr%s.bim',i),quote=F,row.names=F,col.names=F)
  
  # Update the sample information with Pedigree
	data <- read.table(sprintf('chrvcf/test_chr%s.fam',i),h=F,stringsAsFactors=F)
	sampleinfo <- read.table('../../pedigree.txt',h=T,sep='\t',stringsAsFactors=F)
	rownames(sampleinfo) <- as.vector(sampleinfo$SampleID)
	rownames(data) <- as.vector(data$V2)
	temp <- sampleinfo[rownames(data),c(2,3,4,5,6,7)]
	write.table(temp,sprintf('chrvcf/test_chr%s.fam',i),quote=F,row.names=F,col.names=F)
	rm(data,sampleinfo,temp)
}

# Merge PLINK files for all chromsomes into One file
beddat <- list.files('chrvcf',pattern='test_chr[0-9]*.bed',full.names=T)
bimdat <- list.files('chrvcf',pattern='test_chr[0-9]*.bim',full.names=T)
famdat <- list.files('chrvcf',pattern='test_chr[0-9]*.fam',full.names=T)
flist <- cbind(beddat,cbind(bimdat,famdat))
temp <- as.numeric(gsub('^(.*?)test_chr([0-9]*?)\\.bed$','\\2',flist[,1],perl=T))
flist <- flist[order(temp),]
write.table(flist[2:nrow(flist),],'chrvcf/allfiles.txt',quote=F,sep=' ',row.names=F,col.names=F)
system(sprintf('%s --bfile chrvcf/test_chr1 --merge-list chrvcf/allfiles.txt --keep-allele-order --allow-no-sex --make-bed --recode transpose --out chrvcf/testf --snps-only just-acgt',plinkdir))
data <- read.table('chrvcf/testf.fam',h=F,stringsAsFactors=F)
write.table(data[grep('\\.s[0-9]$',data$V2,perl=T),c(1,2,1)],'chrvcf/removesamples.txt',quote=F,row.names=F,col.names=F)

# Clean the Variants and Samples that show duplicates and 
system(sprintf('%s --bfile chrvcf/testf --remove chrvcf/removesamples.txt --missing --me 0.05 0.02 --chr 1-22 --keep-allele-order --allow-no-sex --out chrvcf/disorder_filtered --snps-only just-acgt',plinkdir))
lmiss <- read.table('chrvcf/disorder_filtered.lmiss',h=T,stringsAsFactors=F)
imiss <- read.table('chrvcf/disorder_filtered.imiss',h=T,stringsAsFactors=F)

# Plot the Genotype Rate
pdf("../Plots/IndGenotypeRate.pdf", onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
print(ggplot(imiss,aes(x = F_MISS)) + xlab("IndividualMissingness") + geom_histogram(colour="black",bins=60))
dev.off()
pdf("../Plots/LocusGenotypeRate.pdf", onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
print(ggplot(lmiss,aes(x = F_MISS)) + xlab("SNPMissingness") + geom_histogram(colour="black",bins=100))
dev.off()
pdf("../Plots/LocusGenotypeRate_Zoomed.pdf", onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
print(ggplot(lmiss,aes(x = F_MISS)) + xlab("SNPMissingness") + geom_histogram(colour="black",bins=200) + coord_cartesian(xlim = c(0,0.05)))
dev.off()



# Annotate the SNPs using ANNOVAR per Chromsome and Extract only Non-synonymous, putative splice sites, Stopgain, Stoploss

data <- read.table('chrvcf/disorder_filtered.fam',h=F,stringsAsFactors=F)
write.table(data$V2,'finalsamples.txt',quote=F,row.names=F,col.names=F)
for(i in 1:22){
  
  # Extract the SNPs 
	data <- read.table('chrvcf/disorder_filtered.bim',h=F,stringsAsFactors=F)
	write.table(data[which(data$V1 == i),c(1,4)],'finalsnps.txt',quote=F,sep="\t",row.names=F,col.names=F)
	
	system(sprintf('vcftools --gzvcf chrvcf/test%s.vcf.gz --positions finalsnps.txt --keep finalsamples.txt --remove-indels --recode --recode-INFO-all --stdout | gzip -c > chrvcf/testfiltered%s.vcf.gz',i,i))
	system(sprintf('%s/convert2annovar.pl -format vcf4old chrvcf/testfiltered%s.vcf.gz -withfreq -outfile chrvcf/test%s.avinput',annovar,i,i))
	system(sprintf('perl %s/annotate_variation.pl -out chrvcf/test%s.txt -build hg19 chrvcf/test%s.avinput %s/humandb/',annovar,i,i,annovar))
	edata <- read.table(sprintf('chrvcf/test%s.txt.exonic_variant_function',i),h=F,stringsAsFactors=F,sep='\t')
	edata <- edata[-grep("^syno|unknown",edata$V2,perl=T),]
	fsnps <- as.vector(paste(edata$V4,edata$V5,sep=":"))
	if(i == 1){
		write.table(fsnps,'nonsynonymous_splicesites.txt',quote=F,row.names=F,col.names=F)
	} else {
		write.table(fsnps,'nonsynonymous_splicesites.txt',quote=F,row.names=F,col.names=F,append=T)
	}
	rm(data,edata,fsnps)
  
  # Remove the intermediary files
	system(sprintf('rm chrvcf/testfiltered%s.vcf.gz chrvcf/test%s.txt.* chrvcf/test%s.avinput',i,i,i))
}

# Filter out for Nonsynonymous, Stopgain, Stoploss, Putative splice
system(sprintf('%s --bfile chrvcf/disorder_filtered --extract nonsynonymous_splicesites.txt --keep-allele-order --allow-no-sex --make-bed --chr 1-22 --out chrvcf/disorder_filtered_nonsynonymous',plinkdir))

# Phase the Variant data of Nonsynonymous, Stopgain, Stoploss, Putative splics
seedgen <- 12345
disordername <- 'DISODER_NAME'
threads <- 20 # Number of threads to be used

# Split the Nonsynonymous Variant data PLINK file into Chromosomes and Phase it using Shapeit software
for(i in 1:22){
  system(sprintf('%s --bfile chrvcf/disorder_filtered_nonsynonymous --chr %s --geno 0.1 --allow-no-sex --keep-allele-order --make-bed --out chrvcf/phase/disorder_chr%s',plinkdir,i,i))
  sink(sprintf('../Scripts/%s_chr_%s.sh',disordername,i))
  cat(sprintf('#!/bin/bash\n'))
  cat(sprintf('\n'))
  cat(sprintf('%s --seed %s -B chrvcf/phase/disorder_chr%s -M %s/CEU_chr%s.txt -O chrvcf/phase/disorder_chr%s_phased -T %s\n',shapeit,seedgen,i,ceudir,i,i,threads))
  sink()
  system(sprintf('qsub -l nodes=1:ppn=1,naccesspolicy=singleuser -l walltime=100:00:00 -N Phasing ../Scripts/%s_chr_%s.sh',disordername,i))
}





