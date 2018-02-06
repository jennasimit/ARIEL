args           = (commandArgs(TRUE));
chr            = as.numeric(args[1]);
mafrv          = as.numeric(args[2]);
genotype.file  = args[3];# converted vcf file; output from vcf_format_change_v2.pl
phenotype.file = args[4];
gene.file      = args[5];
output.file    = args[6];
nc				= args[7];

cat("Input parameters:\n","chr              = ",chr,"\n","mafrv            = ",mafrv,"\n","genotype file    = ",genotype.file,"\n","phenotype file   = ",phenotype.file,"\n","gene file        = ",gene.file,"\n","output file name = ",output.file,"\n",sep="");

source("ARIELqQT.R") # ARIEL for SNP quality incorporation, QT data

#fregions <- "gene_info" # 4 columns: gene.name, chr, start, end
regions     <- read.table(gene.file);#,quote="",header=FALSE,sep="\t");
chr.regions <- regions[regions[,2]==chr,];
starts      <- chr.regions[,3];
ends        <- chr.regions[,4];
nregions    <- length(ends);
rm(regions); gc();

#fpheno <- paste("pheno.txt",sep="");

pheno <- read.table(phenotype.file);#,quote="",header=FALSE,sep="\t");
pheno <- unlist(pheno[,6])  # column 6 of PED file is phenotype
}

idlength <- length(pheno)*3;

#genofile <- paste("vcf_converted.txt",sep="");
genodata <- read.table(genotype.file,na.strings=".",colClasses=c('character','numeric',rep('character',3),rep('numeric',idlength)));#,header=FALSE,quote="",sep="\t"); 
genodata <- genodata[,-3]; # rm rsIDs column b/c not needed for analysis
chrgenodata <- genodata[genodata[,1]==chr,];
rm(genodata);
gc();


output <-file(output.file,"w"); 

for (i in 1:nregions)
  {
    data <- ariel.data.fn(pheno,chrgenodata,chr=chr,start=starts[i],end=ends[i],nc=nc);
    nsnps <- length(data$geno[,1]);
    gene <- c(as.character(chr.regions[i,1]),nsnps,chr,starts[i],ends[i]);
    print(c(i,gene,dim(data$geno)[1]));
    if (dim(data$geno)[1]<=1)
      {
        print("Geno is null!");
      }
    else
      {
        probs <- ARIELq.fn(y=data$pheno,g=data$geno,qgeno=data$qgeno,mafrv=mafrv);
        print(probs);
        cat(gene,probs,sep="\t",file=output);
        cat("\n",file=output);
      }
  }
close(output);				

q(save="no");

