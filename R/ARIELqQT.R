
ariel.data.fn <- function(pheno,allgenodata,chr,start,end,nc=1) {
# outputs the data required for runnning amelia function for a given region (start,end)
# chr is the chromosome number
# nc=TRUE for non-consensus qs analysis, othrwise consensus qs analysis

# input: 
# pheno: N-vector of QTs 
# genodata: columns are chr, position,reference allele, alternative allele, SNPquality,
#          set of 3 columns for indiv i (allele 1, allele 2, Gquality)

# G is K by N matrix of genotype scores
# genotype scores: 0,1,2 for number of minor alleles
# Q[k] SNP quality score at SNP k (based on reads from all indivs)


#pheno <- unlist(pheno)
#Nold <- length(pheno)

#indrm <- which(pheno<0)
#if(length(indrm)>0)     {
# pheno <- pheno[-indrm]
#}

N <- length(pheno)
#print(indrm)
#n1 <- sum(pheno==1)
#n2 <- N-n1

print(c("N",N))

names(allgenodata)[1:2] <- c("chr","pos")
chrgeno <- allgenodata[allgenodata$chr==chr,]

genodata <- chrgeno[chrgeno$pos >= start & chrgeno$pos <= end,]

print(dim(genodata))

if (dim(genodata)[1]==0){
#	return(list(pheno=pheno,geno=genodata,qSNP=0,qgeno=0))
out <- list(pheno=pheno,geno=genodata,qgeno=0)
}

if (dim(genodata)[1]>0) {
names(genodata)[1:4] <- c("chr","pos","ref","alt")

if(nc==0) {
# SNP quality scores
Q <- genodata[,5] # SNP scores: 1 for each variant
Q <- 1-10^(-Q/10) # transform phred score to prob of correct

Qrm <- which(Q==0) # remove data with SNP quality = 0
if(length(Qrm)>0) {
Q <- Q[-Qrm]
genodata <- genodata[-Qrm,]
}
q <- Q
}
#####

K <- dim(genodata)[1]
ref <- genodata$ref
print(c("K",K))

if(nc==1) {
# genotype quality scores 
q <- genodata[,(8+3*(0:(N-1)))] # K by N matrix
q <- matrix(as.numeric(unlist(q)),nrow=K,ncol=N)
q <- 1-10^(-q/10) # transform phred scores to probs
}
#q <- q[,-indrm]

# haplotypes
g1 <- genodata[,(6+3*(0:(N-1)))] # K by N,  1=alt allele, 0=ref allele
g2 <- genodata[,(7+3*(0:(N-1)))]
#print(dim(g1))
#g1 <- g1[,-indrm]
#g2 <- g2[,-indrm]



MAF.fn <- function(g1,g2,j)	{
# g1j is vector of 1{1st allele = alt allele} for variant j
 g1j <- g1[j,]
 g2j <- g2[j,]
 mj <- sum(g1j!="NA",na.rm=TRUE) # no. of typed indivs at variant j
 alt.count <- sum(g1j,na.rm=TRUE) + sum(g2j,na.rm=TRUE)
 ref.count <- mj*2 - alt.count
 if(ref.count <= alt.count) {
  maf <- ref.count/(2*mj) # each indiv has 2 alleles at a variant
  a <- 1 # indicator that minor allele = ref allele
			   }					
if(ref.count > alt.count) {
  maf <- alt.count/(2*mj)
  a <- 0 # indicator that minor allele = ref allele
			   }
return(c(maf,a))					
				}

j.mat <- matrix(1:K,nrow=K)

MAF.out <-apply(j.mat,1,MAF.fn,g1=g1,g2=g2)

maf <- MAF.out[1,]
minor.ref <- MAF.out[2,]
# where minor.ref=1 have 1{minor al = ref al} so change g1,g2 coding
indm <- which(minor.ref==1) 
# where minor.ref=0 have minor allele = ref allele so no change to g1,g2 coding

print(summary(maf))

if(length(indm)>0){

g1new <- g1
g1new[indm,] <- 1-g1[indm,]
g1 <- g1new; rm(g1new)  # 1=minor allele

g2new <- g2
g2new[indm,] <- 1-g2[indm,]
g2 <- g2new; rm(g2new) # 1=minor allele
}

ind0 <- which(maf==0) # remove non-polymorphic SNPs
if(length(ind0)>0) 	{
 maf <- maf[-ind0]
 g1 <- g1[-ind0,]
 g2 <- g2[-ind0,]
 if(nc==1) q <- q[-ind0,]
 if(nc==0) {q <- q[-ind0]; q<- matrix(rep(q,N),nrow=length(q),ncol=N) }
					}
G <- g1+g2


out <- list(pheno=pheno,geno=G,qgeno=q)
}


#return(list(pheno=pheno,geno=G,qSNP=Q))
return(out)
}



ARIELq.fn <- function(y,g,qgeno,mafrv) {
# n indivs, m variants
# y= case(1)-control(0) status
# assume g (genotype scores) has m rows and n columns
# assume high quality score is high quality
# mafrv is cut-off for variants included in collapsing method


m <- dim(g)[1]
n <- dim(g)[2]

MAF.fn <- function(g,j)	{
# g is m by n matrix of genotypes=minor allele counts {0,1,2}
 gj <- g[j,]
 nj <- sum(gj!="NA",na.rm=TRUE) # no. of typed indivs at variant j
 maf <- sum(gj,na.rm=TRUE)/(2*nj)
return(maf)					
			}


j.mat <-matrix(1:m,nrow=m)
MAF <-apply(j.mat,1,MAF.fn,g=g) 

ind.rvs <- which(MAF < mafrv)

if(length(ind.rvs)<2)	{
out <- rep("NA",8)
						}

if(length(ind.rvs)>1)	{
 
grv <- g[ind.rvs,] # genotypes of rare variants 
qrv <- qgeno[ind.rvs,]

mrv <- dim(grv)[1]
print(c("mrv=",mrv))

# original collapsing method
r.fn <- function(g) return(sum(g>0,na.rm=TRUE)) # for indiv i, count of number of variants with at least 1 minor allele
len.fn <- function(g) return(sum(g!="NA",na.rm=TRUE)) # number of variants successfully typed

w.fn <- function(i,q,g)	{
 qi <- q[,i]
 gi <- g[,i]
 mr <- length(gi)
# w <- mean(qi)*( sum(gi==0)==mr | sum(gi>0)==mr ) +  ( sum(gi==0)<mr & sum(gi>0)<mr )*(mean(qi[gi==0]) + sum(qi[gi>0]))/(1+sum(gi>0))
if(sum(gi==0)==mr | sum(gi>0)==mr) w <- mean(qi)
if( sum(gi==0)<mr & sum(gi>0)<mr ) w <- (mean(qi[gi==0]) + sum(qi[gi>0]))/(1+sum(gi>0))

return(w)
				}


i.mat <- matrix(1:N,nrow=N)
w <- apply(i.mat,1,w.fn,q=qrv,g=grv) # weight for indiv i

r <- apply(grv,2,r.fn) # n vector

ni <- apply(grv,2,len.fn) 
prop <- r/ni

##
oreg <- glm(y~prop)
opval <- summary(oreg)$coefficients[2,4]
obeta <- summary(oreg)$coefficients[2,1]
se <- summary(oreg)$coefficients[2,2]
oOR <- obeta
oORadj <- obeta/mrv
Loadj <- obeta/mrv-1.96*se/mrv
Uoadj <- obeta/mrv+1.96*se/mrv


wreg <- glm(y~prop,weights=w)
wpval <- summary(wreg)$coefficients[2,4]
wbeta <- summary(wreg)$coefficients[2,1]
se <- summary(wreg)$coefficients[2,2]
wOR <- wbeta
wORadj <- wbeta/mrv
Lwadj <- wbeta/mrv-1.96*se/mrv	
Uwadj <- wbeta/mrv+1.96*se/mrv

#


###########################

							


# output = number of rare snps,original and weighted pvalues, ORs, lower 95% confidence limits, upper 95% confidence limits

#out <- c(o1pval,o2pval,w1pval,w2pval,o1OR,o2OR,w1OR,w2OR,Lo1,Lo2,Lw1,Lw2,Uo1,Uo2,Uw1,Uw2)
out <- c(mrv,opval,wpval,oORadj,wORadj,Loadj,Lwadj,Uoadj,Uwadj,oOR,w2OR)

#write.table(matrix(c(o1pval,o2pval,w1pval,w2pval,o1OR,o2OR,w1OR,w2OR),nrow=1),sep="\t",row.names=FALSE,append=FALSE,file=pfile,col.names=FALSE)

	}
return(out)
}

######


