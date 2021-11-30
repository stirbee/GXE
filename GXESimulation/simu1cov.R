
################################################################################################
###note: check file.exist: plink1.9 mtg2 XXX.fam XXX.bim XXX.bed
###To run the simulation:
###simu_phen(ip[1],ip[2],ip[3],ip[4],ip[5],ip[6],ip[7],ip[8],ip[9],ip[10],ip[11])
###There are 11 parameters used in the simulation, which needs to specify in this order:
###V1 nm: the prefix of genotype data, which is in plink format
###V2 nr.qlt: number of QTL
###V3 alpha0: additive genetic variance of first trait
###V4 beta: additive genetic variance of second trait
###V5 alpha1: variance of interaction 
###V6 cova0beta: covariance between alpha0 and beta
###V7 e0: residual variance of first trait
###V8 epsilon: residual variance of second trait
###V9 e1: variance of interaction
###V10 cove0ep: variance between e0 and epsilon
###V11 mod: h0 or h1 
###the simulation for h0 is as followed:
###h0: y=a0+e0, c=beta + epsilon
###the simulation for h1 is as followed:
###h1: y=a0+a1*c+e0+e1*c, c=beta + epsilon
################################################################################################

options(stringsAsFactors=FALSE)
#nm<-'after_qc'
#alpha0<-1
#cova0beta<-0.5
#beta<-1
#alpha1<-0.25
#
#e0<-1
#cove0ep<-0.3
#epsilon<-1
#e1<-1
#nr.qtl<-100
#mod<-'h0'

ip<-commandArgs(trailingOnly=TRUE)
options(warn=1)

simu_phen<-function(nm,nr.qtl,alpha0,beta,alpha1,cova0beta,
		e0,epsilon,e1,cove0ep,mod){
		
	nr.qtl<-as.numeric(nr.qtl)
	alpha1<-as.numeric(alpha1)
	alpha0<-as.numeric(alpha0)
	beta<-as.numeric(beta)
	cova0beta<-as.numeric(cova0beta)
	e0<-as.numeric(e0)
	epsilon<-as.numeric(epsilon)
	e1<-as.numeric(e1)
	cove0ep<-as.numeric(cove0ep)
		
	cova0a1<-0.05
	cova1beta<-0
	
	cove0e1<-0.05
	cove1ep<-0
	nr.trait<-2
#------------------------
### genetic and residual variance-covariance structure 
	gmat<-matrix(c(alpha0,cova0beta,cova0a1,
				cova0beta,beta,cova1beta,
				cova0a1,cova1beta,alpha1),ncol=3,byrow=TRUE)
	mu_g<-array(0,3)
	
	emat<-matrix(c(e0,cove0ep,cove0e1,
				cove0ep,epsilon,cove1ep,
				cove0e1,cove1ep,e1),ncol=3,byrow=TRUE)
	mu_e<-array(0,3)
	
	if(mod=='h0'){
		gmat<-gmat[1:2,1:2]
		emat<-emat[1:2,1:2]
		
		mu_g<-mu_g[1:2]
		mu_e<-mu_e[1:2]
	}
	
#------------------------	
### causal SNPs	
	bim<-read.table(paste(nm,'.freq',sep=''),header=FALSE)
	fam<-read.table(paste(nm,'.fam',sep=''),header=FALSE)
	plinkbim<-read.table(paste(nm,'.bim',sep=''),header=FALSE)
	colnames(plinkbim)<-c('chr','SNP','pos','bp','A1','A2')
	
	v1<-seq(1:nrow(plinkbim))
	set.seed((nr.qtl+alpha0+beta+alpha1+cova0beta+e0+epsilon+e1+cove0ep)*100+ceiling(nrow(plinkbim)/100)+1)
	v2<-sort(sample(v1,nr.qtl,replace=FALSE))
#------------------------
### effects of causal SNPs
	if(!file.exists(paste(nm,'.freq',sep=''))){
		system(paste('./mtg2 -plink ',nm,' -frq 1 > mtgfrq.log',sep=''))
	}
		
	library(MASS)
	set.seed((nr.qtl+alpha0+beta+alpha1+cova0beta+e0+epsilon+e1+cove0ep)*100+ceiling(nrow(plinkbim)/100)+2)
	gv<-mvrnorm(nr.qtl,t(mu_g),gmat)
	v12<-gv*sqrt(1/(bim$V3[v2]*(1-bim$V3[v2]/2)))
	v3<-cbind(v2,as.character(bim$V2[v2]),v12[,1],bim$V3[v2])
	v4<-cbind(v2,as.character(bim$V2[v2]),v12[,2],bim$V3[v2])
	
	
	sink('snp.lst1')
	write.table(v3,quote=F,col.name=F,row.name=F)
	sink()
	
	sink('snp.lst2')
	write.table(v4,quote=F,col.name=F,row.name=F)
	sink()
	
	if(mod!='h0'){
		v5<-cbind(v2,as.character(bim$V2[v2]),v12[,3],bim$V3[v2])
		sink('snp.lst3')
		write.table(v5,quote=F,col.name=F,row.name=F)
		sink()	
	}
#-----------------------------------		
### to get breeding value for each individual
	system(paste('./mtg2 -plink ',nm, ' -simreal snp.lst1 > tbv1.log',sep=''))
	system(paste('mv ',nm,'.bv tbv1.bv',sep=''))
	
	system(paste('./mtg2 -plink ',nm, ' -simreal snp.lst2 > tbv2.log',sep=''))
	system(paste('mv ',nm,'.bv tbv2.bv',sep=''))
	
	if(mod!='h0'){
	system(paste('./mtg2 -plink ',nm, ' -simreal snp.lst3 > tbv3.log',sep=''))
	system(paste('mv ',nm,'.bv tbv3.bv',sep=''))
	}
#------------------------------------------------------------	
### scale variance of breeding values according to gmat
	yv1<-read.table('tbv1.bv')
	yv1<-scale(yv1[,1])
	yv1=yv1*gmat[1,1]^.5
		
	yv2<-read.table('tbv2.bv')
	yv2<-scale(yv2[,1])
	yv2=yv2*gmat[2,2]^.5
	
	if(mod!='h0'){
	yv3<-read.table('tbv3.bv')
	yv3<-scale(yv3[,1])
	yv3=yv3*gmat[3,3]^.5
	}else{yv3<-0}
		
	yv<-as.data.frame(cbind(yv1,yv2,yv3))
#------------------------------------------------------------	
### simu phenotype
	library(MASS)
	set.seed((nr.qtl+alpha0+beta+alpha1+cova0beta+e0+epsilon+e1+cove0ep)*100+ceiling(nrow(plinkbim)/100)+3)
	v12=mvrnorm(length(yv$V1),t(mu_e),emat)
### phenotype of covariate (eg smk) c=beta +epsilon
	yv$V2=yv$V2+v12[,2]
### phenotype of main response (eg bmi)
	if(mod!='h0'){
	yv$V1=yv$V1+yv$V3*yv$V2+v12[,1]#+v12[,3]*yv$V2 ### h1 y=alpha0 + alpha1 * c + e0 + e1 * c
	}else{
	yv$V1=yv$V1+v12[,1] ### h0, null model y=alpha0 + e 
	}

	out=yv[,1:2]
	v1=cbind(as.character(fam$V1),as.character(fam$V2),out)
	sink(paste('sample_',mod,'.dat',sep=''))
	write.table (v1,row.names=F,col.names=F,quote=F)
	sink()
		
}

simu_phen(ip[1],ip[2],ip[3],ip[4],ip[5],ip[6],ip[7],ip[8],ip[9],ip[10],ip[11])

# Rscript --vanilla function_simu_phen_gcec.r 'after_qc' 1000 1 1 0.25 0.5 1 1 0.25 0.3 'h0'
# Rscript --vanilla function_simu_phen_gcec.r 'after_qc' 1000 1 1 0.25 0.5 1 1 0.25 0.3 'h1'






