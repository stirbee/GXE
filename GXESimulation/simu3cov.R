
################################################################################################
###note: check file.exist: plink1.9 mtg2 XXX.fam XXX.bim XXX.bed
###To run the simulation:
###simu_phen(ip[1],ip[2],ip[3],ip[4],ip[5],ip[6],ip[7],ip[8],ip[9],ip[10],ip[11],ip[12],ip[13])
###There are 11 parameters used in the simulation, which needs to specify in this order:
###V1 nm: the prefix of genotype data, which is in plink format
###V2 nr.qlt: number of QTL
###V3 alpha0: additive genetic variance of first trait
###V4 beta: additive genetic variance of second trait
###V5 beta2: additive genetic variance of third trait
###V6 beta3: additive genetic variance of third trait
###V7 alpha1: variance of interaction 
###V8 alpha2: variance of interaction
###V9 alpha3: variance of interaction
###V10 cova0beta: covariance between alpha0 and beta
###V11 cova0beta2: covariance between alpha0 and beta2
###V12 e0: residual variance of first trait
###V13 epsilon: residual variance of second trait
###V14 epsilon2: residual variance of third trait
###V15 epsilon3: residual variance of four trait
###V16 mod: h0 or h1 
###the simulation for h0 is as followed:
###h0: y=a0+e0, c=beta + epsilon
###the simulation for h1 is as followed:
###h1: y=a0+a1*c1+a2*c2+e0, c1=beta + epsilon, c2=beta2+epsilon2
################################################################################################

options(stringsAsFactors=FALSE)

ip<-commandArgs(trailingOnly=TRUE)
options(warn=1)

simu_phen<-function(nm,nr.qtl,alpha0,beta,beta2,beta3,alpha1,alpha2,alpha3,cova0beta,cova0beta2,
		e0,epsilon,epsilon2,epsilon3,mod){
		
	nr.qtl<-as.numeric(nr.qtl)
	alpha1<-as.numeric(alpha1)
	alpha0<-as.numeric(alpha0)
	alpha2<-as.numeric(alpha2)	
	alpha3<-as.numeric(alpha3)
	beta<-as.numeric(beta)
	beta2<-as.numeric(beta2)	
	beta3<-as.numeric(beta3)
	cova0beta<-as.numeric(cova0beta)
	cova0beta2<-as.numeric(cova0beta2)
	e0<-as.numeric(e0)
	epsilon<-as.numeric(epsilon)
	epsilon2<-as.numeric(epsilon2)
	epsilon3<-as.numeric(epsilon3)
		
	cova0a1<-0.05
	cova0a2<-0.05
	cova1a2<-0
	cova1beta<-0
	cova2beta<-0	
	cova1beta2<-0
	cova2beta2<-0
	covbetabeta2<-0.5

	
	cove0ep<-0
	cove0ep2<-0
	covepep2<-0

	nr.trait<-4
#------------------------
### genetic and residual variance-covariance structure 
	gmat<-matrix(c(alpha0,cova0beta,cova0beta2,0.5,cova0a1,cova0a2,0.05,
				cova0beta,beta,covbetabeta2,0.5,cova1beta,cova2beta,0,				
				cova0beta2,covbetabeta2,beta2,0.5,cova1beta2,cova2beta2,0,
				0.5,0.5,0.5,beta3,0,0,0,
				cova0a1,cova1beta,cova1beta2,0,alpha1,cova1a2,0,
				cova0a2,cova2beta,cova2beta2,0,cova1a2,alpha2,0,
				0.05,0,0,0,0,0,alpha3),ncol=7,byrow=TRUE)
	mu_g<-array(0,7)
	
	emat<-matrix(c(e0,cove0ep,cove0ep2,0,
				cove0ep,epsilon,covepep2,0,
				cove0ep2,covepep2,epsilon2,0,
				0,0,0,epsilon3),ncol=4,byrow=TRUE)
	mu_e<-array(0,4)
	
	if(mod=='h0'){
		gmat<-gmat[1:4,1:4]
		emat<-emat[1,1]
		
		mu_g<-mu_g[1:4]
		mu_e<-mu_e[1]
	}
	
#------------------------	
### causal SNPs	
	bim<-read.table(paste(nm,'.freq',sep=''),header=FALSE)
	fam<-read.table(paste(nm,'.fam',sep=''),header=FALSE)
	plinkbim<-read.table(paste(nm,'.bim',sep=''),header=FALSE)
	colnames(plinkbim)<-c('chr','SNP','pos','bp','A1','A2')
	
	v1<-seq(1:nrow(plinkbim))

	
	set.seed((nr.qtl+alpha0+beta+beta2+beta3+alpha1+alpha2+alpha3+cova0beta+cova0beta2+e0+epsilon+epsilon2+epsilon3)*100+ceiling(nrow(plinkbim)/100)+1)
	v2<-sort(sample(v1,nr.qtl,replace=FALSE))
#------------------------
### effects of causal SNPs
	if(!file.exists(paste(nm,'.freq',sep=''))){
		system(paste('./mtg2 -plink ',nm,' -frq 1 > mtgfrq.log',sep=''))
	}
		
	library(MASS)
	set.seed((nr.qtl+alpha0+beta+beta2+beta3+alpha1+alpha2+alpha3+cova0beta+cova0beta2+e0+epsilon+epsilon2+epsilon3)*100+ceiling(nrow(plinkbim)/100)+2)
	gv<-mvrnorm(nr.qtl,t(mu_g),gmat)
	v12<-gv*sqrt(1/(bim$V3[v2]*(1-bim$V3[v2]/2)))
	v3<-cbind(v2,as.character(bim$V2[v2]),v12[,1],bim$V3[v2])
	v4<-cbind(v2,as.character(bim$V2[v2]),v12[,2],bim$V3[v2])
	v5<-cbind(v2,as.character(bim$V2[v2]),v12[,3],bim$V3[v2])	
	v6<-cbind(v2,as.character(bim$V2[v2]),v12[,4],bim$V3[v2])
	
	sink('snp.lst1')
	write.table(v3,quote=F,col.name=F,row.name=F)
	sink()
	
	sink('snp.lst2')
	write.table(v4,quote=F,col.name=F,row.name=F)
	sink()
	
	sink('snp.lst3')
	write.table(v5,quote=F,col.name=F,row.name=F)
	sink()
	
	sink('snp.lst4')
	write.table(v6,quote=F,col.name=F,row.name=F)
	sink()
	

	if(mod!='h0'){
		v7<-cbind(v2,as.character(bim$V2[v2]),v12[,5],bim$V3[v2])
		sink('snp.lst5')
		write.table(v7,quote=F,col.name=F,row.name=F)
		sink()
		v8<-cbind(v2,as.character(bim$V2[v2]),v12[,6],bim$V3[v2])
		sink('snp.lst6')
		write.table(v8,quote=F,col.name=F,row.name=F)
		sink()	
		v9<-cbind(v2,as.character(bim$V2[v2]),v12[,7],bim$V3[v2])
		sink('snp.lst7')
		write.table(v9,quote=F,col.name=F,row.name=F)
		sink()
		
	}
#-----------------------------------		
### to get breeding value for each individual
	system(paste('./mtg2 -plink ',nm, ' -simreal snp.lst1 > tbv1.log',sep=''))
	system(paste('mv ',nm,'.bv tbv1.bv',sep=''))
	
	system(paste('./mtg2 -plink ',nm, ' -simreal snp.lst2 > tbv2.log',sep=''))
	system(paste('mv ',nm,'.bv tbv2.bv',sep=''))
	
	system(paste('./mtg2 -plink ',nm, ' -simreal snp.lst3 > tbv3.log',sep=''))
	system(paste('mv ',nm,'.bv tbv3.bv',sep=''))	

	system(paste('./mtg2 -plink ',nm, ' -simreal snp.lst4 > tbv4.log',sep=''))
	system(paste('mv ',nm,'.bv tbv4.bv',sep=''))	
	
	if(mod!='h0'){
	system(paste('./mtg2 -plink ',nm, ' -simreal snp.lst5 > tbv5.log',sep=''))
	system(paste('mv ',nm,'.bv tbv5.bv',sep=''))
	system(paste('./mtg2 -plink ',nm, ' -simreal snp.lst6 > tbv6.log',sep=''))
	system(paste('mv ',nm,'.bv tbv6.bv',sep=''))
	system(paste('./mtg2 -plink ',nm, ' -simreal snp.lst7 > tbv7.log',sep=''))
	system(paste('mv ',nm,'.bv tbv7.bv',sep=''))
	}
#------------------------------------------------------------	
### scale variance of breeding values according to gmat
	yv1<-read.table('tbv1.bv')
	yv1<-scale(yv1[,1])
	yv1=yv1*gmat[1,1]^.5
		
	yv2<-read.table('tbv2.bv')
	yv2<-scale(yv2[,1])
	yv2=yv2*gmat[2,2]^.5
	
	yv3<-read.table('tbv3.bv')
	yv3<-scale(yv3[,1])
	yv3=yv3*gmat[3,3]^.5
	
	yv4<-read.table('tbv4.bv')
	yv4<-scale(yv4[,1])
	yv4=yv4*gmat[4,4]^.5	
	
	if(mod!='h0'){
	yv5<-read.table('tbv5.bv')
	yv5<-scale(yv5[,1])
	yv5=yv5*gmat[5,5]^.5
	
	yv6<-read.table('tbv6.bv')
	yv6<-scale(yv6[,1])
	yv6=yv6*gmat[6,6]^.5

	yv7<-read.table('tbv7.bv')
	yv7<-scale(yv7[,1])
	yv7=yv7*gmat[7,7]^.5	
	
	}else{
	yv5<-0
	yv6<-0
	yv7<-0
	}
		
	yv<-as.data.frame(cbind(yv1,yv2,yv3,yv4,yv5,yv6,yv7))
#------------------------------------------------------------	
### simu phenotype
	library(MASS)
	set.seed((nr.qtl+alpha0+beta+beta2+beta3+alpha1+alpha2+alpha3+cova0beta+cova0beta2+e0+epsilon+epsilon2+epsilon3)*100+ceiling(nrow(plinkbim)/100)+3)
	v12=mvrnorm(length(yv$V1),t(mu_e),emat)
### phenotype of covariate (eg smk) c=beta +epsilon
	yv$V2=yv$V2+v12[,2]
	yv$V3=yv$V3+v12[,3]	
	yv$V4=yv$V4+v12[,4]
### phenotype of main response (eg bmi)
	if(mod!='h0'){
	yv$V1=yv$V1+yv$V5*yv$V2+yv$V6*yv$V3+yv$V7*yv$V4+v12[,1]#   +v12[,3]*yv$V2 ### h1 y=alpha0 + alpha1 * c1 + alpha2 * c2 +alpha3 * c3+e0 
	}else{
	yv$V1=yv$V1+v12[,1] ### h0, null model y=alpha0 + e 
	}

	out=yv[,1:4]
	v1=cbind(as.character(fam$V1),as.character(fam$V2),out)
	sink(paste('sample_',mod,'.dat',sep=''))
	write.table (v1,row.names=F,col.names=F,quote=F)
	sink()
		
}

simu_phen(ip[1],ip[2],ip[3],ip[4],ip[5],ip[6],ip[7],ip[8],ip[9],ip[10],ip[11],ip[12],ip[13],ip[14],ip[15],ip[16])

# Rscript  simu2cov.R 'real' 500 1 1 1  0.2 0.2 0.5 0.5 1 1 1 'h1'







