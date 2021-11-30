import sys,os,cStringIO
phe=open(sys.argv[1])
val=open(sys.argv[2])
snp=sys.argv[3]

VAL={}
for i in val:
	f=i.split()
	VAL[f[0]]=f[0]

PHE={}
for i in phe:
	f=i.split()
	PHE[f[1]]='\t'.join(f[1:])


rel=open('rel_id','w')
for i in open(snp+'.fam'):
	f=i.split()
	if f[1] not in VAL:
		rel.write(f[0]+'\t'+f[1]+'\n')
rel.close()

fin=open('rel_id')
p=open('newphe.csv','w')
env=open('env.txt','w')
line1=[]
line2=[]
for i in fin:
	f=i.split()
	if f[1] in PHE:
		line1.append(PHE[f[1]].split()[0])
		line2.append(PHE[f[1]].split()[1])
		env.write('\t'.join(PHE[f[1]].split()[2:3])+'\n')
p.writelines(','+','.join(line1)+'\n')
p.writelines('phe'+','+','.join(line2)+'\n')
p.close()
fin.close()
env.close()
os.system("../../bin/plink --cow --bfile %s  --keep %s --make-bed --out genotype >& /dev/null"%( snp, 'rel_id'))


