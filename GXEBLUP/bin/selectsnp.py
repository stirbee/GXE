import sys,os
file=open(sys.argv[1])
n=float(sys.argv[2])
bigmap=open(sys.argv[3])
bigpen=sys.argv[4]

q=open('gxesnp.txt','w')
snp={}
for i in file.readlines()[1:]:
	f=i.split(',')
	if float(f[7])<n:	
		q.write(f[1]+'\n')
		snp[f[1]]=f[1]
print len(snp)
q.close()
q2=open('remainsnp.txt','w')
for i in bigmap:
	f=i.split()
	if f[1] not in snp:
		q2.write(f[1]+'\n')
q2.close()
os.system("../../bin/plink --cow --ped  %s --map %s --extract  gxesnp.txt --recode  --out gxe >& /dev/null"%(sys.argv[4],sys.argv[3]))
os.system("../../bin/plink --cow  --ped  %s --map %s  --extract  remainsnp.txt --recode  --out remain >& /dev/null"%(sys.argv[4],sys.argv[3]))

def marker(file):
	fi=open(file+'.ped')
	q=open('marker.dat','w')
	for i in fi:
		f=i.split()
		q.write(f[1]+'\t'+'\t'.join(f[6:])+'\n')
	q.close()
	fi.close()
	os.system("python2 gmatrix_input.py >& /dev/null")
	os.system("./G-matrix.sh >& /dev/null")
	os.system("mv igmatrix.dat %s >& /dev/null"%(file+'.iga'))

marker('gxe')
marker('remain')


