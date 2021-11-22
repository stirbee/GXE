import sys
sol=open(sys.argv[1])
val=open(sys.argv[2])
q=open(sys.argv[3],'w')

g1={}
for i in sol:
	f=i.split()
	if f[0]=='3':
		g1[f[4]]=f[7]
for i in val:
	f=i.split()
	if f[0] in g1:
		gebv=float(g1[f[0]])
		q.write(f[0]+'\t'+f[1]+'\t'+str(gebv)+'\n')
sol.close()
val.close()
q.close()