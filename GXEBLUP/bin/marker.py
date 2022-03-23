import sys
ped=open(sys.argv[1])
q=open(sys.argv[2],'w')
for i in ped:
	f=i.split()
	q.write(f[1]+'\t'+'\t'.join(f[6:])+'\n')
q.close()
ped.close()