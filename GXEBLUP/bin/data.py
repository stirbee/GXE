import sys,os,random

y=open(sys.argv[1])
q=open(sys.argv[2],'w')
n=sys.argv[3]

ID=[]

for i in y:
	f=i.split()
	ID.append(f[1])

val=random.sample(ID,len(ID))	
	
for i in range(0,len(val)):
		q.write(val[i]+'\n')	
		
y.close()
q.close()	



k = n    #sys.argv[2]   ###split k parts
y1=open(sys.argv[2],'r')

M=len(y1.readlines())/int(k)
y1.seek(0)


for i in range(1,int(k)+1):
    pfile = open(str(i)+'.val','w')
    for j in range(M):
        pfile.write(y1.readline())        
    pfile.close()
    
    
pfile = open(str(i)+'.val','a')  
for j in range(M):
    pfile.write(y1.readline())
pfile.close()

