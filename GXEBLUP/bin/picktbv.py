#!/usr/bin/python
#-*- coding:utf-8 -*-

import sys
id=open(sys.argv[1],'r')
y=open(sys.argv[2],'r')
q=open(sys.argv[3],'w')

dick={}
for i in y:
	f=i.split()
	dick[f[1]] =f[-1]

for e in id:
	f2=e.split()
	if f2[0]  in dick:
		q.write(f2[0]+'\t'+dick[f2[0]]+'\n')
id.close()
y.close()
q.close()

