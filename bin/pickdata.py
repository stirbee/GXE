#!/usr/bin/python
#-*- coding:utf-8 -*-

import sys
id=open(sys.argv[1],'r')
y=open(sys.argv[2],'r')
q=open(sys.argv[3],'w')

dick={}
for i in id:
	f=i.split()
	dick[f[0]] =f[0]

for e in y:
	f2=e.split()
	if f2[1] not in dick:
		q.write(f2[1]+'\t'+'1'+'\t'+f2[2]+'\n')
id.close()
y.close()
q.close()

