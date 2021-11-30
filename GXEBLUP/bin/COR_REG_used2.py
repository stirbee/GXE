#-*- coding:utf-8 -*-
###########################The file with two line,firt line is Y£¬the second line is X.
from math import sqrt

def multiply(a,b):
    sum_ab=0.0
    for i in range(len(a)):
        temp=a[i]*b[i]
        sum_ab+=temp
    return sum_ab

def cal_pearson(x,y):
    n=len(x)
    sum_x=sum(x)
    sum_y=sum(y)
    sum_xy=multiply(x,y)
    sum_x2 = sum([pow(i,2) for i in x])
    sum_y2 = sum([pow(j,2) for j in y])
    molecular=sum_xy-(float(sum_x)*float(sum_y)/n)
    denominator=sqrt((sum_x2-float(sum_x**2)/n)*(sum_y2-float(sum_y**2)/n))
    return molecular/denominator

import sys
file=sys.argv[1]
result=sys.argv[2]
g=open(result,'a')
f=open(file,'r')
data={}
lines=f.readlines()
for line in lines:
    cols=line.strip('\n').split()
    for i in range(len(cols)):
        data.setdefault(i,[]).append(float(cols[i]))
	
x=data[1]
y=data[2]

if __name__=='__main__':
	print ('corr'+'\t'+str(round(cal_pearson(x,y),3)))
	g.write('corr'+'\t'+str(round(cal_pearson(x,y),3))+'\n')
########################################################  REG(not use)
#import numpy as np
#import pandas as pd
##x = np.array([2,3,4,6])
#xx = pd.DataFrame({"k": y})
#yy = pd.Series(x)   
#res = pd.ols(y=yy, x=xx)  
#print ('regression='+str(res.beta[0]))

######################################################Corr(not used)
#import scipy.stats as stats 
#r, p=stats.pearsonr(x,y)
#print r , p  ########cor   p_value
########################################################REG(used)
from scipy import stats
import numpy as np
slope, intercept, r_value, p_value, std_err = stats.linregress(y,x)  
slope = round(slope,3)  
intercept = round(intercept,3)  
print 'reg'+'\t'+str(slope) #, intercept 
g.write('reg'+'\t'+str(slope)+'\n')      