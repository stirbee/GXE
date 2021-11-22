#!/usr/bin/env python

input_snp=open('marker.dat')
output_map=open('map.dat','w')
output_ped=open('pedigree.dat','w')
output_par=open('par.dat','w')
partxt="""$MINMAF 
0.01
$FREQMETHOD 
1
$SCALEMETHOD 
1
$DIAG_ADD 
0.01
$G_ADD  
0
$DIAG_ONE 
1
$PROP_A_to_G 
0
$CAL_DET 
1
$OUT_GMATRIX 
1
$OUT_IGMATRIX 
1
"""
n=0
count = 0
for aline in input_snp :
    linelst=aline.split()
    if count == 0:
        n = len(linelst)
    output_ped.write(linelst[0].rjust(4) + '  0     0 \n')
    count += 1
print 'The num of ind:', count
n = (n - 1) / 2
print 'The num of snp:', n
for i in range(n) :
    output_map.write(str(i+1).rjust(5) +'  ' + str(i+1).rjust(5) + '\n')
output_par.write(partxt) 
input_snp.close()
output_map.close()
output_ped.close()
output_par.close()
