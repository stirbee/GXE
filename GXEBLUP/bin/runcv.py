import sys
dir=sys.argv[1]
phe=sys.argv[2]
snp=sys.argv[3]
fold=sys.argv[4]


q=open(sys.argv[5],'w')
title='''#!/bin/bash
#PBS -q unlimited
#PBS -l nodes=1:ppn=1

if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
fi
'''
q.write(title+'\n')
q.write('cd '+dir+'\n')
q.write('PHE='+phe+'\n')
q.write('SNP='+snp+'\n')
q.write('fold='+fold+'\n')

lmm='''chr=(`seq 1 ${fold}`)
for i in "${chr[@]}";do
	mkdir cv"$i"
	cd cv"$i"
	python2 ../../bin/LMMphe.py ${PHE} ../"$i".val ${SNP}
	cd ../
done

#run structLMM 
python ../bin/driverLMM.py genotype newphe.csv env.txt result.csv ${fold}
'''
q.write(lmm+'\n')


code='''chr=(`seq 1 ${fold}`)
for i in "${chr[@]}";do
	cd cv"$i"
	cp ../../bin/gmatrix_input.py  ../../bin/gmatrix ../../bin/selectsnp.py ../../bin/G-matrix.sh  ./
	cp ../../bin/ggblup.DIR ./
	cp ../../bin/r_dmuai ./
	cp ../../bin/gblup.DIR ./

	#gblup
	python2 ../../bin/pickdata.py ../"$i".val ${PHE} data.txt
	python2 ../../bin/picktbv.py ../"$i".val ${PHE} val_yc
	./r_dmuai gblup
	python2 ../../bin/gebv.py gblup.SOL val_yc yc_gebv
	python2 ../../bin/COR_REG_used2.py yc_gebv resultG.txt

	chr=(`seq 1 5`)
	for ii in "${chr[@]}";do
		python2 selectsnp.py  result.csv  1.0E-0"$ii" ../../G/marker.map ../../G/marker.ped

		./r_dmuai ggblup
		
		python2 ../../bin/ggebv.py ggblup.SOL val_yc yc_ggebv
		python2 ../../bin/COR_REG_used2.py yc_ggebv resultGXE.txt
	done
	cd ../
done
'''
q.write(code+'\n')
q.write('paste ./cv*/resultGXE.txt >gxeblup.txt'+'\n')
q.write('cat ./cv*/resultG.txt > gblup.txt'+'\n')
q.close()