#!/bin/bash
#PBS -q unlimited
#PBS -l nodes=1:ppn=5
#PBS -N GS_cv_age

if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
fi

###-------./run.sh /media/disk5/songhl/structLMM/mtg2/vgxe0.25/phe0.25.txt /media/disk5/songhl/structLMM/mtg2/real 5 5
##文件为全路径，表型文件五列：FID ID phe cov TBV； 芯片文件为质控后的bed格式
 
PHE=$1
SNP=$2
repeat=$3
fold=$4

#---random pheotype files for repeat
chr=(`seq 1 ${repeat}`)
for i in "${chr[@]}";do
	if [ ! -d repeat"$i" ];then
		mkdir repeat"$i"
		cd repeat"$i"
		python2 ../bin/data.py  ${PHE}  phe.txt  ${fold}
		cd ../
	else
		cd repeat"$i"
		python2 ../bin/data.py  ${PHE}  phe.txt  ${fold}
		cd ../
	fi
done

#--ped\map\G
mkdir G
cd G
../bin/plink --cow  --bfile ${SNP} --recode12 --out marker
cp ../bin/gmatrix_input.py  ../bin/gmatrix  ../bin/G-matrix.sh  ./
python2 ../bin/marker.py marker.ped marker.dat
python2 gmatrix_input.py 
./G-matrix.sh
cd ../

chr=(`seq 1 ${repeat}`)
for i in "${chr[@]}";do
	python2 ./bin/runcv.py repeat"$i" $1 $2 $4 qsub"$i".sh
	chmod 777 qsub"$i".sh
	#nohup ./qsub"$i".sh &
done
