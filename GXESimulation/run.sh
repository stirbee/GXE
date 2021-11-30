./plink --cow --file cow --maf 0.01 --make-bed --out realcow
./mtg2 -plink realcow -frq 1

Rscript simu1cov.R 'realcow' 306  1 1 0.25 0.5 1 1 0.25 0.3 'h1'
mkdir vgxe0.25
mv sample_h1.dat snp.* tbv* ./vgxe0.25


Rscript simu1cov.R 'realcow' 306  1 1 1 0.5 1 1 1 0.3 'h1'
mkdir vgxe1
mv sample_h1.dat snp.* tbv* ./vgxe1

Rscript simu1cov.R 'realcow' 306  1 1 2 0.5 1 1 2 0.3 'h1'
mkdir vgxe2
mv sample_h1.dat snp.* tbv* ./vgxe2

Rscript simu2cov.R 'realcow' 306 1 1 1  0.1 0.15 0.5 0.5 1 1 1 'h1'
mkdir 2cor
mv sample_h1.dat  snp.*  tbv* ./2cor


Rscript simu3cov.R 'realcow' 306 1 1 1 1 0.1 0.1 0.05 0.5 0.5 1 1 1 1 'h1'
mkdir 3cor
mv sample_h1.dat  snp.*  tbv* ./3cor






