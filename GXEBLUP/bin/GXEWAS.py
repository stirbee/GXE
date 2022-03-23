import os,sys
import pandas as pd
import scipy as sp
from limix_core.util.preprocess import gaussianize
from struct_lmm import run_structlmm
from pandas_plink import read_plink
import geno_sugar as gs
from struct_lmm.utils.sugar_utils import norm_env_matrix
from limix.qtl import st_sscan
import numpy as np  
from pysnptools.snpreader import Bed


random = sp.random.RandomState(1)

# import genotype file
bedfile = sys.argv[1]
(bim, fam, G) = read_plink(bedfile, verbose=False)

#load snps
#Isnp = gs.is_in(bim, ("10", 25946396, 26022685))
#G, bim = gs.snp_query(G, bim, Isnp)

snps = G.compute().T
print(snps.shape)

# load phenotype file
phenofile = sys.argv[2]
dfp = pd.read_csv(phenofile, index_col=0)
pheno = gaussianize(dfp.loc["phe"].values[:, None])
# load environment
envfile = sys.argv[3]
E = sp.loadtxt(envfile)
E = norm_env_matrix(E)
# define fixed effect covs
covs = sp.ones((E.shape[0], 1))

M = Bed(bedfile, count_A1 = False).read()
freq = np.sum(M.val,  axis = 0) / (2*M.iid_count)
scale = np.sum(2*freq*(1-freq))
Z = M.val - 2*freq
kinship = np.dot(Z, Z.T)/scale
d = np.diag(kinship) + 0.001
np.fill_diagonal(kinship, d)
K = norm_env_matrix(kinship[0])

print(pheno.shape)
print(E.shape)
print(covs.shape)

# run struct lmm (both interaction and association tests)
#r = st_sscan(snps, pheno, E, M=covs, tests=["inter", "assoc"], verbose=False)
#print(r)

#Genome-wide analysis with Struct-LMM
from sklearn.impute import SimpleImputer
import geno_sugar.preprocess as prep
imputer = SimpleImputer(missing_values=sp.nan, strategy="mean")
preprocess = prep.compose(
    [
        prep.filter_by_missing(max_miss=0.10),
        prep.impute(imputer),
        prep.filter_by_maf(min_maf=0.01),
        prep.standardize(),
    ]
)

bims =[]
res = []
queue = gs.GenoQueue(G, bim, batch_size=1, preprocess=preprocess)
for _G, _bim in queue:
	try:
		r = st_sscan(_G, pheno, E, M=K, tests=["inter", "assoc"], verbose=False)
		bims.append(_bim)
		res.append(r)
	except:
		pass
    
    
# concatenate results
bims = pd.concat(bims).reset_index(drop=True)
res = pd.concat(res).reset_index(drop=True)
result=pd.concat([bims, res], axis=1)


#if not os.path.exists("out"):
#	os.makedirs("out")
result.to_csv(sys.argv[4], index=False)
