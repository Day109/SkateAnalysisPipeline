import os, glob, sys
import pandas as pd
from bioinfokit.analys import norm, get_data

InputDf = sys.argv[1]
OutputDf = InputDf.split('.')[0]+'_RPKM.txt'

CntDf = pd.read_csv(InputDf,sep='\t')
CntDf = CntDf.set_index('Geneid')

CntDf = CntDf.iloc[:,4:]

nm = norm()
nm.rpkm(df=CntDf,gl='Length')
RpkmDf = nm.rpkm_norm
RpkmDf.to_csv(OutputDf,sep='\t')

