import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

methyl = pd.read_csv("/Volumes/Tim/thesis/data/methylation/450k.csv",nrows = 100)
metab = pd.read_csv("/Volumes/Tim/thesis/data/metabolomics/gctof.bc.csv",nrows = 100)
data = pd.concat([methyl,metab],join = "inner")
print(data) # Empty

# for meth in methyl.columns.values:
# 	print(meth + "~ gctof_1")
# 	mod = smf.mixedlm(meth + "~gctof_1",)
# 	print(meth)
