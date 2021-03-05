import multiprocessing as mp
import pandas as pd
import statsmodels.api as sm
# P value cutof for association
pval = 0.01
# Home directory
wd = "/Users/timvigers/Dropbox/School/MS Thesis/data/"
# Import data (from pickle_csvs.py)
psv = pd.read_pickle(wd+"raw_data/psv.pkl")
sv = pd.read_pickle(wd+"raw_data/sv.pkl")
# Get outcome
ia = psv["IAgroup2"]
ia = ia.replace({'case': 1,'control':0})
# Get lists of probes and metabolites to test - convert to strings
probes = pd.read_csv(wd+"raw_data/probes.csv",header=None)
probes = probes.astype(str).values.tolist()
metabolites = pd.read_csv(wd+"raw_data/metabolites.csv",header=None)
metabolites = metabolites.astype(str).values.tolist()
# Check probe association with outcome
probe_candidates = [] # To store results
for i in range(0,len(probes)):
    c = probes[i][0]
    logit_model=sm.Logit(ia,sm.add_constant(psv[[c]]))
    result=logit_model.fit(disp=0)
    if result.pvalues[[1]][0] < pval:
        probe_candidates.append(c)
# Now check metabolite association with outcome
metab_candidates = [] # To store results
for i in range(0,len(metabolites)):
    m = metabolites[i][0]
    logit_model=sm.Logit(ia,sm.add_constant(sv[[m]]))
    result=logit_model.fit(disp=0)
    if result.pvalues[[1]][0] < pval:
        metab_candidates.append(m)
# Check which pairs are associated with each other
methyl_psv_candidates = []
for i in range(0,len(probe_candidates)):
    p = probe_candidates[i]
    for j in range(0,len(metab_candidates)):
        m = metab_candidates[j]
        model=sm.OLS(sv[[m]],sm.add_constant(psv[[p]]))
        result=model.fit(disp=0)
        if result.pvalues[[1]][0] < pval:
            methyl_psv_candidates.append([p,m])
# Write results
pd.DataFrame(methyl_psv_candidates).to_csv(wd+"mediation/methyl_psv_py.csv",index=False)
# Same again, but switch psv and sv
# Check probe association with outcome
probe_candidates = [] # To store results
for i in range(0,len(probes)):
    c = probes[i][0]
    logit_model=sm.Logit(ia,sm.add_constant(sv[[c]]))
    result=logit_model.fit(disp=0)
    if result.pvalues[[1]][0] < pval:
        probe_candidates.append(c)
# Now check metabolite association with outcome
metab_candidates = [] # To store results
for i in range(0,len(metabolites)):
    m = metabolites[i][0]
    logit_model=sm.Logit(ia,sm.add_constant(psv[[m]]))
    result=logit_model.fit(disp=0)
    if result.pvalues[[1]][0] < pval:
        metab_candidates.append(m)
# Check which pairs are associated with each other
metab_psv_candidates = []
for i in range(0,len(probe_candidates)):
    p = probe_candidates[i]
    for j in range(0,len(metab_candidates)):
        m = metab_candidates[j]
        model=sm.OLS(psv[[m]],sm.add_constant(sv[[p]]))
        result=model.fit(disp=0)
        if result.pvalues[[1]][0] < pval:
            metab_psv_candidates.append([p,m])
# Write results
pd.DataFrame(metab_psv_candidates).to_csv(wd+"mediation/metab_psv_py.csv",index=False)