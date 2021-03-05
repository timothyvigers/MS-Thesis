import pandas as pd
df = pd.read_csv("/Users/timvigers/Dropbox/School/MS Thesis/data/raw_data/psv.csv")
df.to_pickle("/Users/timvigers/Dropbox/School/MS Thesis/data/raw_data/psv.pkl")
df = pd.read_csv("/Users/timvigers/Dropbox/School/MS Thesis/data/raw_data/sv.csv")
df.to_pickle("/Users/timvigers/Dropbox/School/MS Thesis/data/raw_data/sv.pkl")