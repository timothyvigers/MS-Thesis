import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
# Working directory
dir = "C:/Users/tim/Documents/GitHub/MS-Thesis/"
# Read in data
data = pd.read_csv(dir+"data/networks/all_data_methyl_scaled.csv")
# Dummy variables
data = pd.get_dummies(data)
# Outcome
labels = np.array(data["T1Dgroup_T1D case"])
# Drop outcome from dataframe
data = data.drop("T1Dgroup_T1D case", axis = 1)
# Save feature list
feature_list = list(data.columns)
# Convert to numpy array
data = np.array(data)
