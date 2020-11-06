import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
# Working directory
dir = "C:/Users/tim/Documents/GitHub/MS-Thesis/"
# Read in data
data = pd.read_csv(dir+"data/networks/all_data_methyl_scaled.csv")
# Dummy variables
data = pd.get_dummies(data)
data = data.dropna(axis = 0, how = 'any')
# Outcome
labels = np.array(data["T1Dgroup_T1D case"])
# Drop outcome from dataframe
data = data.drop("T1Dgroup_T1D case", axis = 1)
# Save feature list
feature_list = list(data.columns)
# Convert to numpy array
data = np.array(data)
# Split into training and test
train_features, test_features, train_labels, test_labels = \
train_test_split(data, labels, test_size = 0.25, random_state = 1017)
# Train
tree = RandomForestClassifier(n_estimators = 1000,random_state=1017)
tree.fit(train_features, train_labels)
# Test
rf_predictions = tree.predict(test_features)
rf_probs = tree.predict_proba(test_features)[:, 1]
# Feature importance
importance = tree.feature_importances_
feat_importance = dict(zip(feature_list,importance))
# Write to CSV
df = pd.DataFrame(list(feat_importance.items()),columns = ['feature','importance'])
df.to_csv(dir+"data/networks/rf_variable_importance.csv",index=False)
