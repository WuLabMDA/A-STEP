import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import pathlib
import os
import fn_CCM as ccm

from sklearn.linear_model import LogisticRegression as lr
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from scipy.special import logit
import matplotlib.pyplot as plt
import seaborn as sns


#---Data prep----
root_dir = pathlib.Path.cwd()
df_data = pd.read_csv('sample_data.csv')
y = df_data['prog_3_mo']
T = df_data['treatment']
columns = ["Gender","Age","Tobacco.Use","Pathology","PD.L1.expression","Line.of.IO_conden"]
X = df_data[columns]

#---build a regression model----
model = Pipeline([('scaler', StandardScaler()),('logistic_classifier', lr())])
model.fit(X, T)
predictions = model.predict_proba(X)
predictions_binary = model.predict(X)
print('Accuracy: {:.4f}\n'.format(metrics.accuracy_score(T, predictions_binary)))

#Convert propability to logit
predictions_logit = np.array([logit(xi) for xi in predictions[:,1]])

#--plot density distribution---
fig, ax = plt.subplots(1,2, figsize=(11,7))
fig.suptitle('Density distribution plots for propensity score and logit(propensity score).')
sns.kdeplot(x = predictions[:,1], hue = T , ax = ax[0])
ax[0].set_title('Propensity Score')
sns.kdeplot(x = predictions_logit, hue = T , ax = ax[1])
ax[1].axvline(-0.4, ls='--')
ax[1].set_title('Logit of Propensity Score')
plt.show()

#---create a separate df---
df_psm = pd.DataFrame(df_data.patient_id)
df_psm.loc[:,'outcome'] = y
df_psm.loc[:,'treatment'] = T
df_psm.loc[:,'propensity_score'] = predictions[:,1]
df_psm.loc[:,'propensity_score_logit'] = predictions_logit

#---Start finding match pairs-----------
sigma = 0.25
interested_pid = df_psm.patient_id
df_all_match = pd.DataFrame()
distances, indexes = ccm.KNNeighbours(sigma,df_psm)
df_psm['results'] = df_psm.reset_index().apply(ccm.perfom_matching_v2, axis = 1, args = (indexes, df_psm,distances))
val = ~df_psm.results.isna()
val = val.replace({True: 1, False: 0})

if sum(val)>0:
    df_psm[['matched_element','distance']] = df_psm['results'].apply(pd.Series)
    df_psm = df_psm.drop(['results'],axis=1)
    matched_cohort = ~df_psm.matched_element.isna()
    matched_cohort = df_psm[matched_cohort][df_psm.columns]
    temp_out = ccm.One2One(matched_cohort,indexes,distances) #One to one matching
    #conversion to patient id
    df_match = ccm.index2pid(temp_out,df_psm)
else:
    
     print('No matches found..')
