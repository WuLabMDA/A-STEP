import pandas as pd
import numpy as np
#import plotly
import shap
import pathlib
from sklearn import linear_model
from sklearn.preprocessing import StandardScaler, MinMaxScaler

import fn_shap as calc


#----------------------------------------------------------------------
root_dir = pathlib.Path.cwd()
all_features = []
'''----------------------------
#       Loss B
----------------------------'''
df = pd.read_csv('SHAP_prep.csv')
df = df.drop(df.columns[0],axis=1)
Y = df['b.scores']
X = df.iloc[:,:-1]
all_features.append(X)

#--- measure shap values from feature matrix-----
model = linear_model.LinearRegression()
model.fit(X, Y) 
explainer = shap.Explainer(model.predict, X)
shap_values_A = explainer(X) 
shap_A = pd.DataFrame(shap_values_A.values,columns=X.columns)

#--- Define feature direction-----
k2 = calc.abs_shap(shap_values_A,X)
k2['Recom'] = 'NA'
for i in range(len(k2)):
    if k2['Magnitude'][i]< 0:
        k2['Recom'][i]='ICI-chemo'
    else:
        k2['Recom'][i]='ICI-mono'
        
        
k2 = k2.sort_values(by='Magnitude',ascending = True)
shap_val =  k2['Magnitude'].reset_index(drop=True)
shap_val = StandardScaler().fit_transform(shap_val.to_numpy().reshape(-1,1))
k2['SHAP_norm'] = shap_val
calc.plot_bar2(k2)


