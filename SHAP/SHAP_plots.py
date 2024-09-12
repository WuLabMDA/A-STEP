import pandas as pd
import numpy as np
#import plotly
import shap
import pathlib
from sklearn import linear_model
import matplotlib.patches as mpatches
from sklearn.preprocessing import StandardScaler, MinMaxScaler

import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

def shap_comp(shap_v,all_features):
    feature_list = all_features.columns
    shap_v.columns = feature_list
    df_v = all_features.copy().reset_index().drop('index',axis=1)
    
    # Determine the correlation in order to plot with different colors
    corr_list = list()
    for i in feature_list:
        b = np.corrcoef(shap_v[i],df_v[i])[1][0]
        corr_list.append(b)
    corr_df = pd.concat([pd.Series(feature_list),pd.Series(corr_list)],axis=1).fillna(0)
    corr_df.columns  = ['Variable','Corr']
    corr_df['Color'] = np.where(corr_df['Corr']>0,'red','blue')
    corr_df['Direction'] = np.where(corr_df['Corr']>0,'High-value','Low-value')
    
    k=pd.DataFrame(shap_v.mean()).reset_index()
    k.columns = ['Variable','Magnitude']
    k2 = k.merge(corr_df,left_on = 'Variable',right_on='Variable',how='inner')
    
    #--- to set reference to all high values--
    for i in range(len(k2)):
        if k2.Direction[i] == 'Low-value':
            k2.Magnitude[i] = k2.Magnitude[i] *-1
            k2.Corr[i] = k2.Corr[i]*-1    
            k2.Color[i] = 'red'   
            k2.Direction[i] = 'High-value'  

    return k2

def abs_shap(explainer,df):
    shap_v = explainer.values
    shap_v = pd.DataFrame(shap_v)
    feature_list = df.columns
    shap_v.columns = feature_list
    df_v = df.copy().reset_index().drop('index',axis=1)
    
    # Determine the correlation in order to plot with different colors
    corr_list = list()
    for i in feature_list:
        b = np.corrcoef(shap_v[i],df_v[i])[1][0]
        corr_list.append(b)
    corr_df = pd.concat([pd.Series(feature_list),pd.Series(corr_list)],axis=1).fillna(0)
    corr_df.columns  = ['Variable','Corr']
    corr_df['Color'] = np.where(corr_df['Corr']>0,'red','blue')
    corr_df['Direction'] = np.where(corr_df['Corr']>0,'High-value','Low-value')
    

    k=pd.DataFrame(shap_v.mean()).reset_index()
    k.columns = ['Variable','Magnitude']
    k2 = k.merge(corr_df,left_on = 'Variable',right_on='Variable',how='inner')
    
    #--- to set reference to all high values--
    for i in range(len(k2)):
        if k2.Direction[i] == 'Low-value':
            k2.Magnitude[i] = k2.Magnitude[i] *-1
            k2.Corr[i] = k2.Corr[i]*-1    
            k2.Color[i] = 'red'   
            k2.Direction[i] = 'High-value'                   

    return k2


def plot_bar2(k2):

    ax = k2.plot.barh(x='Variable',y='SHAP_norm', figsize=(6,5),legend=False,color=k2['SHAP_norm'].apply(lambda x: 'indianred' if x > 0 else 'royalblue'))
    ax.set_xlabel('Different of Average Treatment Effects (ICI-mono - ICI-chemo)')
    red_patch = mpatches.Patch(color='indianred', label='ICI-mono')
    blue_patch = mpatches.Patch(color='royalblue', label='ICI-chemo')
    plt.legend(handles=[red_patch,blue_patch])
    plt.show()

#----------------------------------------------------------------------
root_dir = pathlib.Path.cwd()
all_features = []
'''----------------------------
#       Loss A
----------------------------'''
df = pd.read_csv('SHAP_A.csv')
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
k2 = abs_shap(shap_values_A,X)
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
plot_bar2(k2)


'''----------------------------
#       Loss B
----------------------------'''

df = pd.read_csv('SHAP_B.csv')
df = df.drop(df.columns[0],axis=1)
Y = df['b.scores']
X = df.iloc[:,:-1]
all_features.append(X)

#--- measure shap values from feature matrix-----
model = linear_model.LinearRegression()
model.fit(X, Y) 
explainer = shap.Explainer(model.predict, X)
shap_values_B = explainer(X)
shap_B = pd.DataFrame(shap_values_B.values,columns=X.columns)


#--- Define feature direction-----
k2 = abs_shap(shap_values_B,X)
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
plot_bar2(k2)


'''----------------------------
#       Loss C
----------------------------'''

df = pd.read_csv('SHAP_C.csv')
df = df.drop(df.columns[0],axis=1)
Y = df['b.scores']
X = df.iloc[:,:-1]
all_features.append(X)

#--- measure shap values from feature matrix-----
model = linear_model.LinearRegression()
model.fit(X, Y) 
explainer = shap.Explainer(model.predict, X)
shap_values_C = explainer(X)
shap_C = pd.DataFrame(shap_values_C.values,columns=X.columns)


#--- Define feature direction-----
k2 = abs_shap(shap_values_C,X)
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
plot_bar2(k2)


'''----------------------------
#       Loss D
----------------------------'''

df = pd.read_csv('SHAP_D.csv')
df = df.drop(df.columns[0],axis=1)
Y = df['b.scores']
X = df.iloc[:,:-1]
all_features.append(X)

#--- measure shap values from feature matrix-----
model = linear_model.LinearRegression()
model.fit(X, Y) 
explainer = shap.Explainer(model.predict, X)
shap_values_D = explainer(X)
shap_D = pd.DataFrame(shap_values_D.values,columns=X.columns)


#--- Define feature direction-----
k2 = abs_shap(shap_values_D,X)
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
plot_bar2(k2)


'''----------------------------
#       Loss E
----------------------------'''

df = pd.read_csv('SHAP_E.csv')
df = df.drop(df.columns[0],axis=1)
Y = df['b.scores']
X = df.iloc[:,:-1]
all_features.append(X)

#--- measure shap values from feature matrix-----
model = linear_model.LinearRegression()
model.fit(X, Y) 
explainer = shap.Explainer(model.predict, X)
shap_values_E = explainer(X)
shap_E = pd.DataFrame(shap_values_E.values,columns=X.columns)


#--- Define feature direction-----
k2 = abs_shap(shap_values_E,X)
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
plot_bar2(k2)


'''----------------------------
#       composite
----------------------------'''
w = [0.002,0.1,0.0,0.007,0.65]
shap_A = shap_A*w[0]
shap_B = shap_B*w[0]
shap_C = shap_C*w[0]
shap_D = shap_D*w[0]
shap_E = shap_E*w[0]
all_features = np.hstack(all_features)

shap_all = pd.concat((shap_A,shap_B,shap_C,shap_D,shap_E),axis=1)
df_new = []
col_names =[]

for i in range(shap_all.shape[1]):
    curr_features = shap_all[shap_all.columns[i]]
    if isinstance(curr_features,pd.DataFrame):
        curr_features = curr_features.sum(axis=1)
    df_new.append(curr_features)
    col_names.append(shap_all.columns[i])

df_new = np.vstack(df_new)
df_new = df_new.transpose()
shap_v = pd.DataFrame(df_new,columns=col_names)
all_features = pd.DataFrame(all_features,columns=col_names)


shap_v = shap_v.loc[:,~shap_v.columns.duplicated()].copy()
all_features = all_features.loc[:,~all_features.columns.duplicated()].copy()
k2 = shap_comp(shap_v,all_features)

k2['Recom'] = 'NA'
for i in range(len(k2)):
    if k2['Corr'][i]< 0:
        k2['Recom'][i]='ICI-chemo'
    else:
        k2['Recom'][i]='ICI-mono'
        
        
k2 = k2.sort_values(by='Magnitude',ascending = True)
shap_val =  k2['Magnitude'].reset_index(drop=True)
shap_val = StandardScaler().fit_transform(shap_val.to_numpy().reshape(-1,1))
k2['SHAP_norm'] = shap_val
plot_bar2(k2)
    



