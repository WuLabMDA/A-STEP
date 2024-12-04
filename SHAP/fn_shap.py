import pandas as pd
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

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
