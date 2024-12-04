from sklearn.neighbors import NearestNeighbors

import sys
import pathlib
import os


# Enabled to remove warnings for demo purposes.
import warnings
warnings.filterwarnings('ignore')


import math
import numpy as np
import pandas as pd
import statsmodels.api as sm

#---------------------------functions------------------------
def KNNeighbours(sigma,X):
    caliper = np.std(X.propensity_score) * sigma
    knn = NearestNeighbors(n_neighbors=15 , p = 2, radius=caliper)
    knn.fit(X[['propensity_score_logit']].to_numpy()) 
    
    distances , indexes = knn.kneighbors(X[['propensity_score_logit']].to_numpy(),n_neighbors=15)    
    return distances, indexes


def perfom_matching_v2(row, indexes, df_current,distances):
    current_index = int(row['index']) # Obtain value from index-named column, not the actual DF index.
    prop_score_logit = row['propensity_score_logit']
    #curr_distances = distances[current_index,:]
    #print(current_index)
    num = 0
    for idx in indexes[current_index,:]:  
        if (current_index != idx) and (row.treatment == 1) and (df_current.loc[idx].treatment == 0):
            dist = distances[current_index,num]
            num = num+1
            return int(idx),dist
        else:
            num = num +1
        
        
def One2One(matched_cohort,indexes,distances):
    control_sample =[]
    treated_sample =[]
    pairwise_distance = []
    for treat_idx in (matched_cohort.index):
        num = 0
        match_idx = matched_cohort['matched_element'][treat_idx]
        test = matched_cohort['matched_element'] == match_idx
        test = test.replace({True: 1, False: 0})
        arr = test.to_numpy().nonzero()
        candidate = np.zeros([len(arr[0]),2],dtype=float)
        
        for k in range(len(test)):
            if test.iloc[k] == 1:
                candidate[num,0] = test.index[k]
                candidate[num,1] =  matched_cohort['distance'][matched_cohort.index[k]]
                num =  num +1
        candidate = candidate[candidate[:, 1].argsort()]
    
        if (match_idx not in control_sample) and (candidate[0,0] not in treated_sample):
            control_sample.append(match_idx)
            treated_sample.append(candidate[0,0]) # the first sorted rows (min value)
            pairwise_distance.append(candidate[0,1]) # the first sorted rows (min value)
    
    df_match = pd.DataFrame(control_sample,columns=['Control_idx'])
    df_match['Treated_idx'] = treated_sample
    df_match['P_distance'] = pairwise_distance
    
    return df_match

def index2pid(temp_out,df_current):
    control = df_current.loc[temp_out.Control_idx]
    treated = df_current.loc[temp_out.Treated_idx]
    
    control_pid = control.patient_id
    control_pid = control_pid.reset_index(drop=True)
    treated_pid = treated.patient_id
    treated_pid = treated_pid.reset_index(drop=True)
    
    df_match = pd.DataFrame()
    df_match['Control_PID'] = control_pid
    df_match['Treated_PID'] = treated_pid
    df_match['P_distance'] = temp_out.P_distance
    
    return df_match
        
