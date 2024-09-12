# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 02:41:55 2023

@author: qaal
"""


from sklearn.linear_model import LogisticRegression as lr
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from sklearn import metrics
import sys
import pathlib
import os


c

#---------------------------functions------------------------
def KNNeighbours(sigma,X):
    caliper = np.std(X.propensity_score) * sigma
    knn = NearestNeighbors(n_neighbors=15 , p = 2, radius=caliper)
    knn.fit(X[['propensity_score_logit']].to_numpy()) 
    
    distances , indexes = knn.kneighbors(X[['propensity_score_logit']].to_numpy(),n_neighbors=15)    
    return distances, indexes


def perfom_matching_v2(row, indexes, df_current):
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
        
        
def One2One(cohort,indexes,distances):
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
        
#---------------------------Main section------------------------    
root_dir = pathlib.Path.cwd()
result_dir = os.path.join(root_dir,'match_data')
os.makedirs(result_dir,exist_ok=True)

df_data = pd.read_csv('NewMDACCTrain.csv')
y = df_data[['prog_3_mo']]
df_data = df_data.drop(columns = ['drug'])

#df_psm = df_data[["patient_id","treatment","Gender","Age","Tobacco.Use","Pathology","PD.L1.expression","Line.of.IO_conden","Brain.met","Met.status","Liver.met"]]
df_psm = df_data[["patient_id","treatment","Gender","Age","Tobacco.Use","Pathology","PD.L1.expression","Line.of.IO_conden"]]
T = df_psm.treatment
#X = df_psm.loc[:,df_psm.columns !='treatment' !]
X = df_psm.drop(['treatment','patient_id'],axis=1)

# Design pipeline to build the treatment estimator
model = Pipeline([('scaler', StandardScaler()),('logistic_classifier', lr())])
model.fit(X, T)
predictions = model.predict_proba(X)
predictions_binary = model.predict(X)

print('Accuracy: {:.4f}\n'.format(metrics.accuracy_score(T, predictions_binary)))
print('Confusion matrix:\n{}\n'.format(metrics.confusion_matrix(T, predictions_binary)))
print('F1 score is: {:.4f}'.format(metrics.f1_score(T, predictions_binary)))

#Convert propability to logit
predictions_logit = np.array([logit(xi) for xi in predictions[:,1]])
# Density distribution of propensity score (logic) broken down by treatment status
fig, ax = plt.subplots(1,2, figsize=(11,7))
fig.suptitle('MDA Cohort: Density distribution plots for propensity score and logit(propensity score).')
sns.kdeplot(x = predictions[:,1], hue = T , ax = ax[0])
ax[0].set_title('Propensity Score')
sns.kdeplot(x = predictions_logit, hue = T , ax = ax[1])
ax[1].axvline(-0.4, ls='--')
ax[1].set_title('Logit of Propensity Score')
plt.show()

common_support = (predictions_logit > -10) & (predictions_logit < 10)

df_psm.loc[:,'propensity_score'] = predictions[:,1]
df_psm.loc[:,'propensity_score_logit'] = predictions_logit
df_psm.loc[:,'outcome'] = y.prog_3_mo
df_psm.loc[:,'patient_id'] = df_data.patient_id

X.loc[:,'propensity_score'] = predictions[:,1]
X.loc[:,'propensity_score_logit'] = predictions_logit
X.loc[:,'outcome'] = y.prog_3_mo
X.loc[:,'treatment'] = df_data.treatment
X.loc[:,'patient_id'] = df_data.patient_id

#---Start finding nearest neighbours-----------
sigma = 0.25
interested_pid = df_psm.patient_id
df_all_match = pd.DataFrame()
X_origin = X
for run in range(100):    
   #print('-----------------------')
    print('Resampling iter:',run)
    #print('-----------------------')
    temp =  df_psm.isin(interested_pid)
    temp = temp.patient_id
    df_current = df_psm[temp][df_psm.columns] 
    X = X_origin[temp][X_origin.columns] 

    
    df_current = df_current.reset_index(drop=True)
    X = X.reset_index(drop=True)
    distances, indexes = KNNeighbours(sigma,X)
    
    df_current['results'] = df_current.reset_index().apply(perfom_matching_v2, axis = 1, args = (indexes, df_current))
    val = ~df_current.results.isna()
    val = val.replace({True: 1, False: 0})
    if sum(val)>0:
        df_current[['matched_element','distance']] = df_current['results'].apply(pd.Series)
        df_current = df_current.drop(['results'],axis=1)

        #matched cohort
        matched_cohort = ~df_current.matched_element.isna()
        matched_cohort = df_current[matched_cohort][df_current.columns]

        #One to one matching
        temp_out = One2One(matched_cohort,indexes,distances)
    
        #conversion to patient id
        df_match = index2pid(temp_out,df_current)
        interested_pid = list(set(interested_pid)-set(df_match.Control_PID))
        interested_pid = list(set(interested_pid)-set(df_match.Treated_PID))
        
        df_all_match = pd.concat((df_all_match,df_match),axis=0)
        del indexes,distances
        
    else:
        print('No more matches..')
        #--Take out matched the data-----
        df_all_match = df_all_match[df_all_match['P_distance'] <0.12].reset_index(drop=True)
        control_cohort = df_data['patient_id'].isin(df_all_match.Control_PID)
        control_cohort = df_data[control_cohort][df_data.columns]
        treated_cohort = df_data['patient_id'].isin(df_all_match.Treated_PID)
        treated_cohort = df_data[treated_cohort][df_data.columns]
        matched_cohort = pd.concat((control_cohort,treated_cohort),axis=0)
        matched_cohort.to_csv(os.path.join(result_dir,'Matched_MDA.csv'))
        #--Take out unmatched the data-----
        unmatched_cohort = ~df_data['patient_id'].isin(matched_cohort.patient_id)
        unmatched_cohort = df_data[unmatched_cohort][df_data.columns]
        unmatched_cohort.to_csv(os.path.join(result_dir,'Unmatched_MDA.csv'))

        sys.exit()
        

            
