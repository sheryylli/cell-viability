import zipapp
import streamlit as st
import pandas as pd
import numpy as np
from sklearn import datasets
from sklearn.ensemble import RandomForestClassifier
import time
import random
from datetime import datetime
random.seed(datetime.now())
import random
import math
from scipy.stats import spearmanr
from scipy.stats import pearsonr
import plotly.express as px
st.write("""
# Cell Viability and Fold Change
""")
st.sidebar.header('User Input Parameters')
sample_size = st.sidebar.slider('sample size', 100, 500, 100)
ratio = st.sidebar.slider('ratio (cell viaibility : protein abundance)', 1, 5, 1)
range_of_viability_values = st.sidebar.slider('range of cell viability values', 100, 30000, 25000)
#range_variance = st.sidebar.slider('range of nosie added', 0, 10, 2)

#running a thousand runs- range of correlation (summarizing) //without a button

def user_input_features():

    data = {
        'sample size': sample_size,
        'ratio (cell viability : protein abundance)': ratio,
        'range of cell viability values': range_of_viability_values
        
    }
    features = pd.DataFrame(data, index=[0])
    return features

df = user_input_features()

st.subheader('User Input parameters')
st.write(df)
st.subheader("""
Cell viability values
""")
progress_bar = st.sidebar.progress(0)
status_text = st.sidebar.empty()
last_rows = np.random.randn(1, 1)
chart = st.line_chart(last_rows)


def randomNumber(range_of_viability_values):
    randNum = random.randint(1, range_of_viability_values)%(range_of_viability_values-1+1)+1 #<- set this as a parameter
    return randNum

for i in range(1, sample_size):
    #make the second one a parameter
    new_rows = randomNumber(range_of_viability_values)+ np.random.randn(0, 5)
    #status_text.text("%i%% Complete" % i)
    chart.add_rows(new_rows)
    
    #progress_bar.progress(i)
    last_rows = new_rows
    time.sleep(0.00)

progress_bar.empty()
st.button("Re-run")
r_s = []
normalization_value = []

#Plot 

for _ in range(sample_size):
    protein = cell_viability / ratio #relationship between protein and cell viability but also add noise here
    #select mode / relationship

    #no relationship:
    #protein = random number()

    #slider for range for protein values 




    cell_viability = randomNumber(range_of_viability_values)
    #in which this is the case:
 
    normalization_value.append(cell_viability)
    a  = protein / cell_viability
    r_s.append(a)

control = 10
total_control = 0
for i in range(sample_size):
    random_choice = random.choice(r_s)
    total_control += random_choice

c = (total_control) / control

fold_change_values = []
for i in range(0, len(r_s)):
    fold_change = r_s[i] / c
    fold_change_values.append(fold_change)

  #plot

data = {'Normalization value': normalization_value[:],
        'Fold change': fold_change_values[:]
    }
    
df2 = pd.DataFrame(data,columns=['Normalization value','Fold change'])
df2.plot(x ='Normalization value', y='Fold change', kind = 'scatter')
  


  #return fold_change_values
st.subheader('Fold Change')
plot = px.scatter(df2, normalization_value, fold_change_values)
st.plotly_chart(plot)

#calculate the correlation
st.subheader('Correlation')
coef, p = spearmanr(normalization_value, fold_change_values)
st.write('Spearmans correlation coefficient: %.3f' % coef)
alpha = 0.05
if p > alpha:
    st.write('Samples are uncorrelated (fail to reject H0) p=%.3f' % p)
else:
    st.write('Samples are correlated (reject H0) p=%.3f' % p)
  

  #pearson
pearsoncorr, _ = pearsonr(normalization_value, fold_change_values)
st.write('Pearsons correlation: %.3f' % pearsoncorr)
