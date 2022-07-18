import zipapp
import streamlit as st
import pandas as pd
import numpy as np
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
relationship = st.sidebar.selectbox("relationship between protein abundance and cell viability?",
            ("Yes", "No"))
ratio = st.sidebar.slider('ratio (cell viaibility : protein abundance)', 1, 5, 1)
range_of_viability_values = st.sidebar.slider('range of cell viability values', 100, 30000, 25000)
variance = st.sidebar.slider('choose a range of noise to be added', 0, 100, 10)
run = st.selectbox("how many times this will be run",
        (1,10, 100, 1000))
def user_input_features():

    data = {
        'sample size': sample_size,
        'ratio': ratio,
        'range of cell viability values': range_of_viability_values,
        'range of noise added': variance,
        'relationship' : relationship
        
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
    randNum = random.randint(1, range_of_viability_values)%(range_of_viability_values-1+1)+ variance
    return randNum

normalization_value = []
def test(sample_size):
    for i in range(sample_size):
        a = randomNumber(range_of_viability_values)
        normalization_value.append(a)
        #if i ==0:
        new_rows = a+ np.random.randn(5, 1).cumsum(axis=0)
        
        #new_rows = last_rows[-1, :] + np.random.randn(5, 1).cumsum(axis=0)
        status_text.text("%i%% Complete" % i)
        chart.add_rows(new_rows)
        progress_bar.progress(i)
        last_rows = new_rows
        time.sleep(0.00)

progress_bar.empty()
st.button("Re-run")
#protein values + plot
r_s = []
correlation_spearman = []
correlation_pearson = []
def proteinfunction(sample_size):
#get protein values
    for i in range(sample_size):
        if relationship == "No":
            protein = randomNumber(range_of_viability_values)
        if relationship == "Yes":
            #probably add some noise here
            #protein_value = (normalization_value[i]/ ratio) 
            #noise = np.random.randint(1, 10)
            #protein = protein_value + noise
            protein = (normalization_value[i] + variance) / ratio

        #cell_viability = cell_values[i]
        b  = protein / normalization_value[i]
        r_s.append(b)
        
    #plot protein values?
    total_control = 0
    for i in range(sample_size):
        #random_choice = random.choice(r_s)
        #total_control += random_choice
        total_control += r_s[i]

    c = (total_control) / sample_size

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
    correlation_spearman.append(p)
    alpha = 0.05
    if p > alpha:
        st.write('Samples are uncorrelated (fail to reject H0) p=%.3f' % p)
    else:
        st.write('Samples are correlated (reject H0) p=%.3f' % p)
    

  #pearson
    pearsoncorr, _ = pearsonr(normalization_value, fold_change_values)
    st.write('Pearsons correlation: %.3f' % pearsoncorr)
    correlation_pearson.append(pearsoncorr)

    #clear it
    r_s.clear()
    normalization_value.clear()

def average(l):
    avg = sum(l) / len(l)
    return avg
#running it like 1000 times:
for _ in range(run):
    test(sample_size)
    proteinfunction(sample_size)

st.subheader("All of correlation values")
pearson_average = average(correlation_pearson)
spearman_average = average(correlation_spearman)
st.write('average pearson correlation: ' )
st.write(pearson_average)
st.write('average spearman correlation: ')
st.write(spearman_average)
