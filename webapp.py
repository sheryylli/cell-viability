import zipapp
import streamlit as st
import pandas as pd
import numpy as np
import time
import random
import scipy
from datetime import datetime
import statistics
import random
from scipy.stats import spearmanr
from scipy.stats import pearsonr
import plotly.express as px
import plotly.figure_factory as ff
from numpy import array
from random import randint
#setting up format
st.write("""
# Cell Viability and Fold Change""")
st.write('determining the relationship between cell viability and fold change.')
st.write(' Protein abundance and cell viability are measure for screens. Cell viability is how many cells are still alive at the time of measuring protein ')
 
st.write('To interpret screening results, we consider the protein measurement relative to cell viability. These values are then compared to average of controls to calculate a fold change')
st.sidebar.header('User Input Parameters')
sample_size = st.sidebar.slider('sample size', 100, 500, 100)
relationship = st.sidebar.selectbox("relationship between protein abundance and cell viability?",
            ("Yes", "No"))
ratio = st.sidebar.slider('ratio (cell viaibility : protein abundance)', 1, 5, 1)
range_of_viability_values = st.sidebar.slider('range of cell viability values', 1000, 30000, 15000)
variance = st.sidebar.slider('choose a range of noise to be added', 0, 100, 10)
run = st.selectbox("how many times this will be run",
        (1,10, 100, 1000))

def user_input_features():
#gets the user input features
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

def randomNumber(range_of_viability_values: int)-> int:
#gets a random number within range_of_viability_values set by user
    randNum = random.randint(1, range_of_viability_values)%(range_of_viability_values-1+1)+ variance
    return randNum

normalization_value = []
def test(sample_size: int):
#makes the list of normalization values and graphs 
    for i in range(sample_size):
        a = randomNumber(range_of_viability_values)
        normalization_value.append(a)

    np.random.normal(0, variance, size=(1, sample_size))


#r_s to represent protein/cell viability values
r_s = []
correlation_spearman = []
correlation_pearson = []
mu=0.0
std = 0.1
def gaussian_noise(ratio,mu,std : float) -> float:
#function used to add noise to cell viability and protein abundance values
    noise = np.random.normal(mu, std, size= None)
    x_noisy = ratio + noise
    return x_noisy 

rand_protein = randint(0, run-1)
fold_change_values = []
run_counter = 0
def proteinfunction(sample_size : int): 
#get protein values using previously determined cell viability values
    for i in range(sample_size):
        if relationship == "No":
            protein = randomNumber(range_of_viability_values)
        #if there is a relationship, now use the ratio
        if relationship == "Yes":
            protein = (normalization_value[i]) / (gaussian_noise(ratio, 0.0, variance))
    
        b  = protein / normalization_value[i]
        r_s.append(b)

    #determining C, the sum of all the values in r_s divided by the sample size
    total_control = 0
    for i in range(sample_size):
        total_control += r_s[i]

    c = (total_control) / sample_size

    fold_change_values = []
    for i in range(0, len(r_s)):
        fold_change = r_s[i] / c
        if fold_change / -1 > 0:
            fold_change = -1 * fold_change
        fold_change_values.append(fold_change)
  
    data = {'Normalization value': normalization_value[:],
            'Fold change': fold_change_values[:]
        }

    #if its the first run, or if the user chose only one run, print the fold change graph (since dont want every single graph)
    if (run_counter ==0) or (run==1):
        fig2 = px.scatter(data, x="Normalization value", y="Fold change",
                title="Fold Change ")
        st.plotly_chart(fig2)

    coef, p = spearmanr(normalization_value, fold_change_values)
    correlation_spearman.append(coef)

    pearsoncorr, _ = pearsonr(normalization_value, fold_change_values)
    correlation_pearson.append(pearsoncorr)
    #clear it
    r_s.clear()
    normalization_value.clear()
    
#running it for however times the user has decided
for i in range(run):
    test(sample_size)
    proteinfunction(sample_size)
    run_counter += 1

#summary stats
st.subheader("All of correlation values")
col1, col2 = st.columns(2)
with col1:
    data2 = {
            #'run': _,
            'spearman': correlation_spearman[:],
            'pearson' : correlation_pearson[:],
        }
    correlation_values = pd.DataFrame(data2)
    st.write(correlation_values)
#column on right shows summary stats for the entire run, each consisting of sample_size samples
with col2:
    pearson_average = sum(correlation_pearson) / len(correlation_pearson)
    spearman_average = sum(correlation_spearman) / len(correlation_spearman)
    #pearson
    st.header('summary stats for pearson')
    st.write('average pearson correlation: ' )
    st.write(pearson_average)
    st.write('median pearson correlation: ')
    st.write(statistics.median(correlation_pearson))
    st.write('range pearson correlation: ')
    st.write(max(correlation_pearson) - min(correlation_pearson))
    #spearman
    st.header('summary stats for spearman')
    st.write('average spearman correlation: ')
    st.write(spearman_average)
    st.write('median spearman correlation: ')
    st.write(statistics.median(correlation_spearman))
    st.write('range spearman correlation: ')
    st.write(max(correlation_spearman) - min(correlation_spearman))
st.header('Distribution')
# Add histogram data
print(len(correlation_spearman))
df3 = px.data.tips()
#histogram for correlation coefficient values
fig = px.histogram(data2, x=['spearman','pearson'])
st.plotly_chart(fig)
#scatterplot of correlation coefficient values
fig5 = px.scatter(data2, x=['spearman','pearson'])
st.plotly_chart(fig5)
