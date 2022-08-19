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
import math
st.write("""
# Cell Viability and Fold Change""")
st.write('Determining the relationship between cell viability and fold change.')
st.write(' Protein abundance and cell viability are measure for screens. Cell viability is how many cells are still alive at the time of measuring protein ')
 
st.write('To interpret screening results, we consider the protein measurement relative to cell viability. These values are then compared to average of controls to calculate a fold change')

st.sidebar.header('User Input Parameters')

sample_size = st.sidebar.slider('Sample size', 100, 500, 100)
relationship = st.sidebar.selectbox("Relationship between protein abundance and cell viability?",("Yes", "No"))
if relationship == "Yes":
    ratio = st.sidebar.slider('Ratio (cell viaibility : protein abundance)', 1, 5, 1)
    noise_type = st.sidebar.selectbox("Which type of noise will be added", ("Gaussian",))
range_of_viability_values = st.sidebar.slider('Range of cell viability values', 1000, 30000, 25000)
lower_bound = st.sidebar.slider('Lower end of the cell viability range', 0, 10000, 0)
variance = st.sidebar.slider('Choose a range of noise to be added', 0, 10000, 10)
run = st.sidebar.selectbox("How many times this will be run", (1,10, 100, 1000))
num_controls = st.sidebar.selectbox("How many controls there will be", (1, 10, 100, 1000))

#protein value should never be negative- reject it if its negative and just get a new value thats positive

def user_input_features():
#gets the user input features
    data = {
        'sample size': sample_size,
        'range of cell viability values': range_of_viability_values,
        'lower end of range': lower_bound,
        'range of noise added': variance,
        'relationship' : relationship
        
    }
    features = pd.DataFrame(data, index=[0])
    return features

df = user_input_features()
st.subheader('User Input parameters')
st.write(df)

#generate noice function
mu=0.0
std = math.sqrt(variance)
def noise(mu, variance: float)-> float:
    if relationship == "Yes":
        if noise_type == "Gaussian":
            noise = np.random.normal(mu, std, size= None)
        
    if relationship == "No":
         noise = np.random.normal(mu, std, size= None)
    return noise

def randomNumber(range_of_viability_values, lower_bound):
#gets a random number within range_of_viability_values set by user
    randNum = random.randint(lower_bound, range_of_viability_values)
    return randNum


normalization_value = []
def test(sample_size: int):
#makes the list of normalization values and graphs 
    for _ in range(sample_size):
        a = randomNumber(range_of_viability_values, lower_bound) + noise(0.0, variance)
        while a <= 0:
            a = randomNumber(range_of_viability_values, lower_bound) + noise(0.0, variance)
        normalization_value.append(a)
    
r_s = []
correlation_spearman = []
correlation_pearson = []
rand_protein = randint(0, run-1) #for printing out a random run of fold change values
run_counter = 0
fold_change_values = []
r_s_noise = []
fold_change_values_noise=[]

def getcorrelations(normalization, fold_change_list):
    coef, p = spearmanr(normalization, fold_change_list)
    correlation_spearman.append(coef)

    pearsoncorr, _ = pearsonr(normalization, fold_change_list)
    correlation_pearson.append(pearsoncorr)
   
    data2 = {
                #'run': _,
                'spearman': correlation_spearman[:],
                'pearson' : correlation_pearson[:],
            }
    correlation_values = pd.DataFrame(data2)
    #also add a plot graph of values of cell viability and protein abundance values
    if run_counter == run-1:
        col1, col2 = st.columns(2)
        with col1:
            
            st.header("All of correlation values")
            
            
            st.write(correlation_values)
        #column on right shows summary stats for the entire run, each consisting of sample_size samples
        with col2:
            pearson_average = sum(correlation_pearson) / len(correlation_pearson)
            spearman_average = sum(correlation_spearman) / len(correlation_spearman)
            #pearso\

            st.header('Summary stats for pearson')
            st.write('Average pearson correlation: ' )
            st.write(pearson_average)
            st.write('Median pearson correlation: ')
            st.write(statistics.median(correlation_pearson))
            st.write('Range pearson correlation: ')
            st.write(max(correlation_pearson) - min(correlation_pearson))
            #spearman
            st.header('Summary stats for spearman')
            st.write('Average spearman correlation: ')
            st.write(spearman_average)
            st.write('Median spearman correlation: ')
            st.write(statistics.median(correlation_spearman))
            st.write('Range spearman correlation: ')
            st.write(max(correlation_spearman) - min(correlation_spearman))
        st.header('Distribution')
        # Add histogram data
        #print(len(correlation_spearman))
        df3 = px.data.tips()
        #histogram for correlation coefficient values
        #scatterplot of correlation coefficient values
        fig5 = px.scatter(data2, x=['spearman','pearson'])
        st.plotly_chart(fig5)
   
    
    return


    

    
#describe methodology
#if theres no relationship, then do expect a correlation
#fix the webapp 
#thursday: focus on report + presentation


#as normalization values increase, why are the fold change values so low

#plot showing protein value and normalization value next to fold change plot
protein_values =[]
def no_relationship(sample_size : int): 
  
    for i in range(sample_size):
        protein = randomNumber(range_of_viability_values, lower_bound) + noise(0.0, variance)
        while protein <= 0:
            protein = randomNumber(range_of_viability_values, lower_bound) + noise(0.0, variance)
        protein_values.append(protein)
        a  = protein / normalization_value[i]
        r_s.append(a)
        
    c = (np.mean(r_s[:num_controls]))

    fold_change_values = []

    for i in range(0, sample_size):
        fold_change = r_s[i] / c
        fold_change_values.append(fold_change)

    data = {'Normalization value': normalization_value,
            'Fold change': fold_change_values
        }
    data2 = {'Normalization value': normalization_value,
            'Protein abundance': protein_values }
    
    if (run_counter ==0) or (run==1):
        fig1 = px.scatter(data2, x = "Normalization value", y = "Protein abundance",
                title= "Cell viability and protein abundance values")
        st.plotly_chart(fig1)
        st.header("Fold Change")
        fig2 = px.scatter(data, x="Normalization value", y="Fold change",
                title="Fold Change (noise added to only ratio) with noise (normal distribution)")
        st.plotly_chart(fig2)

    
    getcorrelations(normalization_value, fold_change_values)
    r_s.clear()
    normalization_value.clear()
    

def yes_relationship(sample_size):

    for i in range(sample_size):
        #noise added to overall 
        noise_protein = ( (normalization_value[i]) / ratio) + noise(0.0, variance)
        while noise_protein< 0:
            noise_protein = ( (normalization_value[i]) / ratio) + noise(0.0, variance)
        b = (noise_protein / normalization_value[i]) 
        
        protein_values.append(noise_protein)
        r_s_noise.append(b)



    d = (np.mean(r_s_noise[:num_controls]))
        
    
    fold_change_values_noise = []


    for i in range(0, sample_size):
     
        fold_change_noise = r_s_noise[i] / d
        fold_change_values_noise.append(fold_change_noise)
   
    
    data = {'Normalization value': normalization_value,
            'Protein abundance': protein_values
    }
    data_log = {'Normalization value': normalization_value,
        'Fold change': fold_change_values_noise
    }

    #if its the first run, or if the user chose only one run, print the fold change graph (since dont want every single graph)
    if run_counter ==0 or run==1:
        fig1 = px.scatter(data, x = "Normalization value", y = "Protein abundance",
                title= "Cell viability and protein abundance values")
        st.plotly_chart(fig1)
        st.header("Fold Change")
        fig3 = px.scatter(data_log, x="Normalization value", y="Fold change",
                title="Fold Change (noise added to overall value) )")
        st.plotly_chart(fig3)
    
    getcorrelations(normalization_value, fold_change_values_noise)
    r_s.clear()
    r_s_noise.clear()
    normalization_value.clear()

   
#running it for however times the user has decided
for i in range(run):
    test(sample_size)
    if relationship == "Yes":
        yes_relationship(sample_size)
        
        
    if relationship == "No":
        no_relationship(sample_size)
        
    run_counter = run_counter + 1
 

fold_change_values.clear()
fold_change_values_noise.clear()










