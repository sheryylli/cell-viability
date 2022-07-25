import zipapp
import streamlit as st
import pandas as pd
import numpy as np
import time
import random
import scipy
from datetime import datetime
import statistics
random.seed(datetime.now())
#import seaborn as sns
import random
from scipy.stats import spearmanr
from scipy.stats import pearsonr
import plotly.express as px
import plotly.figure_factory as ff
from numpy import array


st.write("""
# Cell Viability and Fold Change

""")
st.write('determining the relationship between cell viability and fold change.')
st.write(' Protein abundance and cell viability are measure for screens. Cell viability is how many cells are still alive at the time of measuring protein ')
 
st.write('To interpret screening results, we consider the protein measurement relative to cell viability. These values are then compared to average of controls to calculate a fold change')
st.sidebar.header('User Input Parameters')
sample_size = st.sidebar.slider('sample size', 100, 500, 100)
relationship = st.sidebar.selectbox("relationship between protein abundance and cell viability?",
            ("Yes", "No"))
ratio = st.sidebar.slider('ratio (cell viaibility : protein abundance)', 1, 5, 1)
range_of_viability_values = st.sidebar.slider('range of cell viability values', 100, 30000, 25000)
variance = st.sidebar.slider('choose a range of noise to be added', 0, 100, 10)
run = st.selectbox("how many times this will be run",
        (1,10, 100, 1000))
#decide if user wants to see all fold change plots

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
st.subheader("""
Cell viability values
""")
progress_bar = st.sidebar.progress(0)
status_text = st.sidebar.empty()
last_rows = np.random.randn(1, 1)
chart = st.line_chart(last_rows)

def randomNumber(range_of_viability_values: int)-> int:
#gets a random number within the range
    randNum = random.randint(1, range_of_viability_values)%(range_of_viability_values-1+1)+ variance
    return randNum

normalization_value = []
def test(sample_size: int):
#makes the list of normalization values and graphs 
    for i in range(sample_size):
        a = randomNumber(range_of_viability_values)
        normalization_value.append(a)
        #if i ==0:
        new_rows = a+ np.random.randn(5, 1).cumsum(axis=0)
        
        #new_rows = last_rows[-1, :] + np.random.randn(5, 1).cumsum(axis=0)
        #status_text.text("%i%% Complete" % i)
        chart.add_rows(new_rows)
        #progress_bar.progress(i)
        last_rows = new_rows
        time.sleep(0.00)
    #add noise
    np.random.normal(0, variance, size=(1, sample_size))
    #print(normalization_value)

progress_bar.empty()
st.button("Re-run")
#protein values + plot
r_s = []
correlation_spearman = []
correlation_pearson = []
mu=0.0
std = 0.1
def gaussian_noise(ratio,mu,std : float) -> float:
    noise = np.random.normal(mu, std, size= None)
    x_noisy = ratio + noise
    return x_noisy 

show_fold = st.selectbox('see all fold change plots?', ['yes', 'no'])
def proteinfunction(sample_size : int): #type hint #doc for each  function
#get protein values
    #apply noise beforehand
    
    for i in range(sample_size):
        if relationship == "No":
            protein = randomNumber(range_of_viability_values) #add values here
        if relationship == "Yes":
            #np.random.normal adds noise
            #https://numpy.org/doc/stable/reference/random/generated/numpy.random.normal.html
            protein = (normalization_value[i]) / (gaussian_noise(ratio, 0.0, variance))
            #        + np.random.normal(0, variance, size=(1,0))) / ratio

            #write it out mathematically
            #apply noise beforehand
            

        #cell_viability = cell_values[i]
        b  = protein / normalization_value[i]
        r_s.append(b)

        
    #plot protein values?
    total_control = 0
    for i in range(sample_size):
        total_control += r_s[i]

    c = (total_control) / sample_size

    fold_change_values = []
    for i in range(0, len(r_s)):
        fold_change = r_s[i] / c
        fold_change_values.append(fold_change)
    #np.random.normal(0, variance, size=(1, sample_size))

    data = {'Normalization value': normalization_value[:],
            'Fold change': fold_change_values[:]
        }

    #print("FOLD CHANGE")
    print(len(fold_change_values)) #100

    #fig2 = px.scatter(data, x="Normalization value", y="Fold change",
    #            title="Fold Change ")
    #st.plotly_chart(fig2)

    #return fold_change_values
        
    
    if show_fold == 'yes':
        fig2 = px.scatter(data, x="Normalization value", y="Fold change",
                title="Fold Change ")
        st.plotly_chart(fig2)
    #st.subheader('Correlation')
    for _ in range(0, len(r_s)):
        coef, p = spearmanr(normalization_value, fold_change_values)
    #st.write('Spearmans correlation coefficient: %.3f' % coef)
        correlation_spearman.append(coef)
        print("spearman")
        print(coef)
        
    #alpha = 0.05
    #if p > alpha:
    #    st.write('Samples are uncorrelated (fail to reject H0) p=%.3f' % p)
    #else:
    #    st.write('Samples are correlated (reject H0) p=%.3f' % p)
    

  #pearson
        pearsoncorr, _ = pearsonr(normalization_value, fold_change_values)
        print("pearson")
        print(pearsoncorr)
    #st.write('Pearsons correlation: %.3f' % pearsoncorr)
        correlation_pearson.append(pearsoncorr)

    #clear it
    r_s.clear()
    normalization_value.clear()
#getting dataframe of corr values

#running it like 1000 times:
for _ in range(run):
    test(sample_size)
    proteinfunction(sample_size)
    
st.subheader("All of correlation values")
col1, col2 = st.columns(2)
with col1:
    data2 = {
            #'run': _,
            'spearman': correlation_spearman[:],
            'pearson' : correlation_pearson[:],
        }
    correlation_values = pd.DataFrame(data2)
                                                # index=[0])
    st.write(correlation_values)
with col2:
    pearson_average = sum(correlation_pearson) / len(correlation_pearson)
    spearman_average = sum(correlation_spearman) / len(correlation_spearman)
    st.header('summary stats for pearson')
    st.write('average pearson correlation: ' )
    st.write(pearson_average)
    st.write('median pearson correlation: ')
    st.write(statistics.median(correlation_pearson))
    st.write('range pearson correlation: ')
    st.write(max(correlation_pearson) - min(correlation_pearson))
    st.header('summary stats for spearman')
    st.write('average spearman correlation: ')
    st.write(spearman_average)

    #st.write(correlation_pearson.describe())

    st.write('median spearman correlation: ')
    st.write(statistics.median(correlation_spearman))
    st.write('range spearman correlation: ')
    st.write(max(correlation_spearman) - min(correlation_spearman))
    #st.write(correlation_spearman.describe())
st.header('Distribution')
# Add histogram data
# Group data together
#hist_data = [x1, x2]
#hist_data = [data2[0], data2[1]]
    
#df = px.data.tips()
#fig = px.histogram(data_frame = data2
   #                 ,x = 'pearson'
   #         )
#fig.show()
#fig = px.histogram(data2, x="correlation_spearman")
#fig = px.scatter(data2, correlation_spearman, y)
#st.plotly_chart(fig)
# Add histogram data
#x1 = np.random.randn(200) - 2
#x2 = np.random.randn(200)
#x3 = np.random.randn(200) + 2

# Group data together
#hist_data = [x1, x2, x3]

#group_labels = ['Group 1', 'Group 2', 'Group 3']

# Create distplot with custom bin_size
#fig = ff.create_distplot(
         #hist_data, group_labels, bin_size=[.1, .25, .5])

# Plot!

#st.plotly_chart(fig, use_container_width=True)
print(len(correlation_spearman))

#for i in range(run):
#    for j in range(run):
       # print(correlation_spearman[i][j])

df3 = px.data.tips()
fig = px.histogram(data2, x=['spearman','pearson'])

#fig.show()
st.plotly_chart(fig)
