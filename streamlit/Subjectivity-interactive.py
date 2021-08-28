import streamlit as st
import streamlit.components.v1 as components
import os 

# Set configuration
st.set_page_config(layout = 'wide')

# Get image location
img_loc = os.getcwd() + '/streamlit/WebApp_legend.png'

# Write header
st.markdown(''' 
         *Interactive results*  
         # Hotspots of social and ecological impacts from freshwater stress and storage loss 
         **Huggins et al.** (submitted, 2021). Link to preprint to come.
         
         --- 
         
         This site allows you to explore the sensitivitiy of the social-ecological vulnerability hotspot reuslts to subjective aspects of the study's methodology. Select from the alternative methodological configurations below, and see how the hotspot results change in response. 
         
         <font style='font-size:12px' color='orange'> Note to Firefox users: you will need to interact with the dropdown menu below for results to appear. </font>
         
         ''', 
         unsafe_allow_html=True)

# Interactive map
components.iframe('https://xanderhuggins.shinyapps.io/Hotspot-web-app/?_ga=2.195577964.1941895401.1619383215-1450163622.1619383215',
                  height = 500,
                  scrolling = True)         

# Make columns for additional information
col1, col2 = st.beta_columns(2)

# More info
col1.markdown('''
Alternatives in **bold** are those used in the study.               
              
| Options        | Alternatives | 
| ------------- |:-------------:| 
| Scale  | [HydroBASINS](https://www.hydrosheds.org/page/hydrobasins) levels 3 OR **4** OR 5 | 
| Freshwater demand   | **Withdrawal** OR Consumption -- from [Huang et al. (2018)](https://zenodo.org/record/1209296#.YJqs3bVKj4Y) | 
| Streamflow | [**GSCD**](http://www.gloh2o.org/gscd/) (Beck et al., 2013) OR GRUN [(Ghiggi et al. 2019)](https://essd.copernicus.org/articles/11/1655/2019/) |
| Sensitivity aggregation| **Fuzzy sum** OR Mean |   
''')

# Legend
col2.image(img_loc, caption = None, width = 400)
