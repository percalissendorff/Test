## companion_frequency_test.py ##
#################################

import streamlit as st
import numpy as np
from scipy import integrate
from streamlit import session_state

st.set_page_config(
	page_title="CompaionFrequency",
	page_icon="",
	layout="wide",
	)

# Create a text input widget
host_mass = st.number_input("Host star [Solar mass]:", value=1.0)
Jup_min= st.number_input("Min companion mass [Jupiter mass]:", value=0.3)
Jup_max = st.number_input("Max companion mass [Jupiter mass]:", value=3.0)
a_min = st.number_input("a min [AU]:", value=1.0)
a_max = st.number_input("a max [AU]:", value=10.0)


## Do calculations
q_Jupiter = 0.0009545942339693249/host_mass # Msun
beta = -0.36 
# mean_bd = 1.70 is the mean of the log-10-normal in log-AU log(50) Raghavan et al. 2010 (physical separation) for FGK stars
# mean_bd = 1.43 for M dwarfs log-10(27) Winters et al. (2019) - physical separation not projected. Meyer et al. 2024 
# mean_bd = 2.72 for A stars log-10(525) De Rosa et al. (2014) - physical separation not projected. Meyer et al. 2024
mean_bd = 1.43
# sigma_bd = 1.21 M dwarf Winters et al. (2019) aghavan et al. 20
# sigmba_bd = 0.79 A stars De Rosa et al. (2014)
sigma_bd = 1.21
A_bd_ln = -3.78 # with 68 % CI [-0.47,0.7]
A_bd = np.exp(A_bd_ln)
alpha_pl = 1.43 # 68 % CI [-0.13,0.12] 
A_pl_ln = -5.52 # 68 % CI [-0.83,0.93] 
mu_ln = 1.32  # 68 % CI [-0.22,0.18] 
sigma_pl_ln = 0.53 # 68 % CI [-0.07,0.06]
A_pl = np.exp(A_pl_ln)
mu_pl = (mu_ln)/np.log(10)
sigma_natural = np.exp(sigma_pl_ln)
sigma_pl = sigma_natural/np.log(10)

# Defining the functions for mass and orbital separation distibutions for both brown dwarf and planets 
def mass_fctn_bd(q):
	return (q**(beta))#dq
def orbital_dist_bd(a):
	return (A_bd*np.exp((-(np.log10(a)-(mean_bd))**2.)/(2.*sigma_bd**2.)))/(2.0*np.pi*sigma_bd*a)#da
def mass_fctn_pl(q):
	return (q**(-alpha_pl))#dm
def orbital_dist_pl(a):
	 return (A_pl/(2*np.pi*sigma_pl*a))*np.exp((-(np.log10(a)-(mu_pl))**2.)/(2.*sigma_pl**2.))#da

f_bd = (integrate.quad(mass_fctn_bd,(Jup_min*q_Jupiter),(Jup_max*q_Jupiter))[0]*integrate.quad(orbital_dist_bd,a_min,a_max)[0])
f_pl = (integrate.quad(mass_fctn_pl,(Jup_min*q_Jupiter),(Jup_max*q_Jupiter))[0]*integrate.quad(orbital_dist_pl,a_min,a_max)[0])

# Display the entered number
st.write("BD companion frequency:", f_bd)
st.write("Planet companion frequency:", f_pl)


