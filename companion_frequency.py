## companion_frequency.py ##
##############################



import streamlit as st
import pandas as pd
import numpy as np
# ~ from bokeh.plotting import figure
# ~ from bokeh.io import curdoc
# ~ from bokeh.palettes import Magma256, Viridis256

# 
from scipy import integrate
from scipy import interpolate
import tkinter as tk



from streamlit import session_state

st.set_page_config(
	page_title="SLSdb",
	page_icon="",
	layout="wide",
	)
    
    


## From GUI script

def calculate_companion_frequency():
	host_mass = float(e1.get()) # Msun   
	Jup_min = float(e2.get())      
	Jup_max = float(e3.get())      
	a_min = float(e4.get())   
	a_max = float(e5.get()) 
	
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
	result.set([f_bd, f_pl])





## Initial setup
master = tk.Tk()

tk.Label(master, text="Host mass [Solar]").grid(row=0)
tk.Label(master, text="Min companion mass [Jupiter]").grid(row=1)
tk.Label(master, text="Max companion mass [Jupiter]").grid(row=2)
tk.Label(master, text="a min [AU]").grid(row=3)
tk.Label(master, text="a max [AU]").grid(row=4)


tk.Label(master, text="Host mass [Solar]").grid(row=0)
tk.Label(master, text="Min companion mass [Jupiter]").grid(row=1)
tk.Label(master, text="Max companion mass [Jupiter]").grid(row=2)
tk.Label(master, text="a min [AU]").grid(row=3)
tk.Label(master, text="a max [AU]").grid(row=4)

e1 = tk.Entry(master)
e2 = tk.Entry(master)
e3 = tk.Entry(master)
e4 = tk.Entry(master)
e5 = tk.Entry(master)
b1 = tk.Button(master, text="Calculate", command=calculate_companion_frequency)
tk.Label(master, text="f_bd").grid(row=8, column=0)
tk.Label(master, text="f_pl").grid(row=8, column=1)

e1.grid(row=0, column=1)
e2.grid(row=1, column=1)
e3.grid(row=2, column=1)
e4.grid(row=3, column=1)
e5.grid(row=4, column=1)
b1.grid(row=6, column=0, columnspan=2)

result = tk.StringVar()
result_label = tk.Label(master, textvariable=result)
result_label.grid(row=9, column=0, columnspan=2)

# Run the main loop
master.mainloop()


