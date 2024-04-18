import numpy as np
import scipy
from brian2 import *
import matplotlib.pyplot as plt
from scipy import signal

start_scope()

# PARAMETERS
Cm        = 1*uF/cm**2          # Membrane capacitance
gNa       = 100*mS/meter**2     # Conductance of persistent sodium current
gK        = 40*mS/meter**2      # Conductance of potassium current
gAHP      = 0.01*mS/meter**2    # Conductance of afterhyperpolarization current
gKL       = 0.05*mS/meter**2    # Conductance of potassium leak current
gNaL      = 0.0175*mS/meter**2  # Conductance of sodium leak current
gClL      = 0.05*mS/meter**2    # Conductance of chloride leak current
phi       = 3*Hz                # Time constant of gating variables
VCl       = -81.93*mV           # Reversal potential of chloride current
gCa       = 0.1*mS/meter**2     # Calcium conductance
VCa       = 120*mV              # Reversal potential of calcium
beta      = 7.0                 # Ratio of intracellular to extracellular volume of the cell
rho       = 1.25*mM/second      # Pump strength
Gglia     = 66*mM/second        # Strength of glial uptake
epsilon   = 1.2*Hz              # Diffusion constant
k_o_inf   = 4.0*mM              # Steady state extracellular potassium concentration
Cl_i      = 6.0*mM              # Intracellular chloride concentration
Cl_o      = 130.0*mM            # Extracellular chloride concentration

alpha_K     = 0.94
alpha_Na    = 1.0
A_1         = 0.75
B_1         = 0.92
mu_1        = 2.6
lambda_2     = 7.41
sigma_2     = 1.0
mu_2        = 2.6
sigma_3     = 35.7
mu_3        = 1.94
lambda_3    = 24.3
sigma_4     = 0.88
mu_4        = 1.48
lambda_4    = 24.6
A_1Na       = 1.5
A_1K        = 2.6
lambda_1K   = 32.5

eqn = '''
dv/dt = (INa + IK + ICl)/Cm : volt
INa = -gNa*(m_inf**3)*h*(v-VNa)-gNaL*(v-VNa) : amp/meter**2
IK = -(gK*(n**4)+((gAHP*Ca_i)/(1*mM+Ca_i)))*(v-VK)-gKL*(v-VK) : amp/meter**2
ICl = -gClL*(v-VCl) : amp/meter**2

dn/dt = phi*(alpha_n*(1-n) - beta_n*n) : 1
dh/dt = phi*(alpha_h*(1-h) - beta_h*h) : 1
dCa_i/dt = (-0.002*gCa*(v-VCa)/(1+exp(-(v/mV+25)/2.5)))*mmol/(meter*ucoulomb)-(Ca_i/80)/second : mmolar

dK_o/dt = 0.33*IK*mmol/(meter*ucoulomb) - 2*beta*Ipump - Iglia - Idiff : mmolar
Ipump = (rho/(1+exp((25-Na_i/mM)/3)))*(1/(1+exp(5.5-K_o/mM))) :  mmolar/second
Iglia = Gglia/(1+exp((18-K_o/mM)/2.5)) :  mmolar/second
Idiff = epsilon*(K_o-k_o_inf) :  mmolar/second

K_i = 140.0*mM+(18.0*mM-Na_i) : mmolar  
dNa_i/dt = (0.33*INa*mmol/(meter*ucoulomb))/beta-3*Ipump :  mmolar
Na_o = 144.0*mM-beta*(Na_i-18.0*mM) : mmolar  

VNa = 26.64*log(abs(Na_o/Na_i))*mV : volt
VK = 26.64*log(K_o/K_i)*mV : volt
VCl = 26.64*log(Cl_i/Cl_o)*mV : volt

m_inf = alpha_m/(alpha_m+beta_m) : 1
alpha_m = 1/-exprel(-0.1*(v/mV+30)) : 1
beta_m = 4*exp(-(v/mV+55)/18) : 1
alpha_n = 1/-exprel(-0.1*(v/mV+34)) : 1
beta_n = 0.125*exp(-(v/mV+44)/80) : 1 
alpha_h = 0.07*exp(-(v/mV+44)/20) : 1
beta_h = 1/(1+exp(-0.1*(v/mV+4))) : 1

K_oi = K_o/K_i : 1
Na_io = Na_i/Na_o : 1
g1 = 420.0*(1-A_1*(1-B_1*exp(-mu_1*Na_io))**(1/3)) : 1
g2 = exp(sigma_2*(1-lambda_2*K_oi)/(1.0+exp(-mu_2*Na_io))) : 1
g3 = 1/(1+exp(sigma_3*(1+mu_3*Na_io-lambda_3*K_oi)))**5 : 1
g4 = 1/(1+exp(sigma_4*(1+mu_4*Na_io-lambda_4*K_oi)))**5 : 1
g_1K = A_1K*exp(-lambda_1K*K_oi) : 1
g_1Na = A_1Na : 1
'''

prefs.codegen.target = 'numpy'
prefs.codegen.loop_invariant_optimisations = False
np.seterr(all='raise')

neuron = NeuronGroup(1, eqn, method='euler')

variables = ['v', 'h', 'n', 'Na_i', 'K_o', 'Ca_i']

for var in variables:
    setattr(neuron, var, 0)

neuron.v = -55*mV
neuron.h = 0.6
neuron.n = 0.3
neuron.Na_i = 100*mM
neuron.K_o = 20*mM
neuron.Ca_i = 100*mM

p = StateMonitor(neuron, ['v', 'h', 'n'], record=True)
defaultclock.dt = 0.01*ms

run(10*ms)

fig, ax = plt.subplots()
ax.plot(p.t/ms, p.v[0]/mV)
ax.set(xlabel='Time (ms)', ylabel='Voltage (mV)')
title('Q3.1 Voltage vs. Time')
plt.show()