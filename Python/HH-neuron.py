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


# HH EQUATIONS
eqn = '''
dv/dt = (INa + IK + ICl)/Cm : volt
INa = -gNa*(m_inf**3)*h*(v-VNa)-gNaL*(v-VNa) : amp/meter**2
IK = -(gK*(n**4)+((gAHP*Ca_i)/(1*mM+Ca_i)))*(v-VK)-gKL*(v-VK) : amp/meter**2
ICl = -gClL*(v-VCl) : amp/meter**2

dn/dt = phi*(alpha_n*(1-n) - beta_n*n) : 1
dh/dt = phi*(alpha_h*(1-h) - beta_h*h) : 1
dCa_i/dt = (-0.002*gCa*(v-VCa)/(1+exp(-(v/mV+25)/2.5)))*mmol/(meter*ucoulomb)-(Ca_i/80)/second : mmolar

dK_o/dt = 0.33*((mM*cm**2)/ucoulomb)*IK - 2*beta*Ipump - Iglia - Idiff : mmolar
Ipump = (rho/(1+exp((25-Na_i/mM)/3)))*(1/(1+exp(5.5-K_o/mM))) :  mmolar/second
Iglia = Gglia/(1+exp((18-K_o/mM)/2.5)) :  mmolar/second
Idiff = epsilon*(K_o-k_o_inf) :  mmolar/second

K_i = 140.0*mM+(18.0*mM-Na_i) : mmolar  
dNa_i/dt = (0.33*INa*mmol/(meter*ucoulomb))/beta-3*Ipump :  mmolar
Na_o = 144.0*mM-beta*(Na_i-18.0*mM) : mmolar  

VNa = 26.64*log(abs(Na_o/Na_i))*mV : volt
VK = 26.64*log(abs(K_o/K_i))*mV : volt
VCl = 26.64*log(Cl_i/Cl_o)*mV : volt

m_inf = alpha_m/(alpha_m+beta_m) : 1
alpha_m = 1/-exprel(-0.1*(v/mV+30)) : 1
beta_m = 4*exp(-(v/mV+55)/18) : 1
alpha_n = 1/-exprel(-0.1*(v/mV+34)) : 1
beta_n = 0.125*exp(-(v/mV+44)/80) : 1 
alpha_h = 0.07*exp(-(v/mV+44)/20) : 1
beta_h = 1/(1+exp(-0.1*(v/mV+4))) : 1
'''

prefs.codegen.target = 'numpy'
prefs.codegen.loop_invariant_optimisations = False
np.seterr(all='raise')


neuron = NeuronGroup(1, eqn, method='rk4')


variables = ['v', 'h', 'n', 'Na_i', 'K_o', 'Ca_i']

for var in variables:
    setattr(neuron, var, 0)

neuron.v = -50*mV
neuron.h = 0.08553
neuron.n = 0.96859
neuron.Na_i = 15.5*mM
neuron.K_o = 7.8*mM
neuron.Ca_i = 0.0*mM

p = StateMonitor(neuron, ['v', 'h', 'n', 'Na_i', 'K_o', 'Ca_i'], record=True)
defaultclock.dt = 0.01*ms

run(20*second)


fig, ax = plt.subplots()
ax.plot(p.t, p.v[0]/mV)
ax.set(xlabel='Time (s)', ylabel='Voltage (mV)')
title('')
plt.show()

fig, ax = plt.subplots()
ax.plot(p.t, p.Na_i[0]/mM, label = '[Na]_i')
ax.plot(p.t, p.K_o[0]/mM, label = '[K]_o')
# ax.plot(p.t, p.Ca_i[0], label = '[Ca]_i')
ax.set(xlabel='Time (s)', ylabel='Concentration (mM)')
title('')
ax.legend()
plt.show()
