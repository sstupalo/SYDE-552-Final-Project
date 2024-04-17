import numpy as np
import scipy
from brian2 import *
import matplotlib.pyplot as plt
from scipy import signal

start_scope()

# PARAMETERS
Cm        = 1*uF/cm**2          # Membrane capacitance
Rl        = 100*ohm*cm          # Longitudinal intracellular resistivity
gNa       = 100*mS/meter**2     # Conductance of persistent sodium current
gK        = 40*mS/meter**2      # Conductance of potassium current
gAHP      = 0.01*mS/meter**2    # Conductance of afterhyperpolarization current
gKL       = 0.05*mS/meter**2    # Conductance of potassium leak current
gNaL      = 0.0175*mS/meter**2  # Conductance of sodium leak current
gClL      = 0.05*mS/meter**2    # Conductance of chloride leak current
f         = 3*Hz                # Time constant of gating variables
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
# Na_i_rest = 18.0*mM             # Normal resting intracellular sodium concentration
# K_i_rest  = 140.0*mM            # Normal resting intracellular potassium concentration
# Na_o_rest = 144.0*mM             # Normal resting extracellular sodium concentration

Na_i = 10*mM
K_o = 4*mM
Ca_i = 0.1*uM

# DYNAMIC VARIABLES
# V     = 0 * mV             # Membrane potential (millivolts)
# h     = 0                  # Inactivating sodium gate (dimensionless)
# n     = 0                  # Activating potassium gate (dimensionless)
# Na_i  = 0 * mM             # Intracellular sodium concentration (millimolar)
# K_o   = 0 * mM             # Extracellular potassium concentration (millimolar)
# Ca_i   = 0 * mM            # Intracellular calcium concentration (millimolar)


# HH EQUATIONS
# Im = Cm*dv/dt
# Why is the unit of dCa_i/dt seemingly A/m^2 - mM?
eqn = '''
Im = INa + IK + ICl : amp/meter**2
INa = -gNa*(m_inf**3)*h*(v-VNa)-gNaL(v-VNa) : amp/meter**2
IK = -(gK*(n**4)+((gAHP*Ca_i)/(1+Ca_i)))*(v-VK)-gKL(v-VK) : amp/meter**2
ICl = -gClL(v-VCl) : amp/meter**2

dn/dt = rho*(alpha_n*(1-n) - beta_n*n) : Hz
dh/dt = rho*(alpha_h*(1-h) - beta_h*h) : Hz
dCa_i/dt = -0.002*gCa*(v-VCa)/(1+exp(-(v+25)/2.5))-(Ca_i/80) : amp (point current)


dK_o/dt = -0.33*IK - 2*beta*Ipump - Iglia - Idiff : amp (point current)
Ipump = (rho/(1+exp((25-Na_i)/3)))*(1/(1+exp(5.5-K_o))) : amp/meter**2
Iglia = Gglia/(1+exp((18-K_o)/2.5)) : amp/meter**2
Idiff = epsilon*(K_o-k_o_inf) : amp/meter**2

K_i = 140.0+(18.0-Na_i) : mmolar/meter**3
dNa_i/dt = 0.33*INa/beta-3*Ipump : amp (point current)
Na_o = 144.0-beta*(Na_i-18.0) : mmolar/meter**3

VNa = 26.64*log(Na_o/Na_i) : volt
VK = 26.64*log(K_o/K_i) : volt
VCl = 26.64*log(Cl_i/Cl_o) : volt

m_inf = alpha_m/(alpha_m+beta_m) : 1 (constant over dt)
alpha_m = 0.1*(v+30)/(1-exp(-0.1*(v+30))) : Hz (constant over dt)
beta_m = 4*exp(-(v+55)/18) : Hz (constant over dt)
alpha_n = 0.01*(v+34)/(1-exp(-0.1*(v+34))) : Hz (constant over dt)
beta_n = 0.125*exp(-(v+44)/80) : Hz (constant over dt)
alpha_h = 0.07*exp(-(v+44)/20) : Hz (constant over dt)
beta_h = 1/(1+exp(-0.1*(v+4))) : Hz (constant over dt)
'''

soma = Soma(diameter=30*um)
soma.axon = Cylinder(length=100*um, diameter=1*um, n=10)
soma.dendrite = Cylinder(length=50*um, diameter=2*um, n=5)
soma.dendrite.branch1 = Cylinder(length=10*um, diameter=1*um, n=3)
soma.dendrite.branch2 = Cylinder(length=10*um, diameter=1*um, n=3)

# not sure if we should use cable or soma for the morphology
#cable = Cylinder(length=10*cm, diameter=2*238*um, n=1000, type='axon')

neuron = SpatialNeuron(morphology=soma, model=eqn, method="exponential_euler", Cm=Cm, Ri=Rl)

# variables = ['v', 'h', 'n', 'Na_i', 'K_o', 'Ca_i']

#
# for var in variables:
#     setattr(neuron, var, 0)

# neuron.v = 0*mV
# neuron.h = 0.6
# neuron.n = 0.3
# neuron.Na_i = 100*mM
# neuron.K_o = 20*mM
# neuron.Ca_i = 100*mM

# p = StateMonitor(neuron, ['v', 'h', 'n'], record=True)
# defaultclock.dt = 0.01*ms

# run(10*ms)
# neuron.I[0] = 1*uA # current injection at one end of the axon
# run(3*ms)
# neuron.I = 0*uA
# run(10*ms)

# # plots here
# fig, ax = plt.subplots()
# ax.plot(p.t/ms, p.v[0]/mV)
# ax.set(xlabel='Time (ms)', ylabel='Voltage (mV)')
# title('Q3.1 Voltage vs. Time')
# plt.show()

# fig, ax = plt.subplots()
# ax.plot(p.t/ms, (p.n[0]/max(p.n[0])), label = 'n')
# ax.plot(p.t/ms, (p.h[0]/max(p.h[0])), label = 'h')
# ax.set(xlabel='Time (ms)', ylabel='Voltage (mV)')
# title('Q3.1 Normalized m, n, and h vs. Time')
# ax.legend(loc='lower right')
# plt.show()
