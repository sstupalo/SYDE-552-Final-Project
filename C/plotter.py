import numpy as np
import matplotlib.pyplot as plt

# Load data from datafile.dat
data = np.loadtxt('datafile.dat')

# Extract columns
time = data[:, 0]
voltage = data[:, 1]
potassium_concentration = data[:, 2]
sodium_concentration = data[:, 3]

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(time, voltage, label='Voltage')
plt.plot(time, potassium_concentration, label='Potassium Concentration')
plt.plot(time, sodium_concentration, label='Sodium Concentration')
plt.xlabel('Time (s)')
plt.ylabel('Values')
plt.title('Neuron Data')
plt.legend()
plt.grid(True)
plt.show()