import numpy as np
import matplotlib.pyplot as plt

# Load data from datafile.dat
data = np.loadtxt('datafile.dat')

# Extract columns
time = data[:, 0]
voltage = data[:, 1]
potassium_concentration = data[:, 2]
sodium_concentration = data[:, 3]

# Plot the voltage data
plt.figure(figsize=(10, 6))
plt.plot(time, voltage, label='Voltage', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('Voltage')
plt.title('Voltage vs Time')
plt.legend()
plt.grid(True)
plt.show()

# Plot the potassium concentration data
plt.figure(figsize=(10, 6))
plt.plot(time, potassium_concentration, label='Potassium Concentration', color='green')
plt.xlabel('Time (s)')
plt.ylabel('Potassium Concentration')
plt.title('Potassium Concentration vs Time')
plt.legend()
plt.grid(True)
plt.show()

# Plot the sodium concentration data
plt.figure(figsize=(10, 6))
plt.plot(time, sodium_concentration, label='Sodium Concentration', color='red')
plt.xlabel('Time (s)')
plt.ylabel('Sodium Concentration')
plt.title('Sodium Concentration vs Time')
plt.legend()
plt.grid(True)
plt.show()