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
plt.plot(time/1000, voltage, label='Voltage', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('Voltage')
plt.title('Voltage vs Time')
plt.legend()
plt.grid(True)
plt.show()

# Create a figure and axis object
fig, ax1 = plt.subplots(figsize=(10, 6))

# Plot the potassium concentration data
ax1.plot(time/1000, potassium_concentration, label='Potassium Concentration', color='green')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Potassium Concentration (mM)', color='green')
ax1.tick_params(axis='y', labelcolor='green')
ax1.grid(True)

# Create a secondary y-axis for sodium concentration
ax2 = ax1.twinx()
ax2.plot(time/1000, sodium_concentration, label='Sodium Concentration', color='red')
ax2.set_ylabel('Sodium Concentration (mM)', color='red')
ax2.tick_params(axis='y', labelcolor='red')

# Add legends
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2)

plt.title('Concentration vs Time')
plt.show()
