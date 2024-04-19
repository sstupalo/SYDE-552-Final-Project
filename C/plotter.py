import numpy as np
import matplotlib.pyplot as plt
import scipy

pH_vals = [7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7]
ISI_result = []

fig, axs = plt.subplots(4, 2, figsize=(8, 12))
# fig.suptitle('Voltage vs Time for Different pH Values', fontsize=16)

for i, pH in enumerate(pH_vals):

    filename = f"results/datafile_{pH}.dat"
    data = np.loadtxt(filename)

    # Extract columns
    time = data[:, 0]
    voltage = data[:, 1]
    potassium_concentration = data[:, 2]
    sodium_concentration = data[:, 3]


    spike_indicies = scipy.signal.find_peaks(voltage)[0]
    spike_times = time[spike_indicies]

    ISI = np.diff(spike_times)
    ISI = ISI[6:]  # Removing initial outlier spikes

    mean_ISI = np.mean(ISI)
    ISI_result.append(mean_ISI)

    # # Plot the voltage data
    # plt.figure(figsize=(10, 6))
    # plt.plot(time/1000, voltage, label='Voltage', color='blue')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Voltage')
    # plt.title(f'Voltage vs Time at pH {pH}')
    # plt.legend()
    # plt.grid(True)
    # plt.show()

    # Create a figure and axis object
    # fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot the potassium concentration data
    # ax1.plot(time/1000, potassium_concentration, label='Potassium Concentration', color='green')
    # ax1.set_xlabel('Time (s)')
    # ax1.set_ylabel('Potassium Concentration (mM)', color='green')
    # ax1.tick_params(axis='y', labelcolor='green')
    # ax1.grid(True)

    # # Create a secondary y-axis for sodium concentration
    # ax2 = ax1.twinx()
    # ax2.plot(time/1000, sodium_concentration, label='Sodium Concentration', color='red')
    # ax2.set_ylabel('Sodium Concentration (mM)', color='red')
    # ax2.tick_params(axis='y', labelcolor='red')

    # # Add legends
    # lines, labels = ax1.get_legend_handles_labels()
    # lines2, labels2 = ax2.get_legend_handles_labels()
    # ax2.legend(lines + lines2, labels + labels2)

    # plt.title('Concentration vs Time')
    # plt.show()

    # Plot voltage vs time for each pH value
    row = i % 4 
    col = i // 4
    axs[row, col].plot(time, voltage, label=f'pH {pH}', color='blue')
    axs[row, col].set_xlabel('Time (s)', fontsize=10)
    axs[row, col].set_ylabel('Voltage (mV)', fontsize=10)
    axs[row, col].set_title(f'pH {pH}', fontsize=12)
    axs[row, col].grid(True) 

plt.tight_layout()
plt.show()

ISI_result = np.array(ISI_result)

# Plot the ISI data
plt.figure(figsize=(8, 5))
plt.plot(pH_vals, ISI_result, label='ISI', color='blue', marker='o')
plt.xlabel('pH')
plt.ylabel('ISI (ms)')
plt.title(f'ISI vs pH')
plt.legend()
plt.grid(True)
plt.show()
