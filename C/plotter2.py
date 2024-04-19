import numpy as np
import matplotlib.pyplot as plt
import scipy

pH_vals = [7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7]
I_ext_vals = np.arange(0, 201, 10)
threshold = 20
results = []

for i, pH in enumerate(pH_vals):
    for j, I_ext in enumerate(I_ext_vals):
        filename = f"results/{pH}/datafile_{I_ext}.dat"
        data = np.loadtxt(filename)
        voltage = data[:, 1]
        voltage = voltage[int(len(voltage) * 0.03):]  # Skip the first 3% of voltage data

        if any(voltage > threshold):
            results.append(I_ext)
            break

# Plotting
plt.plot(pH_vals, results, marker='o', linestyle='-')
plt.xlabel('pH')
plt.ylabel('I_ext at threshold voltage (uA/cm^2)')
plt.title('I_ext at Threshold Voltage vs pH')
plt.grid(True)
plt.show()



pH = float(input("Enter pH value: "))
I_ext = int(float(input("Enter I_ext value: ")))

filename = f"results/{pH}/datafile_{I_ext}.dat"
print(filename)
data = np.loadtxt(filename)

time = data[:, 0]
voltage = data[:, 1]

plt.figure(figsize=(10, 6))
plt.plot(time/1000, voltage, label='Voltage', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('Voltage')
plt.title(f'Voltage vs Time at pH {pH}')
plt.legend()
plt.grid(True)
plt.show()



