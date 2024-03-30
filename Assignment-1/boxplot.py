import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Define parameters
N = [262144, 4194304]
stencils = [5, 9]
executions = range(1, 4)

# Read timing data from text file
with open("Datafile.txt", "r") as file:
    timing_data = file.readlines()

# Convert timing data to float and reshape into a 3D array
timing_array = np.array(timing_data).astype(float).reshape(len(executions), len(N), len(stencils))

# Create DataFrame
data = pd.DataFrame(columns=['Time', 'N', 'Stencil', 'Execution'])

for execution, exec_timing in zip(executions, timing_array):
    for i, n in enumerate(N):
        for j, stencil in enumerate(stencils):
            time = exec_timing[i][j]  # Get the timing value directly
            data = data.append({'Time': time, 'N': n, 'Stencil': stencil, 'Execution': execution}, ignore_index=True)

# Plotting
plt.figure(figsize=(10, 6))
sns.boxplot(data=data, x='N', y='Time', hue='Stencil')
plt.xlabel('(N, Stencil)')
plt.ylabel('Time (seconds)')
plt.title('Timing Data for Different Data Sizes and Stencil Configurations')
plt.legend(title='Stencil')

# Save the plot as an image file
plt.savefig('boxplot.png')

# Show the plot
plt.show()
