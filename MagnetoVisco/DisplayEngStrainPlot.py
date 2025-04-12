import numpy as np
import matplotlib.pyplot as plt



FP_files = ['simple_compression_BCC_v0.txt', 
            'simple_compression_BCC_v1.txt', 
            'simple_compression_BCC_v2.txt',
            'simple_compression_BCC_v3.txt',
            'simple_compression_BCC_v4.txt']
colors = ['r', 'g', 'b', 'c', 'm']

# Load the data from the saved file
plt.figure(figsize=(8, 6))
for i, file in enumerate(FP_files): 
    data = np.loadtxt(f'./results/homo/{file}')
    ttt = data[:, 0]   # First column: time
    F33 = data[:, 1]   # Second column: F33 (engineering strain component)
    P33 = data[:, 2]   # Third column: P33 (engineering stress)
    # Plot P33 vs F33 in reverse order
    plt.plot(F33, P33, colors[i], label=file) 


plt.xlabel('Engineering Strain Component (F33)')
plt.ylabel('Engineering Stress (P33)')
plt.title('Engineering Stress vs. Engineering Strain for Centered Cube. Max F: 0.1, Strain rate 0.1/s')
plt.grid(True)
plt.legend()
plt.gca().invert_xaxis()  # Invert x-axis (F33)
plt.gca().invert_yaxis()  # Invert y-axis (P33)
plt.show()



FP_files = ['R-model-SC-VF-010-VI_FP_graph.txt', 
            'R-model-SC-VF-020-VI_FP_graph.txt', 
            'R-model-SC-VF-030-VI_FP_graph.txt']
colors = ['r', 'g', 'b']

# Load the data from the saved file
plt.figure(figsize=(8, 6))
for i, file in enumerate(FP_files): 
    data = np.loadtxt(f'./results/Compression/FP_Plot/{file}')
    ttt = data[:, 0]   # First column: time
    F33 = data[:, 1]   # Second column: F33 (engineering strain component)
    P33 = data[:, 2]   # Third column: P33 (engineering stress)
    # Plot P33 vs F33 in reverse order
    plt.plot(F33, P33, colors[i], label=file) 


plt.xlabel('Engineering Strain Component (F33)')
plt.ylabel('Engineering Stress (P33)')
plt.title('Engineering Stress vs. Engineering Strain for SC. Strain rate 0.3/s')
plt.grid(True)
plt.legend()
plt.gca().invert_xaxis()  # Invert x-axis (F33)
plt.gca().invert_yaxis()  # Invert y-axis (P33)


FP_files = ['R-model-AR_15-RI-VF-07-N_2-R_20-VI_FP_graph.txt', 
            'R-model-AR_15-RI-VF-017-N_5-R_20-VI_FP_graph.txt', 
            'R-model-AR_15-RI-VF-027-N_8-R_20-VI_FP_graph.txt']
# Load the data from the saved file
plt.figure(figsize=(8, 6))
for i, file in enumerate(FP_files): 
    data = np.loadtxt(f'./results/Compression/FP_Plot/{file}')
    ttt = data[:, 0]   # First column: time
    F33 = data[:, 1]   # Second column: F33 (engineering strain component)
    P33 = data[:, 2]   # Third column: P33 (engineering stress)
    # Plot P33 vs F33 in reverse order
    plt.plot(F33, P33, colors[i], label=file) 


plt.xlabel('Engineering Strain Component (F33)')
plt.ylabel('Engineering Stress (P33)')
plt.title('Engineering Stress vs. Engineering Strain for RI with aspect ratio 1.5. Strain rate 0.3/s')
plt.grid(True)
plt.legend()
# Invert both axes
plt.gca().invert_xaxis()  # Invert x-axis (F33)
plt.gca().invert_yaxis()  # Invert y-axis (P33)


FP_files = ['R-model-AR_05-RI-VF-07-N_2-R_20-VI_FP_graph.txt', 
            'R-model-AR_05-RI-VF-017-N_5-R_20-VI_FP_graph.txt', 
            'R-model-AR_05-RI-VF-027-N_8-R_20-VI_FP_graph.txt']
# Load the data from the saved file
plt.figure(figsize=(8, 6))
for i, file in enumerate(FP_files): 
    data = np.loadtxt(f'./results/Compression/FP_Plot/{file}')
    ttt = data[:, 0]   # First column: time
    F33 = data[:, 1]   # Second column: F33 (engineering strain component)
    P33 = data[:, 2]   # Third column: P33 (engineering stress)
    # Plot P33 vs F33 in reverse order
    plt.plot(F33, P33, colors[i], label=file) 

plt.xlabel('Engineering Strain Component (F33)')
plt.ylabel('Engineering Stress (P33)')
plt.title('Engineering Stress vs. Engineering Strain for RI with aspect ratio 0.5. Strain rate 0.3/s')
plt.grid(True)
plt.legend()
# Invert both axes
plt.gca().invert_xaxis()  # Invert x-axis (F33)
plt.gca().invert_yaxis()  # Invert y-axis (P33)

plt.show()