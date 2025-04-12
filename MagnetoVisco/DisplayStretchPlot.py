import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------- #
# Plot the stretch against time for different models
# ---------------------------------- #


# ----------------- #
# SIMPLE CUBE
# ----------------- #
FP_files=['R-model-SC-VF-010-VIM_FP_graph.txt', 'R-model-SC-VF-020-VIM_FP_graph.txt', 'R-model-SC-VF-030-VIM_FP_graph.txt']
colors = ['r', 'g', 'b']
# Load the data from the saved file
plt.figure(figsize=(8, 6))
for i, file in enumerate(FP_files): 
    data = np.loadtxt(f'./results/Magnetostriction/B_1mTs/FP_Plot/{file}')
    ttt = data[:, 0]   # First column: time
    F33 = data[:, 1]   # Second column: F33 (engineering strain component)
    P33 = data[:, 2]   # Third column: P33 (engineering stress)
    # Plot P33 vs F33 in reverse order
    plt.plot(ttt/15, F33, c=colors[i], label=file) #/15 to normalize time

plt.xlabel('Normalized time')
plt.xlim(0,3)
plt.ylabel('Stretch (F33)')
plt.title('Stretch against normalized time for b_rate = 1 mT/s for Simple cube')
plt.grid(True)
plt.legend()


FP_files=['R-model-RI-VF-07-N_2-R_20-VIM_FP_graph.txt', 'R-model-RI-VF-017-N_5-R_20-VIM_FP_graph.txt', 'R-model-RI-VF-027-N_8-R_20-VIM_FP_graph.txt',
          'R-model-AR_05-RI-VF-07-N_2-R_20-VIM_FP_graph.txt', 'R-model-AR_05-RI-VF-017-N_5-R_20-VIM_FP_graph.txt', 'R-model-AR_05-RI-VF-027-N_8-R_20-VIM_FP_graph.txt']
colors = ['r', 'g', 'b']
# ----------------- #
# RVI inclusions
# b = 1 mT/s
# Aspect ratio 1.5
# ----------------- #
# Load the data from the saved file
plt.figure(figsize=(8, 6))
for i, file in enumerate(FP_files[:3]): 
    data = np.loadtxt(f'./results/Magnetostriction/B_1mTs/FP_Plot/{file}')
    ttt = data[:, 0]   # First column: time
    F33 = data[:, 1]   # Second column: F33 (engineering strain component)
    P33 = data[:, 2]   # Third column: P33 (engineering stress)
    # Plot P33 vs F33 in reverse order
    plt.plot(ttt/15, F33, c=colors[i], label=file) #/15 to normalize time

plt.xlabel('Normalized time')
plt.xlim(0,5)
plt.ylabel('Stretch (F33)')
plt.title('Stretch against normalized time for b_rate = 1 mT/s for random arangement. Aspect ratio 1.5')
plt.grid(True)
plt.legend()


# ----------------- #
# RVI inclusions
# b = 1 mT/s
# Aspect ratio 0.5
# ----------------- #
# Load the data from the saved file
plt.figure(figsize=(8, 6))
for i, file in enumerate(FP_files[3:]): 
    data = np.loadtxt(f'./results/Magnetostriction/B_1mTs/FP_Plot/{file}')
    ttt = data[:, 0]   # First column: time
    F33 = data[:, 1]   # Second column: F33 (engineering strain component)
    P33 = data[:, 2]   # Third column: P33 (engineering stress)
    # Plot P33 vs F33 in reverse order
    plt.plot(ttt/15, F33, c=colors[i], label=file) #/15 to normalize time

plt.xlabel('Normalized time')
plt.xlim(0,5)
plt.ylabel('Stretch (F33)')
plt.title('Stretch against normalized time for b_rate = 1 mT/s for random arangement')
plt.grid(True)
plt.legend()

# ----------------- #
# RVI inclusions
# b = 10 mT/s
# ----------------- #
FP_files=['R-model-RI-VF-07-N_2-R_20-VIM_FP_graph.txt', 'R-model-RI-VF-017-N_5-R_20-VIM_FP_graph.txt', 'R-model-RI-VF-027-N_8-R_20-VIM_FP_graph.txt',
          'R-model-RI-VF-07-N_2-R_20-VIM_FP_graph.txt', 'R-model-RI-VF-017-N_5-R_20-VIM_FP_graph.txt', 'R-model-RI-VF-027-N_8-R_20-VIM_FP_graph.txt']
plt.figure(figsize=(8, 6))
for i, file in enumerate(FP_files[3:]): 
    data = np.loadtxt(f'./results/Magnetostriction/B_10mTs/FP_Plot/{file}')
    ttt = data[:, 0]   # First column: time
    F33 = data[:, 1]   # Second column: F33 (engineering strain component)
    P33 = data[:, 2]   # Third column: P33 (engineering stress)
    # Plot P33 vs F33 in reverse order
    plt.plot(ttt/1.5, F33, c=colors[i], label=file) #/15 to normalize time


plt.xlabel('Normalized time')
plt.xlim(0,5)
plt.ylabel('Stretch (F33)')
plt.title('Stretch against normalized time for b_rate = 10 mT/s for random arangement')
plt.grid(True)
plt.legend()


# ----------------- #
# RVI inclusions
# b = 100 mT/s
# ----------------- #
plt.figure(figsize=(8, 6))
for i, file in enumerate(FP_files[3:]): 
    data = np.loadtxt(f'./results/Magnetostriction/B_100mTs/FP_Plot/{file}')
    ttt = data[:, 0]   # First column: time
    F33 = data[:, 1]   # Second column: F33 (engineering strain component)
    P33 = data[:, 2]   # Third column: P33 (engineering stress)
    # Plot P33 vs F33 in reverse order
    plt.plot(ttt/0.15, F33, c=colors[i], label=file) #/15 to normalize time


plt.xlabel('Normalized time')
plt.xlim(0,5)
plt.ylabel('Stretch (F33)')
plt.title('Stretch against normalized time for b_rate = 100 mT/s for random arangement')
plt.grid(True)
plt.legend()

# ----------------- #
# RVI inclusions
# b = 1000 mT/s
# ----------------- #
plt.figure(figsize=(8, 6))
for i, file in enumerate(FP_files[3:]): 
    data = np.loadtxt(f'./results/Magnetostriction/B_1000mTs/FP_Plot/{file}')
    ttt = data[:, 0]   # First column: time
    F33 = data[:, 1]   # Second column: F33 (engineering strain component)
    P33 = data[:, 2]   # Third column: P33 (engineering stress)
    # Plot P33 vs F33 in reverse order
    plt.plot(ttt/0.015, F33, c=colors[i], label=file) #/15 to normalize time


plt.xlabel('Normalized time')
plt.xlim(0,5)
plt.ylabel('Stretch (F33)')
plt.title('Stretch against normalized time for b_rate = 1000 mT/s for random arangement')
plt.grid(True)
plt.legend()

plt.show()