import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

base_filenames = ['18515004', '18515005', '18515008', '18515010']

for base_filename in base_filenames:
    # Create the figure and subplots for each file
    fig, axs = plt.subplots(3, 1, figsize=(15, 9), sharex=True)
    
    # Load the .npy data for channel 0 and the corresponding spikes
    data_ch0 = np.load(f'./numpy_data/{base_filename}.abf_ch0_sweep0.npy')
    spikes = np.load(f'./spikes/{base_filename}.abf_ch0_sweep0.npy_spikes.npy')
    time = np.arange(data_ch0.size) / 10000.0  # Assuming a 10 kHz sampling rate

    # Plotting the first subplot
    axs[0].plot(time, data_ch0, label=f'Ch0 Data')
    axs[0].scatter(time[spikes == 1], data_ch0[spikes == 1], color='red', s=10, label=f'Spikes')
    axs[0].set_title(f'Channel 0 Data with Spikes: {base_filename}')
    axs[0].legend()

    # Load the .npy data for channel 1
    data_ch1 = np.load(f'./numpy_data/{base_filename}.abf_ch1_sweep0.npy')

    # Plotting the second subplot
    axs[1].plot(time, data_ch1, label=f'Ch1 Data')
    axs[1].set_title(f'Channel 1 Data: {base_filename}')
    axs[1].legend()

    # Load and plot the currents from the .csv files
    currents_df = pd.read_csv(f'./results/{base_filename}.abf_ch0_sweep0.npy_spikes_currents.csv')
    axs[2].plot(currents_df['Time (ms)'] / 1000.0, currents_df['Current (pA)'], label=f'Simulated Currents')
    axs[2].set_title(f'Simulated Currents: {base_filename}')
    axs[2].set_xlabel('Time (s)')
    axs[2].legend()

    plt.tight_layout()

    # Save each figure with a unique filename
    plt.savefig(f'./plots/{base_filename}_plot.png', dpi=500)
    plt.close(fig)  # Close the figure to free memory


# Additional code for overlapping plots with 3 subplots each
overlap_pairs = [('18515004', '18515005'), ('18515008', '18515010')]  # Pairs of files for KO and WT
labels = ['KO', 'WT']  # Labels for the conditions

for pair in overlap_pairs:
    fig, axs = plt.subplots(3, 1, figsize=(15, 9), sharex=True)
    
    for i, base_filename in enumerate(pair):
        # Load Channel 0 data and spikes
        data_ch0 = np.load(f'./numpy_data/{base_filename}.abf_ch0_sweep0.npy')
        spikes = np.load(f'./spikes/{base_filename}.abf_ch0_sweep0.npy_spikes.npy')
        time = np.arange(data_ch0.size) / 10000.0  # Assuming 10 kHz sampling rate

        # Plot Channel 0 data with spikes
        axs[0].plot(time, data_ch0, label=f'{labels[i]}: {base_filename} ch0')
        axs[0].scatter(time[spikes == 1], data_ch0[spikes == 1], color='red', s=10, label=f'{labels[i]} Spikes')

        # Load Channel 1 data
        data_ch1 = np.load(f'./numpy_data/{base_filename}.abf_ch1_sweep0.npy')

        # Plot Channel 1 data
        axs[1].plot(time, data_ch1, label=f'{labels[i]}: {base_filename} ch1')

        # Load and plot currents from CSV files
        currents_df = pd.read_csv(f'./results/{base_filename}.abf_ch0_sweep0.npy_spikes_currents.csv')
        axs[2].plot(currents_df['Time (ms)'] / 1000.0, currents_df['Current (pA)'], label=f'{labels[i]} Currents')

    # Setting titles for each subplot
    axs[0].set_title('Channel 0 Data with Spikes')
    axs[1].set_title('Channel 1 Data')
    axs[2].set_title('Simulated Currents')
    axs[2].set_xlabel('Time (s)')  # Only the bottom subplot needs an x-label
    
    # Adding legends to each subplot
    for ax in axs:
        ax.legend()

    plt.tight_layout()
    # Save each figure with a unique filename reflecting the pair of conditions
    plt.savefig(f'./plots/overlap_{pair[0]}_{pair[1]}.png', dpi=500)
    plt.close(fig)  # Close the figure to conserve memory
