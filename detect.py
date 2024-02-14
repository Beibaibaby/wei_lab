import numpy as np
import matplotlib.pyplot as plt
import os

def load_data(file_path):
    """Load the time series data from a numpy file."""
    return np.load(file_path)

def detect_spikes(data, threshold=-20, time_bin_ms=30, sampling_rate=10000):
    """
    Detects spikes in the data where there is a drop greater than the threshold within a specified time bin.
    
    Parameters:
    - data: numpy array of the time series data.
    - threshold: the drop in pA to detect as a spike. Default is -20 pA.
    - time_bin_ms: the time bin in milliseconds to check for drops. Default is 15 ms.
    - sampling_rate: number of samples per second. Default is 10000 (10 points per ms).
    
    Returns:
    - A numpy array with 1s at the times of detected spikes and 0s elsewhere.
    """
    samples_per_bin = int(time_bin_ms * (sampling_rate / 1000))
    spikes = np.zeros(data.shape)
    
    i = 0
    while i < len(data) - samples_per_bin:
        window = data[i:i+samples_per_bin]
        if window.min() - window[0] < threshold:
            spike_pos = np.argmin(window) + i
            spikes[spike_pos] = 1
            # Move i to the end of the current window to skip the rest of the drop
            i = spike_pos + 1  # Start looking for the next spike after the current spike position
            
            # Optional: Skip ahead further to ensure we're past the drop
            while i < len(data) - 1 and data[i] < data[i + 1]:
                i += 1
        else:
            i += 1
    
    return spikes


def plot_data_with_spikes(data, spikes, sampling_rate=10000):
    """
    Plot the original time series data along with the detected spikes.
    
    Parameters:
    - data: numpy array of the time series data.
    - spikes: numpy array with 1s at spike times and 0s elsewhere.
    - sampling_rate: number of samples per second. Default is 10000.
    """
    time = np.arange(data.size) / sampling_rate
    plt.figure(figsize=(35, 5))
    plt.plot(time, data, label='Current')
    plt.scatter(time[spikes==1], data[spikes==1], color='red', label='Spikes', zorder=2)
    plt.xlabel('Time (s)')
    plt.ylabel('Current (pA)')
    plt.legend()
    
    output_dir = "./plot"
    os.makedirs(output_dir, exist_ok=True)
    
    # Construct the output file path
    base_name = os.path.basename(file_path)  # Extracts the file name from the path
    output_file = os.path.join(output_dir, base_name.replace('.npy', '_spikes.png'))

    plt.savefig(output_file,dpi=500)
    plt.close()  # Close the plot to free up memory
    print(f"Plot saved as {output_file}")
    
    
    output_dir = "./spikes"
    # Construct output file path
    output_file = os.path.join(output_dir, f"{os.path.basename(file_path)}_spikes.npy")
          
    # Save the data to a NumPy file
    np.save(output_file, spikes)
    print(f"spikes data saved to {output_file}")


# File paths (replace these with the actual paths)
file_paths = ["./numpy_data/18515004.abf_ch0_sweep0.npy", # Example file path
              "./numpy_data/18515005.abf_ch0_sweep0.npy",
              "./numpy_data/18515008.abf_ch0_sweep0.npy",
              "./numpy_data/18515010.abf_ch0_sweep0.npy"
             ]

# Process each file
for file_path in file_paths:
    data = load_data(file_path)
    spikes = detect_spikes(data)
    plot_data_with_spikes(data, spikes)
