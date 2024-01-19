import pyabf
import matplotlib.pyplot as plt
import numpy as np
import os

def plot_and_save_abf(file_path, output_dir):
    abf = pyabf.ABF(file_path)

    # Determine the number of channels
    n_channels = abf.channelCount

    plt.figure(figsize=(35, 6 * n_channels))

    for channel in range(n_channels):
        for sweep in abf.sweepList:
            abf.setSweep(sweep, channel=channel)
            plt.subplot(n_channels, 1, channel+1)
            plt.plot(abf.sweepX, abf.sweepY, label=f"Channel {channel} Sweep {sweep}")
            plt.xlabel("Time (s)")
            plt.ylabel(abf.sweepLabelY)
            plt.title(f"ABF File: {os.path.basename(file_path)}, Channel {channel}")
            plt.legend()

    # Save the plot
    output_file = os.path.join(output_dir, os.path.basename(file_path) + '.png')
    plt.savefig(output_file, dpi=500)
    plt.close()
    print(f"Plot saved as {output_file}")

def save_abf_to_numpy(file_path, output_dir):
    abf = pyabf.ABF(file_path)

    # Determine the number of channels
    n_channels = abf.channelCount

    for channel in range(n_channels):
        for sweep in abf.sweepList:
            abf.setSweep(sweep, channel=channel)
            data = abf.sweepY  # This gets the sweep data

            # Construct output file path
            output_file = os.path.join(output_dir, f"{os.path.basename(file_path)}_ch{channel}_sweep{sweep}.npy")

            # Save the data to a NumPy file
            np.save(output_file, data)
            print(f"Channel {channel} Sweep {sweep} data saved to {output_file}")

# Define the directory where you want to save the plots and numpy arrays
output_directory_plots = "./plots"
output_directory_numpy = "./numpy_data"

# Ensure the output directories exist
os.makedirs(output_directory_plots, exist_ok=True)
os.makedirs(output_directory_numpy, exist_ok=True)

# Replace these with your actual file paths
file_paths = ["/gpfs/data/doiron-lab/draco/wei_lab/data/ko/18515004.abf", "/gpfs/data/doiron-lab/draco/wei_lab/data/ko/18515008.abf",
              "/gpfs/data/doiron-lab/draco/wei_lab/data/wt/18515005.abf", "./data/wt/18515010.abf"]

for path in file_paths:
    plot_and_save_abf(path, output_directory_plots)
    save_abf_to_numpy(path, output_directory_numpy)
