from scipy.io import loadmat
import os

def load_data(folder):
    data = loadmat(os.path.join(folder, 'Spike_Data_Binned.mat'))
    data = data["Spike_Data_Bin"].squeeze()
    return data
