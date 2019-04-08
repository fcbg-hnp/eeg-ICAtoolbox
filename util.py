import numpy as np
import mne
from mne.channels import Montage


def xyz_to_montage(path) :
    """Convert xyz positions to a mne montage type"""

    n = int(open(path).readline().split(' ')[0])
    coord = np.loadtxt(path, skiprows = 1, usecols = (0,1,2), max_rows = n)
    names = np.loadtxt(path, skiprows = 1, usecols = 3, max_rows = n, dtype = np.dtype(str)).tolist()
    return Montage(coord, names, 'standard_1005', selection = [i for i in range(n)])

def plot_overlay(raw, ica):
    """Custom plot overlay given fitted ica and raw"""
    picks = mne.pick_types(raw.info, meg=True, eeg=True, stim=False, eog=False,
                      misc=False, exclude='bads')
    raw_copy = raw.copy()
    ica.apply(raw_copy)
    mean_raw_signal = np.mean(raw.get_data(picks), axis = 0)
    mean_raw_applied_signal = np.mean(raw_copy.get_data(picks), axis = 0)
    sfreq = raw.info["sfreq"]
    ch_types = ['eeg', 'eeg']
    ch_names = ['before', 'after']
    data = [mean_raw_signal, mean_raw_applied_signal]
    info = mne.create_info(ch_names=ch_names, sfreq=sfreq, ch_types=ch_types)
    raw = mne.io.RawArray(data, info)
    raw.info['bads'] = ["before"]
    raw.plot(scalings="auto",n_channels=2, butterfly=True, block=True, bad_color=(1,0,0))
    return()
