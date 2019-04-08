import mne
import matplotlib.pyplot as plt
from mne.datasets import sample
from mne.preprocessing import ICA
import numpy as np

def compute_gfp(raw):
    return(np.mean(raw.get_data()**2,axis=0)**0.5)
    
def plot_overlay(raw, ica):
    picks = mne.pick_types(raw.info, meg=True, eeg=True, stim=False, eog=False,
                      misc=False, exclude='bads')
    raw_copy = raw.copy()
    ica.apply(raw_copy)
    gfp_raw_signal = compute_gfp(raw)
    gfp_raw_applied_signal = compute_gfp(raw_copy)
    sfreq = raw.info["sfreq"]
    ch_types = ['eeg', 'eeg']
    ch_names = ['before', 'after']
    data = [gfp_raw_signal, gfp_raw_applied_signal]
    info = mne.create_info(ch_names=ch_names, sfreq=sfreq, ch_types=ch_types)
    raw = mne.io.RawArray(data, info)
    raw.info['bads'] = ["before"]
    raw.plot(scalings="auto",n_channels=2, butterfly=True, block=True, bad_color=(1,0,0))
    return()

# getting some data ready
data_path = sample.data_path()
raw_fname = r"C:\Users\vferat\Desktop\brainHackTestOpenfilter-raw.fif"
raw = mne.io.read_raw_fif(raw_fname, preload=True)
#raw.filter(1., None, n_jobs=1, fir_design='firwin')
raw.plot(scalings="auto", butterfly=True, block=True)
method = 'fastica'
ica = ICA(n_components=12, method='fastica')
ica.fit(raw, decim=1)
ica.exclude.extend([0])
plot_overlay(raw,ica)
