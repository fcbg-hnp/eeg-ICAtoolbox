import numpy as np
import scipy
import pandas as pd
import mne
import seaborn as sns
import matplotlib.pyplot as plt
from mne.channels import Montage
from mne.viz.topomap import _plot_topomap
from matplotlib.gridspec import GridSpec


def compute_head_pos(montage):
    pos = montage.get_pos2d()
    scale = 0.85 / (pos.max(axis=0) - pos.min(axis=0))
    center = 0.5 * (pos.max(axis=0) + pos.min(axis=0))
    head_pos = {'scale': scale, 'center': center}
    return(head_pos)


def xyz_to_montage(path):
    """Convert xyz positions to a mne montage type"""
    n = int(open(path).readline().split(' ')[0])
    coord = np.loadtxt(path, skiprows=1, usecols=(0, 1, 2), max_rows=n)
    names = np.loadtxt(path, skiprows=1, usecols=3,
                       max_rows=n, dtype=np.dtype(str)).tolist()
    montage = Montage(coord, names, 'standard_1005',
                      selection=[i for i in range(n)])
    return(montage)


def compute_gfp(raw):
    return(np.mean(raw.get_data()**2, axis=0)**0.5)


def plot_overlay(raw, ica):
    print(type(raw))
    """Custom plot overlay given fitted ica and raw"""
    if type(raw) == 'mne.io.fiff.raw.Raw':
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
        raw.plot(scalings="auto", n_channels=2, butterfly=True,
                 block=True, bad_color=(1, 0, 0))
    else:
        pass
    return()


def find_common_channels(ica_a, ica_b):
    """Find ch_names shared between 2 ica objects"""
    ch_names_a = [ch.lower() for ch in ica_a.ch_names]
    ch_names_b = [ch.lower() for ch in ica_b.ch_names]
    common = [x for x in ch_names_a if x in ch_names_b]
    return(common)


def find_index_by_name(names, ica):
    ch_names = [ch.lower() for ch in ica.ch_names]
    Index = []
    for name in names:
        idx = ch_names.index(name)
        Index.append(idx)
    return(Index)


def extract_common_components(ica_a, ica_b):
    common = find_common_channels(ica_a, ica_b)
    idx_a = find_index_by_name(common, ica_a)
    idx_b = find_index_by_name(common, ica_b)
    components_a = ica_a.get_components()[idx_a]
    components_b = ica_b.get_components()[idx_b]
    return(components_a.T, components_b.T)


def construct__template_from_montage(Index, template_values):
    """Extract a template corresponding to the given indexes"""
    match_template = template_values[Index]
    match_template = match_template / np.linalg.norm(match_template)
    return(match_template)


def compute_correlation(match_templates, ics):
    Templates_names = []
    Ics_names = []
    Correlation = []
    for t, temp in enumerate(match_templates):
        for i, ic in enumerate(ics):
            Templates_names.append("template " + str(t).zfill(3))
            Ics_names.append("ic " + str(i).zfill(3))
            pearson = scipy.stats.pearsonr(temp, ic)[0]
            pearson = np.abs(pearson)
            Correlation.append(pearson)
    data = {'Templates_names': Templates_names,
            'Ics_names': Ics_names,
            'correlation': Correlation}
    df = pd.DataFrame(data)
    return(df)


def plot_correlation(df, match_templates, pos, quality, head_pos=None):
    dfp = df.pivot("Templates_names", "Ics_names", "correlation")
    fig = plt.figure()
    gs = GridSpec(11, len(match_templates))
    axes = []
    for t, temp in enumerate(match_templates):
        axes.append(fig.add_subplot(gs[0:5, t]))
        axes[-1].set_title("template " + str(t).zfill(3))
        _plot_topomap(temp, pos, axes=axes[-1], head_pos=head_pos, show=False,
                      vmin=temp.min(), vmax=temp.max(), outlines="head")
    ax_colorbar = fig.add_subplot(gs[5, :])
    ax_matrix = fig.add_subplot(gs[7:11, :])
    sns.heatmap(dfp, linewidths=0.1, annot=False, ax=ax_matrix, cmap="YlGnBu",
                vmin=0, vmax=1, square=False, cbar_ax=ax_colorbar,
                cbar_kws={"orientation" :'horizontal'})
    ax_matrix.set_ylabel('')
    ax_matrix.set_xlabel("Comparaison Quality: " + str(quality*100) + "%")
    plt.subplots_adjust(left=0.17, bottom=0.15,
                        right=None, top=None,
                        wspace=None, hspace=None)
    plt.show(fig)
    return()
