import os
from tqdm import tqdm
import uproot
import numpy as np
import pandas as pd


def load_data(path, shield_surface, files_to_open):
    """
    Load data from files into unique array

    path: absolute path to the files location
    shield_surface: either "panel" or "lid"
    files_to_open: either string "all" to load data from all the files or int with number of files
    """

    files = sorted([file for file in os.listdir(path) if "appliedmap" in file])

    if files_to_open != "all":
        files = files[:files_to_open]

    df = pd.DataFrame()

    for file in tqdm(files):
        try:
            tmp_df = pd.DataFrame()

            tf = uproot.open(os.path.join(path, file))[f"{shield_surface}hits"]

            for var in tf.keys():
                tmp_df[var]  = tf[var].array(library="np")
            
        except:
            print(f"skipping {file}")
        
        df = pd.concat([tmp_df, df])
        
    evt_list = sorted(df.EventID.unique())
    total_nof_events = len(evt_list)

    return df


def get_cryostat_volume():
    
    cryostat_radius = 3.5 # m
    cryostat_height = 7 # m
    cryostat_volume = np.pi * cryostat_radius * cryostat_radius * cryostat_height # m3    
    return cryostat_volume


def get_RT_volume():

    RT_radius = 0.95 # m
    RT_height = 4    # m, 4m inside the cryostat
    RT_volume = np.pi * RT_radius * RT_radius * RT_height
    return RT_volume


def get_moderator_volume():
    n_panels = 12

    panel_height    = 3   # m
    panel_width     = 1   # m
    panel_thickness = 0.1 # m
    mod_radius      = 2   # m

    RT_radius = 0.95 # m

    panel_volume = panel_height * panel_width * panel_thickness # m3
    bottom_cap   = np.pi * (mod_radius + panel_thickness) * (mod_radius + panel_thickness) * panel_thickness # m3
    top_cap      = bottom_cap - np.pi * RT_radius * RT_radius * panel_thickness  # m3

    mod_volume = panel_volume * n_panels + bottom_cap + top_cap # m3
    return mod_volume


def get_AAr_volume():

    cryostat_volume = get_cryostat_volume()
    RT_volume       = get_RT_volume()
    mod_volume      = get_moderator_volume()

    AAr_volume = cryostat_volume - RT_volume - mod_volume
    return AAr_volume


# assuming a volume of AAr and given the activity of Ar39 in natural argon from literature
# compute rate of Ar39 in Hz
def get_Ar39_rate():

    lAr_density = 1396 # kg/m3
    AAr_volume = get_AAr_volume()

    AAr_mass = lAr_density * AAr_volume # kg

    ar39_activity = 1.01 # Bq / kg
    ar39_rate = ar39_activity * AAr_mass # Hz

    return ar39_rate


# compute how many events of Ar39 you expect in a given time window
def get_expected_ar39_in_window(window = 1e-5):

    ar39_rate = get_Ar39_rate()

    return ar39_rate * window


# muon rate is specific of MUSUN sim
# it won't be correct for other MUSUN sims
def get_muon_rate():

    muon_rate = 263 # muons / hour

    return muon_rate / 3600 # Hz


# input
# H: height of the PMMA panel (where light guides will be attached)
# bar_width: light guide width
# n_bar = number of light guides bar to attach
# remember that bar_widt and n_bar should be such that bar_width * n_bar = H
# output 
# light guides position in bins, possible empyt spaces
def get_guides(H, bar_width, n_bar):
    
    if bar_width * n_bar > H:
        raise ValueError('A very specific bad thing happened.')
    # total free space between light guides
    S = H - n_bar*bar_width
    
    # single space between to adjacent light guides
    s = S // (n_bar - 1)

    # residual space (if != 0 you basically split it in 2 and add it at the top and bottom to make to guides fit the panel and be evenly spaced)
    residual_space = H - n_bar*bar_width - s*(n_bar - 1)

    slices = []
    guide_position = []
    
    # if residual space is not zero, skip corresponding zbands
    if residual_space:
        slices.append(0)

        slices.append(residual_space)
        guide_position.append(0)
    
    # start always with a light guide
    slices.append(residual_space+bar_width*2)
    guide_position.append(1)

    # increase by residual_space (if any, otherwise 0) + bar_width * 2 (because zband width = 0.5cm)
    counter = residual_space+bar_width*2
    
    # loop over the entire height H of the panel
    for i in range(n_bar-1):
        slices.append(counter+s*2)
        guide_position.append(0)

        slices.append(counter + s*2 + bar_width*2)
        guide_position.append(1)

        counter += s*2+bar_width*2

    # finish with eventual additional residual space
    if residual_space:
        slices.append(counter + residual_space)
        guide_position.append(0)
        
    return residual_space, slices, guide_position


# take data from histogrammed panel surface and slice it into a given number of light guides according to
# the their size and number; then apply single light guide Photon Detection Efficiency (PDE)
# output is the number of PE detected per light guide per panel per event
def slicing(df, n_bar, detection_efficiency, time_window = 0):
    
    residual_space, slices, guide_position = get_guides(300, 10, n_bar)
    
    if time_window:
        bins = np.arange(0, 0.005, time_window)
        df['time_bin'] = pd.cut(df['time'], bins=bins)
        grouped_zband = df.groupby(["EventID", "panel", "time_bin"]).zband
        grouped_zband_tolist = grouped_zband.apply(lambda d: list(d)).dropna()
        grouped_df = grouped_zband_tolist.apply(np.sum, axis = 0)
    else:
        grouped_df = df.groupby(["EventID", "panel"]).zband.sum()
        sliced_df = grouped_df.apply(lambda d: np.split(d, slices))

    sliced_df = grouped_df.apply(lambda d: np.split(d, slices))

    if residual_space:
        sliced_df = sliced_df.apply(lambda d: [np.sum(array)*pos for array, pos in zip(d[1:-1], guide_position)])
    else:
        sliced_df = sliced_df.apply(lambda d: [np.sum(array)*pos for array, pos in zip(d, guide_position)])

    sliced_df = sliced_df.apply(lambda d: np.random.binomial(d, detection_efficiency))

    return sliced_df


# take output from function "slicing"
# and sum over all PE detected in a single panel/surface
def slice_panel(df, n_bar, detection_efficiency):
        
    sliced_df = slicing(df, n_bar, detection_efficiency)

    photons_per_panel = sliced_df.apply(lambda d: sum(d))
    photons_per_panel = photons_per_panel[photons_per_panel != 0]
    
    return photons_per_panel  


# this function computes the number of contiguous panels being hit
# 0 in the output means that are no events actually detected
# 1 means that only one panel was hit
# > 1 gives exactly the number of contigous panels (between 2 and max number of contiguous surfaces)
# WARNING: at the moment, surfaces outside the moderator are numbered [1,12] and inside [13,24]
# surfaces 12 and 13 might be considered as contigous, 
# please be careful and shift the inner guides (e.g. +100) when passing data to this function 
def compute_contigous_surface(panels_hit):
    
    if panels_hit == 0:
        return 0

    diff = np.diff(panels_hit)
    stop = np.where(diff != 1)[0]

    nof_contiguous_panels = []

    # if stop array has length 0, means either that we have only panel
    # or that all panels in the initial array are contiguous
    
    #print(stop)
    
    if len(stop > 0):
        
        position = 0

        for s in stop:
            
            nof_contiguous_panels.append(len(panels_hit[position : s + 1]))
            position = s + 1

        nof_contiguous_panels.append(len(panels_hit[position : ]))
        
    else:

        nof_contiguous_panels.append(len(panels_hit))

 
    return_list = nof_contiguous_panels
    
    if (panels_hit[0] == 1) & (panels_hit[-1] == 12):

        return_list = []
        
        if nof_contiguous_panels[1:-1]: 
        
            return_list = nof_contiguous_panels[1:-1]
            return_list.append(nof_contiguous_panels[0] + nof_contiguous_panels[-1])
        
        else:
                
            return_list.append(nof_contiguous_panels[0] + nof_contiguous_panels[-1])

    return return_list
