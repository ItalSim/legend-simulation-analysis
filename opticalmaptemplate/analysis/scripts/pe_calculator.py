import os
import math
import uproot
import random
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm, Normalize


# input file
# just substitute it with the file you want to analyze
file = "/lfs/l1/legend/users/cbarton/data/appliedmap/appliedmap.root"

df = pd.DataFrame()
tf = uproot.open(file)["appliedmapinfo"]
df["EventID"]  = tf["EventID"].array(library="np")
df["sector"]   = tf["sector"].array(library="np")
df["panel"]    = tf["panel"].array(library="np")
df["zband"]    = tf["zband"].array(library="np")

total_nof_events = len(df.EventID.unique())

# given height of PMMA panel, light guide width and number of bars
# this function returns the z position of the light guides along the panel
def get_guides(H, bar_width, n_bar):
    
    if bar_width * n_bar > H:
        raise ValueError('A very specific bad thing happened.')
    # total free space between light guides
    S = (H - n_bar*bar_width)
    
    # single space between to adjacent light guides
    s = math.floor(S/(n_bar - 1))

    # residual space (if != 0 you basically split it in 2 and add it at the top and bottom to make to guides fit the panel and be evenly spaced)
    residual_space = H - n_bar*bar_width - s*(n_bar - 1)

    slices = []
    guide_position = []
    
    # if residula space is not zero, skip corresponding zbands
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



plt.figure(figsize = (8,6))

# random number of colors and markers
# different (color, marker) combination for different number of light guides per panel
color = sns.color_palette("tab10")
marker = ["o", "s", "d", "^", "v", "*", "P", "+", "x"]

total_pe_dict = {}

# set number of light guides per panel
# with 10cm high bars and 300cm panel, maximum is 30 
nof_guides_list = [2, 5, 12, 20, 30]

for bar_idx, n_bar in tqdm(enumerate(nof_guides_list)):
    
    residual_space, slices, guide_position = get_guides(300, 10, n_bar)

    # probability for a 128nm photon to be detected at one end of the light guide (from external simulations)
    detection_efficiency = 0.003

    # loop over different PE threshold on single light guide
    for PE_threshold in np.arange(0, 10, 1):

        multplicity_list = [1]

        # loop over different multiplicity condition 
        for M_idx, M in enumerate(multplicity_list):
            efficiency_list = []
            efficiency_counter = 0

            # loop over events
            for e in df.EventID.unique():
                
                # for some reason in the Vancouver analysis this event was making everything crash, had to skip it
                # can be removed (or commented)
                if e == 34724:
                    continue
                
                # select data for specific EventID, group by panel and sum over number of photons in each panel
                evt_df = df[df.EventID == e].groupby("panel").zband.sum()

                # resample df, reduce number of zbands
                resampled_df = evt_df.apply(lambda d: np.split(d, slices))

                # condition on how to treat the extra residual space at top and bottom
                if residual_space:
                    resampled_df = resampled_df.apply(lambda d: [np.sum(array)*pos for array, pos in zip(d[1:-1], guide_position)])
                else:
                    resampled_df = resampled_df.apply(lambda d: [np.sum(array)*pos for array, pos in zip(d, guide_position)])
                
                    
                # apply effificiency by multiplying number of photons arrived in each bar with binomial distribution 
                efficiency_applied = resampled_df.apply(lambda d: np.random.binomial(d, detection_efficiency))
                
                # surviving threshold
                PE_thr_applied = efficiency_applied.apply(lambda d: [a for a in d if a > PE_threshold])
        
                # Majority condition (M = 1 means "single light guide trigger" which corresponds to the case no majority activated)
                M = 1
                
                is_detectable = PE_thr_applied.apply(lambda d: len(d) >= M).any()
                
                if is_detectable:
                    efficiency_counter += 1

                # this is to put the label only once (kind of stupid way)
                label = None
                if PE_threshold == (10 - 2):
                    label = f"light guides per panel: {n_bar}"
            
            plt.plot(PE_threshold, efficiency_counter / total_nof_events * 100, marker = marker[bar_idx], color = color[bar_idx], label = label, linestyle = "")

plt.xlabel("PE")
plt.ylabel("tagging efficiency [%]")
plt.legend()
plt.savefig("../output/tagging_efficiency.png")
plt.savefig("../output/tagging_efficiency.pdf")