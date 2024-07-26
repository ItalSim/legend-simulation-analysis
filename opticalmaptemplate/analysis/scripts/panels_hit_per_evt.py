import os
import math
import pickle
import uproot
import random
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm.notebook import tqdm
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm, Normalize
import opanalysis

new_path = "/lfs/l1/legend/users/cbarton/simulations/campaigns/opmapprocessing/24-07-18-Ar41+H2-350mwlsrcryo-allmaps/output"
files = [file for file in sorted(os.listdir(new_path)) if "appliedmap" in file]

n_files = "all"

panel_df = opanalysis.load_data(new_path, "panel", n_files)
lid_df = opanalysis.load_data(new_path, "lid", n_files)

panel_df = panel_df[panel_df.protonnumber == 18]
lid_df   = lid_df[lid_df.protonnumber == 18]

lid_df.rename(columns={"xband" : "zband"}, inplace = True)
lid_df["panel"] = lid_df.panel + 24

array1 = sorted(panel_df.EventID.unique())
array2 = sorted(lid_df.EventID.unique())

only_in_array1 = np.setdiff1d(array1, array2)
only_in_array2 = np.setdiff1d(array2, array1)
common_elements = np.intersect1d(array1, array2)

evt_list = np.unique(np.concatenate((only_in_array1, only_in_array2, common_elements)))

distances_data = {"eff" : [], "n_bar" : [], "distance_naiv": [], "distance_geom": []} 

detection_efficiency_list = [0.005, 0.01, 0.1, 1]

for eff in detection_efficiency_list:

    nof_guides_list = [2,3,4,5]
    l = len(nof_guides_list)
    fig, ax = plt.subplots(1,l,figsize = (5*l,4))

    for bar_idx, n_bar in enumerate(tqdm(nof_guides_list)):
        # divide panels into outer and inner and then divide single panel
        # into n_bar light guides and apply to each guide PDE eff

        outer_top_lids = (lid_df.panel >= 1 + 24) & (lid_df.panel <= 12 + 24)
        outer_bot_lids = (lid_df.panel >= 37 + 24) & (lid_df.panel <= 48 + 24)
        inner_top_lids = (lid_df.panel >= 13 + 24) & (lid_df.panel <= 24 + 24)
        inner_bot_lids = (lid_df.panel >= 25 + 24) & (lid_df.panel <= 36 + 24)
        
        photons_per_panel_inner_lateral = opanalysis.slice_panel(panel_df[panel_df.panel > 12],  300, 12, 10, eff)
        photons_per_panel_outer_lateral = opanalysis.slice_panel(panel_df[panel_df.panel <= 12], 300, 12, 10, eff)
        photons_per_panel_inner_lids    = opanalysis.slice_panel(lid_df[inner_top_lids | inner_bot_lids], 52, n_bar, 10, eff)
        photons_per_panel_outer_lids    = opanalysis.slice_panel(lid_df[outer_top_lids | outer_bot_lids], 52, n_bar, 10, eff)

        photons_per_panel_inner_lateral = photons_per_panel_inner_lateral.groupby(level=0).size()
        photons_per_panel_inner_lids    = photons_per_panel_inner_lids.groupby(level=0).size()   
        photons_per_panel_outer_lateral = photons_per_panel_outer_lateral.groupby(level=0).size()
        photons_per_panel_outer_lids    = photons_per_panel_outer_lids.groupby(level=0).size()   
        
        photons_per_panel_inner_lateral = photons_per_panel_inner_lateral.reindex(evt_list, fill_value=0)
        photons_per_panel_inner_lids    = photons_per_panel_inner_lids.reindex(evt_list, fill_value=0)
        photons_per_panel_outer_lateral = photons_per_panel_outer_lateral.reindex(evt_list, fill_value=0)
        photons_per_panel_outer_lids    = photons_per_panel_outer_lids.reindex(evt_list, fill_value=0)

        inn = photons_per_panel_inner_lateral + photons_per_panel_inner_lids
        out = photons_per_panel_outer_lateral + photons_per_panel_outer_lids

        n = 36
        
        h = ax[bar_idx].hist2d(inn, out, range = [[0,n], [0,n]], bins = n, vmin = 0)
        plt.colorbar(h[3])
        
        naiv_median = (np.median(inn), np.median(out))
        ax[bar_idx].plot(naiv_median[0], naiv_median[1], color = "red", marker = "*", markersize = 10)

        combined_array = [[a, b] for a, b in zip(inn, out)]
        geom_median = opanalysis.geometrical_median(combined_array)
        ax[bar_idx].plot(geom_median[0], geom_median[1], color = "#07A9FF", marker = "*", markersize = 10)

        line = (1, -1, 0)  # Line equation: x - y = 0
        distance_geom = opanalysis.distance_point_line(geom_median, line)
        distance_naiv = opanalysis.distance_point_line(naiv_median, line)

        distances_data["eff"].append(eff)
        distances_data["n_bar"].append(n_bar)
        distances_data["distance_naiv"].append(distance_naiv)
        distances_data["distance_geom"].append(distance_geom)

        ax[bar_idx].set_title(f"12 (lateral) - {n_bar} (lids) - PDE:{eff*100}%")
        ax[bar_idx].set_xlabel("inner panels hit")
        ax[bar_idx].set_ylabel("outer panels hit")
        ax[bar_idx].set_xticks(np.arange(0,40,5))
        ax[bar_idx].set_yticks(np.arange(0,40,5))
        ax[bar_idx].axis('scaled')
        
        t = np.arange(0, n, 0.1)
        ax[bar_idx].plot(t,t, color = "red")

    plot_name = f"plots/panels_hit_per_evt/ar41_lateral_and_lids_{int(eff*1000)}"

    plt.savefig(plot_name + ".png")
    plt.savefig(plot_name + ".pdf")

    with open(plot_name + ".pkl", 'wb') as file:
        pickle.dump(fig, file)

data = pd.DataFrame(distances_data)
data.to_hdf(path_or_buf = "panels_hit_per_evt.h5", key = "data")