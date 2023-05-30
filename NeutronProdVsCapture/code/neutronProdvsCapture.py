#! /usr/bin/env python3

import uproot
import numpy as np
from matplotlib import pyplot as plt, colors
from matplotlib.ticker import (LinearLocator, MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle

import time
import urllib.parse
import sys

#%config Completer.use_jedi = False

def prodPosVsCap(nCAr_a, nCAr_z, nCAr_ID, nCAr_EventID, NGe77, NeutronID, NeutronEventID, NeutronVolume, Neutronxloc, Neutronyloc, Neutronzloc, verbose = False):
    #First we have to find the ID of the neutrons that have been capture on Ge77
    #with Z = 32 and A = 77.
    
    nIDinGe      = []
    nIDEventinGe = []
    histVolume   = []

    st_nCAr_a       = np.hstack(nCAr_a)
    st_nCAr_z       = np.hstack(nCAr_z)
    st_nCAr_ID      = np.hstack(nCAr_ID)
    st_nCAr_EventID = np.hstack(nCAr_EventID)

    nIDinGe      = st_nCAr_ID[(st_nCAr_a == 77) & (st_nCAr_z == 32)]
    nIDEventinGe = st_nCAr_EventID[(st_nCAr_a == 77) & (st_nCAr_z == 32)]
    #nIDinGe = (st_nCAr_a == 77) & (st_nCAr_z == 32)

    # then we can cross-check this information with the aumount of NGe77

    st_NGe77 = np.hstack(NGe77)
    if sum(st_NGe77) != len(nIDinGe):
        print("Something is odd")

    # Now we have to go the Neutron set of variables and find these IDs there.
    st_NeutronID          = np.hstack(NeutronID)
    st_NeutronEventID     = np.hstack(NeutronEventID)
    st_NeutronVolume      = np.hstack(NeutronVolume)
    st_Neutronxloc        = np.hstack(Neutronxloc)
    st_Neutronyloc        = np.hstack(Neutronyloc)
    st_Neutronzloc        = np.hstack(Neutronzloc)

    if verbose:
        print("N of ge77 (nCap):  ",  len(nIDinGe))
        print("N of ge77 (NGe77): ", sum(st_NGe77))
    #histVolume = []

    for idGe, idEventGe in zip(nIDinGe, nIDEventinGe):
        idBool = ((st_NeutronID == idGe) & (st_NeutronEventID == idEventGe))
        if verbose:
            print("nCar ID:          ", idGe)
            print("Neutron ID:       ", st_NeutronID[idBool])
            print("nCar Event ID:    ", idEventGe)
            print("Neutron Event ID: ", st_NeutronEventID[idBool])
            print("Volume:           ", st_NeutronVolume[idBool])
            print("X pos :           ", st_Neutronxloc[idBool])
            print("Y pos :           ", st_Neutronyloc[idBool])
            print("Z pos :           ", st_Neutronzloc[idBool])
            print("---------------------------------")
        if st_NeutronVolume[idBool] == -4:
            
            Xpos = st_Neutronxloc[idBool]
            Ypos = st_Neutronyloc[idBool]
            Zpos = st_Neutronzloc[idBool]
            
            r = np.sqrt(Xpos**2+Ypos**2)
            
            if ((r<2.0) & (Zpos<0.92) & (Zpos>-2.08)):
                histVolume = histVolume + [-42] ## Internal LAr
            else:
                histVolume = histVolume + [-41] ## External LAr
        else:
            histVolume = histVolume + st_NeutronVolume[idBool].tolist()

    #print(histVolume)
    #histVolume = st_NeutronVolume[nIDinGe].tolist()
    
    return histVolume

def countVolumes(histVolume):
    aux_W    = 0
    aux_C    = 0
    aux_V    = 0
    aux_Le   = 0
    aux_M    = 0
    aux_Li   = 0
    aux_R    = 0
    aux_U    = 0
    aux_G    = 0
    aux_else = 0
    
    for volume in histVolume:
        if volume == -10:
            aux_W += 1
        elif (volume == -9) or (volume == -8) or (volume == -7) or (volume == -5):
            aux_C += 1
        elif volume == -6:
            aux_V += 1
        elif volume == -41:
            aux_Le += 1
        elif volume == -42:
            aux_Li += 1
        elif (volume == -3) or (volume == -2):
            aux_M += 1
        elif volume == -1:
            aux_R += 1
        elif volume == 0:
            aux_U += 1
        elif volume == 1:
            aux_G += 1
        else:
            aux_else += 1
    
    counts = [aux_else, aux_W, aux_C, aux_V, aux_Le, aux_M, aux_Li, aux_R, aux_U, aux_G]
    return counts
    
    
def loadAndCalcProd(final,path,verbose=False):
    histVolume  = []
    n_erro = 0
    for ind in range (0, final):

        if verbose:
            print("Loading and Processing %d out of %d files" % (ind,(final)), end="\r")
        filename = path+"design2_%d_t0.root:Score" % ind
        try:
            tf = uproot.open(filename)
        except:
            n_erro = n_erro + 1
            continue

        nCAr_a             = tf["nCapture_A"].array(library="np")
        nCAr_z             = tf["nCapture_Z"].array(library="np")
        nCAr_ID            = tf["nCapture_ID"].array(library="np")
        nCAr_EventID       = tf["nCapture_EventID"].array(library="np")


        Neutronxloc        = tf["Neutronxloc"].array(library="np")
        Neutronyloc        = tf["Neutronyloc"].array(library="np")
        Neutronzloc        = tf["Neutronzloc"].array(library="np")
        NeutronID          = tf["NeutronID"].array(library="np")
        NeutronEventID     = tf["NeutronEventID"].array(library="np")
        NeutronVolume      = tf["NeutronVolume"].array(library="np")

        NGe77              = tf["NGe77"].array(library="np")

        ## Analysis code
        if verbose:
            print("File %d" % ind)

        histVolume_aux  = prodPosVsCap(nCAr_a, nCAr_z, nCAr_ID, nCAr_EventID, NGe77, NeutronID, NeutronEventID, NeutronVolume, Neutronxloc, Neutronyloc, Neutronzloc, verbose)
        histVolume     = histVolume + histVolume_aux

        del tf


    t1 = time.time()
    if verbose:
        print("Load time for %d files took %.2f seconds (%.2f minutes)" % (ind+1, t1-t0, (t1-t0)/60))
        print("Number of files: %d" % (final-n_erro))
    
    if final-n_erro == 0:
        print("Not a single file found on the given folder, exiting...")
        sys.exit(1)
        
    return histVolume


def plotProdvsCap(histVolume_NoMod, histVolume_Mod, flag_plot):

    import math
    import matplotlib.colors as mcolors

    font = 12

    volume_labels = ['Everything else', 'Water Tank', 'Cryostat', 'Vacuum Gap', 'LAr External Mod', 'Moderator', 'LAr Internal Mod', 'Reentrant Tube', 'Underground LAr', 'Ge Detectors']
    volume_ticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    fig = plt.figure(figsize=(9, 6))
    plt.title("Production volume of the neutron captured on Ge detectors", fontsize=font)
    ax = plt.gca()
    
    if flag_plot == 'noMod':
        count1 = countVolumes(histVolume_NoMod)
         # plot the second histogram with dashed bars
        plt.bar(volume_ticks, count1, width=0.8, bottom=None, align='center', color='white', edgecolor=mcolors.TABLEAU_COLORS, linestyle='--')#, hatch='////')
        ax.bar_label(ax.containers[0], label_type='edge')
    elif flag_plot == 'Mod':
        count2 = countVolumes(histVolume_Mod)
        # plot the first histogram with solid bars
        plt.bar(volume_ticks, count2, width=0.8, bottom=None, align='center', color=mcolors.TABLEAU_COLORS)
        ax.bar_label(ax.containers[0], label_type='edge')
    else:
        count1 = countVolumes(histVolume_NoMod)
        count2 = countVolumes(histVolume_Mod)
        
        # plot the second histogram with dashed bars
        plt.bar(volume_ticks, count1, width=0.8, bottom=None, align='center', color='white', edgecolor=mcolors.TABLEAU_COLORS, linestyle='--')#, hatch='////')
        
        # plot the first histogram with solid bars
        plt.bar(volume_ticks, count2, width=0.8, bottom=None, align='center', color=mcolors.TABLEAU_COLORS)
        
         # add labels to the solid bars
        ax.bar_label(ax.containers[0], label_type='edge')
        ax.bar_label(ax.containers[1], label_type='edge')

    plt.xlabel("Detector volume", fontsize=font, horizontalalignment='right', x=1.0)
    plt.ylabel("Number of Ge77 captures by neutrons", fontsize=font, horizontalalignment='right', y=1.0)
    
    ax.tick_params(axis='both', labelcolor='k', labelsize="large")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='major', direction="in", length=8, width=1)
    ax.tick_params(which='minor', direction="in", length=3, width=1)

    ax.set_xticks(volume_ticks)
    ax.set_xticklabels(volume_labels, rotation=45, ha="right")

    plt.ylim([0, 200])

    if   flag_plot == 'noMod':
        plt.legend(['Without moderators'])
    elif flag_plot == 'Mod':
        plt.legend(['With moderators'])
    else:
        plt.legend(['Without moderators', 'With moderators'])

    plt.savefig('../results/ProductionNeutronCaptured_{flag}.{ext}'.format(flag = flag_plot, ext='pdf'), dpi = 300, bbox_inches='tight', pad_inches=0.1)
    plt.savefig('../results/ProductionNeutronCaptured_{flag}.{ext}'.format(flag = flag_plot, ext='png'), dpi = 300, bbox_inches='tight', pad_inches=0.1)

    
def main(flag_plot, path1, path2, number_files=200, verbose=False):
    
    t0 = time.time()
    import pandas as pd

    # Loop over "n" files to load the variables
    final = number_files # Number of files to be opened at the same time

    histVolume_NoMod = []
    histVolume_Mod   = []

    if flag_plot == 'Mod':
        #path = "/lfs/l1/legend/users/iabritta/Single_550cmInnerCryostat/PMMA_200_latest/"
        path = path1
        histVolume_Mod = loadAndCalcProd(final,path,verbose=False)

    elif flag_plot == 'noMod':
        #path = "/lfs/l1/legend/users/iabritta/Single_550cmInnerCryostat/NoMod/"
        path = path1
        histVolume_NoMod = loadAndCalcProd(final,path,verbose=False)
        
    elif flag_plot == 'both':
        path_Mod = path2
        histVolume_Mod = loadAndCalcProd(final,path_Mod,verbose=False)

        path_NoMod = path1
        histVolume_NoMod = loadAndCalcProd(final,path_NoMod,verbose=False)

    else:
        print("Choose a valid plot option")

    plotProdvsCap(histVolume_NoMod, histVolume_Mod, flag_plot)
    

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage='usage: %prog\t [-ubsv] get/put Localpath Key')
    parser.add_option('-f','--flag-plot', dest='flag_plot', type='string', default='both', 
                      help='Choose between "noMod", "Mod" or "both"');
    parser.add_option('-p','--path1', dest='path1', type='string', default=None, 
                      help='path to your data');
    parser.add_option('-d','--path2', dest='path2', type='string', default=None, 
                      help='path to your data');
    parser.add_option('-n','--number', dest='number_files', type='int', default=200, 
                      help='Number of files to analyse');
    parser.add_option('-v','--verbose', dest='verbose', action="store_true", default=False, help='verbose output;');
    (options, args) = parser.parse_args()
    
    if (options.flag_plot == 'both') and (options.path1 is not None) and (options.path2 is not None):
        
        if options.verbose:
            print("[Warning] Remember to give NoMod path as path1 and Mod as path2")
        main(options.flag_plot, options.path1, options.path2, options.number_files, options.verbose)
        

    elif (options.flag_plot == 'noMod') and (options.path1 is not None):
        if options.verbose:
            print("Calculation Production for NoMod case")
        main(options.flag_plot, options.path1, options.path2, options.number_files, options.verbose)
        
    elif (options.flag_plot == 'Mod') and (options.path1 is not None):
        if options.verbose:
            print("Calculation Production for Mod case")
        main(options.flag_plot, options.path1, options.path2, options.number_files, options.verbose)
        
    else:
        parser.error("incorrect number of arguments")
        
    ##### Example of use:
    
    # python3 neutronProdvsCapture.py -f both --path1 /lfs/l1/legend/users/iabritta/Single_550cmInnerCryostat/NoMod/ --path2 /lfs/l1/legend/users/iabritta/Single_550cmInnerCryostat/PMMA_200_latest/ -n 20 -v
        