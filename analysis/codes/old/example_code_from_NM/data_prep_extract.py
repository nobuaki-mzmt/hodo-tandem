"""
data_preparation.py
N. Mizumoto
This script extracts body part data which is in our interests of the analysis, including
  antennatipr, antennabaser, antennatipl, antennabasel, headtip, pronotumfront,
  pronotumend, abdomentip
Also, measure the body size and antenna length of each termite.
"""

import sys
import os
import glob

import h5py

import numpy as np
from numpy.linalg import norm

import pandas as pd
import pickle

import scipy
from scipy.interpolate import interp1d

import math

#0: antennatipr, 1: antennamiddler, 2: antennabaser
#3: antennatipl, 4: antennamiddlel, 5: antennabasel
#6: headtip, 7: pronotumfront, 8: pronotumend
#9: abdomenfront, 10: abdomentip, 
#11: frontlegl, 12: middlelegr, 13: middlelegl
#14: hindlegr, 15: hindlegl, 16: frontlegr


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
def data_extract(in_dir, species, experiments, treatment, antenna):
  files = glob.glob(in_dir + "/*.h5")
  
  df_bl_al   = pd.DataFrame()
  
  for f_name in files:
    print(f_name)
    pair_name = os.path.basename(f_name.replace(".h5", ""))
    
    ## load data
    with h5py.File(f_name, "r") as f:
      locations = f['locations'][:]
      total_frame = locations.shape[0]
    
    # measure body length and antenna length
    body_length, antenna_r_length, antenna_l_length = [], [], []
    for ind_i in range(locations.shape[3]):
      body_length_temp = ( 
        np.sqrt( (locations[:, 6, 0, ind_i] - locations[:, 7, 0, ind_i])**2 +
              (locations[:, 6, 1, ind_i] - locations[:, 7, 1, ind_i])**2) +
        np.sqrt( (locations[:, 7, 0, ind_i] - locations[:, 8, 0, ind_i])**2 +
              (locations[:, 7, 1, ind_i] - locations[:, 8, 1, ind_i])**2) +
        np.sqrt( (locations[:, 8, 0, ind_i] - locations[:, 9, 0, ind_i])**2 +
              (locations[:, 8, 1, ind_i] - locations[:, 9, 1, ind_i])**2) +
        np.sqrt( (locations[:, 9, 0, ind_i] - locations[:, 10, 0, ind_i])**2 +
              (locations[:, 9, 1, ind_i] - locations[:, 10, 1, ind_i])**2) 
      )
      body_length = np.append(body_length, np.mean(body_length_temp))
      
      antenna_r_length_temp = ( 
        np.sqrt( (locations[:, 0, 0, ind_i] - locations[:, 1, 0, ind_i])**2 +
              (locations[:, 0, 1, ind_i] - locations[:, 1, 1, ind_i])**2) +
        np.sqrt( (locations[:, 1, 0, ind_i] - locations[:, 2, 0, ind_i])**2 +
              (locations[:, 1, 1, ind_i] - locations[:, 2, 1, ind_i])**2)  
      )
      antenna_r_length = np.append(antenna_r_length, np.mean(antenna_r_length_temp))
    
      antenna_l_length_temp = ( 
        np.sqrt( (locations[:, 3, 0, ind_i] - locations[:, 4, 0, ind_i])**2 +
              (locations[:, 3, 1, ind_i] - locations[:, 4, 1, ind_i])**2) +
        np.sqrt( (locations[:, 4, 0, ind_i] - locations[:, 5, 0, ind_i])**2 +
              (locations[:, 4, 1, ind_i] - locations[:, 5, 1, ind_i])**2)  
      )
      antenna_l_length = np.append(antenna_l_length, np.mean(antenna_l_length_temp))
    if treatment == "F":
      body_length      = np.append(body_length, np.nan)
      antenna_r_length = np.append(antenna_r_length, np.nan)
      antenna_l_length = np.append(antenna_l_length, np.nan)
    if treatment == "M":
      body_length      = np.append(np.nan, body_length)
      antenna_r_length = np.append(np.nan, antenna_r_length)
      antenna_l_length = np.append(np.nan, antenna_l_length)

    df_temp = {
        "video":       pair_name,
        "species":     species,
        "antenna":     antenna,
        "experiments": experiments,
        "total_frame": total_frame,
        "body_length_female": body_length[0],
        "body_length_male": body_length[1],
        "antenna_r_length_female": antenna_r_length[0],
        "antenna_r_length_male": antenna_r_length[1],
        "antenna_l_length_female": antenna_l_length[0],
        "antenna_l_length_male": antenna_l_length[1]
        }
    df_bl_al = pd.concat([df_bl_al, pd.DataFrame([df_temp])])
    
    # summarize data
    locations = locations[:, [0, 2, 3, 5, 6, 7, 8, 10], :, :]

    hdf5_file_path = f_name.replace("filter", "extract")
    with h5py.File(hdf5_file_path, 'w') as hdf5_file:
      hdf5_file.create_dataset('locations', data=locations)
    
  return df_bl_al
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
def main():
  df_bl_al = pd.DataFrame()
  data_place_species = glob.glob("data_fmt/data_filter/*")
  
  for data_species_i in data_place_species:
    species = os.path.basename(data_species_i)
    data_place = glob.glob(data_species_i + "/*")
    
    for data_place_i in data_place:
      print(data_place_i)
      experiments = os.path.basename(data_place_i).split("_")[0]
      treatment   = os.path.basename(data_place_i).split("_")[1]
      antenna     = os.path.basename(data_place_i).split("_")[2]
      print(species + ", " + experiments + ", " + treatment + ", antenna:"+antenna)
      df_temp = data_extract(data_place_i, species, experiments, treatment, antenna)
      
      df_bl_al = pd.concat([df_bl_al, pd.DataFrame(df_temp)])
  
  csv_file_path = "data_fmt/df_bl_al.csv"
  df_bl_al.to_csv(csv_file_path, index=False)
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
if __name__ == "__main__":
    main()
#------------------------------------------------------------------------------#
