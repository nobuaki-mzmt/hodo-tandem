"""
data_prep_heatmap.py
N. Mizumoto
Create feather files to plot heatmap using R
"""

import sys
import os
import glob

import h5py

import numpy as np
from numpy.linalg import norm

import pandas as pd

import scipy
from scipy.interpolate import interp1d

import math

import matplotlib.pyplot as plt

#0: antennatipr, 1: antennabaser, 2: antennatipl, 3: antennabasel
#4: headtip, 5: pronotumfront, 6: pronotumend, 7: abdomentip

#------------------------------------------------------------------------------#
# get center and rotation for whole body
#------------------------------------------------------------------------------#
def body_centerize(loc, f_ind, a, b):
  dir_vec = loc[:, a, :, f_ind] - loc[:, b, :, f_ind]
  if a == 7:
    dir_vec = loc[:, b, :, f_ind] - loc[:, a, :, f_ind]
  center_vec = loc[:, a, :, f_ind]
  relative_loc = loc.copy()
  for i_ind in range(loc.shape[3]):
    for i_nodes in range(loc.shape[1]):
      relative_loc[:, i_nodes, :, i_ind] = loc[:, i_nodes, :, i_ind] - center_vec
  rotated_loc = relative_loc.copy()
  
  for i_ind in range(loc.shape[3]):
    for i_nodes in range(loc.shape[1]):
      input_vectors = relative_loc[:, i_nodes, :, i_ind]
  
      angle_rad_input = np.arctan2(input_vectors[:, 1], input_vectors[:, 0])
      angle_rad_dir   = np.arctan2(dir_vec[:, 1], dir_vec[:, 0])
      new_angle = angle_rad_input - angle_rad_dir + math.pi/2
      
      new_angle[new_angle > math.pi] -= 2*math.pi
      new_angle[new_angle < -math.pi] += 2*math.pi
      
      v_lengths = np.linalg.norm(input_vectors, axis=1)
      
      rotated_loc[:, i_nodes, :, i_ind] = np.column_stack((v_lengths*np.cos(new_angle), v_lengths*np.sin(new_angle)))
  return(rotated_loc)
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
def data_preparation(in_dir, species, experiments, treatment, antenna):
  files = glob.glob(in_dir + "/*.h5")
  
  df = pd.DataFrame()
  df2 = pd.DataFrame()
  
  for f_name in files:
    pair_name = os.path.basename(f_name.replace(".h5", ""))
    print(pair_name)
    ## load data
    with h5py.File(f_name, "r") as f:
      locations = f['locations'][:]
      total_frame = locations.shape[0]
    
  
    ## generate dataframes
    
    for i_ind in range(locations.shape[3]):
      
      rotated_locations_head = body_centerize(locations, i_ind, 4, 5) # relative to head front
      rotated_locations_pron = body_centerize(locations, i_ind, 5, 6) # relative to pronotum front
    
      if antenna == 0 and i_ind == 1:
        rotated_locations_head[:, 0, :, i_ind] = np.nan
        rotated_locations_head[:, 2, :, i_ind] = np.nan
      elif antenna == 1 and i_ind == 1:
        rotated_locations_head[:, 0, :, i_ind] = np.nan
      elif antenna == 15 and i_ind == 1:
        rotated_locations_head[:, 1, :, i_ind] = np.nan
      
      relative_antennabaser_x = rotated_locations_head[:, 1, 0, i_ind]
      relative_antennabaser_y = rotated_locations_head[:, 1, 1, i_ind]
      relative_antennatipr_x = rotated_locations_head[:, 0, 0, i_ind]
      relative_antennatipr_y = rotated_locations_head[:, 0, 1, i_ind]
      relative_antennabasel_x = rotated_locations_head[:, 3, 0, i_ind]
      relative_antennabasel_y = rotated_locations_head[:, 3, 1, i_ind]
      relative_antennatipl_x = rotated_locations_head[:, 2, 0, i_ind]
      relative_antennatipl_y = rotated_locations_head[:, 2, 1, i_ind]
      
      relative_head_x = rotated_locations_pron[:, 4, 0, i_ind]
      relative_head_y = rotated_locations_pron[:, 4, 1, i_ind]
      
      sex = treatment
      if treatment == "FM":
        fTip_mHead_dis = np.sqrt( (locations[:, 7, 0, 0] - locations[:, 4, 0, 1])**2 + 
          (locations[:, 7, 1, 0] - locations[:, 4, 1, 1])**2 )
        #fTip_mHead_dis = fTip_mHead_dis.round(2)
        mTip_fHead_dis = np.sqrt( (locations[:, 7, 0, 1] - locations[:, 4, 0, 0])**2 + 
        (locations[:, 7, 1, 1] - locations[:, 4, 1, 0])**2 )
        #mTip_fHead_dis = mTip_fHead_dis.round(2)
        
        # determine if you are behind the partner
        rotated_locations_partnertip = body_centerize(locations, 0, 7, 4)
        m_behind = rotated_locations_partnertip[:, 4, 1, 1] < 0
        rotated_locations_partnertip = body_centerize(locations, 1, 7, 4)
        f_behind = rotated_locations_partnertip[:, 4, 1, 0] < 0
        
        #
        fTip_x = locations[:, 7, 0, 0].round(2)
        mHead_x = locations[:, 4, 0, 1].round(2)
        fTip_y = locations[:, 7, 1, 0].round(2)
        mHead_y = locations[:, 4, 1, 1].round(2)
        
        # 
        fspeed = np.append(np.sqrt(np.diff(locations[:, 5, 0, 0])**2 + np.diff(locations[:, 5, 1, 0])**2), np.nan)
        mspeed = np.append(np.sqrt(np.diff(locations[:, 5, 0, 1])**2 + np.diff(locations[:, 5, 1, 1])**2), np.nan)
        facc   = np.append(np.nan, np.diff(fspeed))
        macc   = np.append(np.nan, np.diff(mspeed))

        
        if i_ind == 0:
          sex = "F"
        else:
          sex = "M"
     
      df_temp = {
        "frame": list(range(0,rotated_locations_head.shape[0],1)),
        "relative_antennabaser_x": rotated_locations_head[:, 1, 0, i_ind].round(2),
        "relative_antennabaser_y": rotated_locations_head[:, 1, 1, i_ind].round(2),
        "relative_antennatipr_x":  rotated_locations_head[:, 0, 0, i_ind].round(2),
        "relative_antennatipr_y":  rotated_locations_head[:, 0, 1, i_ind].round(2),
        "relative_antennabasel_x": rotated_locations_head[:, 3, 0, i_ind].round(2),
        "relative_antennabasel_y": rotated_locations_head[:, 3, 1, i_ind].round(2),
        "relative_antennatipl_x":  rotated_locations_head[:, 2, 0, i_ind].round(2),
        "relative_antennatipl_y":  rotated_locations_head[:, 2, 1, i_ind].round(2),
        "relative_head_x":         rotated_locations_head[:, 4, 0, i_ind].round(2),
        "relative_head_y":         rotated_locations_head[:, 4, 1, i_ind].round(2),
        "fTip_mHead_dis":          fTip_mHead_dis,
        "mTip_fHead_dis":          mTip_fHead_dis,
        "m_behind":                  m_behind,
        "f_behind":                  f_behind,
        "sex": sex,
        "video": pair_name
      }
      df_temp = pd.DataFrame(df_temp)
      df = pd.concat([df, pd.DataFrame(df_temp)])
      
      
      angle_antenna_r =  np.arctan2(relative_antennatipr_y - relative_antennabaser_y,
                                   relative_antennatipr_x - relative_antennabaser_x) - math.pi/2
      angle_antenna_l =  np.arctan2(relative_antennatipl_y - relative_antennabasel_y,
                                   relative_antennatipl_x - relative_antennabasel_x) - math.pi/2
      angle_antenna_r[angle_antenna_r < - math.pi] = angle_antenna_r[angle_antenna_r < - math.pi] + math.pi*2
      angle_antenna_l[angle_antenna_l < - math.pi] = angle_antenna_l[angle_antenna_l < - math.pi] + math.pi*2
      
      
      
      df_temp = {
        "frame": list(range(0,rotated_locations_head.shape[0],1)),
        "fTip_x": fTip_x,
        "fTip_y": fTip_y,
        "mHead_x": mHead_x,
        "mHead_y": mHead_y,
        "angle_antenna_r": angle_antenna_r.round(2),
        "angle_antenna_l": angle_antenna_l.round(2),
        "fTip_mHead_dis":  fTip_mHead_dis,
        "mTip_fHead_dis":  mTip_fHead_dis,
        "m_behind":                  m_behind,
        "f_behind":                  f_behind,
        "fspeed": fspeed,
        "mspeed": mspeed,
        "facc": facc,
        "macc": macc,
        "sex": sex,
        "video": pair_name
      }
      df_temp = pd.DataFrame(df_temp)
      df2 = pd.concat([df2, pd.DataFrame(df_temp)])


      
  filename = "2D_hist_" + species + "_" + os.path.basename(in_dir) + "_df.feather"
  #df.to_hdf("data_fmt/" + filename, key = "df", mode='w',  complevel=9, complib='zlib')
  df.reset_index().to_feather("data_fmt/" + filename)
  
  filename = "antenna_angle_" + species + "_" + os.path.basename(in_dir) + "_df.feather"
  #df.to_hdf("data_fmt/" + filename, key = "df", mode='w',  complevel=9, complib='zlib')
  df2.reset_index().to_feather("data_fmt/" + filename)
    
  return 0
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
def main():
  place = "data_fmt/data_extract/*"
  data_place_species = glob.glob(place)
  
  for data_species_i in data_place_species:
    species = os.path.basename(data_species_i)
    data_place = glob.glob(data_species_i + "/*")
    
    for data_place_i in data_place:
      print(data_place_i)
      experiments = os.path.basename(data_place_i).split("_")[0]
      treatment   = os.path.basename(data_place_i).split("_")[1]
      antenna     = os.path.basename(data_place_i).split("_")[2]
      print(species + ", " + experiments + ", " + treatment + ", antenna:"+antenna)
      
      data_preparation(data_place_i, species, experiments, treatment, antenna)
    
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
def main2(place, species):
  data_place_i = "data_fmt/data_extract/" + place
  print(data_place_i)
  experiments = os.path.basename(data_place_i).split("_")[0]
  treatment   = os.path.basename(data_place_i).split("_")[1]
  antenna     = os.path.basename(data_place_i).split("_")[2]
  print(species + ", " + experiments + ", " + treatment + ", antenna:"+antenna)
  
  data_preparation(data_place_i, species, experiments, treatment, antenna)
    
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
if __name__ == "__main__":
  if len(sys.argv) > 1:
    main2(sys.argv[1], sys.argv[2])
  else:
    # Call the main function without any parameter
    main()
#------------------------------------------------------------------------------#
