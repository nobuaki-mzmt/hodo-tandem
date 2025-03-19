"""
data_prep_filtering.py
N. Mizumoto
This script reads all .h5 results from SLEAP and organize for the further analysis
"""

import glob
import os

import pandas as pd

import h5py

import numpy as np
import scipy
from scipy.interpolate import interp1d

#------------------------------------------------------------------------------#
# interpolate the data
#------------------------------------------------------------------------------#
def fill_missing(Y, kind="linear"):
    initial_shape = Y.shape
    Y = Y.reshape((initial_shape[0], -1))
    # Interpolate along each slice.
    for i in range(Y.shape[-1]):
        y = Y[:, i]
        # Build interpolant.
        x = np.flatnonzero(~np.isnan(y))
        if len(x) > 3:
          f = interp1d(x, y[x], kind=kind, fill_value=np.nan, bounds_error=False)
          # Fill missing
          xq = np.flatnonzero(np.isnan(y))
          y[xq] = f(xq)
          # Fill leading or trailing NaNs with the nearest non-NaN values
          mask = np.isnan(y)
          y[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), y[~mask])
          Y[:, i] = y
          if sum(np.isnan(y)) > 0:
            print("error"+str(i))
            print("error"+(i))
    # Restore to initial shape.
    Y = Y.reshape(initial_shape)
    return Y
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
def data_filter(in_dir, dish_size, species):
  experiments = os.path.basename(in_dir).split("_")[0]
  files = glob.glob(in_dir + "/*.h5")
  df_pair_behavior = pd.DataFrame()
  df_video = pd.DataFrame()
  antenna   = os.path.basename(in_dir).split("_")[2]
  for f_name in files:
    ## load data
    with h5py.File(f_name, "r") as f:
      dset_names = list(f.keys())
      locations = f["tracks"][:].T
      node_names = [n.decode() for n in f["node_names"][:]]
    
    if locations.shape[3] > 2:
      locations = locations[:,:,:,0:2]
    
    # downsample to 30FPS
    if experiments == 'control':
      if species != "Hod-sjo":
        locations = locations[::2, :,:,:]
    
    print(f_name)

    ## processing locations
    # data filling
    locations = fill_missing(locations)
    
    # scaling in mm (2000 pixels = dish_size)
    locations[:, :, :, :] = locations[:, :, :, :] / 2000 * dish_size
    
    # filtering
    for i_ind in range(locations.shape[3]):
      for i_coord in range(locations.shape[2]):
        for i_nodes in range(locations.shape[1]):
          locations[:, i_nodes, i_coord, i_ind] = scipy.signal.medfilt( locations[:, i_nodes, i_coord, i_ind], 5)
    
    if (antenna == "00"):
      locations[:, 0, :, 1] = np.nan
      locations[:, 1, :, 1] = np.nan
      locations[:, 3, :, 1] = np.nan
      locations[:, 4, :, 1] = np.nan
    elif (antenna == "10"):
      locations[:, 0, :, 1] = np.nan
      locations[:, 1, :, 1] = np.nan
    elif (antenna == "15"):
      locations[:, 0, :, 1] = locations[:, 1, :, 1]

    
    hdf5_file_path = f_name.replace("raw", "fmt/data_filter")
    with h5py.File(hdf5_file_path, 'w') as hdf5_file:
      hdf5_file.create_dataset('locations', data=locations)
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#

def main_data_filter(place=None):
  if place is None:
    place = "data_raw/*"
  else:
    place = "data_raw/" + place + "/*"
  data_place_species = glob.glob(place)
  for data_place_species_i in data_place_species:
    print(data_place_species_i)
    species = os.path.basename(data_place_species_i)
    data_place = glob.glob(data_place_species_i + "/*")
    for data_place_i in data_place:
      print(data_place_i)
      if (species == "Ret-spe"):
        dish_size = 55
      else:
        dish_size = 90
      
      data_filter(in_dir = data_place_i, dish_size = dish_size, species = species)
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
if __name__ == "__main__":
    main_data_filter()
#------------------------------------------------------------------------------#
