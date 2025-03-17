import glob
import os
import math

import pandas as pd

import h5py

import numpy as np
import scipy
from scipy.interpolate import interp1d

#------------------------------------------------------------------------------#
def fill_missing(Y, kind="linear"):
    initial_shape = Y.shape
    Y = Y.reshape((initial_shape[0], -1))
    for i in range(Y.shape[-1]):
        y = Y[:, i]
        x = np.flatnonzero(~np.isnan(y))
        if len(x) > 3:
            f = interp1d(x, y[x], kind=kind, fill_value=np.nan, bounds_error=False)
            xq = np.flatnonzero(np.isnan(y))
            y[xq] = f(xq)
            mask = np.isnan(y)
            y[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), y[~mask])
            Y[:, i] = y
    Y = Y.reshape(initial_shape)
    return Y
#------------------------------------------------------------------------------#

def data_filter(in_dir, treatment):
    df = pd.DataFrame()
    df_body = pd.DataFrame()
    files = glob.glob(in_dir + "/*.h5")
    for f_name in files:
        with h5py.File(f_name, "r") as f:
            dset_names = list(f.keys())
            locations = f["tracks"][:].T
            node_names = [n.decode() for n in f["node_names"][:]]
        
        if locations.shape[3] > 2:
            print("too many tracks: " + str(locations.shape[3]))
        
        locations = fill_missing(locations)
        
        for i_ind in range(locations.shape[3]):
            for i_coord in range(locations.shape[2]):
                for i_nodes in range(locations.shape[1]):
                    locations[:, i_nodes, i_coord, i_ind] = scipy.signal.medfilt(locations[:, i_nodes, i_coord, i_ind], 5)
        
        video = os.path.splitext(os.path.basename(f_name))[0]
        
        if 'abdomentip' not in node_names:
            raise ValueError("Error: 'abdomentip' is missing from node_names.")
        
        body_length = [0, 0]
        for i_ind in range(locations.shape[3]):
            head_x = locations[:, node_names.index('headtip'), 0, i_ind]
            head_y = locations[:, node_names.index('headtip'), 1, i_ind]
            tip_x = locations[:, node_names.index('abdomentip'), 0, i_ind]
            tip_y = locations[:, node_names.index('abdomentip'), 1, i_ind]
            center_x = locations[:, node_names.index('abdomenfront'), 0, i_ind]
            center_y = locations[:, node_names.index('abdomenfront'), 1, i_ind]
            body_length[i_ind] = np.median(np.sqrt((head_x - center_x)**2 + (head_y - center_y)**2) +
                                           np.sqrt((tip_x - center_x)**2 + (tip_y - center_y)**2))
        
        df_temp = {
            "video": video,
            "female": body_length[0],
            "male": body_length[1]
        }
        df_temp = pd.DataFrame([df_temp])
        df_body = pd.concat([df_body, df_temp])
        print("locations.shape:", locations.shape)

        locations_abdomen = locations[:, node_names.index('abdomentip'), :, :]
        locations_headtip = locations[:, node_names.index('headtip'), :, :]
        locations = locations[:, node_names.index('abdomenfront'), :, :]
        

        fx = locations[:, 0, 0]
        fy = locations[:, 1, 0]
        mx = locations[:, 0, 1]
        my = locations[:, 1, 1]
        fx_abdomen = locations_abdomen[:, 0, 0]
        fy_abdomen = locations_abdomen[:, 1, 0]
        mx_abdomen = locations_abdomen[:, 0, 1]
        my_abdomen = locations_abdomen[:, 1, 1]
        fx_headtip = locations_headtip[:, 0, 0]
        fy_headtip = locations_headtip[:, 1, 0]
        mx_headtip = locations_headtip[:, 0, 1]
        my_headtip = locations_headtip[:, 1, 1]

        df_temp = {
            "video": video,
            "fx": fx,
            "fy": fy,
            "mx": mx,
            "my": my, 
            "fx_abdomen": fx_abdomen,
            "fy_abdomen": fy_abdomen,
            "mx_abdomen": mx_abdomen,
            "my_abdomen": my_abdomen,
            "fx_headtip": fx_headtip,
            "fy_headtip": fy_headtip,
            "mx_headtip": mx_headtip,
            "my_headtip": my_headtip
        }
        
        df_temp = pd.DataFrame(df_temp)
        df = pd.concat([df, df_temp])
    
    df_body.to_csv('data_fmt/' + treatment + '_bodysize.csv', index=False)
    return df
#------------------------------------------------------------------------------#

def main_data_filter():
    data_place = ["C:/Users/Mizumoto-lab/Desktop/hodo-tandem/analysis/data_raw"]
    for data_place_i in data_place:
        treatment = os.path.basename(data_place_i)
        df = data_filter(in_dir=data_place_i, treatment=treatment)
        filename = treatment + "_df.feather"
        df.reset_index().to_feather("data_fmt/" + filename)
#------------------------------------------------------------------------------#

if __name__ == "__main__":
    main_data_filter()

