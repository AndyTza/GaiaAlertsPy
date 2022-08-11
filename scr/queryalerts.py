import pandas as pd
import warnings
from astropy.time import Time
import numpy as np


def gaia_g_noise(mag):
    """TODO: Source & logistics"""
    return 3.43779 - (mag/1.13759) + (mag/3.44123)**2 - (mag/6.51996)**3 + (mag/11.45922)**4


def query_lightcurve(cand_id):
    """
    TODO: Add docs
    """
    # TODO: asses if file empty 
    _dat = pd.read_csv(f"http://gsaweb.ast.cam.ac.uk/alerts/alert/{cand_id}/lightcurve.csv/", skiprows=1)
    
    #assert _dat:
        #raise ValueError(f"Sorry, {cand_id} was not found.")

    master_frame = np.zeros(shape=(_dat.shape[0], 3))
    
    for i, vals in enumerate(_dat.to_numpy()):
        time_conv = Time(vals[0], format='iso').mjd
        master_frame[i,0] = time_conv

        if vals[-1]=='untrusted' or vals[-1]=='nan':
            master_frame[i,1] = np.nan
        else:
            master_frame[i, 1] = float(vals[-1])

    return master_frame

