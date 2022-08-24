import pandas as pd
import warnings
from astropy.time import Time
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import matplotlib 


base_url = "https://gsaweb.ast.cam.ac.uk/alerts/alert/"

__all__ = [ 'query_lightcurve']


class GaiaAlert:
    def __init__(self, id):
        self.id = id

    def gaia_g_noise_esitmate(self, mag):
        """TODO: Source & logistics"""
        return 3.43779 - (mag/1.13759) + (mag/3.44123)**2 - (mag/6.51996)**3 + (mag/11.45922)**4

    
    def query_lightcurve(self):
        """
        TODO: Add docs
        """

        try:
            _dat = pd.read_csv(f"http://gsaweb.ast.cam.ac.uk/alerts/alert/{self.id}/lightcurve.csv/", skiprows=1)
        except:
            raise ValueError("Sorry, the Gaia alert ID you queried was not found.") 

        master_frame = np.zeros(shape=(_dat.shape[0], 3))
        
        for i, vals in enumerate(_dat.to_numpy()):
            time_conv = Time(vals[0], format='iso').mjd
            master_frame[i,0] = time_conv

            if vals[-1]=='untrusted' or vals[-1]=='nan':
                master_frame[i,1] = np.nan
            else:
                master_frame[i, 1] = float(vals[-1])
        
        indx_nans = ~np.isnan(master_frame[:,1])

        return Table([master_frame[:,0][indx_nans], 
                    master_frame[:,1][indx_nans], 
                    self.gaia_g_noise_esitmate(master_frame[:,1][indx_nans])], 
                    names=("mjd", "mag_G", "mag_G_error")) 