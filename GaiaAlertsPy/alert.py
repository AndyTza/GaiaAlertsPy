from bs4 import BeautifulSoup
import pandas as pd
import warnings
from astropy.time import Time
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import matplotlib 
import re
import astropy.utils.data as dt
import requests
from bs4 import BeautifulSoup


base_url = "https://gsaweb.ast.cam.ac.uk/alerts/alert/"

__all__ = [ 'query_lightcurve']


class GaiaAlert:
    def __init__(self, id):
        self.id = id

    def gaia_g_noise_esitmate(self, mag):
        """TODO: Source & logistics"""
        return 3.43779 - (mag/1.13759) + (mag/3.44123)**2 - (mag/6.51996)**3 + (mag/11.45922)**4

    def query_bprp_history(self):
        """
        Will return the historic BP_RP of each alert!

        This code was used directly from Hogg & Sipőcz 

        TODO: Add docs
        """

        # Scrape images
        htmllines = dt.get_file_contents(base_url + self.id).splitlines()
        
        # Scrape spectra!
        for line in htmllines:
            if "var spectra" in line:
                break
        
        literal = re.sub("^.*= ", "data = ", line)
        literal = re.sub("];.*$", "]", literal)
        exec(literal)


        # Let's extrapolate the spectra index and MJD 
        page = request.get(base_url + self.id)
        soup = BeautifulSoup(page.content, "html.parser")

        mydivs = soup.find_all("div", {"class": "spectra_table scroll"})

        spec_index, spec_date = [], []
        for line in mydivs[0].decode_contents().strip().split('\n'):
            hd = line.split('<tr class="spectrum"')
            fs = line.split('<td')
            
            if len(hd)>1:
                s_indx = hd[1].split('id-spectrum')[1].split('"')[1]
                spec_index.append(int(s_indx))
            
            if len(fs)>1:
                if len(fs[1].split("-"))>2:
                    date = fs[1].split(">")[1].split("</td")[0]
                    spec_date.append(date)

        data_table = Table(data)
        return data_table.add_columns([spec_index, Time(spec_date).mjd],
         names=("spec_index", "mjd"))
        

    def query_lightcurve_alert(self):
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