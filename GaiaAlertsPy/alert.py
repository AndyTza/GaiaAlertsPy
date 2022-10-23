from bs4 import BeautifulSoup
import pandas as pd
import warnings
from astropy.time import Time
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib 
import re
import astropy.utils.data as aud
import requests
import os
import ast
import bs4
import sys
import warnings
from astropy.table import Table, hstack, vstack


base_url = "https://gsaweb.ast.cam.ac.uk/alerts/alert/"

__all__ = [ 'query_lightcurve']


def all_sources():
    """Return astropy.Table of all Gaia Photometric Alerts to-date.

    Returns:
        astropy.Table: all Gaia Photometric Alerts to-date
    """
    return ascii.read("http://gsaweb.ast.cam.ac.uk/alerts/alerts.csv")
    

class GaiaAlertsTable:
    def __init__(self, ra, dec):
        """
        Args:
            ra (float): RA in degs
            dec (float): DEC in degs
        """
        self.ra = ra
        self.dec = dec
        
    def cone_search(self, sep=0.1):
        """Cone search the Gaia Photometric alerts table. Return the closest crossmatch. 

        Args:
            sep (float, optional): Separation in arcseconds. Defaults to 0.1.

        Returns:
            astropy.Table: cone-searched Gaia Photometric Alerts at positio
        """
        master = all_sources()
        
        coord_target = SkyCoord(ra=self.ra*u.deg, dec=self.dec*u.deg, frame='icrs')
        coords = SkyCoord(ra=master['RaDeg']*u.deg,
                         dec=master['DecDeg']*u.deg, frame='icrs')
        
        sep_all = coords.separation(coord_target).arcsec
        
        return master[np.where(sep_all<=sep)]
                

class GaiaAlert:
    def __init__(self, id):
        self.id = id

    def gaia_g_noise_esitmate(self, mag):
        """TODO: Source & logistics"""
        return 3.43779 - (mag/1.13759) + (mag/3.44123)**2 - (mag/6.51996)**3 + (mag/11.45922)**4

    def query_bprp_history(self):
        """ Query BP and RP spectra for each epochal alert.  TODO: Finish docs!

        This function was taken directly from (Sipőcz & Hogg):
            https://github.com/davidwhogg/GaiaAlerts/blob/master/scripts/scrape_spectra.py
            IMPORTANT: If you use this feature, please ensure to cite their initial work. 

        Args:
            name (_type_): _description_

        Returns:
            _type_: _description_
        """

        content = aud.get_file_contents("{}/{}".format(base_url, self.id), cache=True)
        htmldoc = bs4.BeautifulSoup(content, 'html5lib')

        search_text = re.compile('var spectra')
        line = htmldoc.find('script', text=search_text)

        try:
            spectra_data = line.string.split('=')[1].strip()
        except AttributeError:
            warnings.warn("data is not found for {}, check whether it's a "
                        "valid GaiaAlerts object".format(self.id))
            return None

        spectra_data = Table(ast.literal_eval(spectra_data[1:-2]))

        spectra_meta = Table.read(content, format='ascii.html')
        spectra = hstack([spectra_meta, spectra_data], join_type='outer')

        spectra['Name'] = self.id

        return spectra 
        

    def query_lightcurve_alert(self):
        try:
            _dat = pd.read_csv(f"{base_url}{self.id}/lightcurve.csv/", skiprows=1)
        except:
            raise ValueError("Sorry, the Gaia alert ID you queried was not found.") 

        master_frame = np.zeros(shape=(_dat.shape[0], 3))
        
        for i, vals in enumerate(_dat.to_numpy()):
            time_conv = Time(vals[0], format='iso', scale='tcb').jd
            master_frame[i,0] = time_conv

            if vals[-1]=='untrusted' or vals[-1]=='nan':
                master_frame[i,1] = np.nan
            else:
                master_frame[i, 1] = float(vals[-1])
        
        indx_nans = ~np.isnan(master_frame[:,1])

        return Table([master_frame[:,0][indx_nans], 
                    master_frame[:,1][indx_nans], 
                    self.gaia_g_noise_esitmate(master_frame[:,1][indx_nans])], 
                    names=("JD", "mag_G", "mag_G_error")) 


    def query_bprp_mags(self):
        """pseudo- BP/RP magnitudes
        """
        # Fetch BPRP information table
        color_lc = self.query_bprp_history()

        # Take the sum of the ADU 
        bp_adu, rp_adu = [rp_.sum() for rp_ in color_lc['rp']], [bp_.sum() for bp_ in color_lc['bp']]
        zp_bp, zp_rp = 25.3514, 24.7619 # Table 5.2 (https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_calibr_extern.html)

        bp_mag, rp_mag = -2.5*np.log10(bp_adu) + zp_bp , -2.5*np.log10(rp_adu) + zp_rp

        return Table([color_lc['JD'], bp_mag, rp_mag], names=('JD', 'bp_mag', 'rp_mag'))





                    