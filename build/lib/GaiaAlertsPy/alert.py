from bs4 import BeautifulSoup
import pandas as pd
import warnings
from astropy.time import Time
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clip
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

# base URL for Gaia Photometric Alerts
base_url = "https://gsaweb.ast.cam.ac.uk/alerts/alert/"

__all__ = ['query_lightcurve']

def all_sources():
    """Return astropy.Table of all Gaia Photometric Alerts to-date.

    Returns:
        astropy.Table: all Gaia Photometric Alerts to-date
    """
    
    # Try to return the ascii of the alets within the first 20 seconds if not raise an error
    try:
        return ascii.read("http://gsaweb.ast.cam.ac.uk/alerts/alerts.csv")
    except:
        # if the query takes too long, timeout and say that the Gaia alerts server might be temporarily down
        raise ValueError("Sorry, the Gaia alerts server might be temporarily down. Please try again later.")
    

class GaiaAlertsTable:
    # TODO: wip
    """Class to handle Gaia Photometric Alerts table."""
    def __init__(self, ra, dec):
        """
        Args:
            ra (float): RA in degs
            dec (float): DEC in degs
        """
        self.ra = ra
        self.dec = dec
        
    def cone_search(self, sep=0.1):
        """Return astropy.Table of all Gaia Photometric Alerts within a cone search of radius separation.

        Parameters:
            sep (float): separation in arcseconds
        
        Returns:
            astropy.Table: all Gaia Photometric Alerts within a cone search of radius.    
        """
        master = all_sources() # load all sources
        
        coord_target = SkyCoord(ra=self.ra*u.deg, dec=self.dec*u.deg, frame='icrs')
        coords = SkyCoord(ra=master['RaDeg']*u.deg,
                         dec=master['DecDeg']*u.deg, frame='icrs')
        
        sep_all = coords.separation(coord_target).arcsec
        
        return master[np.where(sep_all<=sep)]
                

class GaiaAlert:
    """Class to handle Gaia Photometric Alerts."""
    def __init__(self, id):
        self.id = id

    def gaia_g_noise_esitmate(self, mag):
        """Compute the Gaia G-band noise estimate.
        
        Parameters:
            mag (float): Gaia G-band magnitude.
        
        Returns:
            float: Gaia G-band uncertainty.

        Reference:
            Hodgkin et al. 2021 (https://arxiv.org/abs/2106.00479). Valid for 13<Gmag<21.
        """

        return 3.43779 - (mag/1.13759) + (mag/3.44123)**2 - (mag/6.51996)**3 + (mag/11.45922)**4

    def query_bprp_history(self):
        """ Query BP and RP spectra for each epochal alert.

        Parameters:
        ----------
            id (str): Gaia alert ID.

        Returns:
        -------
            astropy.Table: BP and RP spectra for each epochal alert.

        Notes:
        ------
        This function was taken directly from (Sipőcz & Hogg):
            https://github.com/davidwhogg/GaiaAlerts/blob/master/scripts/scrape_spectra.py
            IMPORTANT: If you use this feature, please ensure to cite their initial work. 
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
        """Query the lightcurve of a Gaia alert.
        
        Parameters:
        ----------
            id (str): Gaia alert ID.
        
        Returns:
        -------
            astropy.Table: lightcurve of a Gaia alert.
        """
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

    def query_bprp_mags(self, sigma_clip=5):
        """ Query the BP and RP magnitudes of a Gaia alrts and does a quick estimate ont he BP and RP magnitudes.

        Parameters:
        ----------
            sigma_clip (int): sigma clipping value
        
        Returns:
        -------
            astropy.Table: BP and RP magnitudes of a Gaia alert.
        """
        # Fetch BPRP information table
        color_lc = self.query_bprp_history()

        # Instrumental zero-point values
        zp_BP, zp_RP = 25.3514, 24.7619 # Table 5.2 (https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_calibr_extern.html)

        bp_mag, rp_mag = [], []
        for _lc in color_lc:
            bp0, rp0 = _lc['bp'], _lc['rp']

            # Count only positive ADU counts & apply 5-sigma clip (see Hodgkin et al. 2021; Section 3.6)
            bp, rp = sigma_clip(bp0[bp0>0], sigma=sigma_clip), sigma_clip(rp0[rp0>0], sigma=sigma_clip)

            # Convert total flux to magnitudes
            bp_mag.append(-2.5*np.log10(bp.sum()) + zp_BP)
            rp_mag.append(-2.5*np.log10(rp.sum()) + zp_RP)

        return Table([color_lc['JD'], bp_mag, rp_mag], names=('JD', 'bp_mag', 'rp_mag'))

def GaiaX_history(year='2021'):
    """Query the GaiaX alert history.

    Parameters:
    ----------
        year (str): year of the GaiaX alert history. Default is '2021'.

    Returns:
    -------
        astropy.Table: GaiaX alert history.
    
    Notes:
    ------
    These are all the alerts raised by the single field-of-view detector (1-FoV detector, see description in Kostrzewa-Rutkowska Z., et al., 2020).
    """
    base_url = "https://gsaweb.ast.cam.ac.uk/alerts/gaiax/"

    return ascii.read(f"{base_url}/{year}")


def GaiaX_alert(alert_id):
    """Fetch the GaiaX alert information.

    Parameters:
    ----------
        alert_id (str): GaiaX alert ID.
    
    Returns:
    -------
        astropy.Table: GaiaX alert information.
    """
    base_url = "https://gsaweb.ast.cam.ac.uk/alerts/gaiax/"
    try:
        return ascii.read(f"{base_url}/{alert_id}.csv")
    except ValueError:
        raise ValueError("Sorry, the GaiaX alert ID you queried was not found.")

