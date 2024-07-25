from GaiaAlertsPy import alert as gaap
import matplotlib
from scipy.interpolate import interp1d
cm = matplotlib.cm.get_cmap("gnuplot")
import numpy as np
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt
from matplotlib import rcParams
from astropy.time import Time
from astropy.io import ascii
import astropy.units as u
from expecto import get_spectrum
from tqdm import tqdm
from PyAstronomy.pyasl import instrBroadGaussFast
from scipy.interpolate import CubicSpline
import extinction

# Define the wavelength-pixel calibration files
file_path = "../WaveCalibFile/"
bp_cal = ascii.read(file_path + "conv_bp.txt")
rp_cal = ascii.read(file_path + "conv_rp.txt")



def pixel_2_nm(px_val, wave='bp'):
    """Convert from gaia pixels to wavelength solution...aproximate"""
    w_bp = np.where((bp_cal['col1']>=10) & (bp_cal['col1']<=50))
    w_rp = np.where((rp_cal['col1']>=10) & (rp_cal['col1']<=50))

    model_bp = interp1d(bp_cal['col1'].data[w_bp], bp_cal['col2'].data[w_bp], kind='quadratic')
    model_rp = interp1d(rp_cal['col1'].data[w_rp], rp_cal['col2'].data[w_rp], kind='quadratic')
    
    if wave=='bp':
        return model_bp(px_val[::-1])
    elif wave=='rp':
        return model_rp(px_val)

def spec_BPRP(bp_rp_spec, interpolation_style='quadratic', size=1_000):

    _bp, _rp = bp_rp_spec
    
    # Spectrum stiching
    pixel = np.arange(0, 60, step=1) # pixel size 
    xc = np.where((pixel>=10) & (pixel<=50)) 
    
    x1, y1 = pixel_2_nm(pixel[xc], wave='bp'), _bp[xc]
    x1c = x1<600

    x2, y2 = pixel_2_nm(pixel[xc], wave='rp'), _rp[xc]
    x2c = np.where((x2>700) & (x2<1000))
    
    X, Y = np.concatenate([x1[x1c], x2[x2c]]), np.concatenate([y1[x1c], y2[x2c]])

    # Interpolate
    sf = interp1d(X, Y, kind=interpolation_style)
    sfx = np.linspace(min(X), max(X), size)

    return sf, sf(sfx)
