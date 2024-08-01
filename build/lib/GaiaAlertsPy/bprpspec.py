from scipy.interpolate import interp1d
import numpy as np
from astropy.io import ascii
import os

class GaiaSpectrumCalibrator:
    """Class to handle the calibration of Gaia BP and RP spectra."""

    def __init__(self):
        calibration_dir = f"{os.getcwd()}/other/WaveCalibFile/"
        self.bp_cal = ascii.read(os.path.join(calibration_dir, "conv_bp.txt")) # TODO: rough calibration. Need to improve this.
        self.rp_cal = ascii.read(os.path.join(calibration_dir, "conv_rp.txt")) # TODO : rough calibration. Need to improve this.

    def pixel_to_nm(self, px_val, wave='bp'):
        """Using an approximate calibration of the Gaia BP and RP spectra, convert pixel values to wavelength values.
        
        Parameters:
        ----------
            px_val (array): pixel values
            wave (str): 'bp' or 'rp'. Default is 'bp'

        Returns:?
        -------
            array: wavelength values in nm.
        """
        if wave == 'bp':
            w_bp = np.where((self.bp_cal['col1'] >= 10) & (self.bp_cal['col1'] <= 50))
            model_bp = interp1d(self.bp_cal['col1'].data[w_bp], self.bp_cal['col2'].data[w_bp], kind='quadratic')
            return model_bp(px_val[::-1])
        elif wave == 'rp':
            w_rp = np.where((self.rp_cal['col1'] >= 10) & (self.rp_cal['col1'] <= 50))
            model_rp = interp1d(self.rp_cal['col1'].data[w_rp], self.rp_cal['col2'].data[w_rp], kind='quadratic')
            return model_rp(px_val)
        else:
            raise ValueError("Invalid wave value. Choose 'bp' or 'rp'.")

class XPStitch:
    """Class to handle the stitching of Gaia BP and RP spectra."""

    def __init__(self, bp_spectra, rp_spectra, interpolation_style='linear'):
        """Initialize the XPStitch object."""
        self.bp_spectra = bp_spectra
        self.rp_spectra = rp_spectra
        self.interpolation_style = interpolation_style
        self.pixel_range = np.arange(0, 60, step=1)
        self.calibrator = GaiaSpectrumCalibrator()

    def stitch_spectra(self, mode='bp'):
        """Stitch the BP and RP spectra together using scipy interp1d.
        
        Parameters:
        ----------
            bp_rp_spec (tuple): BP and RP spectra
            interpolation_style (str): interpolation style. Default is 'quadratic'
        
        Returns:
        -------
            tuple: (spectrum, interpolated spectrum)
        """
        if mode=='both':
            _bp, _rp = self.bp_spectra, self.rp_spectra
        
            # Spectrum stitching
            pixel = self.pixel_range
            xc = np.where((pixel >= 10) & (pixel <= 50))

            x1, y1 = self.calibrator.pixel_to_nm(pixel[xc], wave='bp'), _bp[xc]
            x1c = x1 < 600  # 600 nm is the common wavelength

            x2, y2 = self.calibrator.pixel_to_nm(pixel[xc], wave='rp'), _rp[xc]
            x2c = (x2 > 700) & (x2 < 1000)  # 700 nm is the common wavelength for RP cutoff

            X, Y = np.concatenate([x1[x1c], x2[x2c]]), np.concatenate([y1[x1c], y2[x2c]])

            sf = interp1d(X, Y, kind=self.interpolation_style)
            sfx = np.linspace(min(X), max(X), 60)

            return sfx, sf(sfx)
        
        elif mode=='bp':
            _bp = self.rp_spectra

            # Spectrum stitching
            pixel = self.pixel_range
            xc = np.where((pixel >= 10) & (pixel <= 50))
            sf = interp1d(self.calibrator.pixel_to_nm(pixel[xc], wave='bp'), _bp[xc], kind=self.interpolation_style)
            sfx = np.linspace(350, 600, 60)

            return sfx, sf(sfx)
        elif mode=='rp':
            _rp = self.rp_spectra

            # Spectrum stitching
            pixel = self.pixel_range
            xc = np.where((pixel >= 70) & (pixel <= 110))
            sf = interp1d(self.calibrator.pixel_to_nm(pixel[xc], wave='rp'), _rp[xc], kind=self.interpolation_style)
            sfx = np.linspace(600+1, 998, 60)

            return sfx, sf(sfx)