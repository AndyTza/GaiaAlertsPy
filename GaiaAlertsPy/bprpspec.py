from GaiaAlertsPy import alert as gaap
from scipy.interpolate import interp1d
import numpy as np
from astropy.stats import sigma_clip
from astropy.time import Time
from astropy.io import ascii
import astropy.units as u
from PyAstronomy.pyasl import instrBroadGaussFast
from scipy.interpolate import CubicSpline
import os


class GaiaSpectrumCalibrator:
    """Class to handle the calibration of Gaia BP and RP spectra."""

    def __init__(self, bp_calibration_file, rp_calibration_file):
        self.bp_cal = ascii.read(os.path.join(f"{os.getcwd()}/other/WaveCalibFile/","conv_bp.txt")) # TODO: rough calibration. Need to improve this.
        self.rp_cal = ascii.read(os.path.join(f"{os.getcwd()}/other/WaveCalibFile/","conv_bp.txt")) # TODO : rough calibration. Need to improve this.

    def pixel_to_nm(self, pixel_values, wave='bp'):
        """Convert pixel values to wavelength values using an approximate calibration.

        Parameters:
        ----------
        pixel_values : array-like
            Pixel values to be converted.
        wave : str, optional
            Spectral band to use for conversion ('bp' or 'rp'). Default is 'bp'.

        Returns:
        -------
        array-like
            Wavelength values in nanometers.
        """
        if wave == 'bp':
            calibration_data = self.bp_cal
        elif wave == 'rp':
            calibration_data = self.rp_cal
        else:
            raise ValueError("Invalid wave value. Choose 'bp' or 'rp'.")

        valid_indices = np.where((calibration_data['col1'] >= 10) & (calibration_data['col1'] <= 50))
        model = interp1d(calibration_data['col1'][valid_indices], calibration_data['col2'][valid_indices], kind='quadratic')

        return model(pixel_values[::-1]) if wave == 'bp' else model(pixel_values)

class XPStitch:
    """Class to handle the stitching of Gaia BP and RP spectra."""

    def __init__(self, bp_spectra, rp_spectra, interpolator='linear'):
        """Initialize the XPStitch object."""
        self.bp_spectra = bp_spectra
        self.rp_spectra = rp_spectra
        self.interpolator = interpolator
        self.pixel_range = np.arange(0, 60, step=1)
        self.calibrator = GaiaSpectrumCalibrator("../WaveCalibFile/conv_bp.txt", "../WaveCalibFile/conv_rp.txt")

    def stitch_spectra(self, mode='both'):
        """Stitch the BP and RP spectra together using scipy.interp1d.

        Parameters:
        ----------
        mode : str, optional
            Mode to use for stitching ('both' or 'rp_only'). Default is 'both (i.e BP and RP)

        Returns:
        -------
        tuple
            Interpolated spectrum function and the interpolated spectrum.
        """

        if mode == 'bp':
            stitched_spectra = []
            for bp_spectrum in self.bp_spectra:
                bp_pixels = self.pixel_range[(self.pixel_range >= 10) & (self.pixel_range <= 50)]
                bp_wavelengths = self.calibrator.pixel_to_nm(bp_pixels, wave='bp')

                common_bp_wavelengths = bp_wavelengths < 600

                wavelengths = bp_wavelengths[common_bp_wavelengths]
                intensities = bp_spectrum[common_bp_wavelengths]

                spectrum_function = interp1d(wavelengths, intensities, kind=self.interpolator)
                interpolated_wavelengths = np.linspace(min(wavelengths), max(wavelengths), 1000)

                stitched_spectra.append(spectrum_function(interpolated_wavelengths))

            return stitched_spectra
        elif mode == 'rp':
            stitched_spectra = []
            for rp_spectrum in self.rp_spectra:
                rp_pixels = self.pixel_range[(self.pixel_range >= 10) & (self.pixel_range <= 50)]
                rp_wavelengths = self.calibrator.pixel_to_nm(rp_pixels, wave='rp')

                common_rp_wavelengths = (rp_wavelengths > 700) & (rp_wavelengths < 1000)

                wavelengths = rp_wavelengths[common_rp_wavelengths]
                intensities = rp_spectrum[common_rp_wavelengths]

                spectrum_function = interp1d(wavelengths, intensities, kind=self.interpolator)
                interpolated_wavelengths = np.linspace(min(wavelengths), max(wavelengths), 1000)

                stitched_spectra.append(spectrum_function(interpolated_wavelengths))

            return stitched_spectra
        else:
            stitched_spectra = []
            for bp_spectrum, rp_spectrum in zip(self.bp_spectra, self.rp_spectra):
                bp_pixels = self.pixel_range[(self.pixel_range >= 10) & (self.pixel_range <= 50)]
                rp_pixels = self.pixel_range[(self.pixel_range >= 10) & (self.pixel_range <= 50)]

                bp_wavelengths = self.calibrator.pixel_to_nm(bp_pixels, wave='bp')
                rp_wavelengths = self.calibrator.pixel_to_nm(rp_pixels, wave='rp')

                common_bp_wavelengths = bp_wavelengths < 600
                common_rp_wavelengths = (rp_wavelengths > 700) & (rp_wavelengths < 1000)

                wavelengths = np.concatenate([bp_wavelengths[common_bp_wavelengths], rp_wavelengths[common_rp_wavelengths]])
                intensities = np.concatenate([bp_spectrum[common_bp_wavelengths], rp_spectrum[common_rp_wavelengths]])

                spectrum_function = interp1d(wavelengths, intensities, kind=self.interpolator)
                interpolated_wavelengths = np.linspace(min(wavelengths), max(wavelengths), 1000)

                stitched_spectra.append(spectrum_function(interpolated_wavelengths))

            return stitched_spectra