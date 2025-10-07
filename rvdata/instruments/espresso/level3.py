"""
RVData/rvdata/instruments/espresso/level3.py

UNIGE-ESO - EPRV
Author: Loris JACQUES & Emile FONTANET
Version: 1.0.0

---------------------
Libraries
---------------------
"""

from astropy.io import fits
import os
import pandas as pd
from rvdata.core.models.level3 import RV3
import rvdata.instruments.espresso.config.config as config
from rvdata.instruments.espresso.utils import (
    get_files_names,
    create_PRIMARY,
    validate_fits_file,
    convert_S1D,
)

class ESPRESSORV3(RV3):
    """
    Read ESPRESSO raw and S1D files and convert them into the EPRV
    standard format.

    This class extends the `RV3` base class to handle the reading of ESPRESSO
    (Echelle SPectrograph for Rocky Exoplanets and Stable Spectroscopic
    Observations) raw and S1D files, combining information from both
    sources to produce a standardized EPRV output. It processes various FITS
    extensions and organizes flux, wavelength, variance, and metadata into a
    structured Python object.

    Methods
    -------
        do_conversion(hdul: fits.HDUList) -> None
            Reads the input FITS HDU list, extracts specific extensions related to
            the science data for different chips and fibers, and stores them in a
            standardized format.
            - The method validates the FITS file before conversion to ensure it
            meets the required criteria.
            - Retrieves necessary file paths for additional files required in the
            processing.
            - Creates the ``PRIMARY`` header and necessary metadata.


    Attributes
    ----------
        extensions : dict
            A dictionary containing all the created extensions, where the keys are
            extension names, and the values are the respective data arrays.

        header : dict
            A dictionary containing metadata headers from the FITS files, with
            each extension's metadata stored under its respective key.

    Notes
    -----
        - The ``do_conversion`` method processes and extracts science and
        calibration data.
        - The method ensures the FITS file meets the required criteria before
        conversion.


    Example
    -------
        >>> from core.models.level3 import RV3
        >>> rv3_obj = ESPRESSORV3.from_fits("espresso_raw_file.fits")
        >>> rv3_obj.to_fits("standard_level3.fits")

    """

    def do_conversion(
        self, hdul: fits.HDUList, directory_structure: str = "standard"
    ) -> None:
        """
        Converts FITS files based on certain conditions and configurations.

        This method performs several processing steps:

        1. Validates the FITS file structure before conversion.
        2. Retrieves paths for required additional files.
        3. Converts and stores S1D data.
        4. Creates the ``PRIMARY`` header and integrates necessary metadata.

        Parameters
        ----------
        hdul : fits.HDUList
            The FITS HDU list to be processed.
        directory_structure : str
            Type of database architecture that stores resources. Must be either
            'dace' or 'standard'.

        Raises
        ------
        ValueError
            If the FITS file is invalid and does not meet the required criteria
            for conversion.

        :noindex:
        """
        path = os.path.join(self.dirname, self.filename)
        # Validate the FITS file before conversion. If it does not meet the
        # criteria, raise an error
        try:
            validate_fits_file(path)
        except ValueError as e:
            raise ValueError(e)

        # Retrieve the paths for the necessary files
        names = get_files_names(path, directory_structure, level=3)
        with fits.open(names['s1d_A']) as hdul_ccf:
            self.set_header("INSTRUMENT_HEADER", hdul_ccf["PRIMARY"].header)
        
        nb_trace = 2
        create_PRIMARY(self, names, nb_trace, 1, level = 3)
        convert_S1D(self, names)
        return