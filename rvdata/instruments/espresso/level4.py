"""
RVData/rvdata/instruments/espresso/level4.py

UNIGE-ESO - EPRV
Author: Loris JACQUES & Emile FONTANET
Created: Wed Aug 20 2025
Version: 1.0.0

---------------------
Libraries
---------------------
"""

from astropy.io import fits
import os
import pandas as pd
from rvdata.core.models.level4 import RV4
import rvdata.instruments.espresso.config.config as config
from rvdata.instruments.espresso.utils import (
    get_files_names,
    create_PRIMARY,
    validate_fits_file,
    convert_CCF
)

# ESPRESSO Level4 Reader


class ESPRESSORV4(RV4):
    """
    Read ESPRESSO Level 1 and Level 2 files and convert them into the EPRV
    standard format.

    This class extends the `RV2` base class to handle the reading of ESPRESSO
    (Echelle SPectrograph for Rocky Exoplanets and Stable Spectroscopic
    Observations) Level 1 and Level 2 files, combining information from both
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
            - Converts the spectral blaze functions (``S2D_BLAZE_A``,
            ``S2D_BLAZE_B``, ``BLAZE_A``, ``BLAZE_B``) for the different fibers.
            - Processes and converts the drift file for instrumental calibration.
            - Creates the ``PRIMARY`` header and necessary metadata.
            - For now: Removes unused or redundant extensions such as ``RECEIPT``
            and ``DRP_CONFIG``.

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
        - Blaze correction functions are processed and stored for each fiber.
        - The drift file is processed separately for calibration.
        - Unused extensions (like ``RECEIPT`` and ``DRP_CONFIG``) are removed from the
        final output.

    Example
    -------
        >>> from core.models.level4 import RV4
        >>> rv4_obj = ESPRESSORV4.from_fits("espresso_raw_file.fits")
        >>> rv4_obj.to_fits("standard_level4.fits")

    """

    def do_conversion(
        self, hdul: fits.HDUList, directory_structure: str = "standard"
    ) -> None:
        """
        Converts FITS files based on certain conditions and configurations.

        This method performs several processing steps:

        1. Validates the FITS file structure before conversion.
        2. Retrieves paths for required additional files (e.g., blaze
           functions, drift corrections).
        3. Converts and stores spectral blaze functions for different fibers.
        4. Converts the drift correction data.
        5. Creates the ``PRIMARY`` header and integrates necessary metadata.
        6. Cleans up unused extensions like ``RECEIPT`` and ``DRP_CONFIG``.

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
            print("File is valid for conversion!")
        except ValueError as e:
            raise ValueError(e)

        # Retrieve the paths for the necessary files
        names = get_files_names(path, directory_structure, level=4)
        print(names)
        trace_ind_start = 1
        with fits.open(names['ccf_file']) as hdul_ccf:
            self.set_header("INSTRUMENT_HEADER", hdul_ccf["PRIMARY"].header)
        # -> TODO: Do we want to also store the telcorr CCF?
        # try:
        #     convert_TELLURIC(
        #         self,
        #         names["telluric_file_" + fiber],
        #         trace_ind_start,
        #         config.slice_nb,
        #     )
        #     print("TRACEi_TELLURIC_x extensions " "have been generated.")
        # except Exception:
        #     print(
        #         "No TELLURIC file found, TRACEi_TELLURIC_x extensions "
        #         "will not be generated."
        #     )

        

        
        trace_ind_start += 2
        # Create the PRIMARY header
        nb_trace = config.slice_nb
        create_PRIMARY(self, names, nb_trace, config.slice_nb, level = 4)
        convert_CCF(self, names)
        return

        
    