#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import Angle, SkyCoord, EarthLocation
from astropy.time import Time
import astropy.units as u
from collections import OrderedDict
import os

# import base class
from core.models.level2 import RV2


# MAROONX Level2 Reader
class MAROONXRV2(RV2):
    """
    Read a MAROONX level 1 file and convert it to the EPRV standard format
    Python object.

    This class extends the `RV2` base class to handle the reading of MAROONX
    Level 1 files and converts them into a standardized EPRV format.
    Each extension from the FITS file is read, and relevant data, including
    flux, wavelength, variance, and metadata, are stored as attributes of the
    resulting Python object.

    Methods
    -------
    _read(self, file: str) -> None:
        Reads the input HDF5 file, extracts specific extensions related to the
        science data for different channels and fibers, and stores them in a
        standardized format.

        - The method processes data in different fibers (2, 3, 4, 5, 6) from
          both the BLUE and RED channel.
        - For each channel and fiber, the flux, wavelength, variance, and
          metadata are extracted and stored as a `SpectrumCollection` object.

    Attributes
    ----------
    extensions : dict
        A dictionary containing all the created extensions
        (e.g., `C1_SCI1`, `C1_SKY1`, `C2_CAL1`) where the keys are the
        extension names and the values are `SpectrumCollection` objects
        for each respective dataset.

    header : dict
        A dictionary containing metadata headers from the FITS file, with each
        extension's metadata stored under its respective key.

    Notes
    -----
    - The `_read` method processes science and calibration data from the BLUE
      and RED channels, and it extracts and organizes data for all the
      MAROON-X fibers.
    - The method converts the flux, wavelength, and variance for each extension
      into `SpectrumCollection` objects.
    - Unused extensions are removed from the object.

    Example
    -------
    >>> from core.models.level2 import RV2
    >>> file_MAROONX = 'YYYYMMDDTHHMMSSZ_FFFFF_x_TTTT.hd5'
    >>> rv2_obj = MAROONXRV2()
    >>> rv2_obj._read(file_MAROONX)
    """

    @staticmethod
    def pad_array(arr, target_shape, pad_value=np.nan):
        """
        Pad a 1D array with NaNs to match the target shape.
        If the input array is shorter than the specified target shape, NaNs
        are appended to its end.
        Parameters
        ----------
        arr : np.ndarray
            Input 1D array to pad.
        target_shape : int
            Desired length of the output array.
        pad_value: float
            Value to pad with (default: np.nan)
        Returns
        -------
        padded : np.ndarray
            Padded array of length 'target_shape'.
        """
        arr = np.atleast_2d(arr)
        current_shape = arr.shape
        padded = np.full(target_shape, pad_value, dtype=arr.dtype)
        min_rows = min(current_shape[0], target_shape[0])
        min_cols = min(current_shape[1], target_shape[1])
        padded[:min_rows, :min_cols] = arr[:min_rows, :min_cols]
        return padded

    @staticmethod
    def clean_key(header):
        """
        Standardize a metadata key to ensure compatibility with FITS standards.

        Parameters
        ----------
        header : dict
            Original header dictionary.

        Returns
        -------
        dict
            Cleaned header with standardized keys.
        """
        cleanheader = {}
        for key, value in header.items():
            k = key.replace("MAROONX", "")
            cleanheader[f"HIERARCH {k}" if len(k) > 8 else k] = value
        return cleanheader

    def standardizeMXHeader(self, prime, file, obstype, channel):
        """
        Convert MAROON-X header keywords into the standardized EPRV header
        format.

        Parameters
        ----------
        prime : dict
            Original MAROON-X header dictionary.
        file : str
            Filename of the HDF5 file being processed.
        obstype : str
            Observation type ('SCI' for science, 'CAL' for calibration).

        Returns
        -------
        phead : OrderedDict
            Dictionary containing standardized header keywords.
        """

        hmap_path = os.path.join(os.getcwd(), 'config', 'header_map.csv')
        headmap = pd.read_csv(hmap_path, header=0)

        name = file
        ra_deg = prime['MAROONX TELESCOPE TARGETRA']
        dec_deg = prime['MAROONX TELESCOPE TARGETDEC']
        ra_hms = Angle(ra_deg, unit=u.deg).to_string(unit=u.hour, sep=':',
                                                     precision=3, pad=True)
        dec_dms = Angle(dec_deg, unit=u.deg).to_string(unit=u.degree, sep=':',
                                                       precision=3, pad=True,
                                                       alwayssign=True)
        if channel.upper() == "RED":
            numorder = 28
        else:
            numorder = 34

        header_conv = {
            'OBSTYPE': obstype,
            'NUMORDER': numorder,
            'FILENAME': name,
            'CRA2': ra_hms, 'CRA3': ra_hms, 'CRA4': ra_hms,
            'CDEC2': dec_dms, 'CDEC3': dec_dms, 'CDEC4': dec_dms,
        }

        phead = OrderedDict()
        ihead = prime
        for i, row in headmap.iterrows():
            skey = row['STANDARD']
            if skey in header_conv:
                phead[skey] = header_conv[skey]
            else:
                mxkey = row['INSTRUMENT']
                if pd.notnull(mxkey):
                    mxval = ihead.get(mxkey, None)
                else:
                    mxval = row['DEFAULT']
                phead[skey] = mxval if pd.notnull(mxval) else None

        return phead

    def compute_bjd_from_header(self, hdr):
        """
        Calculates the BJD_TDB value from the MAROON-X header information.

        Parameters
        ----------
        hdr : dict
            Original header dictionary.

        Returns
        -------
        t_bjd : float
            Calculated BJD_TDB value
        """
        lon = -155.469047
        lat = 19.823801
        alt = 4213.0
        jd_utc = float(hdr["JD_UTC_FLUXWEIGHTED_FRD"])
        t = Time(jd_utc, format='jd', scale='utc')

        location = EarthLocation.from_geodetic(lat=lat*u.deg, lon=lon*u.deg,
                                               height=alt*u.m)

        ra = hdr["MAROONX TELESCOPE TARGETRA"]
        dec = hdr["MAROONX TELESCOPE TARGETDEC"]

        target = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")

        ltt_bary = t.light_travel_time(target, location=location)
        t_bjd = t.tdb + ltt_bary

        return t_bjd.value

    def create_wavehead(self, spec, fiber, channel):
        """
        Construct WCS-compliant wavelength header for a given fiber for the
        blue and red cameras.

        Parameters
        ----------
        spec : DataFrame
            Spectrum data, can be blue or red.

        fiber : int
            Fiber number for which to generate the wavelength metadata.

        channel: str
            Flag to select between blue and red camera

        Returns
        -------
        wave_meta : dict
            Dictionary of FITS-compliant wavelength metadata for the fiber.
        """

        if channel.upper() == "RED":
            CRPIX = 2018
            PS_2 = PS_4 = '[1, 4036]'
            naxis1 = 4036
            orders = spec.index.levels[1]   # Red orders 67–94
        else:
            CRPIX = 1977
            PS_2 = PS_4 = '[1, 3954]'
            naxis1 = 3954
            orders = spec.index.levels[1]   # Blue orders 91–124

        wave_meta = {}
        wave_meta["XTENSION"] = "IMAGE"
        wave_meta["BITPIX"] = 16
        wave_meta["NAXIS"] = 2
        wave_meta["NAXIS1"] = naxis1
        wave_meta["NAXIS2"] = len(orders)
        wave_meta["WCSAXES"] = len(orders)
        wave_meta["PCOUNT"] = 0
        wave_meta["GCOUNT"] = 0
        wave_meta["EXTNAME"] = f"{channel}_TRACE{fiber}_WAVE"

        i = 1
        for order in sorted(orders):
            wave = spec['wavelengths'][fiber][order]
            if wave is None:
                continue

            CRVAL = float(wave[CRPIX - 1])
            CDELT = float(wave[CRPIX] - wave[CRPIX - 1])

            wave_meta[f"CTYPE{i}"] = "WAVE"
            wave_meta[f"CPDIS{i}"] = "Non-Parametric"
            wave_meta[f"CUNIT{i}"] = "nanometer"
            wave_meta[f"CRPIX{i}"] = CRPIX
            wave_meta[f"CRVAL{i}"] = CRVAL
            wave_meta[f"CDELT{i}"] = CDELT
            wave_meta[f"PS{i}_0"] = "Spline"
            wave_meta[f"PS{i}_1"] = 0
            wave_meta[f"PS{i}_2"] = PS_2
            wave_meta[f"PS{i}_3"] = order
            wave_meta[f"PS{i}_4"] = PS_4

            i += 1

        return wave_meta

    def make_orderTable(self, spec):
        """
        Generate a combined echelle order table from spectra with
        the starting/ending wavelengths for each order.

        Parameters
        ----------
        spec : DataFrame
            Spectrum data.

        Returns
        -------
        order_table : pandas.DataFrame
            A DataFrame containing echelle order number, internal index, and
            wavelength start/end values for each order.
        """
        orders = spec.index.levels[1]

        wave_start = []
        wave_end = []

        for order in orders:
            w = spec['wavelengths'][3][order]  # fiber 3 is reference
            wave_start.append(np.nanmin(w))
            wave_end.append(np.nanmax(w))

        order_table = pd.DataFrame({
            "echelle_order": np.array(orders, dtype=np.float32),
            "order_index":   np.arange(len(orders), dtype=np.float32),
            "wave_start":    np.array(wave_start, dtype=np.float32),
            "wave_end":      np.array(wave_end, dtype=np.float32),
        })

        return order_table

    def _read(self, file: str) -> None:
        """
        Read the HDF5 file and populate two RV2 objects for red and blue:
          - self.blue_product
          - self.red_product
        """
        # Opening HDF5 file and extracting datasets stored in the file
        store = pd.HDFStore(file, 'r')
        self.spec_blue = store['spec_blue']
        self.header_blue = store['header_blue']
        self.spec_red = store['spec_red']
        self.header_red = store['header_red']
        store.close()

        # read blaze file
        bzfile = os.path.join(os.getcwd(), 
                              '20231012T05_masterflat_backgroundsubtracted_FFFFF_x_blaze.hd5')
        store2 = pd.HDFStore(bzfile, 'r+')
        blaze_blue = store2['blaze_blue']
        blaze_red = store2['blaze_red']
        store2.close()

        obstype, flux_key, fiber_range = (
            ('CAL', 'box_extraction', range(1, 6))
            if any(x in file for x in ['EEEE_', 'LLLL_', 'LLLE_'])
            else ('SCI', 'optimal_extraction', range(1, 7))
        )

        # instantiate two RV2 objects for blue and red
        self.blue_product = RV2()
        self.red_product = RV2()

        # Set up extension description table
        ext_table = {
            "extension_name": [],
            "description": [],
        }

        # Primary header according to EPRV Standard FITS format
        std_primary_blue = self.standardizeMXHeader(self.header_blue, file,
                                                    obstype, channel="BLUE")
        std_primary_red = self.standardizeMXHeader(self.header_red, file,
                                                   obstype, channel="RED")
        self.blue_product.set_header("PRIMARY", std_primary_blue)
        self.red_product.set_header("PRIMARY", std_primary_red)

        # Original Instrument Header
        self.blue_product.set_header("INSTRUMENT_HEADER",
                                     self.clean_key(self.header_blue))
        self.red_product.set_header("INSTRUMENT_HEADER",
                                    self.clean_key(self.header_red))

        # Creating Order tables (separate for blue and red)
        order_table_blue = self.make_orderTable(self.spec_blue)
        order_table_red = self.make_orderTable(self.spec_red)
        self.blue_product.create_extension("ORDER_TABLE", "BinTableHDU",
                                           data=order_table_blue)
        self.red_product.create_extension("ORDER_TABLE", "BinTableHDU",
                                          data=order_table_red)

        ext_table["extension_name"].append("PRIMARY")
        ext_table["description"].append("EPRV Standard Header")
        ext_table["extension_name"].append("INSTRUMENT_HEADER")
        ext_table["description"].append("Inherited from instrument file")
        ext_table["extension_name"].append("RECEIPT")
        ext_table["description"].append("Receipt")
        ext_table["extension_name"].append("DRP_CONFIG")
        ext_table["description"].append("DRP information")
        ext_table["extension_name"].append("EXT_DESCRIPT")
        ext_table["description"].append("Description of each extension")
        ext_table["extension_name"].append("ORDER_TABLE")
        ext_table["description"].append("Table of echelle order information")

        # Extracting flux/wavelength/variance/blaze for blue and red
        blue_data, red_data = {}, {}

        for fiber in fiber_range:
            ref_shape_blue = (34, 3954)
            ref_shape_red = (28, 4036)

            if fiber == 1:
                # Set NaN arrays for fiber 1 as sky fiber is not currently used
                blue_flux = np.full(ref_shape_blue, np.nan)
                red_flux = np.full(ref_shape_red, np.nan)
                blue_wav = np.full(ref_shape_blue, np.nan)
                red_wav = np.full(ref_shape_red, np.nan)
                blue_var = np.full(ref_shape_blue, np.nan)
                red_var = np.full(ref_shape_red, np.nan)
                blue_blz = np.ones(ref_shape_blue)
                red_blz = np.ones(ref_shape_red)
            else:
                if fiber == 5:
                    blue_flux = self.spec_blue['box_extraction'][fiber][:]
                    red_flux = self.spec_red['box_extraction'][fiber][:]
                else:
                    blue_flux = self.spec_blue[flux_key][fiber][:]
                    red_flux = self.spec_red[flux_key][fiber][:]

                blue_wav = self.spec_blue['wavelengths'][fiber][:]
                red_wav = self.spec_red['wavelengths'][fiber][:]

                blue_var = (
                    np.full(ref_shape_blue, np.nan)
                    if obstype == 'CAL'
                    else self.spec_blue['optimal_var'][fiber][:]
                )
                red_var = (
                    np.full(ref_shape_red, np.nan)
                    if obstype == 'CAL'
                    else self.spec_red['optimal_var'][fiber][:]
                )
                blue_blz = blaze_blue['blaze'][fiber][:]
                red_blz = blaze_red['blaze'][fiber][:]

            blue_data[fiber] = {
                'FLUX': np.vstack(blue_flux),
                'WAVE': np.vstack(blue_wav),
                'VAR': np.vstack(blue_var),
                'BLAZE': np.vstack(blue_blz),
            }
            red_data[fiber] = {
                'FLUX': np.vstack(red_flux),
                'WAVE': np.vstack(red_wav),
                'VAR': np.vstack(red_var),
                'BLAZE': np.vstack(red_blz),
            }

        for fiber in fiber_range:
            out_prefix = f"TRACE{fiber}_"
            blue = blue_data[fiber]
            red = red_data[fiber]

            max_cols_blue = blue['FLUX'].shape[1]
            max_cols_red = red['FLUX'].shape[1]

            blue['FLUX'] = self.pad_array(
                blue['FLUX'],
                (blue['FLUX'].shape[0], max_cols_blue),
            )
            blue['WAVE'] = self.pad_array(
                blue['WAVE'],
                (blue['WAVE'].shape[0], max_cols_blue),
            )
            blue['VAR'] = self.pad_array(
                blue['VAR'],
                (blue['VAR'].shape[0], max_cols_blue),
            )
            blue['BLAZE'] = self.pad_array(
                blue['BLAZE'],
                (blue['BLAZE'].shape[0], max_cols_blue),
            )

            red['FLUX'] = self.pad_array(
                red['FLUX'],
                (red['FLUX'].shape[0], max_cols_red),
            )
            red['WAVE'] = self.pad_array(
                red['WAVE'],
                (red['WAVE'].shape[0], max_cols_red),
            )
            red['VAR'] = self.pad_array(
                red['VAR'],
                (red['VAR'].shape[0], max_cols_red),
            )
            red['BLAZE'] = self.pad_array(
                red['BLAZE'],
                (red['BLAZE'].shape[0], max_cols_red),
            )

            for key in ['FLUX', 'WAVE', 'VAR', 'BLAZE']:
                if fiber == 1:
                    self.blue_product.set_data(out_prefix + key, blue[key])
                    self.red_product.set_data(out_prefix + key, red[key])
                else:
                    if key == 'WAVE':
                        # blue wave header:
                        wave_md_b = self.create_wavehead(self.spec_blue, fiber,
                                                         channel="BLUE")
                        wave_hdr_b = fits.Header(wave_md_b)
                        self.blue_product.create_extension(out_prefix + key,
                                                           "ImageHDU",
                                                           data=blue[key],
                                                           header=wave_hdr_b)

                        # red wave header:
                        wave_md_r = self.create_wavehead(self.spec_red, fiber,
                                                         channel="RED")
                        wave_hdr_r = fits.Header(wave_md_r)
                        self.red_product.create_extension(out_prefix + key,
                                                          "ImageHDU",
                                                          data=red[key],
                                                          header=wave_hdr_r)
                    else:
                        self.blue_product.create_extension(out_prefix + key,
                                                           "ImageHDU",
                                                           data=blue[key])
                        self.red_product.create_extension(out_prefix + key,
                                                          "ImageHDU",
                                                          data=red[key])

                ext_table["extension_name"].append(out_prefix + key)
                ext_table["description"].append(f"{key} for fiber {fiber}")

        berv_kms_b = float(self.header_blue['BERV_FLUXWEIGHTED_FRD']) / 1000.0
        berv_kms_r = float(self.header_red['BERV_FLUXWEIGHTED_FRD']) / 1000.0
        berv_z_b = berv_kms_b / 3e5
        berv_z_r = berv_kms_r / 3e5
        bjd_tdb_b = self.compute_bjd_from_header(self.header_blue)
        bjd_tdb_r = self.compute_bjd_from_header(self.header_red)

        drift_b = float(self.header_blue['Relative_Drift'].split()[0])
        drift_r = float(self.header_red['Relative_Drift'].split()[0])

        self.blue_product.set_data("BARYCORR_KMS", np.array([berv_kms_b]))
        self.red_product.set_data("BARYCORR_KMS", np.array([berv_kms_r]))
        self.blue_product.set_data("BARYCORR_Z", np.array([berv_z_b]))
        self.red_product.set_data("BARYCORR_Z", np.array([berv_z_r]))
        self.blue_product.set_data("BJD_TDB", np.array([bjd_tdb_b]))
        self.red_product.set_data("BJD_TDB", np.array([bjd_tdb_r]))
        self.blue_product.set_data("DRIFT", np.array([drift_b]))
        self.red_product.set_data("DRIFT", np.array([drift_r]))

        ext_table["extension_name"].append("BARYCORR_KMS")
        ext_table["description"].append("Barycentric correction velocity (km/s)")
        ext_table["extension_name"].append("BARYCORR_Z")
        ext_table["description"].append("Barycentric correction in redshift")
        ext_table["extension_name"].append("BJD_TDB")
        ext_table["description"].append("Flux weighted midpoint, barycentric dynamical time (JD)")
        ext_table["extension_name"].append("DRIFT")
        ext_table["description"].append("Instrument drift velocity (km/s)")

        self.blue_product.create_extension("EXT_DESCRIPT", "BinTableHDU",
                                           data=pd.DataFrame(ext_table))
        self.red_product.create_extension("EXT_DESCRIPT", "BinTableHDU",
                                          data=pd.DataFrame(ext_table))

    def write_camera_fits(self, input_filename, out_dir='.'):
        """
        Write the two channel-specific FITS files using a timestamp extracted
        from input_filename.

        Example filenames produced:
            MAROONXBLUEL2_YYYYMMDDTHHMMSS.fits
            MAROONXREDL2_YYYYMMDDTHHMMSS.fits
        """
        base = os.path.basename(input_filename)
        timestamp = base.split("_")[0].replace("Z", "")
        blue_name = os.path.join(out_dir, f"MAROONXBLUEL2_{timestamp}.fits")
        red_name = os.path.join(out_dir, f"MAROONXREDL2_{timestamp}.fits")

        self.blue_product.to_fits(blue_name)
        self.red_product.to_fits(red_name)

        return blue_name, red_name
