from astropy.io import fits
import numpy as np
import os
import pandas as pd


def get_neid_instrument_era(obs_jd):
    """Function to determine which NEID instrument era a given
       observation falls under.

    Parameters
    ----------
    obs_jd : float
        The JD of the observation.

    Returns
    -------
    era : int
        The observation's NEID instrument RV era.
    """

    # Read in the map of dates to NEID instrument RV eras
    eramap = pd.read_csv(
        os.path.join(
            os.path.dirname(os.path.dirname(__file__)), "config/neid_inst_eras.csv"
        )
    )

    # Get the era based on the given observation time
    era_time_diffs = obs_jd - eramap["startdate"].values
    era = eramap["era"].values[
        np.argmin(era_time_diffs[np.where(era_time_diffs >= 0)[0]])
    ]

    return era
