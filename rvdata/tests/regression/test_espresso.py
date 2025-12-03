import requests
import os
from rvdata.core.models.level2 import RV2
from rvdata.core.models.level4 import RV4

# from rvdata.instruments.kpf.level2 import KPFRV2
from rvdata.tests.regression.compliance import check_l2_extensions, check_l2_header
from rvdata.tests.regression.compliance import check_l4_extensions, check_l4_header


file_urls = {
    "ESPRESSO": {
        'raw': 'https://dace.unige.ch/downloads/EPRV_standard_data/ESPRE.2017-12-03T02:09:40.348.fits',
        'S2D_BLAZE_A': 'https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T02:09:40.348_S2D_BLAZE_A.fits',
        'S2D_BLAZE_B': 'https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T02:09:40.348_S2D_BLAZE_B.fits',
        'BLAZE_A': 'https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T10:43:59.835_BLAZE_A.fits',
        'BLAZE_B': 'https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T10:43:59.835_BLAZE_B.fits',
        'S1D_A': 'https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T02:09:40.348_S1D_A.fits',
        'S1D_B': 'https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T02:09:40.348_S1D_B.fits',
        'S1D_TELL_CORR_A': 'https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T02:09:40.348_S1D_TELL_CORR_A.fits',
        'DRIFT_MATRIX_B': 'https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T02:09:40.348_DRIFT_MATRIX_B.fits',
        'CCF_A': 'https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T02:09:40.348_CCF_A.fits',
        'CCF_TELL_CORR_A': 'https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T02:09:40.348_CCF_TELL_CORR_A.fits'
    }
}


def download_file(url, filename):
    response = requests.get(url, verify=False)
    response.raise_for_status()  # Check if the request was successful
    with open(filename, "wb") as file:
        file.write(response.content)


def download_files():
    raw_file = "ESPRE.2017-12-03T02:09:40.348.fits"
    s2d_A_file = "r.ESPRE.2017-12-03T02:09:40.348_S2D_BLAZE_A.fits"
    s2d_B_file = "r.ESPRE.2017-12-03T02:09:40.348_S2D_BLAZE_B.fits"
    blaze_A_file = "r.ESPRE.2017-12-03T10:43:59.835_BLAZE_A.fits"
    blaze_B_file = "r.ESPRE.2017-12-03T10:43:59.835_BLAZE_B.fits"
    s1d_A_file = "r.ESPRE.2017-12-03T02:09:40.348_S1D_A.fits"
    s1d_B_file = "r.ESPRE.2017-12-03T02:09:40.348_S1D_B.fits"
    s1d_tell_corr_A_file = "r.ESPRE.2017-12-03T02:09:40.348_S1D_TELL_CORR_A.fits"
    drift_matrix_B_file = "r.ESPRE.2017-12-03T02:09:40.348_DRIFT_MATRIX_B.fits"
    ccf_A_file = "r.ESPRE.2017-12-03T02:09:40.348_CCF_A.fits"
    ccf_tell_corr_A_file = "r.ESPRE.2017-12-03T02:09:40.348_CCF_TELL_CORR_A.fits"

    if not os.path.exists(raw_file):
        download_file(file_urls["ESPRESSO"]['raw'], raw_file)
    if not os.path.exists(s2d_A_file):
        download_file(file_urls["ESPRESSO"]['S2D_BLAZE_A'], s2d_A_file)
    if not os.path.exists(s2d_B_file):
        download_file(file_urls["ESPRESSO"]['S2D_BLAZE_B'], s2d_B_file)
    if not os.path.exists(blaze_A_file):
        download_file(file_urls["ESPRESSO"]['BLAZE_A'], blaze_A_file)
    if not os.path.exists(blaze_B_file):
        download_file(file_urls["ESPRESSO"]['BLAZE_B'], blaze_B_file)
    if not os.path.exists(s1d_A_file):
        download_file(file_urls["ESPRESSO"]['S1D_A'], s1d_A_file)
    if not os.path.exists(s1d_B_file):
        download_file(file_urls["ESPRESSO"]['S1D_B'], s1d_B_file)
    if not os.path.exists(s1d_tell_corr_A_file):
        download_file(file_urls["ESPRESSO"]['S1D_TELL_CORR_A'], s1d_tell_corr_A_file)
    if not os.path.exists(drift_matrix_B_file):
        download_file(file_urls["ESPRESSO"]['DRIFT_MATRIX_B'], drift_matrix_B_file)
    if not os.path.exists(ccf_A_file):
        download_file(file_urls["ESPRESSO"]['CCF_A'], ccf_A_file)
    if not os.path.exists(ccf_tell_corr_A_file):
        download_file(file_urls["ESPRESSO"]['CCF_TELL_CORR_A'], ccf_tell_corr_A_file)

    return raw_file, s2d_A_file, s2d_B_file, blaze_A_file, blaze_B_file, s1d_A_file, s1d_B_file, s1d_tell_corr_A_file, drift_matrix_B_file, ccf_A_file, ccf_tell_corr_A_file


def test_espresso():
    raw_file, s2d_A_file, s2d_B_file, blaze_A_file, blaze_B_file, s1d_A_file, s1d_B_file, s1d_tell_corr_A_file, drift_matrix_B_file, ccf_A_file, ccf_tell_corr_A_file = download_files()
    espr2 = RV2.from_fits(raw_file, instrument="ESPRESSO")
    l2_standard = "./esp_L2_standard.fits"
    espr2.to_fits(l2_standard)
    l2_obj = RV2.from_fits(l2_standard)

    check_l2_extensions(l2_standard)
    check_l2_header(l2_obj.headers['PRIMARY'])

    espr4 = RV4.from_fits(raw_file, instrument="ESPRESSO")
    l4_standard = "./espr_L4_standard.fits"
    espr4.to_fits(l4_standard)
    l4_obj = RV4.from_fits(l4_standard)

    check_l4_extensions(l4_standard)
    check_l4_header(l4_obj.headers['PRIMARY'])


def test_espresso_benchmark(benchmark):
    # run test_kpf() once to download the files
    _ = download_files()
    # now run it again with benchmark
    benchmark(test_espresso)


if __name__ == "__main__":
    test_espresso()
