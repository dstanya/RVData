from astropy.io import fits
import pandas as pd

with fits.open("scrap_l2.fits") as hdul:
    # print HDU list
    print("HDU List:")
    hdul.info()
    print()

    # PRIMARY
    name = "PRIMARY"
    hdu = hdul[name]
    print(name)
    for key, value in hdu.header.items():
        print(f"{key}: {value}")
    print()

    # RECEIPT
    name = "RECEIPT"
    hdu = hdul[name]
    print(name)
    print(pd.DataFrame(hdu.data))
    print()


    # DRP_CONFIG
    name = "DRP_CONFIG"
    hdu = hdul[name]
    print(name)
    print(pd.DataFrame(hdu.data))
    print()

    # EXT_DESCRIPT
    name = "EXT_DESCRIPT"
    hdu = hdul[name]
    print(name)
    print(pd.DataFrame(hdu.data))
    print()

    # ORDER_TABLE
    name = "ORDER_TABLE"
    hdu = hdul[name]
    print(name)
    print(pd.DataFrame(hdu.data))
    print()
