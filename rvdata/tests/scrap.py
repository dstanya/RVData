from rvdata.core.models.level2 import RV2

# TODO make this a CI test
slv2 = RV2()
slv2.to_fits("scrap_l2.fits")

in_slv2 = RV2.from_fits("scrap_l2.fits")

# iterate through and print all PRIMARY header keywords
for key, value in in_slv2.data["PRIMARY"].items():
    print(f"{key}: {value}")