#    if len(cuthdul) != 1:
#        raise AssertionError('something wrong in tesscut! FIXME')
#    else:
#        cuthdul = cuthdul[0]
#    data, data_hdr = cuthdul[1].data, cuthdul[1].header
#    cutout_wcs = wcs.WCS(cuthdul[2].header)
#
#    # flux and flux_err are image cubes of (time x spatial x spatial)
#    quality = data['QUALITY']
#
#    flux = data['FLUX']
#    flux_err = data['FLUX_ERR']
#    time = data['TIME'] # in TJD
#    time += 2457000 # now in BJD
#
#    time = time[quality == 0]
#    flux = flux[quality == 0]
#    flux_err = flux_err[quality == 0]
