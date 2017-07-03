import os
import fnmatch
import time
import math
import sys
import glob
import datetime

import numpy as np
from numpy.lib.recfunctions import append_fields
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import astropy
import astropy.io.fits as fits
from astropy.table import Table
from astropy.time import Time
from astropy import wcs

import atpy
import fitsio
from fitsio import FITS, FITSHDR
sys.path.append('/home/sr525/Python_Code/')
import srpylib as PA
import copy

def file_supercede(dir1, dir2):

    files1 = os.listdir(dir1)
    files2 = os.listdir(dir2)

    print "Number of files in first directory:", len(files1)
    print "Number of files in second directory:", len(files2)

    if len(files1) <> len(files2):
        if len(files1) > len(files2):
            small_files = files2
            long_files = files1
        else:
            small_files = files1
            long_files = files2

    print "Files not in both:", list(set(files1)^set(files2))
    #print set(files1) - set(files2)
    #print set(files2) - set(files1)

    files_both = list(set(files1) & set(files2))

    same = 0
    diff_files = []
    for file in files_both:
        if file[-5:] == ".fits":
            fh1 = fits.open(dir1 + file)
            fh2 = fits.open(dir2 + file)
            if fh1[1].header["DATASUM"] <> fh2[1].header["DATASUM"]:
                diff_files.append(file)
            else:
                same += 1
    print "Number of files that have the same datasum:", same, "out of", len(files1), "/", len(files2)

    return diff_files

def check_prodcatg(path):
    """

    """
    files = fnmatch.filter(os.listdir(path), "*.fits")
    for file in files:
        with fits.open(path + file) as hlist:
            print(path + file)
            print(hlist[0].header["PRODCATG"])

def poly_Area(vertices):
    return np.sum( [-0.5, 0.5] * vertices * np.roll( np.roll(vertices, 1, axis=0), 1, axis=1) )

def rd2sn(ra, dec, racen, deccen):
    beta = ra - racen
    s = np.cos(beta) + np.tan(dec)*np.tan(deccen)
    xi  = np.sin(beta)/np.cos(deccen)/s
    eta = (np.tan(dec) - np.tan(deccen)*np.cos(beta))/s

    return xi, eta

def area(data, racen, deccen, nbins = 1000):

    #xs = np.rad2deg(data["ra"]*np.cos(data["dec"]))
    #ys = np.rad2deg(data["dec"])
    xs = (data["RA"])
    ys = (data["DEC"])

    xs, ys = rd2sn(xs, ys, racen, deccen)
    #xs = np.rad2deg(xs)
    #ys = np.rad2deg(ys)
    plt.plot(xs, ys, "k.")
    #plt.show()

    """
    xwalls = np.linspace( min(xs) - 5.0, max(xs) + 5.0, nbins + 1 )
    ywalls = np.linspace( min(ys) - 5.0, max(ys) + 5.0, nbins + 1 )
    im, xs_bin, ys_bin, ax = plt.hist2d(xs, ys, bins = (xwalls, ywalls) )
    xs_mids = 0.5*(xs_bin[:-1] + xs_bin[1:])
    ys_mids = 0.5*(ys_bin[:-1] + ys_bin[1:])

    plt.close()
    im[im>0] = 1
    conts = plt.contour(xs_mids, ys_mids, im.T, 1)
    plt.plot(xs, ys, "k.")
    plt.show()
    plt.close()
    paths = conts.collections[0]
    paths = paths.get_paths()
    area = poly_Area(paths[0].vertices)
    print area*np.cos(deccen)

    """
    #Uniform grid in RA and Dec
    n = min(xs) - 0.5
    xs1 = []
    ys1 = []

    y = (np.linspace(min(ys)-0.05, max(ys)+0.05, 1000.0))
    x = np.linspace(min(xs)-0.05, max(xs)+0.05, 1000.0)
    y = (y)
    xs1, ys1 = np.meshgrid(x, y)
    xs1 = xs1.flatten()
    ys1 = ys1.flatten()
    area1 = (max(xs)-min(xs)+0.1)*(max(ys)-min(ys)+0.1)

    inside = PA.inside_contour((xs), (ys), xs1, ys1,
                                   c_graph = False, radius = -0.0024,
                                   nbins = 5000)
    area = float(len(xs1[inside]))/float(len(xs1))*area1*(180/np.pi)**2

    return area

def mjd(mdatafile, type = "float"):

    fh = fits.open(mdatafile)
    dqc = fh[1].data
    mjdobs = np.array([dqc['Ymjdobs'],dqc['Jmjdobs'], dqc['Hmjdobs'],dqc['Ksmjdobs']]).flatten()
    mjd_start=min(mjdobs[mjdobs > 0])
    mjd_end=max(mjdobs[mjdobs > 0])
    if type == "string":
        isodate_start=Time(mjd_start, format='mjd').isot
        isodate_end=Time(mjd_end, format='mjd').isot

    if type == "float":
        isodate_start = Time(mjd_start, format='mjd').mjd
        isodate_end = Time(mjd_end, format='mjd').mjd

    return isodate_start, isodate_end

def coords(mdatafile, n):

    fh = fits.open(mdatafile)
    dqc = fh[1].data
    for band in ["Y", "J", "H", "KS"]:
        if dqc[band + "CRPIX1"][n] > 0:
            break
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [dqc[band + "CRPIX1"][n], dqc[band + "CRPIX2"][n]]
    w.wcs.cdelt = [dqc[band + "CD1_1"][n], dqc[band + "CD1_2"][n]]
    w.wcs.crval = [dqc[band + "CRVAL1"][n], dqc[band + "CRVAL2"][n]]
    w.wcs.ctype = [dqc[band + "CTYPE1"][n], dqc[band + "CTYPE2"][n]]
    w.wcs.set_pv([(int(dqc[band + "PV2_1"][n]), int(dqc[band + "PV2_2"][n]), dqc[band + "PV2_3"][n])])

    #pix_coords = np.array([[0,0], [0,dqc[band + "AXIS2LENGTH"][n]], [dqc[band + "AXIS1LENGTH"][n], dqc[band + "AXIS2LENGTH"][n]], [dqc[band + "AXIS1LENGTH"][n], 0]])
    pix_coords = np.array([[0,0], [0,dqc[band + "NAXIS2"][n]], [dqc[band + "NAXIS1"][n], dqc[band + "NAXIS2"][n]], [dqc[band + "NAXIS1"][n], 0]])
    world = w.wcs_pix2world(pix_coords, 1)
    world = [[dqc[band + "MINRA"][n], dqc[band + "MINDEC"][n]], [dqc[band + "MINRA"][n], dqc[band + "MAXDEC"][n]], [dqc[band + "MAXRA"][n], dqc[band + "MAXDEC"][n]], [dqc[band + "MAXRA"][n], dqc[band + "MINDEC"][n]]]

    return world

def mag_limit(band, dqc):

    lims = dqc[band + "ZP"] - 2.5*np.log10((5.0*dqc[band + "SKYNOISE"]*np.sqrt(1.4*3.141593)/(dqc[band + "XPIXSIZE"]*dqc[band + "EXPTIME"])))-dqc[band + "APERCOR3"]
    ids = np.where( (lims == lims) )[0]
    #mag_lim = np.median(dqc[band + "_DEPTH_DYE2006"][dqc[band + "_DEPTH_DYE2006"] > -10])
    mag_lim = np.median(lims[ids])
    return mag_lim

def rename_file(file_path, file, dqc, id, cat_header):

    obsids = sorted([int(dqc["Yobsid"][id]), int(dqc["Jobsid"][id]), int(dqc["Hobsid"][id]), int(dqc["Ksobsid"][id])])
    obsid = list(set(obsids))
    obsid = sorted(obsid)
    obsid = obsid[-1]

    filename = str(dqc["jfilename"][id]).strip("'[]'")

    basename_casu = filename[:-10]

    out_name = file[:-5] + "_" + basename_casu + "_" + str(obsid) + file[-5:]
    return out_name


def batch_ph3cat(mdatafile, catpath, outpath,
    verbose = False, overwrite = False,
    nfilesmax = False, debug = False, redo = False):

    """
    loop through all fits files in the catpath and convert to ESO compatible
    Phase3 format
    """

    files = fnmatch.filter(os.listdir(catpath), "*.fits")

    if files == []:
        print("No files in directory, exiting")
        return

    print('Number of files in input directory:', len(files))

    skipped_files = []
    for (n, file) in enumerate(files):
        print(n, file, 'being processed')
        if nfilesmax and n > nfilesmax - 1:
            print("Maximum number (", nfilesmax, ") of files exceeded.")
            print("Skipped files:", skipped_files)
            return

        outfile = outpath + "/" + file
        if not os.path.isdir(outpath):
            os.mkdir(outpath)

        if (not redo and not os.path.exists(outpath + file)) or redo:
            # Change the -2 back
            # if (not redo and not glob.glob(outpath[:-2] + "/" + file[:-5] + "*")) or redo:
                print(n, "out of", len(files))
                (outfile, skipped_files) = mk_ph3cat(
                    mdatafile, catpath, file, outpath, outfile,
                    verbose=False, debug=debug, mcatalog=False,
                    skipped_files=skipped_files)
                if n == 0:
                    make_mcatalog(mdatafile, outfile, outpath, verbose = False)

    print("Skipped files:", skipped_files)

    return

def change_null(fdata):

    for name in fdata.dtype.names:
        try:
            fdata[name][fdata[name] < -1e8] = np.nan
        except ValueError as e:
            print(e)
    return fdata

def make_mcatalog(mdatafile, file, outpath, verbose = False):

    fhlist = fits.open(file)
    #fhlist = fitsio.read(file)
    fdata = fhlist[1].data

    mhlist = fits.open(mdatafile)
    mdata = mhlist[1].data
    ras = mdata["ra"]
    decs = mdata["dec"]
    if verbose: print "RA range:", min(ras), max(ras)
    if verbose: print "Dec range:", min(decs), max(decs)
    if verbose: print "Number of data rows:", len(ras)

    for (n,col) in enumerate(fhlist[1].header):
        if fhlist[1].header[n] == "phot.mag":
            i = fhlist[1].header.keys()[n][4:]
            band = fhlist[1].header["TTYPE" + i][0]
            if band in ["J", "H", "K"]:
                fhlist[1].header[n] = "em.IR." + band
            elif band == "Y":
                fhlist[1].header[n] = "em.IR.NIR"
            elif band == "Z":
                fhlist[1].header[n] = "em.opt.Z"

    """
    del_list = ["RA2000", "DEC2000"]
    n_list = []
    for key in fhlist[1].header.keys():
        if fhlist[1].header[key] in del_list:
            n_list.append(key[5:])
    for n in n_list:
        for key in ["TTYPE", "TCOMM", "TUCD", "TFORM"]:
            del fhlist[1].header[key + n]

    for n in range(max(np.array(n_list, dtype = int))+1, int(fhlist[1].header.keys()[-3][5:])+1):
        for key in ["TTYPE", "TCOMM", "TUCD", "TFORM"]:
            fhlist[1].header[key + str(n-len(n_list))] = fhlist[1].header[key + str(n)]
            del fhlist[1].header[key + str(n)]
    """
    for n in ["1", "2", "3", "4"]:
        del fhlist[0].header["FPDE" + n]
        del fhlist[0].header["FPRA" + n]
        try:
            del fhlist[0].header["PROV" + n]
        except KeyError as e:
            print e

    mjd_start_s, mjd_end_s = mjd(mdatafile, type = "string")
    mjd_start, mjd_end = mjd(mdatafile)
    fhlist[0].header['PRODCATG'] = ('SCIENCE.MCATALOG', 'Data product category')
    fhlist[0].header['DATE-OBS'] = mjd_start_s
    fhlist[0].header['DATE-END'] = mjd_end_s
    fhlist[0].header['MJD-OBS'] = mjd_start
    fhlist[0].header['MJD-END'] = mjd_end
    fhlist[0].header['TELESCOP'] = 'ESO-VISTA'
    fhlist[0].header['DATE'] = (datetime.datetime.utcnow().isoformat(), 'UTC datetime of header update')
    fhlist[0].header['FILTER'] = 'MULTI'
    fhlist[1].data = fhlist[1].data[0:0]
    fhlist.writeto(outpath + '/ph3_mcatalog.fits', clobber = True, checksum = True)

def mk_ph3cat(mdatafile, catpath, catfile, outpath, outfile, verbose = False, debug = False, cache = False, mcatalog = False, skipped_files = []):

    """
    Convert VDFS CASU or VSA bandmerged catalogue into an ESO compliant Phase 3 catalogue file
    mdatafile is file containing dqc metadata from CASU DQC database or VSA
    catfile is the CASU or VSA bandmerged catalogue
    VSA: FITS file for a single frameset
    """

    print "Astropy version:", astropy.__version__
    print "Reading catfile:", catfile
    print "Reading catpath:", catpath

    t0=time.time()
    mhlist = fits.open(mdatafile)
    mdata = mhlist[1].data

    fh_in = fitsio.FITS(catpath + catfile)
    print "Number of extensions:", len(fh_in)
    if verbose: print fh_in

    try:
        indata = fh_in[1].read()
    except ValueError as e:
        print "Corrupted file:", catfile
        skipped_files.append((catfile, e))
        return catfile, skipped_files

    id = np.where( (mdata["FRAMESETID"] == indata["FRAMESETID"][0]) )[0][0]

    out_name = rename_file(catpath, catfile, mdata, id, fh_in[1].read_header())

    if os.path.exists(outpath + "/" + out_name):
        os.remove(outpath + "/" + out_name)

    fh_cat = fitsio.FITS(outpath + "/" + out_name, "rw")
    fh_cat.write(indata, header = fh_in[1].read_header())
    catdata = fh_cat[1].read()

    hdu0 = fh_cat[0].read_header()
    hdu1 = fh_cat[1].read_header()

    for key in hdu0.keys():
        if key <> "COMMENT":
            hdu0[key] = fh_in[0].read_header()[key]

    if verbose: print hdu0

    ras = catdata["RA"]
    decs = catdata["DEC"]

    if len(ras) == 0:
        print "No data in:", catfile
        skipped_files.append((catfile, "Empty file"))
        return catfile, skipped_files

    if verbose: print "RA range:", min(ras), max(ras)
    if verbose: print "Dec range:", min(decs), max(decs)
    if verbose: print "Number of data rows:", len(ras)

    framesetids = list(set(catdata["FRAMESETID"]))
    print "Number of unique framesetids:", len(framesetids)
    if verbose: print framesetids

    mhlist = fits.open(mdatafile)
    mdata = mhlist[1].data
    id = np.where( (mdata["FRAMESETID"] == catdata["FRAMESETID"][0]) )[0][0]

    year = str(datetime.datetime.now().year)
    month = str(datetime.datetime.now().month)
    day = str(datetime.datetime.now().day)
    if len(month) == 1:
        month = "0" + month
    if len(day) == 1:
        day = "0" + day
    date = year + "-" + month + "-" + day
    prog_id='179.A-2010'
    tile_area = area(catdata, mdata["RA"][id], mdata["DEC"][id])
    [[ra1, dec1], [ra2, dec2], [ra3, dec3], [ra4, dec4]] = coords(mdatafile, id)
    """
    print [[ra1, dec1], [ra2, dec2], [ra3, dec3], [ra4, dec4]]
    plt.plot(np.rad2deg(fh_cat[1].data["RA"]), np.rad2deg(fh_cat[1].data["DEC"]), "k.")
    plt.scatter((mdata["JMAXRA"][id]), (mdata["JMINDEC"][id]))
    plt.scatter((mdata["JMINRA"][id]), (mdata["JMINDEC"][id]))
    plt.scatter((mdata["JMAXRA"][id]), (mdata["JMAXDEC"][id]))
    plt.scatter((mdata["JMINRA"][id]), (mdata["JMAXDEC"][id]))
    plt.scatter( ra1, dec1, color = "red")
    plt.scatter( ra2, dec2, color = "red")
    plt.scatter( ra3, dec3)
    plt.scatter( ra4, dec4)
    plt.show()
    """

    cards_to_add = [('PRODCATG', 'SCIENCE.CATALOGTILE', 'Data product category'), \
                    ('ORIGIN', 'ESO-PARANAL', 'European Southern Observatory'), \
                    ('ORIGIN2', 'VDFS-CASU', 'Cambridge Astronomical Survey Unit'), \
                    ('ORIGIN3', 'VDFS-WFAU', 'Wide Field Astronomy Unit'), \
                    ('ORIGIN4', 'VHS-DM', 'VHS DM team'), \
                    ('DATE', date, ' Date the file was written'), \
                    ('TELESCOP', 'VISTA', 'ESO Telescope designation'), \
                    ('INSTRUME', 'VIRCAM', 'Instrument name'), \
                    ('PROG_ID', prog_id, ' ESO programme identification'), \
                    ('OBSTECH', 'IMAGE,JITTER', ' Technique of observation'), \
                    ('OBJECT', 'TBD', ' Target designation'), \
                    ('RA', np.rad2deg(mdata["RA"][id]), ' Field centre (J2000.0)'), \
                    ('DEC', np.rad2deg(mdata["DEC"][id]), ' Field centre (J2000.0)'), \
                    ('REFERENC', '2013Msngr.154...35M', 'Bibliographic reference'), \
                    ('FPRA1', ra1, ' Footprint (J2000.0)'), \
                    ('FPRA2', ra2, ' Footprint (J2000.0)'), \
                    ('FPRA3', ra3, ' Footprint (J2000.0)'), \
                    ('FPRA4', ra4, ' Footprint (J2000.0)'), \
                    ('FPDE1', dec1, ' Footprint (J2000.0)'), \
                    ('FPDE2', dec2, ' Footprint (J2000.0)'), \
                    ('FPDE3', dec3, ' Footprint (J2000.0)'), \
                    ('FPDE4', dec4, ' Footprint (J2000.0)'), \
                    ('SKYSQDEG', tile_area, ' Sky coverage in units of square degrees'), \
                    ('PHOTSYS', 'Vega', ' Photometric system VEGA or AB'), \
                    ('EPS_REG', 'VHS', 'ESO public survey region name'), \
                    ('PROCSOFT', 'VHSDR4/CASU_Version_1.3/VSA_ VHSv20160507', 'Software Version')]

    for (name, value, details) in cards_to_add:
        #fh_cat[0].header.append(card)
        fh_cat[0].write_key(name, value, comment = details)

    id = np.where( (mdata["FRAMESETID"] == catdata["FRAMESETID"][0]) )[0][0]

    del_list = ["CX", "CY", "CZ", "HTMID"]
    n_list = []
    for key in hdu1.keys():
        if str(hdu1[key]).strip(" ") in del_list:
            n_list.append(key[5:])
    for n in n_list:
        for key in ["TTYPE", "TDIM", "TCOMM", "TUCD", "TFORM"]:
            try:
                hdu1.delete(key + n)
            except KeyError as e:
                print e

    hdu_copy = copy.deepcopy(hdu1)
    n = 1
    m = 1
    while n < 155:
        if m == 7:
            n = 11
        for key in ["TTYPE", "TCOMM", "TUCD", "TFORM"]:
            hdu1[key + str(m)] = hdu_copy[key + str(n)]
        n += 1
        m += 1

    for num in ["151", "152", "153", "154"]:
        for key in ["TTYPE", "TDIM", "TCOMM", "TUCD", "TFORM"]:
            hdu1.delete(key + num)

    for num in np.arange(1, 151):

        hdu1["TUNIT" + str(num)] = ""
        if hdu1["TUCD" + str(num)] == "phot.mag":
            band = hdu1["TTYPE" + str(num)][0]
            if band in ["J", "H", "K"]:
                hdu1["TUCD" + str(num)] = "em.IR." + band
            elif band == "Y":
                hdu1["TUCD" + str(num)] = "em.IR.NIR"
            elif band == "Z":
                hdu1["TUCD" + str(num)] = "em.opt.Z"
            hdu1["TUNIT" + str(num)] = "mag"

    n = 1
    convs = [0.618, 0.937, 1.384, 1.839]
    for band in ["Y", "J", "H", "Ks"]:
        filename = str(mdata[band.lower() + "filename"][id]).strip("'[]'")
        if verbose: print band + " filename:", filename
        if filename <> "NONE":
            filename = "158/" + filename[:-4] + "_cat.fits"
            keyword_prov = 'PROV' + str(n)
            keyword_filter = 'FILTER' + str(n)
            fh_cat[0].write_key(keyword_filter, band, comment='Filter name')
            fh_cat[0].write_key(keyword_prov, filename, comment='Image filename')
            hdu0["MAGLIM" + str(n)] = (mag_limit(band.upper(), mdata) + convs[["Y", "J", "H", "Ks"].index(band)], "VHS MEDIAN Depth (AB)")
            if verbose: print keyword_filter +': ', fh_cat[0].header[keyword_filter]
            if verbose: print keyword_prov +': ', fh_cat[0].header[keyword_prov]
            n += 1


    """
    prihdu = fits.Header()
    #This sometimes needs to be from 4, not sure why.
    for key in hdu0.keys()[4:]:
        print key
        prihdu[key] = fh_cat[0].header[key]
    prihdu = fits.PrimaryHDU(header = prihdu)
    fh_cat = fits.HDUList([prihdu, fh_cat[1]])
    """

    if debug: raw_input("Press ENTER to continue: ")

    fh_cat[1].write_key('EXTNAME', 'PHASE3CATALOG', comment = 'FITS Extension name')
    hdu1['TUCD1']='meta.id'
    hdu1['TCOMM1']='Source name in IAU convention'
    hdu1['TTYPE5']='RAJ2000'
    hdu1['TTYPE6']='DECJ2000'
    hdu1['TUNIT5']='degrees'
    hdu1['TUNIT6']='degrees'
    hdu1["TCOMM151"] = 'Flag indicates primary sources (=1), otherwise (=0)'
    hdu1["TUCD151"] = "meta.code.class"
    catdata["RA"] = np.rad2deg(catdata["RA"])
    catdata["DEC"] = np.rad2deg(catdata["DEC"])

    tdata = catdata
    prims = np.zeros(len(tdata))
    ids = np.where( (tdata["PRIORSEC"] <= 0 ) | (tdata["PRIORSEC"] == tdata["FRAMESETID"]) )[0]
    prims[ids] = 1
    #Not put in order yet
    #tdata = np.lib.recfunctions.append_fields(tdata, ['PRIMARY_SOURCE'], [prims], dtypes = int, asrecarray = False)
    fh_cat[1].insert_column("PRIMARY_SOURCE", prims, colnum = 7)
    fh_cat[1].insert_column("RA2000", catdata["RA"], colnum = 5)
    fh_cat[1].insert_column("DE2000", catdata["DEC"], colnum = 6)
    #for 
    #hdu1.delete("RA")
    #fh_cat[1].rename_column("RA", "RAJ2000")
    new_name_list = list(name_list)
    new_name_list[4] = "RA2000"
    new_name_list[5] = "DEC2000"
    tdata.dtype.names = new_name_list
    del_list = ["CX", "CY", "CZ", "HTMID"]
    tdata = mlab.rec_drop_fields(tdata, del_list)

    tdata = change_null(tdata)

    out_name = rename_file(catpath, catfile, mdata, id, fh_cat[0].header)
    print "Writing:", outpath + "/" + out_name

    fh_cat[1].data = tdata
    fh_cat.writeto(outpath + "/" + out_name, clobber=True, checksum=True)

    return outpath + "/" + out_name, skipped_files

if __name__ == "__main__":

    # VSA metadata file
    mdatafile = "/data/vhs/dqc/vsa/2016/VHSv20160507/vhs_vsa_dqc_tiles_fs_metadata.fits"

    # catpath = "/data/vhs/vsa/VHSv20120417/phase3_dr2_v1/"
    # catpath = "/data/vhsardata/Phase3/vsa/VHSv20160507_v1/"
    catpath = "/data/vhsardata/Phase3/test1/"

    #outpath = "/data/vhsardata/Phase3/dr4_v1/"
    outpath = "/data/vhsardata/Phase3/test1/"
    if not os.path.exists(outpath):
        os.mkdir(outpath)

    print('outpath: ', outpath)
    nfilesmax = 10
    batch_ph3cat(mdatafile, catpath, outpath,
                     verbose=False, overwrite=False, debug=False,
                     redo=True, nfilesmax=nfilesmax)
