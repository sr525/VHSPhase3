def UploadFile(table_name):

    """
    Code for uploading fits tables to the SQL database
    Quite temperamental - don't press control C
    The problem is the code it calls fits2sqlite.py which was copied from someone
    on the internet and will delete input data if upset - be warned and copy 
    input data first
    """

    import subprocess
    import numpy as np
    from astropy.table import Table, vstack, join, Column, hstack

    table_date = table_name[-8:]
    if table_name[0] == "S":
        table_prog = "SExtractor"
    elif table_name[0] == "A":
        table_prog == "Aegean"
    
    filename = "/data/sr525/AlmaArchive/" + table_name + ".fits"

    t = Table.read(filename)
    nums = np.arange(1, len(t)+1)
    t["NUMBER"] = nums
    t["CLASS"] = [0.0]*len(t)
    t.write(filename, overwrite = True)

    subprocess.call(["python2.7", "/home/sr525/Python_Code/Web_Stuff/fits2sqlite.py", \
                    "-t", table_name, "/home/sr525/DB/" + table_date + "/" + \
                    table_prog + ".db", filename, "-d", "0", "-u", "NUMBER"])

    dbs = os.listdir("/home/sr525/DB/" + table_date + "/")

    if table_prog + ".db" not in dbs:

        subprocess.call(["sqlite3", "/home/sr525/DB/" + table_date + "/" + \
                        table_prog + ".db", "CREATE TABLE \
                       Comments(COADD_OBJECTS_ID VARCHAR(10), NAME VARCHAR(100), \
                       TIME VARCHAR(40), COMMENT VARCHAR(1000));"]

def ImageDownload(filename):

    """
    Grabs the images and stores them locally
    """

    from astropy.table import Table
    import astropy.io.fits as fits
    import subprocess
    from astropy import coordinates as coord
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    t = Table.read(filename)

    im_files = t["Image"]
    im_files = set(im_files)
    bash_script = "/home/sr525/bash_scripts/obelics_wget.bash"
    dir = "/data/desardata4/ASTERICS/"

    for filename in im_files:
    
        subprocess.call(["bash", bash_script, filename, "temp.fits"])
    
        temp_file = dir + "temp.fits"
        with fits.open(temp_file) as fhlist:
            #print fhlist
            btype = fhlist[0].header["BTYPE"]
            obsra = fhlist[0].header["OBSRA"]
            obsdec = fhlist[0].header["OBSDEC"]
            c = 299792458.0*1e6
            wav = str(int(c/fhlist[0].header["CRVAL3"]))

            #c = SkyCoord(ra=obsra*u.degree, dec=obsdec*u.degree, frame='icrs')
            #ra = c.ra.to_string(u.degree).replace(".", "p")
            #dec = c.dec.to_string(u.degree, alwayssign = True).replace(".", "p")

            ra = str(obsra).replace(".", "p")
            dec = str(obsdec).replace(".", "p")

            if obsdec > 0:
                dec = "+" + dec

            out_file = dir + ra + dec + "_" + wav + "_" + btype + ".fits"

            subprocess.call(["mv", temp_file, out_file])

def MakeCornerFile():

    """
    Takes a folder full of images, reads each and finds the corners
    Creates a fits file of the minimun and maximum ra and dec
    """

    import os
    from astropy import wcs
    import astropy.io.fits as fits
    from astropy.table import Table

    dir = "/data/desardata4/ASTERICS/"
    files = os.listdir(dir)

    min_ras = []
    max_ras = []
    min_decs = []
    max_decs = []
    filenames = []

    for filename in files:
        if filename[-4:] == "fits":
            with fits.open(dir + filename) as fhlist:
                hdr = fhlist[0].header
                w = wcs.WCS(hdr)
                corners = [[0,0,0,0], [hdr["NAXIS1"], hdr["NAXIS2"], 0, 0]]
                corners = w.wcs_pix2world(corners, 1)
                corners_ra = [corners[0][0], corners[1][0]]
                corners_dec = [corners[0][1], corners[1][1]]
                min_ras.append(min(corners_ra))
                max_ras.append(max(corners_ra))
                min_decs.append(min(corners_dec))
                max_decs.append(max(corners_dec))
                filenames.append(filename)

    t = Table(data = [filenames, min_ras, max_ras, min_decs, max_decs], names = ["FILENAME", "MIN_RA", "MAX_RA", "MIN_DEC", "MAX_DEC"])

    t.write(dir + "Corners_File.fits")

def MakeCutouts(filename, width = 30.0, units = "arcsecs"):

    """
    Takes a file of input ra and dec - filename
    Uses the corner file made in MakeCornerFile to identify relevant images
    makes a cutout image of width defined by width and units
    """


    from astropy.table import Table
    import numpy as np
    import matplotlib.pyplot as plt
    import srpylib
    import astropy.io.fits as fits
    from astropy import wcs
    import srpylib
    import matplotlib.gridspec as gridspec

    dir = "/data/desardata4/ASTERICS/"
    #dir = "/home/sreed/ASTERICS/"
    outdir = "/home/sr525/public_html/Research/Alma/"

    t = Table.read(filename)
    info = Table.read(dir + "Corners_File.fits")

    if units == "arcsecs":
        width = width/3600.0

    try:
        ra = t["ra"][0]
    except KeyError:
        t["ra"] = t["ALPHA_J2000"]
        t["dec"] = t["DELTA_J2000"]

    n = 0
    while n < len(t):
        #while n < 1:
        ra = t["ra"][n]
        dec = t["dec"][n]
        #ra = 334.2600
        #dec = 0.2663

        ids = np.where( (ra > info["MIN_RA"]) & (ra < info["MAX_RA"]) & (dec > info["MIN_DEC"]) & (dec < info["MAX_DEC"]) )[0]
        filenames = info["FILENAME"][ids]
        #filenames = ["334p259125+0p266277777778_1139_Intensity.fits"]

        if len(filenames) > 4:
            num_rows = int(np.ceil(len(filenames)/4.0))
            print "Number of rows:", int(num_rows)
            print "Number of cols:", 4
            num_cols = 4
        else:
            num_rows = 1
            num_cols = len(filenames)

        gs = gridspec.GridSpec(num_rows, num_cols)
        fig = plt.figure()
        j = 0
        k = 0

        for (i, im_file) in enumerate(filenames):
            if i % 4 == 0 and i <> 0:
                j += 1
                k = 0
            
            with fits.open(dir + im_file) as fhlist:
                hdr = fhlist[0].header
                im_data = fhlist[0].data[0,0,:,:]

                w = wcs.WCS(hdr)
                
                max_corners = [[0,0,0,0], [hdr["NAXIS1"], hdr["NAXIS2"],0,0]]
                max_coords = w.wcs_pix2world(max_corners, 1)
                max_width = np.fabs(max_coords[0][0]-max_coords[1][0])

                if width > max_width:
                    im = im_data
                    print "Specififed width larger than image, returning whole image"

                else:
                
                    corner1 = [ra-width/2.0/np.cos(np.deg2rad(dec)), dec-width/2.0,0,0]
                    corner4 = [ra+width/2.0/np.cos(np.deg2rad(dec)), dec+width/2.0,0,0]
                    coords = [corner1, corner4]
                    pix = w.wcs_world2pix(coords, 1)
                    [[pix1, pix2, _, _], [pix3, pix4, _, _,]] = pix
                    #print pix1, pix2, pix3, pix4

                    if pix2 < 0.0 or pix3 < 0.0 or pix4 > hdr["NAXIS1"] or pix1 > hdr["NAXIS2"]:
                        [[cra, cdec,_,_]] = w.wcs_pix2world([[0.0,0.0,0.0,0.0]], 1)
                        c2 = w.wcs_world2pix([[cra-width, cdec+width, 0,0]], 1)
                        #print c2
                        im_size = int(np.ceil(c2[0][1]))
                        im = np.zeros((im_size, im_size))
                        off1 = im_size
                        off2 = 0
                        off3 = 0
                        off4 = im_size
                        if pix1 > hdr["NAXIS2"]:
                            off1 = im_size - int(np.ceil(pix1 - hdr["NAXIS2"]))
                            pix1 = hdr["NAXIS2"]
                        if pix4 > hdr["NAXIS1"]:
                            off4 = im_size - int(np.ceil(pix4 - hdr["NAXIS1"]))
                            pix4 = hdr["NAXIS1"]
                        if pix2 < 0.0:
                            off2 = int(np.ceil(0.0 - pix2))
                            pix2 = 0
                        if pix3 < 0.0:
                            off3 = int(np.ceil(0.0 - pix3))
                            pix3 = 0
                        pix1 = int(np.floor(pix1))
                        pix2 = int(np.ceil(pix2))
                        pix3 = int(np.ceil(pix3))
                        pix4 = int(np.floor(pix4))
                        #print pix1, pix2, pix3, pix4
                        #print off1, off2, off3, off4

                        im[off2:(pix4-pix2+off2),off3:(pix1-pix3+off3)] = im_data[int(pix2):int(pix4), int(pix3):int(pix1)]

                    else:
                        im = im_data[int(pix2):int(pix4),int(pix3):int(pix1)]

                med = np.median(im.flatten())
                MAD = srpylib.MAD(im.flatten(), med)
                vmin = med - 5.0*MAD
                vmax = med + 5.0*MAD

                ax = fig.add_subplot(gs[j, k])
                ax.axes.get_xaxis().set_visible(False)
                ax.axes.get_yaxis().set_visible(False)
                for i in ax.spines.itervalues():
                    i.set_linewidth(3.0)

                c = 299792458.0*1e6
                wav = str(int(c/hdr["CRVAL3"]))
                cmap = "Blues"

                ax.set_title(wav + " - " + str(width*3600.0) + '"')
                ax.imshow(im, vmin = vmin, vmax = vmax, cmap = cmap)
                k += 1

        print outdir + filename[24:-5] + "_" + str(t["NUMBER"][n]) + "_AlmaIms.png"
        plt.savefig(outdir + filename[24:-5] + "_" + str(t["NUMBER"][n]) + "_AlmaIms.png", transparent = True)
        plt.close()


        n += 1




filename = "/data/sr525/AlmaArchive/AObjectTable20170714.fits"
table_name = "AObjectTable20170714"
MakeCutouts(filename, width = 20.0)
