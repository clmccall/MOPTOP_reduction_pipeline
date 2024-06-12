import moptop_settings as mopset
import moptop_dicts as mopinf
import math
import numpy as np
import pandas as pd
from astropy.io import fits
import scipy.stats as ss
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as plt
from photutils import aperture_photometry, CircularAperture, CircularAnnulus
# from photutils.segmentation import SegmentationImage
from astroquery.simbad import Simbad 
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval,ImageNormalize,SqrtStretch
from photutils.utils import calc_total_error
from photutils import find_peaks
import math
from pylab import *
import csv
import glob
import os
from astropy.time import Time
from datetime import date,timedelta,datetime
from pathlib import Path
import requests
import re
from astroplan import Observer
import shutil
from photutils import DAOStarFinder
from astropy.modeling import models, fitting
import cv2

"""
All functions used for the reduction of MOPTOP data.

down_recentdata         - Downloading of specificed data from recent data.
sort_files              - Sorts the fits files in dir1 into source name folders.
delete_row              - When replacing reduced data entries, deletes the old version.
no_wcs                  - Finds the pixel coordinates for frames without WCS coordinates.
optimum_aperture        - Uses manual coordinates to find the true centres and optimum radii for the objects aperture and annuli.
one_cam_photometry      - For when MOPTOP is operating with only one camera. Photometry of the source and calibration star, using the created stacked fits files.
one_cam_polarimetry     - For when MOPTOP is operating with only one camera. Polarimetry of the source and calibration star.
two_cam_photometry      - For when MOPTOP is operating normally with both cameras. Photometry of the source and calibration star.
two_cam_polarimetry     - For when MOPTOP is operating normally with both cameras. Polarimetry of the source and calibration star.
remove_fits_files       - Once reduced, removes the fits files from all folders within dir1 to clear up space on you machine.
novel_pol_error         - Calculation of the values required to find degree of polarization.
EVPA                    - Calculation of theta and its error, used to find evpa.
evpa_unwrap             - Correction of the evpa 180 degree ambiguity.

THIS SECTION SHOULD NOT NEED EDITING.
"""

def down_recentdata(path,days):
    sources = []
    for day in days:
        # try:
            for proposals_arc in mopset.proposals:
                # try:
                    #setting data to todays data in format YYYYMMDD
                    t = date.today()
                    TIME = datetime.now().strftime("%H:%M:%S")
                    TODAY = t.strftime("%Y%m%d") #todays date in YYYYMMDD format

                    #setting correct date of observation
                    year,month,dayy = t.strftime("%Y"),t.strftime("%m"),t.strftime("%d")
                    today = Time(f'{year}-{month}-{dayy} {TIME}')
                    apo = Observer.at_site("lapalma")
                    sun_set_today = apo.sun_set_time(today, which="nearest") #sun set on day of observing
                    time_suns_today = "{0.iso}".format(sun_set_today)[-12:]
                    sun_set_tomorrow = apo.sun_set_time(today,which="next")
                    time_suns_tomorrow = "{0.iso}".format(sun_set_tomorrow)[-12:]

                    if time_suns_today<TIME<'23:59:59' and day == '':
                        DAY = TODAY
                    if '00:00:00'<TIME<time_suns_tomorrow and day == '':
                        DAY = re.sub('-','',str(t-timedelta(days=1)))
                        print(DAY)

                    if day !='':
                        DAY = str(day)

                    for proposal in proposals_arc.keys():

                        #opening recentdata homepage
                        url1 = f'https://telescope.livjm.ac.uk/DataProd/RecentData/{proposal}'
                        response1 = requests.get(url1,auth=(proposal,proposals_arc[proposal][0]))
                        avail = {} #dictionary containing all observations available on recentdata
                        #searching for string of current year in the response and appending the full YYYYMMDD if on recentdata
                        for i in range(len(response1.text)):
                            if year in response1.text[i:i+4]:
                                url = f'https://telescope.livjm.ac.uk/DataProd/RecentData/{proposal}/{response1.text[i:i+8]}'
                                avail[f'{response1.text[i:i+8]}'] = url

                        ordered_avail = {} #ordered dictionary containing all observations available on recentdata
                        for m,n in sorted(avail.items(),reverse=True):   #reversing the order to get newest first
                            ordered_avail[f'{m}'] = n

                        #opening recentdata homepage for specififc proposal
                        url2 = ordered_avail[DAY]
                        response2 = requests.get(url2,auth=(proposal,proposals_arc[proposal][0]))
                        print(f'{proposal} observations found on RecentData for ({DAY})')

                        #collecting gzipped names of fits images
                        gzip_fits = []
                        for j in range(len(response2.text)):
                            if 'fits' in response2.text[j:j+4]:
                                fits_name = re.sub('"','',response2.text[j-24:j+4])
                                if fits_name[0]!='v' or fits_name[0]!='h':
                                    fits_name = re.sub('"','',response2.text[j-24:j+4])
                                # fits_name = re.sub('>','',fits_name)
                                # fits_name = re.sub('F','',fits_name)
                                # fits_name = re.sub('^s','',fits_name)
                                # fits_name = re.sub(' ','',fits_name)
                                # fits_name = re.sub('=','',fits_name)
                                # fits_name = re.sub('\(','',fits_name)
                                # fits_name = re.sub('^[a-zA-Z]', '', fits_name)
                                fits_name = re.sub(r'^[^0-9_]{1,2}', '', fits_name)
                                
                                if fits_name not in gzip_fits:
                                    gzip_fits.append(fits_name)
                        
                        #creating a txt file with all fits.gz files and appending to it everytime it checks
                        if os.path.exists(f'{path}{proposal}_{DAY}_fits_log.txt'):
                            os.remove(f'{path}{proposal}_{DAY}_fits_log.txt')
                        print(f'Fits log for {DAY} does not yet exist, creating file and writing')
                        fits_name_txt = open(f'{path}{proposal}_{DAY}_fits_log.txt', 'w+')
                        tot = []
                        for name in gzip_fits:
                            if name not in os.listdir(f'{path}'):                        
                                fits_name_txt.write(f'{name}\n')
                                tot.append(name)
                        fits_name_txt.close()

                        #reading in and downloading fits files to data/{TODAY} directory
                        down_fits = open(f'{path}{proposal}_{DAY}_fits_log.txt').readlines()
                        fil = []
                        for Fits in down_fits:
                            try:

                                file = Fits[:-1]
                                fil.append(file)

                                if os.path.exists(f'{path}{file}')==False and os.path.exists(f'{path}{file}.gz')==False and os.path.exists(f"{path}{re.sub('.gz','',file)}")==False:# and os.path.exists(f'data/recentdata/{DAY}/{fits}')==False:
                                    url3 = f'https://telescope.livjm.ac.uk/DataProd/RecentData/{proposal}/{DAY}/{file}'
                                    url3 = re.sub('	','',url3)
                                    # print(url3)

                                    response3 = requests.get(url3,auth=(proposal,proposals_arc[proposal][0]))
                                    open(f'{path}{file}', 'wb').write(response3.content)

                                    reverse_index = {}
                                    for master_name, entry in mopinf.source_info.items():
                                        alt_names = entry.get('alt_names', [])
                                        for alt_name in alt_names:
                                            reverse_index[alt_name] = master_name

                                    # Lookup the master name using the reverse index
                                    hdulist = fits.open(f'{path}{file}')
                                    data = hdulist[0].data.astype(float)
                                    header = hdulist[0].header
                                    source = header['OBJECT']
                                    # print(source,type(source))
                                    name = reverse_index.get(source)
                                    # print(name,type(name))

                                    # Path(mopset.dir1+name+'/single-cam_data/').mkdir(parents=True, exist_ok=True)
                                    # Path(mopset.dir1+name+'/rotation-stacked').mkdir(parents=True, exist_ok=True)
                                    Path(mopset.dir2+name+'/').mkdir(parents=True, exist_ok=True)
                                    # Path(mopset.dir3+name+'/Joined/').mkdir(parents=True, exist_ok=True)

                                    if os.path.exists(mopset.dir1+name+'/'+file):
                                        os.remove(mopset.dir1+name+'/'+file)
                                    shutil.move(mopset.dir1+file,mopset.dir1+name+'/')

                                    print(f"Download progress: {(len(fil)/len(tot))*100:.2f} % {name}")
                                    sources.append(name)
                            except Exception as e:
                                print(e)
                                continue
                    os.remove(f'{path}{proposal}_{DAY}_fits_log.txt')
        #         except Exception as e:
        #             print(e)
        #             continue
        # except Exception as e:
        #     print(e)
        #     continue
    sources = pd.Series(sources,dtype=str).unique()
    return sources

def sort_files(path):
    os.chdir(path)
    files = glob.glob('*.fits')
    tot = 0
    for file in files:
        try:
            reverse_index = {}
            for master_name, entry in mopinf.source_info.items():
                alt_names = entry.get('alt_names', [])
                for alt_name in alt_names:
                    reverse_index[alt_name] = master_name

            hdulist = fits.open(f'{path}{file}')
            data = hdulist[0].data.astype(float)
            header = hdulist[0].header
            source = header['OBJECT']
            name = reverse_index.get(source)

            if os.path.exists(mopset.dir1+name+'/'+file):
                os.remove(mopset.dir1+name+'/'+file)
            shutil.move(mopset.dir1+file,mopset.dir1+name+'/')
            tot += 1
            print(f"sorting files: {(tot/(len(files)))*100:.2f} %")
        except:
            continue  

def delete_row(csv_file, datetime_value):
    temp_file = csv_file + '.tmp'  # Create a temporary file to write updated data

    with open(csv_file, mode='r', newline='') as infile, \
            open(temp_file, mode='w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            if row['datetime'] != str(datetime_value):
                writer.writerow(row)  # Write rows except the one to be deleted

    # Replace the original file with the updated one
    os.replace(temp_file, csv_file)

def no_wcs(path,source):
    # files = []
    # os.chdir(path)
    # files.append(glob.glob('*_e_*_*_*_1_1.fits'))
    # files = [element for sublist in files for element in sublist]

    os.chdir(path)
    files = glob.glob('*_e_'+mopset.time_frame+'_*_1_1.fits')
    coord_frames = []
    for file in files:
        hdulist = fits.open(file)
        header = hdulist[0].header
        wave = header['FILTER1'][4:]
        if wave == 'B':
            coord_frames.append(file)
        if wave == 'V':
            coord_frames.append(file)
        if wave == 'R':
            coord_frames.append(file)
        if wave == 'I':
            coord_frames.append(file)
        if wave == 'L':
            coord_frames.append(file)

    runs = []
    for file in coord_frames:
        run = file[-26:-11]
        if run.startswith('/'):
            run = run[1:]
        if run.endswith('_'):
            run = run[:-1]
        runs.append(run)

    for i in range(len(coord_frames)):
        print(coord_frames[i])
        cam = runs[i][:1]
        date_run = runs[i][4:15]
        if date_run.endswith('_'):
            date_run = date_run[:-1]
        cam_date_run = cam+'_'+date_run
        hdulist = fits.open(coord_frames[i])
        data = hdulist[0].data.astype(float)
        header = hdulist[0].header
        wcs = WCS(header)
        wcs_err  = header['WCS_ERR']

        ra = mopinf.source_info[source]['ra']
        dec = mopinf.source_info[source]['dec']
        cal_ra = mopinf.source_info[source]['cal_ra']
        cal_dec = mopinf.source_info[source]['cal_dec']

        # mask = make_source_mask(data, nsigma=2, npixels=5, dilate_size=11)
        mean, median, std = sigma_clipped_stats(data,sigma=4.0)
        data -= median

        # if wcs_err == 0:
        #     star_coord = SkyCoord(ra,dec, unit=(u.hourangle, u.deg))
        #     xy = SkyCoord.to_pixel(star_coord, wcs=wcs, origin=1)
        #     x = xy[0].item(0)
        #     y = xy[1].item(0) 

        #     star_coord_cal = SkyCoord(cal_ra,cal_dec, unit=(u.hourangle, u.deg))
        #     xy_cal = SkyCoord.to_pixel(star_coord_cal, wcs=wcs, origin=1)
        #     x_cal = xy_cal[0].item(0)
        #     y_cal = xy_cal[1].item(0)
        #     globals()['pix_coords_'+cam_date_run] =[(x_cal,y_cal),(x,y)]

        # else:
        plt.rcParams['figure.figsize'] = (10,6)
        def on_double_click(event):
            if event.dblclick:
                plt.gcf().canvas.mpl_disconnect(cid)
                plt.close()

        def select_target_center(data):
            norm = ImageNormalize(data, interval=ZScaleInterval(), stretch=SqrtStretch())
            plt.imshow(data, cmap='gray', origin='lower', norm=norm)

            print('Right-click on the source center')

            click_coords = plt.ginput(1, timeout=-1, mouse_add=3, mouse_pop=2)
            x, y = click_coords[0]
            print(f'x={x:.2f}, y={y:.2f}')

            # print('Right-click on the cal star center')

            # click_coords_cal = plt.ginput(1, timeout=-1, mouse_add=3, mouse_pop=2)
            # x_cal, y_cal = click_coords_cal[0]
            # print(f'x_cal={x_cal:.2f}, y_cal={y_cal:.2f}')

            x_cal, y_cal = 0,0

            positions = [(x_cal, y_cal), (x, y)]
            return positions

        with open(path+'pix_coords.csv', 'a') as csvfile:
            fieldnames = ['date_run','x_src','y_src','x_cal','y_cal']
            csvfile = csv.DictWriter(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL, fieldnames=fieldnames)
            fileEmpty = os.stat(path+'pix_coords.csv').st_size == 0
            if fileEmpty:
                csvfile.writeheader()
            plt.close()

        done_frames = pd.read_csv(path+'pix_coords.csv')    
        if cam_date_run in done_frames['date_run'].values:
            print(f'{cam_date_run} already done')
            globals()['pix_coords_'+cam_date_run] = [done_frames.loc[done_frames['date_run'] == cam_date_run]['x_cal'].values[0],done_frames.loc[done_frames['date_run'] == cam_date_run]['y_cal'].values[0]],[done_frames.loc[done_frames['date_run'] == cam_date_run]['x_src'].values[0],done_frames.loc[done_frames['date_run'] == cam_date_run]['y_src'].values[0]]
            continue
        else:
            # Connect the double-click event to the function
            cid = plt.gcf().canvas.mpl_connect('pick_event', on_double_click)
            # Should be double click but is right click for some reason
            globals()['pix_coords_'+cam_date_run] = select_target_center(data)
            #save the pixel coordinates to a csv 
            with open(path+'pix_coords.csv', 'a') as csvfile:
                fieldnames = ['date_run','x_src','y_src','x_cal','y_cal']
                csvfile = csv.DictWriter(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL, fieldnames=fieldnames)
                fileEmpty = os.stat(path+'pix_coords.csv').st_size == 0
                if fileEmpty:
                    csvfile.writeheader()
                csvfile.writerow({'date_run':cam_date_run,'x_src':globals()['pix_coords_'+cam_date_run][1][0],'y_src':globals()['pix_coords_'+cam_date_run][1][1],'x_cal':globals()['pix_coords_'+cam_date_run][0][0],'y_cal':globals()['pix_coords_'+cam_date_run][0][1]})
                plt.close()
    return

def optimum_aperture(data, rough_position, obj, camera, std, gain):

    rough_x, rough_y = rough_position
    rough_x = int(round(rough_x))
    rough_y = int(round(rough_y))
    search_size = 10

    y_min = max(0, rough_y - search_size)
    y_max = min(data.shape[0], rough_y + search_size + 1)
    x_min = max(0, rough_x - search_size)
    x_max = min(data.shape[1], rough_x + search_size + 1)
    
    window = data[y_min:y_max, x_min:x_max]
    y_indices, x_indices = np.indices(window.shape)
    
    centroid_y = np.sum(y_indices * window) / np.sum(window)
    centroid_x = np.sum(x_indices * window) / np.sum(window)
    
    pos_center = (x_min + centroid_x, y_min + centroid_y)

    # Determine the optimum aperture size based on SNR
    if (obj == 'src') & ((camera == 1) | (camera == 3)):
        min_radius = 1
        max_radius = 15
        step = 0.5
        snr_values = []
        radii = []
        for radius in np.arange(min_radius, max_radius+step, step):
            aperture = CircularAperture(pos_center, radius)
            error = calc_total_error(data,std,gain)
            phot_table = aperture_photometry(data, aperture, error=error)
            flux = phot_table['aperture_sum'][0]
            flux_err = phot_table['aperture_sum_err'][0]
            snr = flux / flux_err
            snr_values.append(snr)
            radii.append(radius)

        max_snr_index = np.argmax(snr_values)
        optimum_radius = radii[max_snr_index]

        # plt.plot(np.arange(min_radius, max_radius + step, step), snr_values)
        # plt.axvline(x=optimum_radius, color='r', linestyle='--', label='Optimum Radius = ' + str(optimum_radius))
        # plt.xlabel('Aperture Radius [pixels]')
        # plt.ylabel('Signal-to-Noise Ratio (SNR)')
        # plt.title('SNR vs. Aperture Size')
        # plt.legend()
        # plt.grid()
        # plt.show()

        min_factor = 1
        max_factor = 10
        step = 0.5
        snr_values = []
        inner = []
        outer = []
        factors = np.arange(min_factor, max_factor+step, step)
        for factor1 in range(len(factors)-1):
            for factor2 in range(len(factors)-1):
                if factor2 > factor1:
                    aperture = CircularAperture(pos_center, radius)
                    annulus_aperture = CircularAnnulus(pos_center, r_in=optimum_radius*factors[factor1], r_out=optimum_radius*factors[factor2])
                    apertures = [aperture, annulus_aperture]
                    error = calc_total_error(data,std,gain)
                    phot_table = aperture_photometry(data, apertures, error=error)
                    bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
                    bkg_sum = bkg_mean * aperture.area
                    final_sum = phot_table['aperture_sum_0'] - bkg_sum
                    phot_table['residual_aperture_sum'] = final_sum
                    phot_table['residual_aperture_sum'].info.format = '%.8g'
                    flux = phot_table['residual_aperture_sum'][0]
                    flux_err = phot_table['aperture_sum_err_0'][0]
                    snr = flux / flux_err
                    snr_values.append(snr)
                    inner.append(factors[factor1])
                    outer.append(factors[factor2])
        
        max_snr_index = np.argmax(snr_values)
        optimum_annulus_radius = (inner[max_snr_index], outer[max_snr_index])

        # plt.scatter(inner, outer, c=snr_values, cmap='viridis')
        # plt.colorbar(label='SNR')
        # plt.xlabel('Inner Factor')
        # plt.ylabel('Outer Factor')
        # plt.title('SNR vs. Inner and Outer Factors')
        # plt.grid()
        # plt.scatter(optimum_annulus_radius[0], optimum_annulus_radius[1], color='red', marker='*', label='Optimum Factors = ' + str(optimum_annulus_radius))
        # plt.legend()
        # plt.show()

    else:
        optimum_radius = np.nan
        optimum_annulus_radius = np.nan

    return pos_center, optimum_radius, optimum_annulus_radius

def one_cam_photometry(filename,run,source):
    pix = []
    counts_cal_1=[]
    counts_cal_1_err=[]
    counts_cal_2=[]
    counts_cal_2_err=[]
    counts_src_1=[]
    counts_src_1_err=[]
    counts_src_2=[]
    counts_src_2_err=[]
    rotor = []
    bkg = []
    date_time = []
    #Extracts relevant fits header data
    for i in range(1,5):
        for n in range(1,17):
            file=str(i)+filename+str(n)+'_1.fits'
            camera = i
            try:
                hdulist = fits.open(file)
            except:
                continue
            print(file)
            hdulist = fits.open(file)
            data = hdulist[0].data.astype(float)
            header = hdulist[0].header
            gain = header['GAIN']
            wcs = WCS(header)
            rotang = header['ROTANGLE']
            rotskypa = header['ROTSKYPA']
            wave = header['FILTER1'][4:]
            date = header['DATE']
            wcs_err  = header['WCS_ERR']
            mjd = header['MJD']
            time = header['UTSTART']
            seeing = header['L1SEESEC']
            bkg.append(header['BACKGRD'])
            try:
                pho = header['SCHEDPHT']
            except:
                pho = np.nan
                pass

            datetime_string = f"{date} {time}"
            date_time.append(datetime.strptime(datetime_string, "%Y-%m-%d %H:%M:%S.%f"))

            #Extracts relevant coordinate information from dictionary
            ra = mopinf.source_info[source]['ra']
            dec = mopinf.source_info[source]['dec']
            cal_ra = mopinf.source_info[source]['cal_ra']
            cal_dec = mopinf.source_info[source]['cal_dec']
            
            #Background subtraction
            # mask = make_source_mask(data, nsigma=2, npixels=5, dilate_size=11)
            mean, median, std = sigma_clipped_stats(data,sigma=4.0)
            data -=median

            # astroquery_results = Simbad.query_object('OJ287')
            # TARGET_RA = str(astroquery_results[0]['RA'])
            # TARGET_DEC = str(astroquery_results[0]['DEC']).replace('+','').replace('-','')
            # star_coord = SkyCoord(TARGET_RA,TARGET_DEC, unit=(u.hourangle, u.deg))

            #Plot highlighting apertures on both source and cal star to check size and position on relevant objects
            star_coord = SkyCoord(ra,dec, unit=(u.hourangle, u.deg))
            xy = SkyCoord.to_pixel(star_coord, wcs=wcs, origin=1)
            x = xy[0].item(0)
            y = xy[1].item(0) 

            star_coord_cal = SkyCoord(cal_ra,cal_dec, unit=(u.hourangle, u.deg))
            xy_cal = SkyCoord.to_pixel(star_coord_cal, wcs=wcs, origin=1)
            x_cal = xy_cal[0].item(0)
            y_cal = xy_cal[1].item(0)
            positions =[(x_cal,y_cal),(x,y)]

            r = mopinf.source_info[source]['aperture_size']

            pixel_scale = np.round(420/2048,3)
            aperture_size = r*pixel_scale
        
            aperture = CircularAperture(positions, r)
            annulus_aperture = CircularAnnulus(positions, r_in=1.5*r, r_out=2.5*r)
            apertures = [aperture, annulus_aperture]
            if mopset.plot == 'y':
                if n == 1:
                    norm=ImageNormalize(data,interval=ZScaleInterval(),stretch=SqrtStretch())
                    plt.imshow((data),cmap='gray',origin='lower',norm=norm)
                    aperture.plot(color='red', lw=1.5, alpha=0.5)
                    annulus_aperture.plot(color='blue',lw=1.5)
                    plt.show()

            #Photometry performance
            error = calc_total_error(data,std,gain)
            phot_table = aperture_photometry(data,apertures,error=error)
            bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
            bkg_sum = bkg_mean * aperture.area
            final_sum = phot_table['aperture_sum_0'] - bkg_sum
            phot_table['residual_aperture_sum'] = final_sum
            phot_table['residual_aperture_sum'].info.format = '%.8g'

            #Consistent table output
            for col in phot_table.colnames:
                phot_table[col].info.format = '%.8g'

            if (camera == 1):
                counts_cal_1.append(phot_table['residual_aperture_sum'][0])
                counts_cal_1_err.append(phot_table['aperture_sum_err_0'][0])
            if (camera == 2):
                counts_cal_2.append(phot_table['residual_aperture_sum'][0])
                counts_cal_2_err.append(phot_table['aperture_sum_err_0'][0])
            if (camera == 3):
                counts_cal_1.append(phot_table['residual_aperture_sum'][0])
                counts_cal_1_err.append(phot_table['aperture_sum_err_0'][0])
            if (camera == 4):
                counts_cal_2.append(phot_table['residual_aperture_sum'][0])
                counts_cal_2_err.append(phot_table['aperture_sum_err_0'][0])
            
            if (camera == 1):
                counts_src_1.append(phot_table['residual_aperture_sum'][1])
                counts_src_1_err.append(phot_table['aperture_sum_err_0'][1])
            if (camera == 2):
                counts_src_2.append(phot_table['residual_aperture_sum'][1])
                counts_src_2_err.append(phot_table['aperture_sum_err_0'][1])
            if (camera == 3):
                counts_src_1.append(phot_table['residual_aperture_sum'][1])
                counts_src_1_err.append(phot_table['aperture_sum_err_0'][1])
            if (camera == 4):
                counts_src_2.append(phot_table['residual_aperture_sum'][1])
                counts_src_2_err.append(phot_table['aperture_sum_err_0'][1])

            rotor.append(n)

    return date,mjd,date_time[0],wave,aperture_size,rotskypa,rotang,counts_src_1,counts_src_1_err,counts_cal_1,counts_cal_1_err,gain,wcs_err,seeing,pho,np.mean(bkg),x,y

def one_cam_polarimetry(filename,run,source):
    date,mjd,date_time,wave,aperture_size,rotskypa,rotang,m1_src,m1_src_err,m1_cal,m1_cal_err,gain,wcs_err,seeing,pho,bkg,x,y = one_cam_photometry(filename,run,source)

    s1_1_1_src = m1_src[0]+m1_src[1]+m1_src[2]+m1_src[3]+m1_src[4]+m1_src[5]+m1_src[6]+m1_src[7]
    s2_1_1_src = m1_src[0]+m1_src[4]+m1_src[1]+m1_src[5]
    s3_1_1_src = m1_src[1]+m1_src[5]+m1_src[2]+m1_src[6]
    
    s1_1_1_src_err = math.sqrt(m1_src_err[0]*m1_src_err[0]+m1_src_err[1]*m1_src_err[1]+m1_src_err[2]*m1_src_err[2]+m1_src_err[3]*m1_src_err[3]+m1_src_err[4]*m1_src_err[4]+m1_src_err[5]*m1_src_err[5]+m1_src_err[6]*m1_src_err[6]+m1_src_err[7]*m1_src_err[7])
    s2_1_1_src_err = math.sqrt(m1_src_err[0]*m1_src_err[0]+m1_src_err[1]*m1_src_err[1]+m1_src_err[4]*m1_src_err[4]+m1_src_err[5]*m1_src_err[5])
    s3_1_1_src_err = math.sqrt(m1_src_err[1]*m1_src_err[1]+m1_src_err[5]*m1_src_err[5]+m1_src_err[2]*m1_src_err[2]+m1_src_err[6]*m1_src_err[6])

    s1_1_1_cal = m1_cal[0]+m1_cal[1]+m1_cal[2]+m1_cal[3]+m1_cal[4]+m1_cal[5]+m1_cal[6]+m1_cal[7]

    s1_1_1_cal_err = math.sqrt(m1_cal_err[0]*m1_cal_err[0]+m1_cal_err[1]*m1_cal_err[1]+m1_cal_err[2]*m1_cal_err[2]+m1_cal_err[3]*m1_cal_err[3]+ m1_cal_err[4]*m1_cal_err[4]+m1_cal_err[5]*m1_cal_err[5]+m1_cal_err[6]*m1_cal_err[6]+m1_cal_err[7]*m1_cal_err[7])

    q_avg_1 = np.pi*(0.5-s3_1_1_src/s1_1_1_src)
    u_avg_1 = np.pi*(s2_1_1_src/s1_1_1_src-0.5)

    q_err_1 = np.pi*np.sqrt((s3_1_1_src_err/s1_1_1_src)**2+((s3_1_1_src*s1_1_1_src_err)/s1_1_1_src**2)**2)
    u_err_1 = np.pi*np.sqrt((s2_1_1_src_err/s1_1_1_src)**2+((s2_1_1_src*s1_1_1_src_err)/s1_1_1_src**2)**2)

    s1_2_1_src = m1_src[8]+m1_src[9]+m1_src[10]+m1_src[11]+m1_src[12]+m1_src[13]+m1_src[14]+m1_src[15]
    s2_2_1_src = m1_src[8]+m1_src[12]+m1_src[9]+m1_src[13]
    s3_2_1_src = m1_src[9]+m1_src[13]+m1_src[10]+m1_src[14]

    s1_2_1_src_err = math.sqrt(m1_src_err[8]*m1_src_err[8]+m1_src_err[9]*m1_src_err[9]+m1_src_err[10]*m1_src_err[10]+m1_src_err[11]*m1_src_err[11]+m1_src_err[12]*m1_src_err[12]+m1_src_err[13]*m1_src_err[13]+m1_src_err[14]*m1_src_err[14]+m1_src_err[15]*m1_src_err[15])
    s2_2_1_src_err = math.sqrt(m1_src_err[8]*m1_src_err[8]+m1_src_err[9]*m1_src_err[9]+m1_src_err[12]*m1_src_err[12]+m1_src_err[13]*m1_src_err[13])
    s3_2_1_src_err = math.sqrt(m1_src_err[9]*m1_src_err[9]+m1_src_err[13]*m1_src_err[13]+m1_src_err[10]*m1_src_err[10]+m1_src_err[14]*m1_src_err[14])

    s1_2_1_cal = m1_cal[8]+m1_cal[9]+m1_cal[10]+m1_cal[11]+m1_cal[12]+m1_cal[13]+m1_cal[14]+m1_cal[15]

    s1_2_1_cal_err = math.sqrt(m1_cal_err[8]*m1_cal_err[8]+m1_cal_err[9]*m1_cal_err[9]+m1_cal_err[10]*m1_cal_err[10]+m1_cal_err[11]*m1_cal_err[11]+m1_cal_err[12]*m1_cal_err[12]+m1_cal_err[13]*m1_cal_err[13]+m1_cal_err[14]*m1_cal_err[14]+m1_cal_err[15]*m1_cal_err[15])

    q_avg_2 = np.pi*(0.5-s3_2_1_src/s1_2_1_src)
    u_avg_2 = np.pi*(s2_2_1_src/s1_2_1_src-0.5)

    q_err_2 = np.pi*np.sqrt((s3_2_1_src_err/s1_2_1_src)**2+((s3_2_1_src*s1_2_1_src_err)/s1_2_1_src**2)**2)
    u_err_2 = np.pi*np.sqrt((s2_2_1_src_err/s1_2_1_src)**2+((s2_2_1_src*s1_2_1_src_err)/s1_2_1_src**2)**2)

    s1_src = (s1_1_1_src+s1_2_1_src)/2
    s1_src_err = 0.5*np.sqrt(s1_1_1_src_err**2+s1_2_1_src_err**2)

    s1_cal = (s1_1_1_cal+s1_2_1_cal)/2
    s1_cal_err = 0.5*np.sqrt(s1_1_1_cal_err**2+s1_2_1_cal_err**2)

    inst_src = -2.5*np.log10(s1_src)
    inst_src_err = (2.5*0.434*np.array(s1_src_err))/np.array(s1_src)

    inst_cal = -2.5*np.log10(s1_cal)
    inst_cal_err = (2.5*0.434*np.array(s1_cal_err))/np.array(s1_cal)

    q_avg = (q_avg_1+q_avg_2)/2
    u_avg = (u_avg_1+u_avg_2)/2

    q_err = 0.5*np.sqrt(q_err_1**2+q_err_2**2)
    u_err = 0.5*np.sqrt(u_err_1**2+u_err_2**2)

    return date,mjd,date_time,wave,aperture_size,gain,rotskypa,rotang,wcs_err,seeing,pho,s1_src,s1_src_err,inst_src,inst_src_err,s1_cal,s1_cal_err,inst_cal,inst_cal_err,q_avg,q_err,u_avg,u_err,bkg,x,y

def two_cam_photometry(filename,run,source):
    counts_cal_1=[]
    counts_cal_1_err=[]
    counts_cal_2=[]
    counts_cal_2_err=[]
    counts_src_1=[]
    counts_src_1_err=[]
    counts_src_2=[]
    counts_src_2_err=[]
    rotor = []
    bkg = []
    date_time = []

    #Extracts relevant fits header data
    for i in range(1,5):
        for n in range(1,17):
            file=str(i)+filename+str(n)+'_1.fits'
            camera = i
            try:
                hdulist = fits.open(file)
            except:
                continue
            print(file)
            data = hdulist[0].data.astype(float)
            header = hdulist[0].header
            gain = header['GAIN']
            wcs = WCS(header)
            rotang = header['ROTANGLE']
            rotskypa = header['ROTSKYPA']
            wave = header['FILTER1'][4:]
            date = header['DATE']
            wcs_err  = header['WCS_ERR']
            mjd = header['MJD']
            time = header['UTSTART']
            seeing = header['L1SEESEC']
            bkg.append(header['BACKGRD'])
            pix_scale = header['CCDSCALE']
            try:
                pho = header['SCHEDPHT']
            except:
                pho = np.nan
                pass
            
            datetime_string = f"{date} {time}"
            if "." not in datetime_string:
                datetime_string += ".000000"
            dt = datetime.strptime(datetime_string, "%Y-%m-%d %H:%M:%S.%f")
            date_time.append(dt.strftime("%Y-%m-%d %H:%M:%S.%f"))

            #Extracts relevant coordinate information from dictionary
            ra = mopinf.source_info[source]['ra']
            dec = mopinf.source_info[source]['dec']
            cal_ra = mopinf.source_info[source]['cal_ra']
            cal_dec = mopinf.source_info[source]['cal_dec']
            
            # #Background subtraction
            # mask = make_source_mask(data, nsigma=2, npixels=5, dilate_size=11)
            # mean, median, std = sigma_clipped_stats(data,sigma=4.0, mask=mask)
            # data -= median

            #Background subtraction
            def create_circular_mask(h, w, center=None, radius=None):
                if center is None:
                    center = (int(w/2), int(h/2))
                if radius is None:
                    radius = min(center[0], center[1], w-center[0], h-center[1])

                Y, X = np.ogrid[:h, :w]
                dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

                mask = dist_from_center <= radius
                return mask

            def create_source_mask(h, w, sources):
                source_mask = np.ones((h, w), dtype=bool)

                for i in range(len(sources)):
                    center = (int(sources['xcentroid'][i]), int(sources['ycentroid'][i]))
                    radius = 5
                    circular_mask = create_circular_mask(h, w, center=center, radius=radius)
                    source_mask &= ~circular_mask

                return source_mask
            
            try:
                h, w = 1024, 1024

                circular_mask = ~create_circular_mask(h, w, radius=512*0.98)

                mean, median, std = sigma_clipped_stats(data, sigma=4.0, mask=None)

                daofind = DAOStarFinder(fwhm=3, threshold=5.0*std,sigma_radius=3)
                sources = daofind(data - median, mask=circular_mask)

                source_mask = create_source_mask(h, w, sources)

                mask = circular_mask | ~source_mask
                mean, median, std = sigma_clipped_stats(data, sigma=4.0, mask=mask)
                data -= mean
                background = mean

            # masked_data = np.ma.masked_array(data, mask=mask)
            # norm=ImageNormalize(masked_data,interval=ZScaleInterval(),stretch=SqrtStretch())
            # plt.imshow((masked_data),cmap='gray',origin='lower',norm=norm)
            # plt.show()
            
            except:
                mean, median, std = sigma_clipped_stats(data,sigma=4.0)
                data -= median
                background = median

            # astroquery_results = Simbad.query_object('OJ287')
            # TARGET_RA = str(astroquery_results[0]['RA'])
            # TARGET_DEC = str(astroquery_results[0]['DEC']).replace('+','').replace('-','')
            # star_coord = SkyCoord(TARGET_RA,TARGET_DEC, unit=(u.hourangle, u.deg))

            #Converts ra and dec to pixel coordinates for source
            star_coord = SkyCoord(ra,dec, unit=(u.hourangle, u.deg))
            xy = SkyCoord.to_pixel(star_coord, wcs=wcs, origin=1)
            x = xy[0].item(0)
            y = xy[1].item(0) 
            star_coord_cal = SkyCoord(cal_ra,cal_dec, unit=(u.hourangle, u.deg))
            xy_cal = SkyCoord.to_pixel(star_coord_cal, wcs=wcs, origin=1)
            x_cal = xy_cal[0].item(0)
            y_cal = xy_cal[1].item(0)
            positions =[(x_cal,y_cal),(x,y)]

            if source in mopset.no_wcs_sources:
                if wave == 'L':
                    factor = 0
                if wave == 'I':
                    factor = 0
                if wave == 'R':
                    factor = 0
                if wave == 'V':
                    factor = 0
                if wave == 'B':
                    factor = 0
                cam = file[:1]
                date_run = file[4:15]
                if date_run.endswith('_'):
                    date_run = date_run[:-1]
                date, run = date_run.split('_')
                date = int(date)
                run = int(run)+factor
                positions = globals()[f'pix_coords_{cam}_{date}_{run}']

            # print("old positions",positions)
            # print("old aperture size",r)
            # print("old inner annulus",inner_annulus)
            # print("old outer annulus",outer_annulus)

            if mopset.optimum_aperture == 'y':

                if (n == 1):
                    if (i == 1) | (i == 3):
                        src_pos, r, annulus_r = optimum_aperture(data, positions[1],'src',i , std, gain)
                    elif (i == 2) | (i == 4):
                        src_pos, _, _ = optimum_aperture(data, positions[1],'src',i , std, gain)

                    cal_pos, _, _ = optimum_aperture(data, positions[0],'cal', None, std, gain)

                positions = [cal_pos, src_pos]

                inner_annulus = r*annulus_r[0]
                outer_annulus = r*annulus_r[1]

            else:
                r = mopinf.source_info[source]['aperture_size']
                inner_annulus = r*mopinf.source_info[source]['inner_annulus']
                outer_annulus = r*mopinf.source_info[source]['outer_annulus']

            aperture_size = f"{r*pix_scale:.4f}"
            inner_annulus_size = f"{inner_annulus*pix_scale:.4f}"
            outer_annulus_size = f"{outer_annulus*pix_scale:.4f}"

            # print(f'Aperture size: {r} pixels')
            # print(f'Inner annulus size: {annulus_r[0]} pixels')
            # print(f'Outer annulus size: {annulus_r[1]} pixels')

            # print(f'Aperture size: {aperture_size} arcsec')
            # print(f'Inner annulus size: {inner_annulus_size} arcsec')
            # print(f'Outer annulus size: {outer_annulus_size} arcsec')

            # print("new positions",positions)
            # print("new aperture size",r)
            # print("new inner annulus",annulus_r[0])
            # print("new outer annulus",annulus_r[1])
                
            #Creation of apertures
            aperture = CircularAperture(positions, r)
            annulus_aperture = CircularAnnulus(positions, r_in=inner_annulus, r_out=outer_annulus)
            apertures = [aperture, annulus_aperture]

            #Plot highlighting apertures on both source and cal star to check size and position on relevant objects
            if mopset.plot == 'y':
                if n == 1:
                    norm=ImageNormalize(data,interval=ZScaleInterval(),stretch=SqrtStretch())
                    plt.imshow((data),cmap='gray',origin='lower',norm=norm)
                    aperture.plot(color='red', lw=1.5, alpha=0.5)
                    annulus_aperture.plot(color='blue',lw=1.5)
                    plt.show()

            #Photometry performance
            error = calc_total_error(data,std,gain)
            phot_table = aperture_photometry(data,apertures,error=error)
            bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
            bkg_sum = bkg_mean * aperture.area
            final_sum = phot_table['aperture_sum_0'] - bkg_sum
            phot_table['residual_aperture_sum'] = final_sum
            phot_table['residual_aperture_sum'].info.format = '%.8g'

            #Consistent table output
            for col in phot_table.colnames:
                phot_table[col].info.format = '%.8g'

            if (camera == 1):
                counts_cal_1.append(phot_table['residual_aperture_sum'][0])
                counts_cal_1_err.append(phot_table['aperture_sum_err_0'][0])
            if (camera == 2):
                counts_cal_2.append(phot_table['residual_aperture_sum'][0])
                counts_cal_2_err.append(phot_table['aperture_sum_err_0'][0])
            if (camera == 3):
                counts_cal_1.append(phot_table['residual_aperture_sum'][0])
                counts_cal_1_err.append(phot_table['aperture_sum_err_0'][0])
            if (camera == 4):
                counts_cal_2.append(phot_table['residual_aperture_sum'][0])
                counts_cal_2_err.append(phot_table['aperture_sum_err_0'][0])
            
            if (camera == 1):
                counts_src_1.append(phot_table['residual_aperture_sum'][1])
                counts_src_1_err.append(phot_table['aperture_sum_err_0'][1])
            if (camera == 2):
                counts_src_2.append(phot_table['residual_aperture_sum'][1])
                counts_src_2_err.append(phot_table['aperture_sum_err_0'][1])
            if (camera == 3):
                counts_src_1.append(phot_table['residual_aperture_sum'][1])
                counts_src_1_err.append(phot_table['aperture_sum_err_0'][1])
            if (camera == 4):
                counts_src_2.append(phot_table['residual_aperture_sum'][1])
                counts_src_2_err.append(phot_table['aperture_sum_err_0'][1])

            rotor.append(n)

    return date,mjd,date_time[0],wave,aperture_size,rotskypa,rotang,counts_src_1,counts_src_1_err,counts_src_2,counts_src_2_err,counts_cal_1,counts_cal_1_err,counts_cal_2,counts_cal_2_err,gain,wcs_err,seeing,pho,background,positions[1][0],positions[1][1]

def two_cam_polarimetry(filename,run,source):
    #Extracts relavent photometry data 
    date,mjd,date_time,wave,aperture_size,rotskypa,rotang,m1_src,m1_src_err,m2_src,m2_src_err,m1_cal,m1_cal_err,m2_cal,m2_cal_err,gain,wcs_err,seeing,pho,bkg,x,y = two_cam_photometry(filename,run,source)

    s1_1_1_src = m1_src[0]+m1_src[1]+m1_src[2]+m1_src[3]+m1_src[4]+m1_src[5]+m1_src[6]+m1_src[7]
    s1_2_1_src = m1_src[8]+m1_src[9]+m1_src[10]+m1_src[11]+m1_src[12]+m1_src[13]+m1_src[14]+m1_src[15]
    s1_1_2_src = m2_src[0]+m2_src[1]+m2_src[2]+m2_src[3]+m2_src[4]+m2_src[5]+m2_src[6]+m2_src[7]
    s1_2_2_src = m2_src[8]+m2_src[9]+m2_src[10]+m2_src[11]+m2_src[12]+m2_src[13]+m2_src[14]+m2_src[15]

    s1_1_1_src_err = math.sqrt(m1_src_err[0]**2+m1_src_err[1]**2+m1_src_err[2]**2+m1_src_err[3]**2+m1_src_err[4]**2+m1_src_err[5]**2+m1_src_err[6]**2+m1_src_err[7]**2)
    s1_2_1_src_err = math.sqrt(m1_src_err[8]**2+m1_src_err[9]**2+m1_src_err[10]**2+m1_src_err[11]**2+m1_src_err[12]**2+m1_src_err[13]**2+m1_src_err[14]**2+m1_src_err[15]**2)
    s1_1_2_src_err = math.sqrt(m2_src_err[0]**2+m2_src_err[1]**2+m2_src_err[2]**2+m2_src_err[3]**2+m2_src_err[4]**2+m2_src_err[5]**2+m2_src_err[6]**2+m2_src_err[7]**2)
    s1_2_2_src_err = math.sqrt(m2_src_err[8]**2+m2_src_err[9]**2+m2_src_err[10]**2+m2_src_err[11]**2+m2_src_err[12]**2+m2_src_err[13]**2+m2_src_err[14]**2+m2_src_err[15]**2)

    counts_src_array = np.array([s1_1_1_src,s1_2_1_src,s1_1_2_src,s1_2_2_src])
    counts_src = np.mean(counts_src_array)

    counts_src_err_array = np.array([s1_1_1_src_err,s1_2_1_src_err,s1_1_2_src_err,s1_2_2_src_err])
    counts_src_err_array = counts_src_err_array[~np.isnan(counts_src_err_array)]
    if len(counts_src_array) > 1:
        counts_src_err = ss.sem(counts_src_array,axis=0,nan_policy='omit')
    else:
        counts_src_err = nanmean(counts_src_err_array)
    
    inst_src = -2.5*np.log10(counts_src)
    inst_src_err = (2.5*0.434*counts_src_err)/counts_src

    s1_1_1_cal = m1_cal[0]+m1_cal[1]+m1_cal[2]+m1_cal[3]+m1_cal[4]+m1_cal[5]+m1_cal[6]+m1_cal[7]
    s1_2_1_cal = m1_cal[8]+m1_cal[9]+m1_cal[10]+m1_cal[11]+m1_cal[12]+m1_cal[13]+m1_cal[14]+m1_cal[15]
    s1_1_2_cal = m2_cal[0]+m2_cal[1]+m2_cal[2]+m2_cal[3]+m2_cal[4]+m2_cal[5]+m2_cal[6]+m2_cal[7]
    s1_2_2_cal = m2_cal[8]+m2_cal[9]+m2_cal[10]+m2_cal[11]+m2_cal[12]+m2_cal[13]+m2_cal[14]+m2_cal[15]

    s1_1_1_cal_err = math.sqrt(m1_cal_err[0]**2+m1_cal_err[1]**2+m1_cal_err[2]**2+m1_cal_err[3]**2+m1_cal_err[4]**2+m1_cal_err[5]**2+m1_cal_err[6]**2+m1_cal_err[7]**2)
    s1_2_1_cal_err = math.sqrt(m1_cal_err[8]**2+m1_cal_err[9]**2+m1_cal_err[10]**2+m1_cal_err[11]**2+m1_cal_err[12]**2+m1_cal_err[13]**2+m1_cal_err[14]**2+m1_cal_err[15]**2)
    s1_1_2_cal_err = math.sqrt(m2_cal_err[0]**2+m2_cal_err[1]**2+m2_cal_err[2]**2+m2_cal_err[3]**2+m2_cal_err[4]**2+m2_cal_err[5]**2+m2_cal_err[6]**2+m2_cal_err[7]**2)
    s1_2_2_cal_err = math.sqrt(m2_cal_err[8]**2+m2_cal_err[9]**2+m2_cal_err[10]**2+m2_cal_err[11]**2+m2_cal_err[12]**2+m2_cal_err[13]**2+m2_cal_err[14]**2+m2_cal_err[15]**2)

    counts_cal_array = np.array([s1_1_1_cal,s1_2_1_cal,s1_1_2_cal,s1_2_2_cal])
    counts_cal = np.mean(counts_cal_array)

    counts_cal_err_array = np.array([s1_1_1_cal_err,s1_2_1_cal_err,s1_1_2_cal_err,s1_2_2_cal_err])
    counts_cal_err_array = counts_cal_err_array[~np.isnan(counts_cal_err_array)]
    if len(counts_cal_array) > 1:
        counts_cal_err = ss.sem(counts_cal_array,axis=0,nan_policy='omit')
    else:
        counts_cal_err = nanmean(counts_cal_err_array)

    inst_cal = -2.5*np.log10(counts_cal)
    inst_cal_err = (2.5*0.434*counts_cal_err)/counts_cal

    # Alt inst mag calculations
    # counts_1 = m1_src[0]+m1_src[1]+m2_src[0]+m2_src[1]
    # counts_2 = m1_src[2]+m1_src[3]+m2_src[2]+m2_src[3]
    # counts_3 = m1_src[4]+m1_src[5]+m2_src[4]+m2_src[5]
    # counts_4 = m1_src[6]+m1_src[7]+m2_src[6]+m2_src[7]
    # counts_5 = m1_src[8]+m1_src[9]+m2_src[8]+m2_src[9]
    # counts_6 = m1_src[10]+m1_src[11]+m2_src[10]+m2_src[11]
    # counts_7 = m1_src[12]+m1_src[13]+m2_src[12]+m2_src[13]
    # counts_8 = m1_src[14]+m1_src[15]+m2_src[14]+m2_src[15]

    # counts_src = np.nanmean([counts_1,counts_2,counts_3,counts_4,counts_5,counts_6,counts_7,counts_8])
    # counts_src_err = ss.sem([counts_1,counts_2,counts_3,counts_4,counts_5,counts_6,counts_7,counts_8],nan_policy='omit')

    # inst_src = -2.5*np.log10(counts_src)
    # inst_src_err = (2.5*0.434*np.array(counts_src_err))/np.array(counts_src)

    # counts_1 = m1_cal[0]+m1_cal[1]+m2_cal[0]+m2_cal[1]
    # counts_2 = m1_cal[2]+m1_cal[3]+m2_cal[2]+m2_cal[3]
    # counts_3 = m1_cal[4]+m1_cal[5]+m2_cal[4]+m2_cal[5]
    # counts_4 = m1_cal[6]+m1_cal[7]+m2_cal[6]+m2_cal[7]
    # counts_5 = m1_cal[8]+m1_cal[9]+m2_cal[8]+m2_cal[9]
    # counts_6 = m1_cal[10]+m1_cal[11]+m2_cal[10]+m2_cal[11]
    # counts_7 = m1_cal[12]+m1_cal[13]+m2_cal[12]+m2_cal[13]
    # counts_8 = m1_cal[14]+m1_cal[15]+m2_cal[14]+m2_cal[15]

    # counts_cal = np.mean([counts_1,counts_2,counts_3,counts_4,counts_5,counts_6,counts_7,counts_8])
    # counts_cal_err = ss.sem([counts_1,counts_2,counts_3,counts_4,counts_5,counts_6,counts_7,counts_8])

    # inst_cal = -2.5*np.log10(counts_cal)
    # inst_cal_err = (2.5*0.434*np.array(counts_cal_err))/np.array(counts_cal)

    pos_q = 0
    pos_u = 1

    q_all = []
    u_all = []

    q_err_all = []
    u_err_all = []

    for i in range(4):
        try:
            f_q = np.sqrt((m2_src[pos_q]*m2_src[pos_q+2])/(m1_src[pos_q]*m1_src[pos_q+2]))
            f_u = np.sqrt((m2_src[pos_u]*m2_src[pos_u+2])/(m1_src[pos_u]*m1_src[pos_u+2]))
            
            c1 = m1_src[pos_q]
            c2 = m1_src[pos_u]
            c3 = m1_src[pos_q+2]
            c4 = m1_src[pos_u+2]
            
            d1 = m2_src[pos_q]/f_q
            d2 = m2_src[pos_u]/f_u
            d3 = m2_src[pos_q+2]/f_q
            d4 = m2_src[pos_u+2]/f_u
            
            Q1 = c1 - d1
            U1 = c2 - d2
            I1 = (c1+d1+c2+d2)/2
            q1 = Q1/I1
            u1 = U1/I1
            
            Q2 = -(c3 - d3)
            U2 = -(c4 - d4)
            I2 = (c3+d3+c4+d4)/2
            q2 = Q2/I2
            u2 = U2/I2
            
            q = np.nanmean([q1,q2])
            u = np.nanmean([u1,u2])
            
            A_q = np.sqrt(m2_src_err[pos_q]**2+(((m1_src[pos_q+2]/(2*m2_src[pos_q]*m2_src[pos_q+2]*m1_src[pos_q]))**2)*(1/(m2_src[pos_q]*m2_src[pos_q+2]*m1_src[pos_q])**2)*((m2_src_err[pos_q]/m2_src[pos_q])**2+(m2_src_err[pos_q+2]/m2_src[pos_q+2])**2+(m1_src_err[pos_q]/m1_src[pos_q])**2))+(m1_src_err[pos_q+2]/m1_src[pos_q+2])**2)
            q_err = q*np.sqrt((A_q/(f_q*m1_src[pos_q]-m2_src[pos_q]))**2+(A_q/(f_q*m1_src[pos_q]+m2_src[pos_q]))**2)

            A_u = np.sqrt(m2_src_err[pos_u]**2+(((m1_src[pos_u+2]/(2*m2_src[pos_u]*m2_src[pos_u+2]*m1_src[pos_u]))**2)*(1/(m2_src[pos_u]*m2_src[pos_u+2]*m1_src[pos_u])**2)*((m2_src_err[pos_u]/m2_src[pos_u])**2+(m2_src_err[pos_u+2]/m2_src[pos_u+2])**2+(m1_src_err[pos_u]/m1_src[pos_u])**2))+(m1_src_err[pos_u+2]/m1_src[pos_u+2])**2)
            u_err = u*np.sqrt((A_u/(f_u*m1_src[pos_u]-m2_src[pos_u]))**2+(A_u/(f_u*m1_src[pos_u]+m2_src[pos_u]))**2)

            q_all.append(q)
            q_err_all.append(q_err)
            u_all.append(u)
            u_err_all.append(u_err)

        except:
            q_all.append(np.nan)
            q_err_all.append(np.nan)
            u_all.append(np.nan)
            u_err_all.append(np.nan)
            continue

        pos_q += 4
        pos_u += 4

    q_avg = np.nanmean(q_all)
    u_avg = np.nanmean(u_all)

    try:
        q_err_all_no_nan = np.nan_to_num(q_err_all, nan=0.0)
        q_err = np.sqrt(np.sum(q_err_all_no_nan**2)) / np.count_nonzero(~np.isnan(q_err_all))

        u_err_all_no_nan = np.nan_to_num(u_err_all, nan=0.0)
        u_err = np.sqrt(np.sum(u_err_all_no_nan**2)) / np.count_nonzero(~np.isnan(u_err_all))
    except:
        q_err = np.nan
        u_err = np.nan
    
    return date,mjd,date_time,wave,aperture_size,gain,rotskypa,rotang,wcs_err,seeing,pho,counts_src,counts_src_err,inst_src,inst_src_err,counts_cal,counts_cal_err,inst_cal,inst_cal_err,q_avg,q_err,u_avg,u_err,bkg,x,y

def remove_fits_files(root_folder, file_extension='.fits'):
    for folder_path, _, file_names in os.walk(root_folder):
        for file_name in file_names:
            if file_name.endswith(file_extension):
                file_path = os.path.join(folder_path, file_name)
                try:
                    os.remove(file_path)
                except Exception as e:
                    print(f"Error removing {file_path}: {e}")

def novel_pol_error(q_i, q_err, u_i, u_err):
    
    sigma_a_sq = 0.5*(q_err*q_err + u_err*u_err)

    p_i = np.sqrt(q_i*q_i + u_i*u_i)

    b_i_sq = (q_i*q_i*u_err*u_err + u_i*u_i*q_err*q_err) / (q_i*q_i+u_i*u_i)

    p_mas = p_i - (b_i_sq / (2.0*p_i)) * (1.0-np.exp(-p_i*p_i/b_i_sq))
    
    sigma = np.sqrt(sigma_a_sq)
    
    P_alpha = 1.0
    beta = 0.97
    gamma = 2.01

    p_max_over_sigma = p_mas/sigma + P_alpha * (1.0 - beta*np.exp(-gamma*p_mas/sigma))
    
    p_max = p_max_over_sigma * sigma
    
    P_alpha = 1.0
    beta = 0.72
    gamma = 0.60
    omega = -0.83
    phi = 4.41

    p_min = p_mas - sigma * P_alpha * ( 1.0 + beta*np.exp(-gamma*p_mas/sigma) * np.sin(omega*p_mas/sigma + phi))

    # p_err = np.sqrt(q_i**2*q_err**2+u_i**2*u_err**2)/p_i

    return p_mas,p_min,p_max

def EVPA(q_i, q_err, u_i, u_err):
    #Calculation of theta for evpa using q and u.
    theta = np.degrees(0.5*np.arctan2(u_i,q_i)) % 180
    # theta = (0.5*np.arctan2(u_i,q_i))*180/np.pi() % 180
    # theta_err = np.sqrt((1/(2*q_i*(1+(u_i**2/q_i**2)))**2*u_err**2+(-u_i/(2*((q_i**2)+u_i**2)))**2*q_err**2))*180/np.pi
    theta_err = np.degrees(np.sqrt(u_i**2*q_err**2+q_i**2*u_err**2)/(2*(q_i**2+u_i**2)))#*180/np.pi
    return theta,theta_err

def evpa_unwrap(evpa,evpa_err,n):
    for i in range(n,len(evpa)):

        #Added errors of the previous n points
        m = 1
        comb_err = 0
        while m <= n:
            comb_err += evpa_err[i-m]
            m += 1
            
        #Weights of the previous n points
        m = 1
        while m <= n:
            globals()['weight'+str(m)] = 1-(evpa_err[i-m]/comb_err)
            m += 1

        #Weight calculation of the previous n points
        m = 1
        weights = 0
        while m <= n:
            weights += evpa[i-m] * globals()['weight'+str(m)]
            m += 1
            
        #Added weights of the previous n points 
        m = 1
        comb_weights = 0
        while m <= n:
            comb_weights += globals()['weight'+str(m)]
            m += 1
            
        #Finds the difference between a point and the weighted average of the previous n points
        diff = evpa[i] - (weights/comb_weights)
        
        #Added errors of all the points
        m = 0
        comb_err = 0
        while m <= n:
            comb_err += evpa_err[i-m]
            m += 1

        #Weight calculation of all the points
        m = 0
        while m <= n:
            globals()['weight'+str(m)] = 1-(evpa_err[i-m]/comb_err)
            m += 1

        #Weight calculation of all the points
        m = 0
        weights = 0
        while m <= n:
            weights += evpa_err[i-m] * globals()['weight'+str(m)]
            m += 1
        
        #Added weights of all the points
        m = 0
        comb_weights = 0
        while m <= n:
            comb_weights += globals()['weight'+str(m)]
            m += 1
        
        #Finds the combined error factor of all points involved
        err = weights/comb_weights

        #The number of 180 degree shifts required
        N = np.round((abs(diff)-err)/180)

        #Performs the necessary shift
        if abs(diff) - err > 90 and diff < 0:
            evpa[i] += N*180
        if abs(diff) - err > 90 and diff > 0:
            evpa[i] -= N*180

    evpa_cor = evpa.reset_index(drop=True)
    return evpa_cor
