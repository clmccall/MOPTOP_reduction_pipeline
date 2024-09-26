import moptop_settings as mopset
import moptop_dicts as mopinf
import math
import numpy as np
import pandas as pd
from astropy.io import fits
import scipy.stats as ss
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS, FITSFixedWarning
import astropy.units as u
import matplotlib.pyplot as plt
from photutils.aperture import aperture_photometry, CircularAperture, CircularAnnulus
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval,ImageNormalize,SqrtStretch
from photutils.utils import calc_total_error
import math
from pylab import *
import csv
import glob
import os
from datetime import datetime
import requests
import re
import shutil
from photutils.detection import DAOStarFinder
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", category=FITSFixedWarning)

def down_recentdata(path,days):
    """
    Downloads and organises data from the Liverpool Telescope Recent Data webpage.

    Parameters:
    path (str): The directory path where the downloaded files will be stored.
    days (list): A list of specific days (in YYYYMMDD format) for which the data should be downloaded. 

    Returns:
    list: A list of unique source names corresponding to the downloaded FITS files.

    Description:
    This function performs the following steps:
    1. Iterates through the specified days to download data.
    2. Constructs URLs to access recent data from the telescope's archive for each proposal.
    3. Searches for FITS files in the archive for the specified day.
    4. Writes the names of the found FITS files to a log file.
    5. Downloads the FITS files listed in the log file to the specified directory.
    6. Moves the downloaded FITS files to the corresponding source directories.
    7. Removes the log file after processing.

    Notes:
    - It uses authentication to access the telescope data, which should be provided in the `proposals_arc` dictionary in moptop_settings.py.
    - The `mopinf` and `mopset` modules are assumed to be pre-configured with relevant information such as directory paths and source info.

    """

    sources = []
    for day in days:
        for proposals_arc in mopset.proposals:
            DAY = str(day)
            for proposal in proposals_arc.keys():

                #opening recentdata homepage
                url1 = f'https://telescope.livjm.ac.uk/DataProd/RecentData/{proposal}'
                response1 = requests.get(url1,auth=(proposal,proposals_arc[proposal][0]))
                avail = {} #dictionary containing all observations available on recentdata
                #searching for string of current year in the response and appending the full YYYYMMDD if on recentdata
                for i in range(len(response1.text)):
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

                #reading in and downloading fits files to directory
                down_fits = open(f'{path}{proposal}_{DAY}_fits_log.txt').readlines()
                fil = []
                for Fits in down_fits:
                    try:

                        file = Fits[:-1]
                        fil.append(file)

                        if os.path.exists(f'{path}{file}')==False and os.path.exists(f'{path}{file}.gz')==False and os.path.exists(f"{path}{re.sub('.gz','',file)}")==False:# and os.path.exists(f'data/recentdata/{DAY}/{fits}')==False:
                            url3 = f'https://telescope.livjm.ac.uk/DataProd/RecentData/{proposal}/{DAY}/{file}'
                            url3 = re.sub('	','',url3)

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
                            name = reverse_index.get(source)

                            os.makedirs(mopset.dir1+name, exist_ok=True)
                            if os.path.exists(mopset.dir1+name+'/'+file):
                                os.remove(mopset.dir1+name+'/'+file)
                            shutil.move(mopset.dir1+file,mopset.dir1+name+'/')

                            print(f"Download progress: {(len(fil)/len(tot))*100:.2f} % {name}")
                            sources.append(name)
                    except Exception as e:
                        print(e)
                        continue
            os.remove(f'{path}{proposal}_{DAY}_fits_log.txt')
    sources = pd.Series(sources,dtype=str).unique()
    return sources

def sort_files(path):
    """
    Sorts and moves FITS files from the specified directory to source-specific subdirectories.

    Parameters:
    path (str): The directory path containing the FITS files to be sorted.

    Description:
    This function performs the following steps:
    1. Changes the current working directory to the specified path.
    2. Finds all FITS files in the directory.
    3. For each FITS file, reads the header to determine the source name.
    4. Uses a reverse index to map alternate source names to master names.
    5. Moves each FITS file to a subdirectory named after its corresponding source.
    6. Prints the sorting progress as a percentage.

    Notes:
    - The function relies on the `mopinf` and `mopset` modules for source information and directory paths.
    - Assumes that FITS files contain a header keyword 'OBJECT' that identifies the source name.
    - If a file already exists in the target directory, it will be overwritten.

    """

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
            os.makedirs(mopset.dir1+name, exist_ok=True)
            shutil.move(mopset.dir1+file,mopset.dir1+name+'/')
            tot += 1
            print(f"sorting files: {(tot/(len(files)))*100:.2f} %")
        except:
            continue  

def delete_row(csv_file, datetime_value):
    """
    Deletes a row from a CSV file based on a specific datetime value.

    Parameters:
    csv_file (str): The path to the CSV file.
    datetime_value (str or datetime): The datetime value identifying the row to be deleted.

    Description:
    This function performs the following steps:
    1. Creates a temporary file to store the updated data.
    2. Reads the original CSV file.
    3. Writes all rows to the temporary file except the row with the specified datetime value.
    4. Replaces the original CSV file with the updated temporary file.

    Notes:
    - The CSV file must have a column named 'datetime'.
    - The datetime value in the CSV file should be comparable to the string representation of `datetime_value`.

    """

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

def no_wcs(path, source):
    """
    Processes FITS files in the given directory to identify and select target centers without WCS information.

    Parameters:
    path (str): The directory path containing the FITS files.
    source (str): The source name to fetch RA and DEC coordinates from the `mopinf.source_info` dictionary.

    Description:
    This function performs the following steps:
    1. Changes the current directory to the specified path.
    2. Identifies FITS files matching a specific pattern based on the `mopset.time_frame`.
    3. Filters FITS files by the filter wavelength ('B', 'V', 'R', 'I', 'L').
    4. Extracts run identifiers from the filenames.
    5. For each filtered FITS file, reads the image data and header information.
    6. If WCS error is zero, calculates pixel coordinates for the source and calibration star based on WCS data.
    7. If WCS error is non-zero, displays the image and prompts the user to select the source center manually.
    8. Saves the selected pixel coordinates to a CSV file for future reference.

    Notes:
    - Requires `astropy`, `matplotlib`, and `numpy` for handling FITS files, WCS, and image display.
    - Assumes that the `mopinf` and `mopset` modules provide necessary source information and directory paths.
    - Will filter out observations with no WCS error and will reduce as normal.
    - If the pointing between filters is accurate, it is possible to only manually select the source centers once and share coords.
    The user must comment out the code which appends the files based on filter to a list.
    In the photometry functions, the user must also input which filter was used to get the coords (more details in the photometry functions).

    """

    os.chdir(path)
    files = glob.glob('*_*_'+mopset.time_frame+'_*_1_1.fits')
    coord_frames = []
    for file in files:
        hdulist = fits.open(file)
        header = hdulist[0].header
        wave = header['FILTER1'][4:]
        # Comment out appending certain filters if sharing the same coords.
        # Requires accurate pointing.
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
        try: 
            cal_ra = mopinf.source_info[source]['cal_ra']
            cal_dec = mopinf.source_info[source]['cal_dec']
        except:
            cal_ra = ra
            cal_dec = dec

        mean, median, std = sigma_clipped_stats(data,sigma=4.0)
        data -= median

        if wcs_err == 0:
            star_coord = SkyCoord(ra,dec, unit=(u.hourangle, u.deg))
            xy = SkyCoord.to_pixel(star_coord, wcs=wcs, origin=1)
            x = xy[0].item(0)
            y = xy[1].item(0) 

            star_coord_cal = SkyCoord(cal_ra,cal_dec, unit=(u.hourangle, u.deg))
            xy_cal = SkyCoord.to_pixel(star_coord_cal, wcs=wcs, origin=1)
            x_cal = xy_cal[0].item(0)
            y_cal = xy_cal[1].item(0)
            globals()['pix_coords_'+cam_date_run] =[(x_cal,y_cal),(x,y)]

        else:
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

                print('Right-click on the cal star center')

                click_coords_cal = plt.ginput(1, timeout=-1, mouse_add=3, mouse_pop=2)
                x_cal, y_cal = click_coords_cal[0]
                print(f'x_cal={x_cal:.2f}, y_cal={y_cal:.2f}')

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
    """
    Determines the optimum aperture size for photometry based on the Signal-to-Noise Ratio (SNR).

    Parameters:
    data (numpy.ndarray): The 2D array of the image data.
    rough_position (tuple): The approximate (x, y) position of the target in the image.
    obj (str): The type of object ('src' for source, 'cal' for calibration star).
    camera (int): The camera number (1 or 3 indicates special handling for sources).
    std (float): The standard deviation of the background.
    gain (float): The gain of the camera.

    Returns:
    tuple: A tuple containing:
        - pos_center (tuple): The (x, y) coordinates of the center of the target.
        - optimum_radius (float): The optimum aperture radius for photometry.
        - optimum_annulus_radius (tuple): The optimum inner and outer radii factors for the annulus.

    Description:
    This function performs the following steps:
    1. Identifies the approximate position of the target in the image and refines it using centroid calculation.
    2. If the object is a source (`obj == 'src'`) and the camera is 1 or 3, it calculates the optimum aperture size for photometry by maximizing the Signal-to-Noise Ratio (SNR).
    3. Determines the best annulus for background subtraction by evaluating different inner and outer radius factors and selecting the combination that maximizes SNR.
    4. Returns the refined position of the target, the optimum aperture radius, and the optimum annulus radii factors.

    Notes:
    - The function utilizes the `CircularAperture` and `CircularAnnulus` classes from the `photutils` library for photometry.
    - Plots showing the SNR vs. aperture size and SNR vs. inner and outer factors are commented out but can be enabled for visualization.

    """

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

def subtract_background(data, radius_factor=0.98, sigma=4.0, fwhm=3, threshold_factor=5.0, source_radius=5):
    """
    Subtracts background from image data using sigma-clipped statistics and star detection.

    Parameters:
    data (numpy.ndarray): The 2D array of the image data.
    radius_factor (float, optional): Factor determining the radius of the circular mask for background subtraction to remove vingetting effects. Defaults to 0.98.
    sigma (float, optional): Sigma value for sigma-clipping during background estimation. Defaults to 4.0.
    fwhm (float, optional): Full-Width at Half-Maximum (FWHM) for star detection. Defaults to 3.
    threshold_factor (float, optional): Factor multiplied by standard deviation for star detection threshold. Defaults to 5.0.
    source_radius (int, optional): Radius of the circular mask around detected sources to exclude from background calculation. Defaults to 5 pixels.

    Returns:
    tuple: A tuple containing:
        - data (numpy.ndarray): The image data with background subtracted.
        - background (float): The estimated background value subtracted from the data.
        - std (float): The standard deviation of the background.

    Description:
    This function performs the following steps:
    1. Creates a circular mask to exclude central regions from background estimation.
    2. Uses sigma-clipped statistics to estimate the background and subtract it from the data.
    3. Detects sources in the image using DAOStarFinder with specified parameters (FWHM and threshold).
    4. Creates a source mask to exclude detected sources from background calculation.
    5. Updates the circular mask to include the source mask and performs sigma-clipped statistics again to refine background estimation.
    6. Subtracts the refined background from the data.

    Notes:
    - Utilizes numpy and astropy libraries for array manipulation and statistical calculations.
    - If an exception occurs during background subtraction or source detection, it defaults to sigma-clipped background subtraction using the median.

    """

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
            radius = source_radius
            circular_mask = create_circular_mask(h, w, center=center, radius=radius)
            source_mask &= ~circular_mask

        return source_mask

    h, w = data.shape
    circular_mask = ~create_circular_mask(h, w, radius=int(512 * radius_factor))

    try:
        mean, median, std = sigma_clipped_stats(data, sigma=sigma, mask=None)
        daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold_factor * std, sigma_radius=sigma)
        sources = daofind(data - median, mask=circular_mask)

        source_mask = create_source_mask(h, w, sources)
        mask = circular_mask | ~source_mask

        mean, median, std = sigma_clipped_stats(data, sigma=sigma, mask=mask)
        data -= mean
        background = mean

    except Exception as e:
        print(f"Exception occurred: {e}")
        mean, median, std = sigma_clipped_stats(data, sigma=sigma)
        data -= median
        background = median

    return data, background, std

def one_cam_photometry(filename, source):
    """
    Perform photometry on FITS files for a single camera setup.

    Parameters:
    filename (str): Base filename pattern for FITS files.
    source (str): Source identifier for extracting coordinate information.

    Returns:
    tuple: A tuple containing:
        - date (str): Date extracted from FITS headers.
        - mjd (float): Modified Julian Date extracted from FITS headers.
        - date_time (datetime.datetime): Date and time of the first FITS file.
        - wave (str): Filter waveband extracted from FITS headers.
        - aperture_size (float): Size of the photometric aperture in arcseconds.
        - rotskypa (float): Rotation sky position angle (ROTSKYPA) from FITS headers.
        - rotang (float): Rotation angle (ROTANGLE) from FITS headers.
        - counts_src_1 (list): Photometric counts of the source in camera 1.
        - counts_src_1_err (list): Photometric errors of the source in camera 1.
        - counts_cal_1 (list): Photometric counts of the calibration star in camera 1.
        - counts_cal_1_err (list): Photometric errors of the calibration star in camera 1.
        - gain (float): Gain value from FITS headers.
        - wcs_err (float): WCS error value from FITS headers.
        - seeing (float): Seeing value (L1SEESEC) from FITS headers.
        - pho (float): Scheduling photometric value from FITS headers.
        - background (float): Estimated background value from background subtraction.
        - x (float): X-coordinate of the source position.
        - y (float): Y-coordinate of the source position.

    Description:
    This function iterates over a series of FITS files generated for a single camera setup.
    It performs the following steps for each file:
    1. Opens the FITS file and extracts necessary header information.
    2. Performs background subtraction on the image data using `subtract_background` function.
    3. Calculates and plots apertures around the source and calibration star positions.
    4. Performs photometry on both the source and calibration star using circular apertures.
    5. Appends photometric results (counts and errors) to respective lists based on camera type.
    6. Collects other relevant header values such as gain, WCS error, seeing, and scheduling photometry.
    7. Returns various extracted and calculated values including date, MJD, time, waveband, aperture size,
       sky position angle, rotation angle, background value, and source coordinates.

    Notes:
    - Requires `numpy`, `matplotlib`, `astropy`, and other related libraries.

    """

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
            try:
                cal_ra = mopinf.source_info[source]['cal_ra']
                cal_dec = mopinf.source_info[source]['cal_dec']
            except:
                cal_ra = ra
                cal_dec = dec
            
            #Background subtraction
            data, background, std = subtract_background(data)

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

    return date,mjd,date_time[0],wave,aperture_size,rotskypa,rotang,counts_src_1,counts_src_1_err,counts_cal_1,counts_cal_1_err,gain,wcs_err,seeing,pho,background,x,y

def one_cam_polarimetry(filename,source):
    """
    Perform polarimetry calculations based on photometric measurements.

    Parameters:
    filename (str): Base filename pattern for FITS files.
    source (str): Source identifier for extracting coordinate information.

    Returns:
    tuple: A tuple containing:
        - date (str): Date extracted from FITS headers.
        - mjd (float): Modified Julian Date extracted from FITS headers.
        - date_time (datetime.datetime): Date and time of the first FITS file.
        - wave (str): Filter waveband extracted from FITS headers.
        - aperture_size (float): Size of the photometric aperture in arcseconds.
        - gain (float): Gain value from FITS headers.
        - rotskypa (float): Rotation sky position angle (ROTSKYPA) from FITS headers.
        - rotang (float): Rotation angle (ROTANGLE) from FITS headers.
        - wcs_err (float): WCS error value from FITS headers.
        - seeing (float): Seeing value (L1SEESEC) from FITS headers.
        - pho (float): Scheduling photometric value from FITS headers.
        - s1_src (float): Average photometric count of the source.
        - s1_src_err (float): Error in the average photometric count of the source.
        - inst_src (float): Instrumental magnitude of the source.
        - inst_src_err (float): Error in the instrumental magnitude of the source.
        - s1_cal (float): Average photometric count of the calibration star.
        - s1_cal_err (float): Error in the average photometric count of the calibration star.
        - inst_cal (float): Instrumental magnitude of the calibration star.
        - inst_cal_err (float): Error in the instrumental magnitude of the calibration star.
        - q_avg (float): Average Stokes Q parameter.
        - q_err (float): Error in the average Stokes Q parameter.
        - u_avg (float): Average Stokes U parameter.
        - u_err (float): Error in the average Stokes U parameter.
        - bkg (float): Background value from background subtraction.
        - x (float): X-coordinate of the source position.
        - y (float): Y-coordinate of the source position.

    Description:
    This function computes polarimetric parameters (Stokes Q and U) based on
    photometric measurements of a source and a calibration star from FITS files.
    It retrieves photometric counts and errors from `one_cam_photometry` function output
    and performs necessary calculations for polarization measurements.

    Notes:
    - Requires `numpy`, `math`, and other related libraries.
    - Relies on `one_cam_photometry` function for initial photometric measurements.
    - Calculates instrumental magnitudes and polarimetric parameters based on the provided formulas.
    """

    date,mjd,date_time,wave,aperture_size,rotskypa,rotang,m1_src,m1_src_err,m1_cal,m1_cal_err,gain,wcs_err,seeing,pho,bkg,x,y = one_cam_photometry(filename,source)

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

def two_cam_photometry(filename, source):
    """
    Perform aperture photometry on FITS image data for MOPTOP two-camera setup.

    Parameters:
    filename (str): Base filename of the FITS images (excluding camera and sequence numbers).
    source (str): Name of the astronomical source for which photometry is being performed.

    Returns:
    tuple: A tuple containing the following elements:
        - date (str): Date of observation.
        - mjd (float): Modified Julian Date of observation.
        - date_time (str): Date and time of observation in ISO format.
        - wave (str): Filter used during observation.
        - aperture_size (str): Size of the aperture used for photometry.
        - rotsky_pa (float): Position angle of the celestial north pole on the sky.
        - rotangle (float): Rotational angle of the instrument.
        - counts_src_1 (list): List of counts for the source in camera 1.
        - counts_src_1_err (list): List of errors in counts for the source in camera 1.
        - counts_src_2 (list): List of counts for the source in camera 2.
        - counts_src_2_err (list): List of errors in counts for the source in camera 2.
        - counts_cal_1 (list): List of counts for the calibrator in camera 1.
        - counts_cal_1_err (list): List of errors in counts for the calibrator in camera 1.
        - counts_cal_2 (list): List of counts for the calibrator in camera 2.
        - counts_cal_2_err (list): List of errors in counts for the calibrator in camera 2.
        - gain (float): Gain value from FITS header.
        - wcs_err (float): Error in the World Coordinate System (WCS) from FITS header.
        - seeing (float): Estimated seeing conditions during observation.
        - pho (float): Photometric zero point from FITS header.
        - background (float): Background level subtracted from the data.
        - src_x (float): X-coordinate of the source position in camera 1.
        - src_y (float): Y-coordinate of the source position in camera 1.

    Description:
    This function performs aperture photometry on FITS image data from two cameras for a specific astronomical source.
    It iterates over multiple images (four cameras, each with multiple sequence numbers) to perform the following steps:
    1. Opens each FITS file and extracts relevant header information (gain, WCS, rotational angles, filter, date, time, etc.).
    2. Subtracts background from the data using a built-in function.
    3. Converts celestial coordinates (RA, Dec) of the source and calibrator to pixel coordinates for both cameras.
    4. Determines optimal aperture size and annulus for photometry based on source and instrument characteristics.
    5. Performs aperture photometry on both source and calibrator using circular apertures and annuli.
    6. Collects and stores photometry results (counts and errors) for both cameras and both source and calibrator.
    7. Returns a tuple containing all collected data for further analysis.

    Notes:
    - Error handling is implemented to skip files that cannot be opened or do not contain required header information.
    - Plotting functionality (`mopset.plot`) is optionally available to visualize aperture positions and sizes for validation. Enable in moptop_settings.py.

    """

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
            try:
                cal_ra = mopinf.source_info[source]['cal_ra']
                cal_dec = mopinf.source_info[source]['cal_dec']
            except:
                cal_ra = ra
                cal_dec = dec

            #Background subtraction
            data, background, std = subtract_background(data)

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

            # "factor" allows the user to use the coordinates obtained with a different filters data.
            # MOPTOP observations are always carried out in BVRIL order, so the run number is a successive sequence.
            # Use factor to adjust the run number to the correct filter.
            # Examples:
            # To use R band coordinates for B band data, set B factor = 2
            # To use B band coordinates for I band data, set I factor = -3
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

def two_cam_polarimetry(filename, source):
    """
    Extracts relevant photometry and polarimetry data from a two-camera system.

    Parameters:
    filename (str): File name or path containing photometry data.
    source (str): Source identifier for which data is extracted.

    Returns:
    tuple: A tuple containing various extracted data fields:
        - date (str): Date of observation.
        - mjd (float): Modified Julian Date of observation.
        - date_time (str): Date and time of observation.
        - wave (float): Wavelength of observation.
        - aperture_size (float): Aperture size used in photometry.
        - gain (float): Gain used in photometry.
        - rotskypa (float): Rotational sky position angle.
        - rotang (float): Rotational angle.
        - wcs_err (float): Error in World Coordinate System.
        - seeing (float): Seeing conditions during observation.
        - pho (float): Photometric value.
        - counts_src (float): Average counts from source photometry.
        - counts_src_err (float): Error in average counts from source photometry.
        - inst_src (float): Instrumental magnitude of the source.
        - inst_src_err (float): Error in instrumental magnitude of the source.
        - counts_cal (float): Average counts from calibration photometry.
        - counts_cal_err (float): Error in average counts from calibration photometry.
        - inst_cal (float): Instrumental magnitude of the calibration source.
        - inst_cal_err (float): Error in instrumental magnitude of the calibration source.
        - q_avg (float): Average Stokes Q parameter.
        - q_err (float): Error in Stokes Q parameter.
        - u_avg (float): Average Stokes U parameter.
        - u_err (float): Error in Stokes U parameter.
        - bkg (float): Background level.
        - x (float): X-coordinate of the source.
        - y (float): Y-coordinate of the source.

    Description:
    This function extracts photometry and polarimetry data from a two-camera system file.
    It calculates average counts, instrumental magnitudes, and errors for both source and calibration data.
    It also computes average Stokes Q and U parameters along with their errors.

    Notes:
    - The function relies on the `two_cam_photometry` function to retrieve photometry data.
    - It performs calculations for average counts, instrumental magnitudes, and Stokes parameters.
    - Error handling ensures that calculations proceed even if individual data points are missing or invalid.
    """

    #Extracts relavent photometry data 
    date,mjd,date_time,wave,aperture_size,rotskypa,rotang,m1_src,m1_src_err,m2_src,m2_src_err,m1_cal,m1_cal_err,m2_cal,m2_cal_err,gain,wcs_err,seeing,pho,bkg,x,y = two_cam_photometry(filename,source)

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
    """
    Removes all files with a specified file extension within a given root folder and its subdirectories.

    Parameters:
    root_folder (str): The root directory path from which to start searching for files.
    file_extension (str, optional): File extension to filter files for deletion (default is '.fits').

    Description:
    This function iterates through all subdirectories starting from `root_folder`.
    For each file found with the specified `file_extension`, it attempts to remove the file.
    If deletion fails due to permission issues or other exceptions, it prints an error message.

    Notes:
    - The function only removes files that exactly match the specified file extension.
    - It recursively searches through all subdirectories of `root_folder`.
    - Errors encountered during file deletion are caught and printed to the console.
    - Be cautious when using this function as it irreversibly deletes files.

    """

    for folder_path, _, file_names in os.walk(root_folder):
        for file_name in file_names:
            if file_name.endswith(file_extension):
                file_path = os.path.join(folder_path, file_name)
                try:
                    os.remove(file_path)
                except Exception as e:
                    print(f"Error removing {file_path}: {e}")

def novel_pol_error(q_i, q_err, u_i, u_err):
    """
    Calculates polarization properties and error bounds according to Plaszczynski et al., 2014, MNRAS, 439, 4048, doi:10.1093/mnras/stu270.

    Parameters:
    q_i (float): Stokes Q parameter.
    q_err (float): Error in Stokes Q parameter.
    u_i (float): Stokes U parameter.
    u_err (float): Error in Stokes U parameter.

    Returns:
    tuple: Tuple containing polarization magnitude (p_mas), minimum polarization estimate (p_min),
           and maximum polarization estimate (p_max).

    Description:
    This function computes the polarization magnitude and associated error bounds using a novel statistical approach.
    It performs the following steps:
    1. Calculates the polarization magnitude (p_i) from Stokes Q and U parameters.
    2. Computes the maximum polarization estimate (p_max) incorporating statistical enhancement factors.
    3. Computes the minimum polarization estimate (p_min) considering both statistical and systematic errors.

    Notes:
    - The function uses mathematical operations and exponential functions to adjust polarization estimates.

    """

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

    return p_mas,p_min,p_max

def EVPA(q_i, q_err, u_i, u_err):
    """
    Calculates the Electric Vector Position Angle (EVPA) from Stokes parameters Q and U.

    Parameters:
    q_i (float): Stokes Q parameter.
    q_err (float): Error in Stokes Q parameter.
    u_i (float): Stokes U parameter.
    u_err (float): Error in Stokes U parameter.

    Returns:
    tuple: A tuple containing EVPA in degrees and its associated error.

    Description:
    This function calculates the Electric Vector Position Angle (EVPA) using Stokes parameters Q and U:
    - EVPA (theta) is calculated as 0.5 * arctan2(U, Q) converted to degrees, modulo 180.
    - Error in EVPA (theta_err) is computed using error propagation formula:
      theta_err = degrees(sqrt((U**2 * q_err**2 + Q**2 * u_err**2) / (2 * (Q**2 + U**2))))

    Notes:
    - The function uses numpy functions for arctan2, square root, and degrees conversion.
    - Modulo 180 ensures the EVPA is returned within the range [0, 180) degrees.
    - Errors in Q and U parameters are required to compute the error in EVPA.

    """

    #Calculation of theta for evpa using q and u.
    theta = np.degrees(0.5*np.arctan2(u_i,q_i)) % 180
    theta_err = np.degrees(np.sqrt(u_i**2*q_err**2+q_i**2*u_err**2)/(2*(q_i**2+u_i**2)))
    return theta,theta_err

def evpa_unwrap(evpa, evpa_err, n):
    """
    Unwraps the EVPA (Electric Vector Position Angle) values based on their errors and neighboring points.

    Parameters:
    evpa (list or pandas Series): List of EVPA values.
    evpa_err (list or pandas Series): List of errors associated with EVPA values.
    n (int): Number of neighboring points to consider for unwrapping.

    Returns:
    evpa_cor (pandas Series): EVPA values after unwrapping, reset to index 0.

    Description:
    This function iterates through the EVPA values starting from index `n` and performs the following steps:
    1. Calculates a weighted average of EVPA values from the previous `n` points, considering errors in EVPA values.
    2. Determines the difference between the current EVPA value and the weighted average.
    3. Adjusts the current EVPA value by multiples of 180 degrees based on the difference and associated error.
    4. Returns the corrected EVPA values as a pandas Series, reset to index 0.

    Notes:
    - The function adjusts EVPA values to avoid abrupt jumps greater than 90 degrees.
    - It uses a weighted average approach considering errors in EVPA values for smoothing.

    """

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
            evpa.loc[i] += N*180
        if abs(diff) - err > 90 and diff > 0:
            evpa.loc[i] -= N*180

    evpa_cor = evpa.reset_index(drop=True)
    return evpa_cor
