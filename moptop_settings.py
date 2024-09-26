"""
CHECK APERTURE SIZE FITS WELL BEOFORE REDUCING FULLY.
"""

dir1 = '/Insert/File/Path/Here/' 
#Path to where fits files are stored or where they are to be downloaded to. 
#Files will be organised based on the "OBJECT" header
#Final and intermediate csv files are saved in the crerated source directories.

change_pol_constants = [59477, 59654, 60561] #List of MJD's when the polarisation constants changed
single_camera_start = [59477] #List of MJD's when MOPTOP started operating with a single camera
single_camera_end = [59654] #List of MJD's when MOPTOP resumed operating with both cameras
                            #If single camera mode is ongoing, enter a future MJD for the end date

sources = ['BL Lac','OJ287','3C 454.3','4C 11.69','TXS 0506+056','PKS1510-089','PG1553+113','3C 345']
            #List of sources to reduce. Only applicable if download_data = 'n'
            #Source names as they appear in source info dict

no_wcs_sources = ['MRK421','SN2023ixf'] #List of sources that frequently do not have WCS coordinates in the header.
                            #User will be asked to manually select the target and calibration star using R band image.

download_data = 'y' #Download data from "Recent Data" (y/n)
download_dates = ['yyyymmdd'] #Dates to download data for (YYYYMMDD)
proposals = [{'ID':['PASSWORD']}] #Proposal iDs and passwords

time_frame = '*'    #If you want to reduce specific files, enter the fits file date and run number here.
                    #Format should be 'yyyymmdd_runnumber'. '*' can be used as a wildcard.
                    #If you want to reduce all files, enter '*'

stacking = 'y'  #This will check all files for any that should be rotation stacked (same run number but different rotation number in filenames).
                #Will also check if any files have data in one camera but not the other, and reduce in single camera mode after stacking.

reduce = 'y' #Performs the data reduction (y/n)
plot = 'n' #Whether or not to display a finding chart plot in the reduction stage (y/n)
optimum_aperture = 'y' #Whether or not to find the optimum aperture/annulus size for the source (y/n)
replace = 'y' #If a file has already been reduced, whether or not to replace the old reduction with a new one (y/n)
clear_dir1 = 'n' #Whether or not to clear the directory where the fits files are stored after they have been reduced (y/n)
calculations = 'y' #Whether or not to perform the calculations to obtain calibrated values (y/n)
evpa_steps = 3 #Number of previous points to be used in the weighted average un unwrapping the EVPA


