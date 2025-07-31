<h1 align="center">
  <br>
MOPTOP reduction pipeline
  <br>
</h1>

---

<p align="center">
  <a href="#introduction">Introduction</a> |
  <a href="#file-summary">File summary</a> |
  <a href="#instructions">Instructions</a> |
  <a href="#acknowledgments">Acknowledgments</a> |
  <a href="#contact">Contact</a>
</p>

---

## Introduction

The code provides a Python pipeline for downloading, calibrating, and reducing photopolarimetric data taken with the MOPTOP polarimeter on the Liverpool Telescope.  

Please cite McCall (2025) (https://researchonline.ljmu.ac.uk/id/eprint/26800/) if using this pipeline or the work within it. 

### How to use

Clone or download the repo. Both can be done by clicking the "code" button found at the top of the repo. To clone the repo using Git, run the command "git clone URL", where URL is the repo HTTPS URL, on your machine. To download a local version, click "Download ZIP". 

The pipeline was developed using Visual Studio Code (https://code.visualstudio.com) so this is the recommended environment.

## File summary

The main pipeline is made up of four main scripts: moptop_settings.py, moptop_reduction.py, moptop_functions.py, and moptop_dicts.py. There is also an additional script, moptop_constants.py. The contents and purpose of these files can be summarised as follows:

### moptop_settings.py

This is the main configuration file for the pipeline and, once set up, should generally be the only file that needs editing run-to-run. This file is where key variables are defined or optional functions are toggled on/off. Descriptions of what each variable or toggle is for can be found as a comment within the file. 

### moptop_reduction.py

This is the main pipeline script, reading in required functions and dictionary parameters based on the configured settings. Details will be given in the Instructions section, but in summary, this script works as follows assuming the required toggles are on:

- Download and organise fits files on your machine.

- Rotation stack frames.

- Reduce frames, determining if single or dual camera methods are required.

- Calibrate the photopolarimetric data, including EVPA unwrapping.

### moptop_functions.py

This file contains all the functions used in the pipeline. These have been taken from the main pipeline body for clarity and aesthetics. Each function is documented in detail.

To note, these functions should not generally need editing and any errors likely arise from setting and dictionary configuration. The following functions may require edits:

- The no_wcs function can have filters commented out if manually determined (using the code) pixel coordinates are to be shared. This is not recommended unless tracking is perfect.

- The optimum_aperture functions contains commented out figures (ap_size vs. SNR and SNR vs. inner/outer annulus radii factor) which can be uncommented if needed. 

- In the photometry functions (one_cam_photometry and two_cam_photometry), the first line starts a for loop in the range of 1-4 (1-5 in Python syntax). These correspond to the possible MOPTOP camera names. If any new cameras are installed in MOPTOP, as in March 2022, these new cameras may be called cameras 5, 6 etc. The upper range limit will need need to be changed to match the name of any new camera. 

### moptop_dicts.py

This is where key parameters are stored in a series of dictionaries. Information relating to sources, including source and calibration star coordinates, calibration star magnitudes, aperture sizes, and polarisation constants. To note, polarisation constants are only required for polarimetric calibration sources when reducing their data. 

The other two dictionaries in this file contain epoch/waveband specific polarimetric calibration values and waveband specific flux calibration values. 

### moptop_constants.py

While this script uses the moptop_settings.py and moptop_dicts.py files, it is not generally required when using the main reduction pipeline. The purpose of this script is to calculate the MOPTOP polarisation coefficients after MOPTOP-altering telescope maintenance (primary/secondary/tertiary mirror cleans, camera changes, re-installing the instrument etc). This script produces the epoch specific q_0, u_0, EVPA k, and instrumental depolarisation values with associated standard error.

The input for this file is the final output of the main pipeline, so one is required to do a "dummy" run of the main pipeline for polarisation calibration data. This run obtains the required uncorrected q and u parameters for the polarisation calibration pipeline and potentially non-sensical "calibrated" values in the case where true coefficients have not yet been determined. These can be ignored.

The main pipeline does not take into account the uncertainties on these constants, but it is recommended to produce values at the start of each epoch until the respective outputs are consistent within uncertainties. This may take several weeks of standard star observations which are obtained through the Liverpool Telescope operations MOPStand proposal. The values published on the Liverpool Telescope website (https://telescope.astro.ljmu.ac.uk/TelInst/Inst/MOPTOP/) originate from this work. 

---

## Instructions 

The only files that need configuring are the moptop_settings.py and moptop_dicts.py files. 

### moptop_dicts.py 

The following steps are compulsory to making the pipeline work. 

Enter any new object parameters into the nested source_info dictionary. The framework can be copied from an existing object, such as BL Lac, and this source will be the example in this description. 

The first entry is an alt_names list. This list should comprise any and all names this object takes in the OBJECT header of the fits files. Different astronomers will call objects by different names, so the purpose of this list to tell the script that frames of different object names are the same source (i.e. BL Lac = bl lac = BL LAC etc). The name of nested dictionary, BL Lac in this case, will be the "master name" which the folders the script creates will be named after.

The next three fields are estimates for the aperture size and inner/outer annulus radii. The user will have the option to have a function refine these values so they need not be accurate.

The remaining fields are source and calibration star RA/Dec coordinates, as well as the BVRIL magnitudes of the calibration star. If calibrated magnitude is not of interest, enter 0 for the magnitude values. The script will use the source coordinates if the calibration star coordinates are left blank.

If new MOPTOP polarisation constants have been calculated following MOPTOP-altering telescope maintenance, these must be added to the polarimetric_constants dictionary. From the constants in the dictionary already, one can see the general naming convention pattern which must be followed. This convention is coefficient_filter_epoch where coefficient is one of q_zero, u_zero, k, or inst_depol, filter is one of B,V,R,I,L, and epoch is an increasing integer. 

### moptop_settings.py

Each setting in this file must be completed carefully in order to produce expected results.

The first setting, dir1, is main directory path to where the individual fits files will be or are already stored. This pipeline is designed to download and open, or just open, fits files located in this directory and move them to a newly created folder of name matching that in the source_info dictionary. Subsequent CSV files will also be saved in this location. It is advised to allow the pipeline to create the filing structure in order to save fits files and CSVs in places the code is set up to find.  

change_pol_constants is a list of MJDs which tells the pipeline which set of polarisations coefficients to use. If new ones are added to the polarimetric_constants dictionary, the MJD of when they take effect must be added to the end of this list. 

Similarly, single_camera_start and single_camera_end are two lists of MJDs which tell the pipeline when to use the single and dual camera reduction methods. Should data be collected when MOPTOP only have one working camera, an MJD must be added to both lists. It is advised to set the end MJD to an arbitrarily distance value and correct it when the second camera is back online. 

The sources list is a list of source names to reduce. These names must be as they appear in the source_info dictionary. This list is only needed when NOT downloading new data. It will be overwritten if downloading new data from the archive.

no_wcs_sources is a list of sources which likely do not have a successful WCS fit, triggering the manual coordinate function. This function will check if a WCS fit has been unsuccessful before displaying a GUI allowing one to right click the source, followed by the calibration star, obtaining pixel coordinates. Only the first frame from each (rotation stacked) camera/filter will be shown as it is assumed the sources do not drift across the frame during an observation. The determined pixel coordinates are saved in an outputted CSV and read again to avoid repeating this process.

The download_data toggle takes a "y" or "n" input and tells the pipeline whether to download new data from the Liverpool Telescope recent data webpage. If set to "y", at least one date (string; YYYYMMDD) must be given in the download_dates list, as well as at least one set of proposal IDs and passwords in the proposal dictionary. Multiple of each can be specified.

The time_frame parameter allows one to specify specific observations to reduce. For example, if observations for 20250101 and 20150102 existed in the directory, one could set this parameter to "20250101_*" to reduce all the data from 20250101 and not 20250102 using the wildcard, *, character. Furthermore, if two observation runs took place on 20250101 - run 14 and 15 (e.g. different filters) - but only one needed to be reduced (run 14), one could specify this using "20250101_14".

The stacking toggle will allowing the user to switch on and off the rotation stacking function should this not be needed (i.e. with GRB observations a higher time resolution may be required). The stacking works by finding frames with the same date and run number but different rotation numbers. The median of these frames in corresponding waveplate positions is calculated and saved as a new file. The unstacked files are removed. This is because the new file has the same name as those frames of rotation 1 and so one does not want to accidentally stack them again.

The reduce toggle switches on and off the main reduction process. The pipeline will automatically choose the single camera vs. dual camera method. This function performs the raw photometry and baseline polarimetry, resulting in instrumental magnitudes and un-corrected q-u values. One may want to switch this off if changes to subsequent functions are introduced (e.g. new polarisation constants for an already reduced epoch). The results from this process are outputted to a CSV called "reduced_data.csv".

The plot toggle switches on and off a display of the first waveplate position frame of each camera and filter with aperture and background annulus displayed. This is to check the positioning is correct.

The optimum_aperture toggle switches on and off the optimum aperture function. This function performs object detection and refines the user inputted RA/Dec coordinates by looking for the closest source. It then tries a series of aperture/annulus sizes, selecting the configuration with the lowest SNR. The function does this for the first waveplate position frame for each camera/filter/rotation number. This function will work with those manually selected coordinates for frames with poor WCS fitting.

The replace toggle switches on and off the function that allows the user to specify whether they want to re-reduce data that has already been done before. 

The clear_dir1 toggle switches on and off the function that removes all fits files from dir1 and all child directories within.  

The calculations toggle switches on and off the calibration of the data within "reduced_data.csv". This process calibrates magnitude and q-u values, as well as calculating flux, linear polarisation degree, and unwrapped EVPA. The results of this process are outputted to a new CSV called "all_data.csv". This is the final output of the pipeline.

One can specify the number of EVPA values to include in the EVPA weighted average using the evpa_steps parameter. This is set to 3, but can take the form of any integer. This pipeline uses a novel approach to unwrap EVPA data, details of which can be found in McCall (2025).

---

## Acknowledgments 

Thank you to Manisha Shrestha, Helen Jermak, Iain Steele for assisting in the starting and continuous development of this pipeline.

Thank you to K-Ryan Hinds for sharing a web-scraping script to automate the data download process from the Liverpool Telescope recent data webpage.

Thank you to James Pryce for helping to QA and troubleshoot the early version of this pipeline.

---

## Contact

Email:

c.mccall@2017.ljmu.ac.uk

callumlmccall@gmail.com
