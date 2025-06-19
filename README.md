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

The first setting, dir1, is main directory path to where the individual fits files will be or are already stored. This pipeline is designed to download and open, or just open, fits files located in this directory and move them to a newly created folder of name matching that in the source_info dictionary. It is advised you allow the pipeline to create the filing structure in order to not save fits files or CSVs in places the code cannot find.  

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