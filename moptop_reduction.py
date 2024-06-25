import moptop_settings as mopset
import moptop_functions as mopfunc
import moptop_dicts as mopinf
import glob
import pandas as pd
import numpy as np
import csv 
import os
from pylab import *
from IPython.display import clear_output
from astropy.io import fits
import re

"""
Reduction loops, csv creation and vetting.
Tables of the final amounts of data will also be displayed.

THIS SECTION SHOULD NOT NEED EDITING.
"""

if mopset.download_data == 'y':
    sources = mopfunc.down_recentdata(mopset.dir1,mopset.download_dates)
else:
    sources = mopset.sources

mopfunc.sort_files(mopset.dir1)

if mopset.stacking == 'y':
    for source in sources:
        try:
            # Rotation stacking
            os.chdir(mopset.dir1+source+'/')
            num_all_files = len(glob.glob('*'+mopset.time_frame+'*.fits'))
            file_list = glob.glob('*'+mopset.time_frame+'*_1_1.fits')

            rot_dates = []
            for file in file_list:
                fil = file[1:-11]
                rot_dates.append(fil)
            rot_dates = pd.Series(rot_dates,dtype=str).unique()
            
            #Stacking
            tot = 0
            globals()['rot_stack_dates_'+source] = []
            clear_output(wait=True)
            for i in range(len(rot_dates)):
                for cam in range(1,5):
                    for rotor_pos in range(1,17):
                        files = []
                        files = glob.glob(str(cam)+rot_dates[i]+'_*_'+str(rotor_pos)+'_1.fits')
                        tot += len(files)
                        if len(files) == 0:
                            continue
                        if len(files) == 1:
                            print(f"Rotation Stacking: {(tot/(num_all_files))*100:.2f} % {source}")
                            continue
                        else:
                            date = re.sub('_','',rot_dates[i][3:-2])
                            globals()['rot_stack_dates_'+source].append(date)
                            for n in range(len(files)):
                                globals()['file'+str(n)] = files[n]
                            image_concat = []
                            for n in range(len(files)):
                                globals()['file'+str(n)+'_path'] = mopset.dir1+source+'/'+globals()['file'+str(n)]
                                globals()['data'+str(n)] = fits.getdata(globals()['file'+str(n)],0)
                                image_concat.append(globals()['data'+str(n)])
                            final_image_med = np.median(image_concat,axis=0)

                            #Assigns fits data from first rotation 
                            hdulist = fits.open(globals()['file0'])
                            header = hdulist[0].header

                            #Saves new stacked image
                            for file in files:
                                os.remove(mopset.dir1+source+'/'+file)
                            with open(mopset.dir1+source+'/'+str(cam)+rot_dates[i]+'_1_'+str(rotor_pos)+'_1.fits', 'wb') as output:
                                fits.writeto(output,final_image_med,header=header,overwrite=True)
                            print(f"Rotation Stacking: {(tot/(num_all_files))*100:.2f} % {source}")
        except:
            continue

if mopset.reduce == 'y':
    for source in sources:
        if source in mopset.no_wcs_sources:
            mopfunc.no_wcs(mopset.dir1+source+'/',source)
        try:
            os.chdir(mopset.dir1+source+'/')

            processed_datetimes = []
            if os.path.exists('reduced_data.csv'):
                with open('reduced_data.csv', mode='r') as existing_data_file:
                    existing_data = csv.DictReader(existing_data_file)
                    for row in existing_data:
                        existing_datetime = row['datetime']
                        existing_dtime = datetime.datetime.strptime(existing_datetime, "%Y-%m-%d %H:%M:%S.%f")
                        processed_datetimes.append(existing_dtime)

            tot = 0
            clear_output(wait=True)
            file_list = glob.glob('[1,3]_*'+mopset.time_frame+'*_1_1.fits')
            for file in file_list:
                hdulist = fits.open(file)
                header = hdulist[0].header
                datetime_string = f"{header['DATE']} {header['UTSTART']}"
                dtime = datetime.datetime.strptime(datetime_string, "%Y-%m-%d %H:%M:%S.%f")
                if mopset.replace == 'y':
                    if os.path.exists(mopset.dir1+source+'/reduced_data.csv'):
                        processed_datetimes = [dt for dt in processed_datetimes if dt != dtime]
                        mopfunc.delete_row(mopset.dir1+source+'/reduced_data.csv', dtime)
                if dtime in processed_datetimes:
                    tot += 1
                    print('Already Reduced:',np.round((tot)/len(file_list)*100,2),'%',source)
                    pass
                else:
                    mjd = header['MJD']
                    method = 'dual_camera'
                    for i in range(len(mopset.single_camera_start)):
                        if (mjd > mopset.single_camera_start[i]) & (mjd < mopset.single_camera_end[i]):
                            method = 'single_camera'
                    tot += 1
                    with open(mopset.dir1+source+'/reduced_data.csv', mode='a') as data_file:
                        fieldnames = ['date','mjd','datetime','wave','ap_size','gain','rotskypa','rotang','wcs_err','seeing','photometric','background','x_pix','y_pix','counts_src','counts_src_err','inst_src','inst_src_err','counts_cal','counts_cal_err','inst_cal','inst_cal_err','q_avg','q_err','u_avg','u_err']
                        data_file = csv.DictWriter(data_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL, fieldnames=fieldnames)
                        fileEmpty = os.stat(mopset.dir1+source+'/reduced_data.csv').st_size == 0
                        if fileEmpty:
                            data_file.writeheader()
                        try:
                            fil = file[1:-8]

                            if method == 'single_camera':
                                date,mjd,date_time,wave,aperture_size,gain,rotskypa,rotang,wcs_err,seeing,pho,s1_src,s1_src_err,inst_src,inst_src_err,s1_cal,s1_cal_err,inst_cal,inst_cal_err,q_avg,q_err,u_avg,u_err,bkg,x,y = mopfunc.one_cam_polarimetry(fil,source)
                                data_file.writerow({'date':date,'mjd':mjd,'datetime':date_time,'wave':wave,'ap_size':aperture_size,'gain':gain,'rotskypa':rotskypa,'rotang':rotang,'wcs_err':wcs_err,'seeing':seeing,'photometric':pho,'background':bkg,'x_pix':x,'y_pix':y,'counts_src':s1_src,'counts_src_err':s1_src_err,'inst_src':inst_src,'inst_src_err':inst_src_err,'counts_cal':s1_cal,'counts_cal_err':s1_cal_err,'inst_cal':inst_cal,'inst_cal_err':inst_cal_err,'q_avg':q_avg,'q_err':q_err,'u_avg':u_avg,'u_err':u_err})
                            if method == 'dual_camera':
                                date,mjd,date_time,wave,aperture_size,gain,rotskypa,rotang,wcs_err,seeing,pho,counts_src,counts_src_err,inst_src,inst_src_err,counts_cal,counts_cal_err,inst_cal,inst_cal_err,q_avg,q_err,u_avg,u_err,bkg,x,y = mopfunc.two_cam_polarimetry(fil,source)
                                data_file.writerow({'date':date,'mjd':mjd,'datetime':date_time,'wave':wave,'ap_size':aperture_size,'gain':gain,'rotskypa':rotskypa,'rotang':rotang,'wcs_err':wcs_err,'seeing':seeing,'photometric':pho,'background':bkg,'x_pix':x,'y_pix':y,'counts_src':counts_src,'counts_src_err':counts_src_err,'inst_src':inst_src,'inst_src_err':inst_src_err,'counts_cal':counts_cal,'counts_cal_err':counts_cal_err,'inst_cal':inst_cal,'inst_cal_err':inst_cal_err,'q_avg':q_avg,'q_err':q_err,'u_avg':u_avg,'u_err':u_err})
                        
                        except Exception as e:
                            print(e)
                            continue
                    print('Reduction Progress:',np.round((tot)/len(file_list)*100,2),'%',source)

            data = pd.read_csv(mopset.dir1+source+'/reduced_data.csv',header=0)
            data = data.sort_values(by=['mjd'])
            data = data.reset_index(drop=True)
            data.to_csv(mopset.dir1+source+'/reduced_data.csv',index=False)
        except Exception as e:
            print(e)
            continue
        
if mopset.clear_dir1 == 'y':
    mopfunc.remove_fits_files(mopset.dir1)

if mopset.calculations == 'y':
    for source in sources:
        try:
            data = pd.read_csv(mopset.dir1+source+'/reduced_data.csv',header=0)
            data.dropna(subset=['q_avg', 'u_avg'],how='any',inplace=True)
            data = data.sort_values(by=['mjd'])
            data = data.reset_index(drop=True)

            data['mag_src'] = data.apply(lambda row: mopinf.source_info[source][row['wave']+'mag_cal'] + row['inst_src'] - row['inst_cal'], axis=1)
            data['mag_src_err'] = np.sqrt(data['inst_src_err']**2+data['inst_cal_err']**2)

            data['mag_cal'] = data.apply(lambda row: mopinf.source_info[source][row['wave']+'mag_cal'], axis=1)

            data['flux_src_mJy'] = data.apply(lambda row: 10**23 * mopinf.flux_constants[row['wave']]['f'] * 10**(-row['mag_src']/2.5) * 1000, axis=1)
            data['flux_src_mJy_err'] = np.log(10)/2.5 * data['flux_src_mJy'] * data['mag_src_err']

            for i in range(1, len(mopset.single_camera_start) + 1):
                mask_1 = data['mjd'] < mopset.single_camera_start[i - 1]
                mask_2 = (data['mjd'] >= mopset.single_camera_start[i - 1]) & (data['mjd'] < mopset.single_camera_end[i - 1])
                mask_3 = data['mjd'] >= mopset.single_camera_end[i - 1]

                key_1 = f'2cam_{i}'
                key_2 = f'1cam_{i}'
                key_3 = f'2cam_{i + 1}'

                data.loc[mask_1, 'q_zero'] = data.loc[mask_1].apply(lambda row: mopinf.polarimetric_constants[f'q_zero_{row["wave"]}_{key_1}'], axis=1)
                data.loc[mask_1, 'u_zero'] = data.loc[mask_1].apply(lambda row: mopinf.polarimetric_constants[f'u_zero_{row["wave"]}_{key_1}'], axis=1)
                data.loc[mask_1, 'k'] = data.loc[mask_1].apply(lambda row: mopinf.polarimetric_constants[f'k_{row["wave"]}_{key_1}'], axis=1)
                data.loc[mask_1, 'inst_depol'] = data.loc[mask_1].apply(lambda row: mopinf.polarimetric_constants[f'inst_depol_{row["wave"]}_{key_1}'], axis=1)

                data.loc[mask_2, 'q_zero'] = data.loc[mask_2].apply(lambda row: mopinf.polarimetric_constants[f'q_zero_{row["wave"]}_{key_2}'], axis=1)
                data.loc[mask_2, 'u_zero'] = data.loc[mask_2].apply(lambda row: mopinf.polarimetric_constants[f'u_zero_{row["wave"]}_{key_2}'], axis=1)
                data.loc[mask_2, 'k'] = data.loc[mask_2].apply(lambda row: mopinf.polarimetric_constants[f'k_{row["wave"]}_{key_2}'], axis=1)
                data.loc[mask_2, 'inst_depol'] = data.loc[mask_2].apply(lambda row: mopinf.polarimetric_constants[f'inst_depol_{row["wave"]}_{key_2}'], axis=1)

                data.loc[mask_3, 'q_zero'] = data.loc[mask_3].apply(lambda row: mopinf.polarimetric_constants[f'q_zero_{row["wave"]}_{key_3}'], axis=1)
                data.loc[mask_3, 'u_zero'] = data.loc[mask_3].apply(lambda row: mopinf.polarimetric_constants[f'u_zero_{row["wave"]}_{key_3}'], axis=1)
                data.loc[mask_3, 'k'] = data.loc[mask_3].apply(lambda row: mopinf.polarimetric_constants[f'k_{row["wave"]}_{key_3}'], axis=1)
                data.loc[mask_3, 'inst_depol'] = data.loc[mask_3].apply(lambda row: mopinf.polarimetric_constants[f'inst_depol_{row["wave"]}_{key_3}'], axis=1)

            data['q'] = data['q_avg'] - data['q_zero']
            data['u'] = data['u_avg'] - data['u_zero']

            data['pol'] = np.nan
            data['pol_err'] = np.nan
            data['evpa'] = np.nan
            data['evpa_err'] = np.nan

            for i in range(len(data)):
                p_mas,p_min,p_max= mopfunc.novel_pol_error(data['q'][i],data['q_err'][i], data['u'][i],data['u_err'][i])
                data.loc[i, 'pol'] = p_mas/data['inst_depol'][i]
                data.loc[i, 'pol_err'] = ((p_max - p_min)/2)/data['inst_depol'][i]

                theta, theta_err = mopfunc.EVPA(data['q'][i],data['q_err'][i], data['u'][i],data['u_err'][i])
                data.loc[i, 'evpa'] = (data['rotskypa'][i]+theta+data['k'][i]) % 180.0
                data.loc[i, 'evpa_err'] = theta_err

            data['q_cor'] = data['pol']*np.cos(np.radians(data['evpa'])*2)
            data['u_cor'] = data['pol']*np.sin(np.radians(data['evpa'])*2)

            data['q_cor_err'] = np.sqrt(np.cos(2*np.radians(data['evpa']))**2*data['pol_err']**2+4*data['pol']**2*np.sin(2*np.radians(data['evpa']))**2*np.radians(data['evpa_err'])**2)
            data['u_cor_err'] = np.sqrt(np.sin(2*np.radians(data['evpa']))**2*data['pol_err']**2+4*data['pol']**2*np.cos(2*np.radians(data['evpa']))**2*np.radians(data['evpa_err'])**2)

            data['pol'] = data['pol']*100
            data['pol_err'] = data['pol_err']*100

            filters = data['wave'].unique()
            for filt in filters:
                globals()[filt+'data'] = data[data['wave'] == filt].reset_index(drop=True)
                globals()[filt+'data']['evpa_cor'] = mopfunc.evpa_unwrap(globals()[filt+'data']['evpa'],globals()[filt+'data']['evpa_err'],mopset.evpa_steps)

            data = pd.concat([globals()[filt+'data'] for filt in filters], ignore_index=True)
            data = data.sort_values(by=['mjd'])
            data = data.reset_index(drop=True)
            data.to_csv(mopset.dir1+source+'/all_data.csv',index=False)

        except Exception as e:
            print(e)
            continue
