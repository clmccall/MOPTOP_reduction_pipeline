from numpy import average

"""
Dictionaries containing important constant information for a source, its calibration star and for MOPTOP itself. 

source_info             -   Gives coordinate infornmation for a source and its associated calibration star. The calibration star magnitudes are also given.
polarimetric_constants  -   Gives q_0, u_0 and k values for all filters at various stages over the lifetime of MOPTOP. 
flux_constants          -   Gives the mid-wavelength and flux calibration constant for each filter on MOPTOP.
"""

source_info = {
    "BL Lac" : {
        "alt_names" : ['bllac','BLLac','BL Lac'],
        "aperture_size" : 4,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "22 02 43.291",
        "dec" : "+42 16 39.979",
        "cal_ra" : "22 02 45.418",
        "cal_dec" : "+42 16 35.473",
        "Bmag_cal" : 14.52,
        "Vmag_cal" : 12.78,
        "Rmag_cal" : 11.93,
        "Imag_cal" : 11.09,
        "Lmag_cal" : 0,
    },
    "OJ287" : {
        "alt_names" : ['oj287','OJ287'],
        "aperture_size" : 5.5,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "08 54 48.875",
        "dec" : "+20 06 30.643",
        "cal_ra" : "08 54 54.44",
        "cal_dec" : "+20 06 13.69",
        "Bmag_cal" : 15.27,
        "Vmag_cal" : 14.56,
        "Rmag_cal" : 14.26,
        "Imag_cal" : 13.54,
        "Lmag_cal" : 0,
    },
    "TXS 0506+056" : {
        "alt_names" : ['txs0506+056','TXS0506+056','TXS 0506+056'],
        "aperture_size" : 5,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "05 09 25.965",
        "dec" : "+05 41 35.334",
        "cal_ra" : "05 09 27.787",
        "cal_dec" : "+05 41 48.39",
        "Bmag_cal" : 16.713,
        "Vmag_cal" : 15.798,
        "Rmag_cal" : 15.260,
        "Imag_cal" : 0,
        "Lmag_cal" : 0,
    },
    "3C 454.3" : {
        "alt_names" : ['3c454.3','3C454.3','3C 454.3'],
        "aperture_size" : 5,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "22 53 57.748",
        "dec" : "+16 08 53.561",
        "cal_ra" : "22 53 59.928",
        "cal_dec" : "+16 09 53.93",
        "Bmag_cal" : 16.85,
        "Vmag_cal" : 15.86,
        "Rmag_cal" : 15.32,
        "Imag_cal" : 14.80,
        "Lmag_cal" : 0,
    },
    "4C 11.69" : {
        "alt_names" : ['4c11.69','4C11.69','4C 11.69'],
        "aperture_size" : 5,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "22 32 36.41",
        "dec" : "+11 43 50.90",
        "cal_ra" : "22 32 32.834",
        "cal_dec" : "+11 42 43.83",
        "Bmag_cal" : 17.18,
        "Vmag_cal" : 16.609,
        "Rmag_cal" : 16.31,
        "Imag_cal" : 0,
        "Lmag_cal" : 0,
    },
    "PG1553+113" : {
        "alt_names" : ['pg1553+113','PG1553+113','PG 1553+113'],
        "aperture_size" : 5,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "15 55 43.044",
        "dec" : "+11 11 24.366",
        "cal_ra" : "15 55 46.076",
        "cal_dec" : "+11 11 19.553",
        "Bmag_cal" : 14.543,
        "Vmag_cal" : 13.923,
        "Rmag_cal" : 13.582,
        "Imag_cal" : 13.230,
        "Lmag_cal" : 0,
    },
    "PKS1510-089" : {
        "alt_names" : ['pks1510-089','PKS1510-089','PKS 1510-089'],
        "aperture_size" : 5,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "15 12 50.533",
        "dec" : "-09 05 59.830",
        "cal_ra" : "15 12 51.665",
        "cal_dec" : "-09 05 23.501",
        "Bmag_cal" : 15.352,
        "Vmag_cal" : 14.687,
        "Rmag_cal" : 14.524,
        "Imag_cal" : 14.347,
        "Lmag_cal" : 0,
    },
    "j0710+591" : {
        "alt_names" : ['j0710+591'],
        "aperture_size" : 4,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "07 10 30.068",
        "dec" : "+59 08 20.372",
        "cal_ra" : "07 10 30.068",
        "cal_dec" : "+59 08 20.372",
        "Bmag_cal" : 0,
        "Vmag_cal" : 0,
        "Rmag_cal" : 0,
        "Imag_cal" : 0,
        "Lmag_cal" : 0,
    },
    "TON599" : {
        "alt_names" : ['ton599','TON599'],
        "aperture_size" : 4,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "11 59 31.834",
        "dec" : "+29 14 43.827",
        "cal_ra" : "11 59 33.591",
        "cal_dec" : "+29 15 00.415",
        "Bmag_cal" : 0,
        "Vmag_cal" : 0,
        "Rmag_cal" : 0,
        "Imag_cal" : 0,
        "Lmag_cal" : 0,
    },
    "OJ248" : {
        "alt_names" : ['oj248','OJ248'],
        "aperture_size" : 4,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "08 30 52.086",
        "dec" : "+24 10 59.820",
        "cal_ra" : "08 30 50.194",
        "cal_dec" : "+24 10 26.16",
        "Bmag_cal" : 14.956,
        "Vmag_cal" : 14.662,
        "Rmag_cal" : 14.44,
        "Imag_cal" : 0,
        "Lmag_cal" : 0,
    },
    "PKS0735+178" : {
        "alt_names" : ['pks0735+178','PKS0735+178','PKS 0735+178'],
        "aperture_size" : 4,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "07 38 07.394",
        "dec" : "+17 42 18.998",
        "cal_ra" : "07 38 02.489",
        "cal_dec" : "+17 41 21.968",
        "Bmag_cal" : 15.48,
        "Vmag_cal" : 14.40,
        "Rmag_cal" : 13.85,
        "Imag_cal" : 0,
        "Lmag_cal" : 0,
    },
    "S50716+714" : {
        "alt_names" : ['s50716+714','S50716+714','S5 0716+714','s5 0716+74'],
        "aperture_size" : 4,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "07 21 53.448",
        "dec" : "+71 20 36.364",
        "cal_ra" : "07 21 54.363",
        "cal_dec" : "+71 19 20.917",
        "Bmag_cal" : 14.145,
        "Vmag_cal" : 13.47,
        "Rmag_cal" : 13.18,
        "Imag_cal" : 0,
        "Lmag_cal" : 0,
    },
    "MRK421" : {
        "alt_names" : ['mrk421','MRK421','mrk421_offset','MRK421_offset'],
        "aperture_size" : 6,
        "inner_annulus" : 2.5,
        "outer_annulus" : 3.5,
        "ra" : "11 04 27.314",
        "dec" : "+38 12 31.798",
        "cal_ra" : "11 04 29.28",
        "cal_dec" : "+38 09 8.13",
        "Bmag_cal" : 18.57,
        "Vmag_cal" : 17.57,
        "Rmag_cal" : 16.67,
        # "Bmag_cal" : 18.1,
        # "Vmag_cal" : 17.55,
        # "Rmag_cal" : 16.49,
        # "cal_ra" : "11 04 15.024",
        # "cal_dec" : "+38 08 27.94",
        # "Bmag_cal" : 18.27,
        # "Vmag_cal" : 16.59,
        # "Rmag_cal" : 16.27,
        # "Imag_cal" : 0,
        # "Lmag_cal" : 0,
    },
    "PKS2155-304" : {
        "alt_names" : ['pks2155-304','PKS 2155-304','pks 2155-304','PKS2155-304'],
        "aperture_size" : 8,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "21 58 52.065",
        "dec" : "-30 13 32.118",
        # "cal_ra" : "21 58 58.645",
        # "cal_dec" : "-30 11 51.940",
        # "Bmag_cal" : 15.530,
        # "Vmag_cal" : 14.660,
        # "Rmag_cal" : 14.550,
        # "Imag_cal" : 13.42,
        # "Lmag_cal" : 0,
        "cal_ra" : "21 58 42.30",
        "cal_dec" : "-30 10 27.00",
        "Bmag_cal" : 14.930,
        "Vmag_cal" : 14.280,
        "Rmag_cal" : 13.920,
        "Imag_cal" : 13.560,
        "Lmag_cal" : 0,
    },
    "3C 345" : {
        "alt_names" : ['3c345','3C345','3C 345','3c 345'],
        "aperture_size" : 8,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "16 42 58.810",
        "dec" : "+39 48 36.994",
        "cal_ra" : "16 42 52.879",
        "cal_dec" : "+39 48 33.610",
        "Bmag_cal" : 16.43,
        "Vmag_cal" : 15.17,
        "Rmag_cal" : 14.50,
        "Imag_cal" : 0,
        "Lmag_cal" : 0,
    },
    "MRK501" : {
        "alt_names" : ['MRK501','mrk501','MRK 501','mrk 501','Mrk501','Mrk 501','NEW_MKR_501_2000'],
        "aperture_size" : 0,
        "inner_annulus" : 0,
        "outer_annulus" : 0, 
        "ra" : "16 53 52.217",
        "dec" : "+39 45 36.609",
        # "cal_ra" : "16 53 57.050",
        # "cal_dec" : "+39 45 35.434",
        "cal_ra" : "16 53 45.866",
        "cal_dec" : "+39 44 9.932",
        # "Bmag_cal" : 16.44,
        # "Vmag_cal" : 15.61,
        # "Rmag_cal" : 15.10,
        # "Imag_cal" : 14.48,
        # "Lmag_cal" : 0,
        "Bmag_cal" : 13.245,
        "Vmag_cal" : 12.534,
        "Rmag_cal" : 12.447,
        "Imag_cal" : 11.28,
        "Lmag_cal" : 0,
    },
    "J004354+404634" : {
        "alt_names" : ['J004354+404634','j004354+404634'],
        "aperture_size" : 8,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "00:43:54.34",
        "dec" : "+40:46:34.80",
        "cal_ra" : "00:43:45.261",
        "cal_dec" : "+40:46:59.152",
        "Bmag_cal" : 17.3,
        "Vmag_cal" : 16.96,
        "Rmag_cal" : 16.5,
        "Imag_cal" : 0,
        "Lmag_cal" : 0,
    },
    "J063520.9-262039" : {
        "alt_names" : ['J063520.9-262039','j063520.9-262039'],
        "aperture_size" : 8,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "06 35 20.909",
        "dec" : "-26 20 39.869",
        "cal_ra" : "06 35 18.412",
        "cal_dec" : "-26 19 34.714",
        "Bmag_cal" : 15.651,
        "Vmag_cal" : 15.140,
        "Rmag_cal" : 14.868,
        "Imag_cal" : 0,
        "Lmag_cal" : 0,
    },
    "AT2022fpx" : {
        "alt_names" : ['at2022fpx'],
        "aperture_size" : 4,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "15 31 03.702",
        "dec" : "+53 24 19.65",
        "cal_ra" : "15 31 16.419",
        "cal_dec" : "+53 22 59.60",
        "Bmag_cal" : 16.253,
        "Vmag_cal" : 15.49,
        "Rmag_cal" : 15.255,
        "Imag_cal" : 14.943,
        "Lmag_cal" : 0,
    },
    "GRB221009A" : {
        "aperture_size" : 4,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "19 13 03.151",
        "dec" : "19 46 22.28",
        "cal_ra" : "19 13 06.89",
        "cal_dec" : "19 46 52.50",
        "Bmag_cal" : 0,
        "Vmag_cal" : 0,
        "Rmag_cal" : 0,
        "Imag_cal" : 0,
        "Lmag_cal" : 0,
    },
    "GRB 230818A" : {
        "alt_names" : ['SWIFT:GRB:1186032:0'],
        "aperture_size" : 4,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "19:03:33.1415",
        "dec" : "+40:54:09.513",
        "cal_ra" : "19 03 40.351",
        "cal_dec" : "+40 51 47.117",
        "Bmag_cal" : 0,
        "Vmag_cal" : 0,
        "Rmag_cal" : 0,
        "Imag_cal" : 0,
        "Lmag_cal" : 0,
    },
    "SN2023ixf" : {
        "alt_names" : ['SN2023ixf'],
        "aperture_size" : 8,
        "inner_annulus" : 2.5,
        "outer_annulus" : 4.5,
        "ra" : "14 03 38.564",
        "dec" : "+54 18 41.94",
        "cal_ra" : "14 03 38.562",
        "cal_dec" : "+54 18 42.02",
        "Bmag_cal" : 0,
        "Vmag_cal" : 0,
        "Rmag_cal" : 0,
        "Imag_cal" : 0,
        "Lmag_cal" : 0,
    },
    "HD14069" : {
        "alt_names" : ['HD14069_zpol'],
        "aperture_size" : 8,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "02 16 45.194",
        "dec" : "+07 41 10.637",
    },
    "BD+32 3739" : {
        "alt_names" : ['BDp32_3739_zpol'],
        "aperture_size" : 11.5,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "20 12 02.149",
        "dec" : "+32 47 43.693",
    },
    "GD319" : {
        "alt_names" : ['GD319_zpol'],
        "aperture_size" : 10.5,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "12 50 04.526",
        "dec" : "+55 06 01.800",
    },
    "HD212311" : {
        "alt_names" : ['HD212311_zpol'],
        "aperture_size" : 10.5,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "22 21 58.583",
        "dec" : "+56 31 52.718",
    },
    "G191B2B" : {
        "alt_names" : ['G191B2B_zpol'],
        "aperture_size" : 7,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "05 05 30.618",
        "dec" : "+52 49 51.919",
    },
    "BD+28 4211" : {
        "alt_names" : ['BDp28_4211_zpol'],
        "aperture_size" : 10,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "21 51 11.022",
        "dec" : "+28 51 50.368",
    },
    "HD251204" : {
        "alt_names" : ['HD251204_pol'],
        "aperture_size" : 4.5,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "06 05 05.667",
        "dec" : "+23 23 38.533",
        "theta_cor_act_B" : 154.83,
        "theta_cor_act_V" : 153.26,
        "theta_cor_act_R" : 152.87,
        "theta_cor_act_I" : 153.70,
        "theta_cor_act_L" : average([154.83,153.26,152.87,153.70]),
        "per_p_act_B" : 4.46,
        "per_p_act_V" : 4.72,
        "per_p_act_R" : 4.78,
        "per_p_act_I" : 4.37,
        "per_p_act_L" : average([4.46,4.72,4.78,4.37]),
    },
    "BD+64 106" : {
        "alt_names" : ['BDp64_106_pol'],
        "aperture_size" : 7.5,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "00 57 36.689",
        "dec" : "+64 51 34.907",
        "theta_cor_act_B" : 97.15,
        "theta_cor_act_V" : 96.63,
        "theta_cor_act_R" : 96.74,
        "theta_cor_act_I" : 96.89,
        "theta_cor_act_L" : average([97.15,96.63,96.74,96.89]),
        "per_p_act_B" : 5.506,
        "per_p_act_V" : 5.687,
        "per_p_act_R" : 5.150,
        "per_p_act_I" : 4.696,
        "per_p_act_L" : average([5.506,5.687,5.150,4.696]),
    },
    "VICyg12" : {
        "alt_names" : ['VICyg12_pol'],
        "aperture_size" : 5,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "20 32 40.957",
        "dec" : "+41 14 29.279",
        "theta_cor_act_B" : 119,
        "theta_cor_act_V" : 115.03,
        "theta_cor_act_R" : 116.23,
        "theta_cor_act_I" : 117,
        "theta_cor_act_L" : average([119,115.03,116.23,117]),
        "per_p_act_B" : 9.67,
        "per_p_act_V" : 8.947,
        "per_p_act_R" : 7.893,
        "per_p_act_I" : 7.06,
        "per_p_act_L" : average([9.67,8.947,7.893,7.06]),
    },
    "HD155197" : {
        "alt_names" : ['HD155197_pol'],
        "aperture_size" : 8,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "17 10 15.747",
        "dec" : "-04 50 03.657",
        "theta_cor_act_B" : 103.06,
        "theta_cor_act_V" : 102.84,
        "theta_cor_act_R" : 102.88,
        "theta_cor_act_I" : 103.18,
        "theta_cor_act_L" : average([103.06,102.84,102.88,103.18]),
        "per_p_act_B" : 4.112,
        "per_p_act_V" : 4.320,
        "per_p_act_R" : 4.274,
        "per_p_act_I" : 3.906,
        "per_p_act_L" : average([4.112,4.320,4.274,3.906]),
    },
    "HILT960" : {
        "alt_names" : ['HILT960_pol','HILT_960_pol'],
        "aperture_size" : 9,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "20 23 28.531",
        "dec" : "+39 20 59.037",
        "theta_cor_act_B" : 55.06,
        "theta_cor_act_V" : 54.79,
        "theta_cor_act_R" : 54.54,
        "theta_cor_act_I" : 53.96,
        "theta_cor_act_L" : average([55.06,54.76,54.54,53,96]),
        "per_p_act_B" : 5.720,
        "per_p_act_V" : 5.663,
        "per_p_act_R" : 5.210,
        "per_p_act_I" : 4.455,
        "per_p_act_L" : average([5.720,5.663,5.210,4.455]),
    },
    "BD+57 2615" : {
        "alt_names" : ['BDp57_2615_pol'],
        "aperture_size" : 4.5,
        "inner_annulus" : 1.5,
        "outer_annulus" : 2.5,
        "ra" : "22 47 49.571",
        "dec" : "+58 08 49.516",
        "theta_cor_act_B" : 42,
        "theta_cor_act_V" : 42,
        "theta_cor_act_R" : 41,
        "theta_cor_act_I" : 41,
        "theta_cor_act_L" : average([42,42,41,41]),
        "per_p_act_B" : 1.91,
        "per_p_act_V" : 2.00,
        "per_p_act_R" : 2.02,
        "per_p_act_I" : 1.71,
        "per_p_act_L" : average([1.91,2.00,2.02,1.71]),
    }
}

polarimetric_constants = {
    "q_zero_B_1" : 0.00189,
    "u_zero_B_1" : -0.01147,
    "q_zero_V_1" : 0.00581,
    "u_zero_V_1" : -0.02141,
    "q_zero_R_1" : 0.01071,
    "u_zero_R_1" : -0.02911,
    "q_zero_I_1" : 0.01147,
    "u_zero_I_1" : -0.03162,
    "q_zero_L_1" : 0.00727,
    "u_zero_L_1" : -0.01818,

    "k_B_1" : 31.88,
    "k_V_1" : 32.69,
    "k_R_1" : 31.54,
    "k_I_1" : 32.38,
    "k_L_1" : 31.85,

    "inst_depol_B_1" : 0.856,
    "inst_depol_V_1" : 0.836,
    "inst_depol_R_1" : 0.873,
    "inst_depol_I_1" : 0.743,
    "inst_depol_L_1" : 0.880,

    "q_zero_B_2" : 0.00834,
    "u_zero_B_2" : -0.00857,
    "q_zero_V_2" : 0.01946,
    "u_zero_V_2" : -0.01136,
    "q_zero_R_2" : 0.02794,
    "u_zero_R_2" : -0.01558,
    "q_zero_I_2" : 0.03056,
    "u_zero_I_2" : -0.01536,
    "q_zero_L_2" : 0.01731,
    "u_zero_L_2" : -0.00884,

    "k_B_2" : 12.53,
    "k_V_2" : 9.26,
    "k_R_2" : 11.94,
    "k_I_2" : 9.48,
    "k_L_2" : 11.50,

    "inst_depol_B_2" : 0.959,
    "inst_depol_V_2" : 0.967,
    "inst_depol_R_2" : 0.989,
    "inst_depol_I_2" : 0.836,
    "inst_depol_L_2" : 0.972,

    "q_zero_B_3" : 0.00127,
    "u_zero_B_3" : -0.01219,
    "q_zero_V_3" : 0.00601,
    "u_zero_V_3" : -0.02326,
    "q_zero_R_3" : 0.01063,
    "u_zero_R_3" : -0.03140,
    "q_zero_I_3" : 0.01135,
    "u_zero_I_3" : -0.03401,
    "q_zero_L_3" : 0.00787,
    "u_zero_L_3" : -0.01979,

    "k_B_3" : 123.77,
    "k_V_3" : 122.08,
    "k_R_3" : 123.55,
    "k_I_3" : 123.68,
    "k_L_3" : 124.85,

    "inst_depol_B_3" : 0.913,
    "inst_depol_V_3" : 0.888,
    "inst_depol_R_3" : 0.901,
    "inst_depol_I_3" : 0.786,
    "inst_depol_L_3" : 0.920,

    "q_zero_B_4" : 0.0007,
    "u_zero_B_4" : -0.0134,
    "q_zero_V_4" : 0.0064,
    "u_zero_V_4" : -0.0250,
    "q_zero_R_4" : 0.0102,
    "u_zero_R_4" : -0.0324,
    "q_zero_I_4" : 0.0109,
    "u_zero_I_4" : -0.0348,
    "q_zero_L_4" : 0.0074,
    "u_zero_L_4" : -0.0200,

    "k_B_4" : 123.85,
    "k_V_4" : 121.55,
    "k_R_4" : 123.89,
    "k_I_4" : 123.97,
    "k_L_4" : 124.33,

    "inst_depol_B_4" : 0.80,
    "inst_depol_V_4" : 0.88,
    "inst_depol_R_4" : 0.93,
    "inst_depol_I_4" : 0.79,
    "inst_depol_L_4" : 0.87,

}

flux_constants = {
    "B" : {
        "wavelength" : 4500,
        "f" : 4.28*10**(-20),
    },
    "V" : {
        "wavelength" : 5300,
        "f" : 3.54*10**(-20),
    },
    "R" : {
        "wavelength" : 6375,
        "f" : 2.94*10**(-20),
    },
    "I" : {
        "wavelength" : 8000,
        "f" : 2.34*10**(-20),
    },
    "L" : {
        "wavelength" : 5500,
        "f" : 3.41*10**(-20),
    }
}