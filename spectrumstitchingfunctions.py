#STANDARD LIBRARY IMPORTS
import numpy as np
import math
import os
import sys
#import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize as sciopt
import scipy.ndimage as ndimage
from scipy.interpolate import griddata
#RELATED THIRD PARTY IMPORTS
from scipy.io.idl import readsav
# import pidly

# import astropy.io.ascii as ascii
# import astropy.io.fits as fits
# from astropy.table import Table
#LOCAL APPLICATION/LIBRARY SPECIFIC IMPORTS WOULD GO HERE

# idl = pidly.IDL('/Applications/exelis/idl85/bin/idl') #Call to IDL functions
# idl('.compile /Users/Nicole/Documents/CU/Research/RY_Lup/Scripts/H2_PPD_modeling_code/ccm_unred.pro')

#WAVELENGTHS OF LINE FREE SECTIONS THAT KEVIN IDENTIFIED FOR CONTINUUM EXTRAPOLATION (SHOULD BE MODIFIABLE BY USER)
wcall_dict = dict({'0': [1110., 1115., 1120., 1125., 1130.],
                   '1': [1142.,1147.5,1150.5,1153.5,1157.,1163.,1168.6],
                   '2': [1181.5,1184.2,1187.,1189.],
                   '3': [1245.,1246.5,1249.4,1252.8,1256.2,1256.8],
                   '4': [1261.,1263.6,1268.,1272.9,1275.6,1285.15,1289.5,1291.7,1294.8,1297.8,1300.],
                   '5': [1312.,1313.9,1321.,1322.5,1324.3,1326.3,1328.5],
                   '6': [1420.5,1421.5,1424.5,1426.1,1428.3,1432.,1440.,1444.1, 1447.7,1448.8,1449.4,1451.5,1454.2,1456.6,1458.2],
                   '7': [1462.7,1466.1,1469.4,1471.1,1477.9,1478.7,1480.7,1482.2,1483.8,1485.3,1490.9,1493.4,1496.,1497.9,1499.4,1503.2,1506.1,1508.85,1510.,1511.,1515.2,1518.8,1519.5],
                   '8': [1512.3,1520.6,1523.8,1528.4],
                   '9': [1664.2,1669.,1675.,1676.,1678.,1681.,1684.5,1687.,1689.3,1691.6,1695.,1698.,1699.7],
                   '10': [1703.,1706.,1709.,1711.,1720.8,1722.1,1726.,1729.,1732.5,1735.5,1739.5,1744.,1746.4,1756.5,1759,1762.5]
                   }) #Dropped the wavelengths that Kevin's code identified as bad to cut out the looping

def read_spec_fromsav(filename):
    """Read in spectrum from .sav file
    :param filename: str, name of file containing spectrum
    :return: dataframe (wavelength, flux, flux uncertainties)
    """
    data_table = readsav(filename)
    wavelength = np.array(data_table["wave"]).astype(np.float64) #[Angstroms] (probably, for UV data)
    flux_observed = np.array(data_table["flux"]).astype(np.float64)
    flux_err = np.array(data_table["err"]).astype(np.float64)

    df = pd.DataFrame({'Wavelength [Angstroms]': wavelength,
                       'Fluxes [erg/s/cm^2/Angstrom]': flux_observed,
                       'Flux Uncertainties [erg/s/cm^2/Angstrom]': flux_err})

    return df

def read_spec_fromfits(filename):
    """Read in spectrum from .fits file
    :param filename: str, name of file containing spectrum
    :return: dataframe (wavelength, flux, flux uncertainties)
    """
    hdulist = fits.open(filename)
    scidata = hdulist[1].data
    wavelength = np.array(scidata['Wavelength']).flatten().astype(np.float64)
    flux_observed = np.array(scidata['Flux']).flatten().astype(np.float64)
    flux_err = np.array(scidata['Error']).flatten().astype(np.float64)

    df = pd.DataFrame({'Wavelength [Angstroms]': wavelength,
                       'Fluxes [erg/s/cm^2/Angstrom]': flux_observed,
                       'Flux Uncertainties [erg/s/cm^2/Angstrom]': flux_err})
    return df

#STILL WORKING ON THIS - NEED TO FIGURE OUT HOW TO CHANGE COLUMN NAMES
def read_spec_fromcsv(filename):
    """Read in spectrum from .csv file
    :param filename: str, name of file containing spectrum
    :return: dataframe (wavelength, flux, flux uncertainties)
    """
    df = pd.read_csv(filename, sep=' ')
    keys = list(df)

    pass

def read_spec_fromtxt(filename):
    """Read in spectrum from .txt file
    :param filename: str, name of file containing spectrum
    :return: dataframe (wavelength, flux, flux uncertainties)
    """
    data_table = ascii.read(filename)
    wavelength = np.array(data_table["col1"]).astype(np.float64)
    flux_observed = np.array(data_table["col2"]).astype(np.float64)
    flux_err = np.array(data_table["col3"]).astype(np.float64)

    df = pd.DataFrame({'Wavelength [Angstroms]': wavelength,
                       'Fluxes [erg/s/cm^2/Angstrom]': flux_observed,
                       'Flux Uncertainties [erg/s/cm^2/Angstrom]': flux_err})
    return df

def read_spec(filename):
    """Determine which function to use to read in spectrum based on end of filename
    :param filename: str, name of file containing spectrum
    :return: dataframe (wavelength, flux, flux uncertainties)
    """
    if filename.endswith('.csv') == True:
        return None #Still working on this function
    else:
        if filename.endswith('.sav') == True:
            df = read_spec_fromsav(filename)
        elif filename.endswith('.fits') == True:
            df = read_spec_fromfits(filename)
        elif filename.endswith('.txt') == True:
            df = read_spec_fromtxt(filename)
        return df

def python_deredden(wavelengths, fluxes, EBV): #THIS REALLY NEEDS TO BE RE-WRITTEN IN PYTHON...
    """Python wrapper around IDL ccm_unred function
    :param wavelengths: array, wavelengths of spectrum
    :param fluxes: array, fluxes of spectrum to be dereddened
    :param EBV: float, color excess
    :return: array, dereddened fluxes
    """
    num_points = len(wavelengths)
    #This is a really stupid way to do this - maybe see if Courtney has a smarter way (when she's done with her defense)
    idl('wavelength_forCCM = DBLARR(' + str(num_points) + ')')
    idl('flux_forCCM = DBLARR(' + str(num_points) + ')')
    for wave_idx, wavelength in enumerate(wavelengths):
        idl('wavelength_forCCM[' + str(wave_idx) + '] = ' + str(wavelength))
        idl('flux_forCCM[' + str(wave_idx)+ '] = ' + str(fluxes[wave_idx]))
    idl('ccm_unred, wavelength_forCCM, flux_forCCM, ' + str(EBV) + ', flux_dereddened')
    #GET PYTHON VARIABLE FROM IDL
    flux_dereddened = idl.flux_dereddened

    return flux_dereddened

def read_deredden_write(infile, Av, R, outfile):
    """Dereddens spectrum, using Python wrapper around IDL ccm_unred.pro
    :param infile: string, filename of original (as-observed) data
    :param Av: float, visual extinction toward target
    :param R: float, ratio of total to selective extinction (e.g. R = 3.1 in MW); related to dust grain size
    :outfile: string, filename to which dereddened spectrum will be saved
    :return: None (output saved to file)
    """
    df = pd.read_csv(infile, sep=' ')
    #CALCULATE COLOR EXCESS
    EBV = Av / R
    #Wavelength range of Cardelli extinction curve is 0.1-3.3 microns, so split spectrum at 1000 Angstroms (0.1 microns)
    todered_df = df.loc[df['Wavelength [Angstroms]'] >= 1000.]
    nodered_df = df.loc[df['Wavelength [Angstroms]'] < 1000.]
    #USE PYTHON WRAPPER AROUND IDL DEREDDEN FUNCTION
    flux_dereddened = python_deredden(np.array(todered_df['Wavelength [Angstroms]']),
                                      np.array(todered_df['Fluxes [erg/s/cm^2/Angstrom]']),
                                      EBV) #This will only return one value (instead of the dereddened flux array if EBV < 0)
    todered_df['Dereddened Fluxes'] = flux_dereddened
    nodered_df['Dereddened Fluxes'] = nodered_df['Fluxes [erg/s/cm^2/Angstrom]']

    final_dereddened_df = pd.concat([nodered_df, todered_df]).drop_duplicates().reset_index(drop=True)
    final_dereddened_df.to_csv(outfile, index=None, sep=' ')

def smooth(wavelengths, fluxes, FWHM, kernel_type):
    """Smooths spectrum with specified kernel (from Allison Youngblood)
    :param wavelengths: array, wavelengths of spectrum
    :param fluxes: array, fluxes of spectrum to smooth
    :param FWHM: float, width of smoothing kernel
    :param kernel_type: str, type of smoothing kernel (e.g. 'Boxcar')
    :return: array, convolution of original fluxes with smoothing kernel
    """
    kernel = make_kernel(wavelengths, FWHM, kernel_type)
    smoothed_flux = np.convolve(fluxes,kernel,mode='same')
    return smoothed_flux

def make_kernel(grid, FWHM, kernel_type):
    """Creates a kernel to convolve with spectral data (from Allison Youngblood)
    :param grid: array, wavelengths ("x"-values) to calculate kernel for
    :param FWHM: float, width of smoothing kernel
    :param kernel_type: str, type of smoothing kernel (e.g. 'Boxcar')
    :return: array, values of smoothing kernel calculated at each wavelength in "grid"
    """
    nfwhm = 4.  ## No idea what this parameter is
    ngrd    = len(grid)
    spacing = (grid[ngrd-1] - grid[0]) / (ngrd - 1.)
    nkpts   = round(nfwhm*FWHM / spacing)

    if (nkpts % 2) != 0:
      nkpts += 1

    kernel = spacing*(np.arange(nkpts) - (nkpts / 2.))

    if kernel_type == "Gaussian":
        kernel = np.exp(-np.log(2.) / (FWHM / 2.)**2*(kernel)**2)

    elif kernel_type == "Boxcar":
        out = np.argwhere(np.absolute(kernel) > FWHM/2.)
        inner = np.argwhere(np.absolute(kernel) <= FWHM/2.)
        kernel[out] = 0.
        kernel[inner] = 1.

    #NORMALIZE KERNEL
    kernel_norm = kernel/np.sum(kernel)
    #MAKE SURE KERNEL IS SYMMETRIC
    kernel_norm=np.append(kernel_norm,kernel_norm[0])

    return kernel_norm

def read_smooth_write(infile, outfile, res_to_smooth, kernel_type):
    """Uses Allison's smoothing function on data then saves smoothed spectrum
    :param infile: str, filename of original (as-observed) data
    :param outfile: str, filename to which dereddened spectrum will be saved
    :param res_to_smooth: float, resolution per pixel of smoothing kernel (e.g. 7*0.01 Angstroms/pixel = 7 pixel smoothing)
    :param kernel_type: str, type of smoothing kernel (e.g. "Boxcar", "Gaussian")
    :return: None (smoothed spectrum will be saved to file)
    """
    df = read_spec(infile)
    df['Smoothed Fluxes'] = smooth(df['Wavelength [Angstroms]'], df['Fluxes [erg/s/cm^2/Angstrom]'],
                                   res_to_smooth, kernel_type)
    df.to_csv(outfile, index=None, sep=' ')

def get_G130_data(filename):
    """Prepare the G130M/G160M data for the panchromatic spectrum
    :param filename: str, name of file with data
    :return: dataframe, re-binned spectrum with LyA feature removed
    """
    G130_df = read_spec(filename)
    #NEED TO REBIN G130M DATA (HIGH-RES) TO GET RID OF SOME OF THE NOISE
    red_factor = 1. / 7. #Smoothed by 7 pixels for spectrum sent to Carlo et al.
    #CAN'T OPERATE IN PLACE BECAUSE THE REBINNING REDUCES THE LENGTH OF THE ARRAYS
    G130_df_smoothed = {}
    for col_head in list(G130_df):
        smoothed_col = ndimage.interpolation.zoom(G130_df[col_head], red_factor)
        G130_df_smoothed[col_head] = smoothed_col
    G130_df_smoothed = pd.DataFrame(G130_df_smoothed)
    return G130_df_smoothed

def get_G140_data(filename):
    """Prepare the G140L data for the panchromatic spectrum
    :param filename: str, name of file with data
    :return: dataframe, re-binned spectrum with LyA feature removed
    """
    G140_df = read_spec(filename)
    #NEED TO REBIN TO 5 ANGSTROMS/PIXEL AT WAVELENGTHS SHORTER THAN 1100 ANGSTROMS
    red_factor = 1. / (5. / 0.0803) #0.0803 Angstroms is original resolution per pixel
    #CAN'T OPERATE IN PLACE BECAUSE THE REBINNING REDUCES THE LENGTH OF THE ARRAYS
    nosmooth_df = G140_df[G140_df['Wavelength [Angstroms]'] >= 1100.]
    tosmooth_df = G140_df[G140_df['Wavelength [Angstroms]'] < 1100.]
    smoothed_df = {}
    for col_head in list(tosmooth_df):
        smoothed_col = ndimage.interpolation.zoom(tosmooth_df[col_head], red_factor)
        smoothed_df[col_head] = smoothed_col
    smoothed_df = pd.DataFrame(smoothed_df)
    G140_df_smoothed = pd.concat([smoothed_df, nosmooth_df]).drop_duplicates().reset_index(drop=True)
    #DROP OBSERVED LYA PROFILE (KEEP THE ONE IN THE G130 DATA)
    G140_df_smoothed = G140_df_smoothed.drop(G140_df_smoothed[(G140_df_smoothed['Wavelength [Angstroms]'] > 1200.)
                                                              & (G140_df_smoothed['Wavelength [Angstroms]'] < 1230.)].index)
    return G140_df

def concat_full_spec(data_dict):
    """Read in data from all five observing modes and stick them together in one dataframe
    :param data_dict: dictionary, filenames containing data for each mode
    :return: dataframe, spectra from all five modes
    """
    #READ IN ALL THE DATA, STORE IN DATAFRAMES
    #G130 DATA NEEDS EXTRA SMOOTHING AND LYA REMOVAL
    G130_df = get_G130_data(data_dict['G130/G160'])
    #G140 DATA ALSO NEEDS EXTRA SMOOTHING AND LYA REMOVAL, BUT DONE DIFFERENTLY FROM G130
    G140_df = get_G140_data(data_dict['G140L'])
    #WANT TO DELETE ROWS WITH WAVELENGTHS THAT ARE ALSO PRESENT IN THE G130 DATA, SINCE THOSE ARE HIGHER RESOLUTION
    G140_df = G140_df.drop(G140_df[G140_df['Wavelength [Angstroms]'] >= G130_df['Wavelength [Angstroms]'].iloc[0]].index)
    G230L_df = read_spec(data_dict['G230L'])
    #AGAIN, KEEPING THE HIGHER RESOLUTION DATA
    G230L_df = G230L_df.drop(G230L_df[G230L_df['Wavelength [Angstroms]'] <= G130_df['Wavelength [Angstroms]'].iloc[0]].index)
    G430L_df = read_spec(data_dict['G430L'])
    #STILL TRYING TO KEEP THE HIGHER RESOLUTION DATA
    G430L_df = G430L_df.drop(G430L_df[(G430L_df['Wavelength [Angstroms]'] >= 2890.)
                                      & (G430L_df['Wavelength [Angstroms]'] <= 3150.)].index)
    #COMBINE ALL THE DIFFERENT DATASETS
    full_spectrum_df = pd.concat([G140_df, G130_df, G230L_df, G430L_df]).drop_duplicates().reset_index(drop=True)
    full_spectrum_df.sort_values('Wavelength [Angstroms]', inplace=True)

    #INTERPOLATE OVER GEOCORONAL LYA, KEEPING OBSERVED WINGS (TURN THIS INTO A SEPARATE FUNCTION AT SOME POINT)
    LyA_indices = full_spectrum_df.index[(full_spectrum_df['Wavelength [Angstroms]'] > 1214.5)
                                         & (full_spectrum_df['Wavelength [Angstroms]'] < 1217.)].tolist()
    replacement_value = np.average(np.concatenate((full_spectrum_df['Fluxes [erg/s/cm^2/Angstrom]'][LyA_indices[0]-10:LyA_indices[0]],
                                                   full_spectrum_df['Fluxes [erg/s/cm^2/Angstrom]'][LyA_indices[-1]:LyA_indices[-1]+10]),
                                                   axis=0))
    full_spectrum_df.loc[(full_spectrum_df['Wavelength [Angstroms]'] > 1214.5)
                         & (full_spectrum_df['Wavelength [Angstroms]'] < 1217.),
                         'Fluxes [erg/s/cm^2/Angstrom]'] = replacement_value

    #INTERPOLATE OVER O I LINES AT 1302, 1304, and 1306 ANGSTROMS
    OI_indices = full_spectrum_df.index[(full_spectrum_df['Wavelength [Angstroms]'] > 1300.)
                                        & (full_spectrum_df['Wavelength [Angstroms]'] < 1308.)].tolist()
    replacement_value = np.average(np.concatenate((full_spectrum_df['Fluxes [erg/s/cm^2/Angstrom]'][OI_indices[0]-10:OI_indices[0]],
                                                   full_spectrum_df['Fluxes [erg/s/cm^2/Angstrom]'][OI_indices[-1]:OI_indices[-1]+10]),
                                                   axis=0))
    full_spectrum_df.loc[(full_spectrum_df['Wavelength [Angstroms]'] > 1300.)
                         & (full_spectrum_df['Wavelength [Angstroms]'] < 1308.),
                         'Fluxes [erg/s/cm^2/Angstrom]'] = replacement_value

    return full_spectrum_df

def interp_spec(full_spectrum_df, wavelength_grid):
    """Interpolate panchromatic spectrum onto finer resolution grid and save to file
    :param full_spectrum_df: dataframe, contains spectra for all five modes
    :param wavelength_grid: array, "x"-values at resolution of new grid
    :return: dataframe, spectrum interpolated onto finer resolution grid
    """
    #INTERPOLATE OBSERVED SPECTRUM ONTO WAVELENGTH GRID
    interp_flux_stitched = griddata(full_spectrum_df['Wavelength [Angstroms]'],
                                    full_spectrum_df['Fluxes [erg/s/cm^2/Angstrom]'],
                                    wavelength_grid,
                                    method='linear') #Cubic interpolation leaves more weird spikes in the data
    interp_fluxerr_stitched = griddata(full_spectrum_df['Wavelength [Angstroms]'],
                                       full_spectrum_df['Flux Uncertainties [erg/s/cm^2/Angstrom]'],
                                       wavelength_grid,
                                       method='linear') #
    full_spectrum_interp_df = pd.DataFrame({'Wavelength [Angstroms]': wavelength_grid,
                                            'Fluxes [erg/s/cm^2/Angstrom]': interp_flux_stitched,
                                            'Flux Uncertainties [erg/s/cm^2/Angstrom]': interp_fluxerr_stitched})
    return full_spectrum_interp_df

def load_interp_save_panchrom(data_dict, spec_res):
    """Read in data, interpolate onto new wavelength grid, save panchromatic spectrum to file
    :param data_dict: dictionary, filenames containing data for each mode
    :param spec_res: float, resolution of interpolated spectrum in Angstroms
    :param outfile: str, filename to save panchromatic spectrum to
    """
    full_spectrum_df = concat_full_spec(data_dict)
    wavelength_grid = np.arange(900., full_spectrum_df['Wavelength [Angstroms]'].iloc[-1], spec_res)
    full_spectrum_interp_df = interp_spec(full_spectrum_df, wavelength_grid)
    #SAVE OUTPUT TO .CSV FILE (WILL NEED TO CHANGE FORMAT TO .TXT TO SEND COLLABORATORS)
    full_spectrum_interp_df.to_csv(data_dict['asobserved_outfile'], index=None, sep=' ',
                                   columns=['Wavelength [Angstroms]',
                                            'Fluxes [erg/s/cm^2/Angstrom]',
                                            'Flux Uncertainties [erg/s/cm^2/Angstrom]'])

def insert_model_LyA(data_df, LyA_df):
    """Add model LyA profile to dereddened panchromatic spectrum
    :param data_df: dataframe, contains panchromatic, dereddened spectrum
    :param LyA_df: dataframe, contains model LyA profile
    :param cont_adj: float, flux shift to get LyA profile to match data
    :return: dataframe, with inserted LyA profile in "Dereddened Fluxes" column
    """
    #CUT OBSERVED SPECTRUM WHERE LYA PROFILE WOULD GO
    LyA_indices = data_df.index[(data_df['Wavelength [Angstroms]'] > 1211.)
                                & (data_df['Wavelength [Angstroms]'] < 1221.)].tolist()
    #THIS IS THE INDEX OF THE LAST DATA POINT BEFORE THE MODEL LYA PROFILE WILL BE INSERTED
    index_to_check = LyA_indices[0]
    #NEGATIVE -> HIGHER LYA, POSITIVE -> HIGHER DATA
    cont_adj = np.average(data_df['Dereddened Fluxes'][index_to_check-10:index_to_check-1]) - LyA_df['Dereddened Fluxes'].iloc[0]
    #SHIFT LYA PROFILE TO MATCH FLUX LEVEL OF OBSERVED DATA
    LyA_df['Dereddened Fluxes'] = LyA_df['Dereddened Fluxes'] + cont_adj

    for LyA_idx in LyA_indices:
        wave_to_match = data_df['Wavelength [Angstroms]'].iloc[LyA_idx]
        #FIND VALUE IN MODEL LYA PROFILE AT WAVELENGTH THAT MOST CLOSELY MATCHES THE WAVELENGTH OF THE DATA
        matched_LyAflux = LyA_df.iloc[(LyA_df['Wavelength [Angstroms]']-wave_to_match).abs().argsort()[0]]['Dereddened Fluxes']
        data_df.set_value(LyA_idx, 'Dereddened Fluxes', matched_LyAflux)

    return data_df

def extrap_FUV_cont(object_name, data_dict, dereddened_df, LyA_df):
    """Use dereddened spectrum to estimate the FUV continuum (converted to Python from Kevin's IDL code, re-factored by Nicole)
    :param object_name: str, name of stellar system
    :param dereddened_df: dataframe, contains panchromatic, dereddened spectrum
    :param LyA_df: dataframe, contains model LyA spectrum (need to know wavelengths to cut out from extrapolation)
    :return: dataframe, with calculated continuum values in "Calculated Continuum Fluxes" column
    """
    #MAKE COPY OF DEREDDENED FLUXES TO USE FOR EXTRAPOLATION
    dereddened_df['Continuum Extrapolated Fluxes'] = dereddened_df['Dereddened Fluxes']
    #SET FLUX ACROSS LYA PROFILE EQUAL TO ZERO
    dereddened_df.loc[(dereddened_df['Wavelength [Angstroms]'] >= min(LyA_df['Wavelength [Angstroms]']))
                      & (dereddened_df['Wavelength [Angstroms]'] <= max(LyA_df['Wavelength [Angstroms]'])),
                      'Continuum Extrapolated Fluxes'] = 0.
    #GET RID OF WEIRD INTERPOLATION EFFECTS FROM NEGATIVE FLUX VALUES
    min_flux = 1.e-17 #Minimum allowed flux value (don't want any < 0)
    dereddened_df.loc[dereddened_df['Continuum Extrapolated Fluxes'] < 0., 'Continuum Extrapolated Fluxes'] = min_flux

    #FORCE POLYNOMIAL TO GO TO ZERO AT LAMBDA < 912 ANGSTROMS AND LAMBDA > 1765 ANGSTROMS
    dereddened_df.loc[(dereddened_df['Wavelength [Angstroms]'] <= 912.)
                      & (dereddened_df['Wavelength [Angstroms]'] >= 1765.), 'Continuum Extrapolated Fluxes'] = 0.

    wavelengths_tofit = [912., 913.]
    fluxes_tofit = [0., min_flux]

    binsize = 0.75 #ANGSTROMS
    for wcall_key, wcall_waves in wcall_dict.iteritems():
        for wave_idx, wave in enumerate(wcall_waves):
            df_to_average = dereddened_df.loc[(dereddened_df['Wavelength [Angstroms]'] >= wave - (binsize / 2.)) &
                                              (dereddened_df['Wavelength [Angstroms]'] <= wave + (binsize / 2.))]
            fluxes_to_average = df_to_average['Continuum Extrapolated Fluxes']
            if np.isnan(np.average(fluxes_to_average)) == False:
                wavelengths_tofit.append(wave)
                fluxes_tofit.append(np.average(fluxes_to_average))

    fluxes_tofit_sorted = [flux for _, flux in sorted(zip(wavelengths_tofit, fluxes_tofit))]
    wavelengths_tofit.sort()

    #WRITE 0.75 ANGSTROM BINNED SPECTRUM TO FILE
    continuum_table = Table([wavelengths_tofit, fluxes_tofit_sorted], names = ("wavelength", "flux"))
    continuum_table.write(data_dict['contfit_outfile'], format = "ascii")

    p, v = np.polyfit(wavelengths_tofit, fluxes_tofit_sorted, 2, cov=True) #This produces a fit, but want flux equal to zero at 912 Angstroms
    dereddened_df['Calculated Continuum Fluxes'] = np.square(dereddened_df['Wavelength [Angstroms]'])*p[0] + dereddened_df['Wavelength [Angstroms]']*p[1] + p[2]
    dereddened_df.loc[dereddened_df['Calculated Continuum Fluxes'] < 0., 'Calculated Continuum Fluxes'] = 0.

    return dereddened_df

def addFUVcont_LyA_obs(dereddened_df, LyA_df):
    """Add new series to dataframe with FUV continuum, LyA profile, and observed data all stitched together
    :param dereddened_df: dataframe, contains dereddened spectrum and calculated FUV continuum values
    :param LyA_df: dataframe, contains model LyA profile
    :return: dataframe, with new 'Continuum + LyA + Data' column
    """
    #START WITH COPY OF COLUMN WITH CALCULATED FUV CONTINUUM FLUXES
    dereddened_df['Continuum + LyA + Data'] = dereddened_df['Calculated Continuum Fluxes']
    #ADD LYA PROFILE BACK TO CONTINUUM (WAVELENGTH RANGE OF PROFILE FROM KEVIN)
    LyA_indices = dereddened_df.index[(dereddened_df['Wavelength [Angstroms]'] >= 1211.)
                                      & (dereddened_df['Wavelength [Angstroms]'] <= 1221.)].tolist()
    #WANT TO DETERMINE HOW MUCH TO SHIFT PROFILE IN FLUX SPACE SO THAT IT MATCHES CALCULATED CONTINUUM VALUES
    #THIS WORKS WELL ENOUGH, ALTHOUGH IT'S NOT A PERFECT MATCH TO THE RED SIDE OF THE LINE
    index_to_check = LyA_indices[0]
    cont_adj = np.average(dereddened_df['Continuum + LyA + Data'][index_to_check-10:index_to_check-1]) - dereddened_df['Dereddened Fluxes'].iloc[index_to_check]
    dereddened_df['Continuum + LyA + Data'][LyA_indices] = dereddened_df['Dereddened Fluxes'][LyA_indices] + cont_adj
    #AND ADD REST OF SPECTRUM (ROS) BACK IN
    ros_indices = dereddened_df.index[dereddened_df['Wavelength [Angstroms]'] >= 1700.].tolist()
    dereddened_df['Continuum + LyA + Data'][ros_indices] = dereddened_df['Dereddened Fluxes'][ros_indices]
    dereddened_df.sort_values('Wavelength [Angstroms]', inplace=True)

    return dereddened_df

def smooth_rebin_highres(dereddened_df, desired_res, old_res):
    """Apply convolution to smooth spectrum, then re-bin to 0.25 Angstroms/pixel
    :param df: dataframe, contains dereddened spectrum, FUV continuum calculations, etc.
    :param desired_res: float, resolution to rebin input to
    :param old_res: float, resolution of original input
    :return: dataframe, contains low-res version of panchromatic spectrum
    """
    #SMOOTH AS-OBSERVED FLUXES
    low_res_obsflux = smooth(dereddened_df['Wavelength [Angstroms]'],
                             dereddened_df['Fluxes [erg/s/cm^2/Angstrom]'],
                             1., 'Gaussian')
    #SMOOTH SPECTRUM WITH FUV CONTINUUM + LYA + OBSERVED DATA
    low_res_fluxwithcontinuum = smooth(dereddened_df['Wavelength [Angstroms]'],
                                       dereddened_df['Continuum + LyA + Data'],
                                       1., "Gaussian")
    #SMOOTH DEREDDENED FLUXES
    low_res_fluxdereddened = smooth(dereddened_df['Wavelength [Angstroms]'],
                                    dereddened_df['Dereddened Fluxes'],
                                    1., "Gaussian")
    #SMOOTH CALCULATED CONTINUUM FLUXES
    low_res_continuumonly = smooth(dereddened_df['Wavelength [Angstroms]'],
                                   dereddened_df['Calculated Continuum Fluxes'],
                                   1., "Gaussian")
    #DECREASE RESOLUTION OF WAVELENGTH GRID TO 0.25 ANGSTROMS
    reduction_factor = 1. / (desired_res / old_res)
    wavelength_final_rebinned = ndimage.interpolation.zoom(dereddened_df['Wavelength [Angstroms]'], reduction_factor)
    fluxobs_rebinned = ndimage.interpolation.zoom(low_res_obsflux, reduction_factor)
    fluxdereddened_rebinned = ndimage.interpolation.zoom(low_res_fluxdereddened, reduction_factor)
    fluxwithcontinuum_rebinned = ndimage.interpolation.zoom(low_res_fluxwithcontinuum, reduction_factor)
    continuumonly_rebinned = ndimage.interpolation.zoom(low_res_continuumonly, reduction_factor)

    lowres_df = pd.DataFrame({"wavelength": wavelength_final_rebinned,
                              "flux_observed": fluxobs_rebinned,
                              "flux_red_corrected": fluxdereddened_rebinned,
                              "flux_red_corrected_continuum": fluxwithcontinuum_rebinned,
                              "continuum_only": continuumonly_rebinned})
    return lowres_df

def stitch_spectra(object_name, Av, Rv, spec_res, low_spec_res, LyA_df, data_dict):
    """Function to carry out full spectrum stitching procedure
    :param object_name: str, name of stellar system
    :param Av: float, extinction along line of sight
    :param Rv: float, ratio of total to selective extinction (e.g. R = 3.1 in MW); related to dust grain size
    :param spec_res: float, desired resolution of interpolated spectrum
    :param low_spec_res: float, desired resolution of low-res interpolated spectrum
    :param LyA_df: dataframe, contains model LyA profile for the stellar system
    :param data_dict: dict, filenames for data from each observing mode of COS and STIS
    :return: None (function output will be saved to text files)
    """
    #CONCATENATE OBSERVED SPECTRA FROM ALL MODES, INTERPOLATE ONTO GRID WITH FINER RESOLUTION, SAVE TO TEXT FILE
    load_interp_save_panchrom(data_dict, spec_res)
    #DEREDDEN SPECTRUM AND WRITE TO FILE
    read_deredden_write(data_dict['asobserved_outfile'], Av, Rv, data_dict['dereddened_outfile'])
    #RE-READ IN DEREDDENED SPECTRUM
    dereddened_df = pd.read_csv(data_dict['dereddened_outfile'], sep=' ')
    #CALCULATE TRANSMISSION (amount of flux observed compared to total flux based on extinction correction)
    dereddened_df['Transmission'] = dereddened_df['Fluxes [erg/s/cm^2/Angstrom]'] / dereddened_df['Dereddened Fluxes']
    #ADD INTRINSIC LYA PROFILE TO PANCHROMATIC SPECTRUM
    dereddened_df = insert_model_LyA(dereddened_df, LyA_df)
    #EXTRAPOLATE FUV CONTINUUM
    dereddened_df = extrap_FUV_cont(object_name, data_dict, dereddened_df, LyA_df)
    #STITCH TOGETHER CALCULATED FUV CONTINUUM, MODEL LYA, AND OBSERVED DATA BEYOND 1700 ANGSTROMS
    dereddened_df = addFUVcont_LyA_obs(dereddened_df, LyA_df)
    #INTERPOLATE OVER CHUNK NEAR 3100 ANGSTROMS (OVERLAP BETWEEN GRATINGS) - THIS SHOULD ALSO BE ITS OWN FUNCTION
    overlap_indices = dereddened_df.index[(dereddened_df['Wavelength [Angstroms]'] > 3120.)
                                          & (dereddened_df['Wavelength [Angstroms]'] < 3288.)].tolist()
    wavelengths_tointerp = np.concatenate((dereddened_df['Wavelength [Angstroms]'][overlap_indices[0]-10:overlap_indices[0]],
                                           dereddened_df['Wavelength [Angstroms]'][overlap_indices[-1]:overlap_indices[-1]+10]),
                                           axis=0)
    #AS OBSERVED FLUXES FIRST
    obs_fluxes_tointerp = np.concatenate((dereddened_df['Fluxes [erg/s/cm^2/Angstrom]'][overlap_indices[0]-10:overlap_indices[0]],
                                          dereddened_df['Fluxes [erg/s/cm^2/Angstrom]'][overlap_indices[-1]:overlap_indices[-1]+10]),
                                          axis=0)
    obs_replacement_values = np.interp(dereddened_df['Wavelength [Angstroms]'][overlap_indices], wavelengths_tointerp, obs_fluxes_tointerp)
    dereddened_df['Fluxes [erg/s/cm^2/Angstrom]'][overlap_indices] = obs_replacement_values

    #NOW DEREDDENED FLUXES
    dered_fluxes_tointerp = np.concatenate((dereddened_df['Dereddened Fluxes'][overlap_indices[0]-10:overlap_indices[0]],
                                            dereddened_df['Dereddened Fluxes'][overlap_indices[-1]:overlap_indices[-1]+10]),
                                            axis=0)
    dered_replacement_values = np.interp(dereddened_df['Wavelength [Angstroms]'][overlap_indices], wavelengths_tointerp, dered_fluxes_tointerp)
    dereddened_df['Dereddened Fluxes'][overlap_indices] = dered_replacement_values

    #Write high-res spectrum to file (as-observed fluxes have been written out elsewhere, don't include LyA model profile)
    highres_table = Table([dereddened_df['Wavelength [Angstroms]'],
                           dereddened_df['Fluxes [erg/s/cm^2/Angstrom]'],
                           dereddened_df['Dereddened Fluxes'],
                           dereddened_df['Continuum + LyA + Data'],
                           dereddened_df['Calculated Continuum Fluxes']],
                          names = ("wavelength",
                                   "flux_observed",
                                   "flux_red_corrected",
                                   "flux_red_corrected_continuum",
                                   "continuum_only"))
    highres_table.write(data_dict['highres_outfile'], format = "ascii")

    #REBIN SPECTRUM TO 0.25 ANGSTROM RESOLUTION (FOR COLLABORATORS)
    lowres_dereddened_df = smooth_rebin_highres(dereddened_df, low_spec_res, spec_res)
    #WRITE LOW RES SPECTRUM TO FILE
    lowres_table = Table([lowres_dereddened_df['wavelength'],
                          lowres_dereddened_df['flux_observed'],
                          lowres_dereddened_df['flux_red_corrected'],
                          lowres_dereddened_df['flux_red_corrected_continuum'],
                          lowres_dereddened_df['continuum_only']],
                         names = ('wavelength',
                                  'flux_observed',
                                  'flux_red_corrected',
                                  'flux_red_corrected_continuum',
                                  'continuum_only'
                                  ))
    lowres_table.write(data_dict['lowres_outfile'], format = "ascii")
