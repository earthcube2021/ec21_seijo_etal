#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Data manipulation
import requests
import numpy as np
import pandas as pd
from scipy import interpolate
from itertools import compress
from datetime import datetime
from datetime import timedelta  

# Visualizations
import matplotlib
import matplotlib.pylab as plt
from matplotlib import cm
import matplotlib.dates as mdates
import cartopy.crs as ccrs
import cartopy.feature as cft
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from svgpath2mpl import parse_path

# Progress bar
from tqdm.notebook import tqdm, trange
import time  

#prevent warnings from showing on screen
import warnings
warnings.filterwarnings('ignore')

# **get_TCs_byNameYear**
# 
# Query Tropical Cyclone data by name ('tc_name') and year ('tc_year').
# 
# name format: e.g. 'maria', date format: 'yyyy-mm-dd'.
def get_TCs_byNameYear(tc_name,tc_year):
    url = 'https://argovis.colorado.edu/tc/findByNameYear?name='+tc_name+'&year='+str(tc_year) #2018-07-15
    print(url)
    resp = requests.get(url)
    # Consider any status other than 2xx an error
    if not resp.status_code // 100 == 2:
        return "Error: Unexpected response {}".format(resp)
    data = resp.json()
    return data


# **get_TCs_byDate**
# 
# Query Tropical Cyclones data by date ('startDate','endDate').
# 
# date format: 'yyyy-mm-dd'.
def get_TCs_byDate(startDate,endDate):
    url = 'https://argovis.colorado.edu/tc/findByDateRange?startDate='+startDate+'T00:00:00&endDate='+endDate+'T00:00:00' #2018-07-15
    print(url)
    resp = requests.get(url)
    # Consider any status other than 2xx an error
    if not resp.status_code // 100 == 2:
        return "Error: Unexpected response {}".format(resp)
    data = resp.json()
    return data

# **get_track_for_storm**

# Load the track for the storm of interest as described by two strings, tc_name (lower case) and tc_year
def get_track_for_storm(tc_name='maria',tc_year='2017'):
    tc_star = get_TCs_byNameYear(tc_name,tc_year)
    df  = pd.DataFrame(tc_star[0]['traj_data'])
    df['_id']  = tc_star[0]['_id']
    return df


# **get_hurricane_marker**
# 
# Generates a hurricane marker for plotting.
#
# create a marker for tropical cyclones
def get_hurricane_marker():
    hurricane = parse_path("""M 188.79857,492.55018 L 180.09663,484.17671 L 188.57725,474.97463 
     C 212.44187,449.07984 230.37031,423.34927 246.04359,392.5 
     C 254.14781,376.5487 265.0005,350.78866 265.0005,347.50373 
     C 265.0005,346.53236 263.21255,347.40666 259.7505,350.07094 
     C 251.67361,356.28665 233.85001,364.64767 222.5005,367.54485 
     C 204.24051,372.20605 178.92084,371.97166 159.45635,366.96123 
     C 147.77122,363.95331 130.93184,355.3283 122.0005,347.77659 
     C 95.11018,325.04006 81.65749,291.36529 81.74139,247 
     C 81.78993,221.33674 85.91479,197.1747 94.55247,171.95714 
     C 111.06665,123.74428 136.98179,82.210848 180.29075,34.54693 
     L 185.6999,28.59386 L 189.6002,31.718323 
     C 191.74536,33.436777 195.9159,37.308373 198.86805,40.32187 
     L 204.23561,45.800955 L 193.66355,57.483549 
     C 168.69038,85.080007 151.53704,109.91644 136.8182,139.79028 
     C 130.67851,152.2516 118.91503,180.17836 119.52809,180.83739 
     C 119.70071,181.02295 122.91512,178.62979 126.67122,175.51926 
     C 144.84799,160.46658 171.06913,152.9127 200.0005,154.39429 
     C 227.96505,155.82638 249.78837,164.40176 267.15103,180.78081 
     C 291.49094,203.74185 302.41509,234.21538 302.36063,279 
     C 302.33536,299.77768 300.97355,312.12979 296.41891,332.89349 
     C 286.70405,377.18157 262.85893,424.36347 228.55502,467.17452 
     C 219.26505,478.76833 199.25099,501.02345 198.17004,500.96183 
     C 197.80179,500.94084 193.58463,497.15559 188.79857,492.55018 
     z M 212.92994,343.99452 C 242.28307,336.85605 266.31414,312.68729 
     273.9846,282.59004 C 276.76052,271.6979 276.75301,253.72727 273.96762,242 
     C 266.78666,211.76606 241.98871,187.12253 211.5005,179.92186 
     C 203.8953,178.12567 200.40831,177.86988 189.0005,178.27134 
     C 173.93019,178.80168 167.30498,180.26871 156.08925,185.55888 
     C 132.8924,196.50023 116.23621,216.81521 109.90648,241.88639 
     C 108.09535,249.06004 107.84969,252.38603 108.2077,264.88639 
     C 108.58615,278.10034 108.93262,280.39476 111.82513,288.842 
     C 113.58452,293.98009 116.23139,300.28009 117.70707,302.842 
     C 137.50495,337.21285 174.70639,353.29022 212.92994,343.99452 z""")
    return hurricane
hurricane = get_hurricane_marker()
hurricane.vertices -= hurricane.vertices.mean()


# #### Sea-ice data functions
# ---

# **get_SOSE_sea_ice**
# 
# This function queries SOSE data given specified region and date.
# 
# xreg, yreg: e.g. [-60.0, -55.0].
# 
# date format: 'yyyy-mm-dd'.
# 
# If printUrl=True, the function print Url for the data query.
def get_SOSE_sea_ice(xreg,yreg,date,printUrl=True):
    # yreg, xreg, date should be lists
    url  = 'https://argovis.colorado.edu/griddedProducts/nonUniformGrid/window?'
    url += 'gridName=sose_si_area_1_day_sparse&presLevel=0&'
    url += 'latRange={}'.format(yreg)
    url += '&lonRange={}'.format(xreg)
    url += '&date={}'.format(date)
    if printUrl:
        print(url)
    resp = requests.get(url)
    # Consider any status other than 2xx an error
    if not resp.status_code // 100 == 2:
        return "Error: Unexpected response {}".format(resp)
    selectionSeaIce = resp.json()    
    return selectionSeaIce


# **parse_into_df_SeaIce**
# 
# This function parses the output of the function get_SOSE_sea_ice.
def parse_into_df_SeaIce(selectionSeaIce):
    meas_keys = selectionSeaIce[0]['data'][0].keys()
    df = pd.DataFrame(columns=meas_keys)
    for data in selectionSeaIce[0]['data']:
        bfrDf = pd.DataFrame(data={'lon':  [data['lon']], 'lat':  [data['lat']], 'value':  [data['value']]}) #pd.DataFrame(data={'lon':  data['lon'], 'lat':  data['lat']})
        df = pd.concat([df, bfrDf], sort=False)
    return df


# #### Argo float data functions
# ---
# **get_selection_profiles**
# 
# This function is from [Tucker, Giglio, Scanderbeg 2020](https://www.essoar.org/doi/10.1002/essoar.10504304.1) and gets profiles in any shape of interest.
# 
# startDate, endDate: 'yyyy-mm-dd'.
# 
# Shape is a list of lists containing [lon, lat] coordinates, e.g. for a squared region:[[[min_longitude,min_latitude],[min_longitude,max_latitude],[max_longitude,max_latitude],[max_longitude,min_latitude],[min_longitude,min_latitude]]].
# 
# For a custom polygon, the user can draw a region using the select region feature in the main Argovis map at https://argovis.colorado.edu: once the shape appears on the map, the corresponding vertices for the polygon appear in URL and can be copied from there to define 'shape' as input.
def get_selection_profiles(startDate, endDate, shape, presRange=None, printUrl=True):
    url = 'https://argovis.colorado.edu/selection/profiles'
    url += '?startDate={}'.format(startDate)
    url += '&endDate={}'.format(endDate)
    url += '&shape={}'.format(shape)
    if presRange:
        pressRangeQuery = '&presRange=' + presRange
        url += pressRangeQuery
    if printUrl:
        print(url)
    resp = requests.get(url)
    # Consider any status other than 2xx an error
    if not resp.status_code // 100 == 2:
        return "Error: Unexpected response {}".format(resp)
    selectionProfiles = resp.json()    
    return selectionProfiles


# **get_profile**
# 
# This function is from [Tucker, Giglio, Scanderbeg 2020](https://www.essoar.org/doi/10.1002/essoar.10504304.1) and gets an Argo float profile of interest.
# 
# profileID format: '5904912_239'.
def get_profile(profileID):
    url = 'https://argovis.colorado.edu/catalog/profiles/{}'.format(profileID)
    resp = requests.get(url)
    # Consider any status other than 2xx an error
    if not resp.status_code // 100 == 2:
        return "Error: Unexpected response {}".format(resp)
    profile = resp.json()
    return profile


# **get_platform_profiles**
# 
# This function is from [Tucker, Giglio, Scanderbeg 2020](https://www.essoar.org/doi/10.1002/essoar.10504304.1) and gets profiles an Argo float of interest.
# 
# platform_number format: '7900379'.
def get_platform_profiles(platform_number):
    url = 'https://argovis.colorado.edu/catalog/platforms/{}'.format(platform_number)
    resp = requests.get(url)
    # Consider any status other than 2xx an error
    if not resp.status_code // 100 == 2:
        return "Error: Unexpected response {}".format(resp)
    platformProfiles = resp.json()
    return platformProfiles


# **parse_into_df**
# 
# This function is from [Tucker, Giglio, Scanderbeg 2020](https://www.essoar.org/doi/10.1002/essoar.10504304.1) and parses profiles from e.g. get_platform_profiles output ('platformProfiles') and get_selection_profiles output ('selectionProfiles') and returns a data frame.
def parse_into_df(profiles):
    meas_keys = profiles[0]['measurements'][0].keys()
    df = pd.DataFrame(columns=meas_keys)
    for profile in profiles:
        profileDf = pd.DataFrame(profile['measurements'])
        profileDf['cycle_number'] = profile['cycle_number']
        profileDf['profile_id'] = profile['_id']
        profileDf['lat'] = profile['lat']
        profileDf['lon'] = profile['lon']
        profileDf['date'] = profile['date']
        profileDf['position_qc'] = profile['position_qc']
        if 'containsBGC' in profile.keys():
            profileDf['containsBGC'] = profile['containsBGC']
        df = pd.concat([df, profileDf], sort=False)
    return df


# **parse_into_df_plev**
# 
# This function is from [Tucker, Giglio, Scanderbeg 2020](https://www.essoar.org/doi/10.1002/essoar.10504304.1) and parses profiles from e.g. get_platform_profiles output ('platformProfiles') and get_selection_profiles output ('selectionProfiles'), and returns a data frame after interpolating profile onto defined pressure levels (e.g. 'plev = np.arange(5,505,5)').
def parse_into_df_plev(profiles, plev):
    plevProfileList = []
    for profile in profiles:
        profileDf_bfr = pd.DataFrame(profile['measurements'])
        plevProfile = profile
        fT = interpolate.interp1d(profileDf_bfr['pres'], profileDf_bfr['temp'], bounds_error=False)
        plevProfile['temp'] = fT(plev)
        # some of the profiles in Argovis may not have salinity 
        # (either because there is no salinity value in the original Argo file or the quality is bad)
        try:
            fS = interpolate.interp1d(profileDf_bfr['pres'], profileDf_bfr['psal'], bounds_error=False)
            plevProfile['psal'] = fS(plev)
        except:
            plevProfile['psal'] = np.nan #  No salinity found in profile
        plevProfile['pres'] = plev
        plevProfileList.append(plevProfile)
    df = pd.DataFrame(plevProfileList)
    df = df.sort_values(by=['cycle_number'])
    df = df.reset_index(drop=True)
    # print all that is available (only some of the info will be stored in the output dataframe, yet users can add more if interested)
    print(df.keys())
    df = df[['cycle_number','_id','date','lon','lat','pres','temp','psal','position_qc','date_qc']]
    return df


# ### Data visualization functions
# ---
# #### Tropical cyclone visualization functions
# ---

# **plot_tracks_time_in_col**
# 
# This function plots Tropical Cyclone tracks ('TCs_Dict') that are output by get_TCs_byDate and get_TCs_byNameYear.
# 
# df_ctag is a tag for the variable plotted in the colorbar: e.g. df_ctag='wind' plots maximum sustained winds for the TC at each position along the track. 
# 
# df_title is a tag for the figure title. If df_title='', the default title is used: 'Tropical Cyclone tracks'.
# 
# tag_TC_or_SH_FILT = 'TC' to map TCs only.
# 
# This function has the capability to map tracks for Southern Hemisphere storms using tag_TC_or_SH_FILT = 'SH_FILT'. Nevertheless, the database for Southern Hemisphere storms is still in development. See Section 1.9.
# 
def plot_tracks_time_in_col(TCs_Dict,df_ctag='wind',df_title='',tag_TC_or_SH_FILT='TC'):
    for i in range(0,len(TCs_Dict)):
        df  = pd.DataFrame(TCs_Dict[i]['traj_data'])
        if (('SH_FILT' == tag_TC_or_SH_FILT) and (tag_TC_or_SH_FILT in TCs_Dict[i]['_id'])) or         (('TC' == tag_TC_or_SH_FILT) and ('SH_FILT' not in TCs_Dict[i]['_id'])):
            plt.scatter(df['lon'],df['lat'],transform=ccrs.PlateCarree(),s=5,c=df[df_ctag],cmap='viridis')#c=mdates.date2num(dti))
    cb = plt.colorbar(orientation='vertical',fraction=0.03,pad=0.02)
    cb.ax.tick_params(labelsize=15)
    cb.set_label('maximum sustained winds, knots', fontsize=16)
    tt = plt.title(df_title,fontsize=24)
    plt.show()


# **TC_and_storms_view**
# 
# Function to map all tropical cyclones for selected time window ('startDate' - 'endDate').
# 
# startDate, endDate: 'yyyy-mm-dd'.
# 
# tag_TC_or_SH_FILT = 'TC' to map TCs only.
# 
# This function has the capability to map tracks for Southern Hemisphere storms using tag_TC_or_SH_FILT = 'SH_FILT'. Nevertheless, the database for Southern Hemisphere storms is still in development. See Section 1.9.
def TC_and_storms_view(startDate,endDate,tag_TC_or_SH_FILT='TC'):
    TCs_Dict = get_TCs_byDate(startDate,endDate=endDate)
    fig = plt.figure(figsize=(15,15))
    ax = plt.axes(projection=ccrs.Mollweide())
    gl = ax.gridlines(draw_labels=True,color='black')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.stock_img()
    bool_list = [] 
    for x in TCs_Dict:
        if ('SH_FILT' in tag_TC_or_SH_FILT):
            df_title = ' Southern Hemisphere storms intensity '
            df_ctag = 'intensity'
            bool_list.append('SH_FILT' in x['_id'])
        else:
            df_title = 'Tropical Cyclone tracks'
            df_ctag = 'wind'
            bool_list.append(~('SH_FILT' in x['_id']))      
    output_select = list(compress(TCs_Dict, bool_list))
    plot_tracks_time_in_col(list(compress(TCs_Dict, bool_list)),df_ctag,df_title,tag_TC_or_SH_FILT)
    return output_select


# **map_TC_and_Argo**

# Co-locate Argo profiles along TC track and map location of profiles and TC track. TC track info is stored in the dataframe 'df' (which is output of get_track_for_storm)
def map_TC_and_Argo(df, delta_days, dx, dy, presRange,printing=False,printing_flag=''):
    prof_beforeTC = []
    prof_afterTC  = []
    col='magenta'
    dx_buffer = 5
    dy_buffer = 5
    fig = plt.figure(figsize=(15,15))
    ax = plt.axes(projection=ccrs.PlateCarree()) #Mollweide
    gl = ax.gridlines(draw_labels=True,color='black')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.coastlines()
    ax.add_feature(cft.LAND)
    ax.add_feature(cft.OCEAN)
    im = plt.scatter(df['lon'],df['lat'],transform=ccrs.PlateCarree(),s=2000,marker=hurricane,
                c=df['wind'], facecolors='none', linewidth=3.5)
    cb = plt.colorbar(orientation='vertical',fraction=0.03,pad=0.02)
    cb.ax.tick_params(labelsize=15)
    cb.set_label('maximum sustained winds, knots', fontsize=16)
    dti = pd.to_datetime(df['timestamp'])
    ax.set_extent([min(df['lon'])-dx_buffer, max(df['lon'])+dx_buffer,
                   min(df['lat'])-dy_buffer, max(df['lat'])+dy_buffer,], crs=ccrs.PlateCarree())
    for i in np.arange(0,len(df['lon']),1):
        for ii in np.arange(0,2,1):
            if ii == 0:# look for "before" profiles
                startDate=str(dti[i]-timedelta(days=delta_days))[0:10]
                endDate=str(dti[i])[0:10]
                mrkr = '*'
                #col = 'k'
            else: #look for "after" profiles
                startDate=str(dti[i])[0:10]
                endDate=str(dti[i]+timedelta(days=delta_days))[0:10]
                mrkr = '*' 
                #col = 'red'
            shape = [[[df['lon'][i]-(dx/2),df['lat'][i]-(dy/2)],[df['lon'][i]-(dx/2),df['lat'][i]+(dy/2)],[df['lon'][i]+(dx/2),df['lat'][i]+(dy/2)],[df['lon'][i]+(dx/2),df['lat'][i]-(dy/2)],[df['lon'][i]-(dx/2),df['lat'][i]-(dy/2)]]]

            strShape = str(shape).replace(' ', '')
            selectionProfiles = get_selection_profiles(startDate, endDate, strShape, str(presRange), printUrl=False)
            selectionProfiles_raw = selectionProfiles
            if len(selectionProfiles) > 0 and not isinstance(selectionProfiles,str):
                selectionDf = parse_into_df(selectionProfiles)
                selectionDf.replace(-999, np.nan, inplace=True)

                gb     = selectionDf.groupby(by='profile_id')
                groups = dict(list(gb))
                gb_list = gb.groups.keys()

                if ii == 0:
                    prof_beforeTC.append(groups)
                else:
                    prof_afterTC.append(groups)
                for tag_id in gb_list:

                    plt.plot(groups[tag_id]['lon'][0],groups[tag_id]['lat'][0],mrkr,transform=ccrs.Geodetic(),markersize=15,linewidth=4,color=col)


                    ax.text(groups[tag_id]['lon'][0]+.25, groups[tag_id]['lat'][0], tag_id,transform=ccrs.Geodetic(),color=col)

            else:
                if ii == 0:
                    prof_beforeTC.append([])
                else:
                    prof_afterTC.append([])
    plt.title('Tropical Cyclone track and location of Argo profiles (magenta)',fontsize=20)
    plt.show()
                
    return prof_beforeTC,prof_afterTC


# #### Argo profiles visualization functions
# ---

# **plot_prof**
# 
# This function plots Argo float profiles.
# 
# dataX is an array vector containing the profile variable. E.g. temperature (df['temp'])  or salinity (df['sal']) at each pressure level from the output of parse_into_df_plev ('df').
# 
# dataY is an array vector containing the profile pressure levels (df['pres']) from the output of  parse_into_df_plev ('df').
# 
# xlab, ylab are the labels for each axis (e.g. xlab = 'Temperature, degC' and ylab=''Pressure, dbar').
# 
# xlim and ylim define the x and y axes limits. If xlim=[], the axis is adjusted automatically to fit the data range, otherwise specified as xlim=[min value, max value]. E.g. xlim=[22,30] for temperature. 
# 
# ylim =[presRange], where presRange=[min value, max_vale]. E.g. presRange=[0,100].
# 
# label sets the label for the plot legend (e.g. label='before' for before TC profiles).
# 
# col sets the color for the profile based on wether the profile was recorded before or after the TC's passage. E.g. col='k' , for black if the profile was recorded before the TC.
def plot_prof(dataX,dataY,xlab,ylab,xlim,ylim,label,col):
    plt.plot(dataX,dataY,label=label,color=col,linewidth=5)
    plt.gca().set_xlabel(xlab,fontsize=24)
    plt.gca().set_ylabel(ylab,fontsize=24)
    plt.ylim(ylim)
    if xlim:
        plt.xlim(xlim)
    plt.gca().invert_yaxis()
    for tick in plt.gca().xaxis.get_majorticklabels():  # example for xaxis
        tick.set_fontsize(24) 
    for tick in plt.gca().yaxis.get_majorticklabels():  # example for xaxis
        tick.set_fontsize(24) 
    return


# **plot_prof_pairs**

# Print ID (i.e. platformNumber_cycleNumber) of oceanic profiles in each item of two lists (prof_beforeTC,prof_afterTC), when the item (i.e. the item that corresponds to a certain index) has profiles in both lists. This function was build to plot profiles before the TC in red and after in black, i.e. when the two lists are indeed for profiles before/after the TC, the plot is done only for locations along the tropical cyclone track of interest where co-located oceanic profiles (stored in prof_beforeTC,prof_afterTC for the TC of interest) are available both before and after the cyclone.
def plot_prof_pairs(prof_beforeTC,prof_afterTC,presRange=[0,100]):
    for x,y in zip(tqdm(prof_beforeTC),prof_afterTC):
        time.sleep(0.01)
        if any(x) and any(y):
            print('-------------')
            # temperature
            fig = plt.figure(figsize=(30,10))
            plt.subplot(121)
            for d in x.keys():
                print('Temperature, before (black): ' + d)
                plot_prof(dataX=x[d]['temp'],dataY=x[d]['pres'],xlab='Temperature, degC',ylab='Pressure, dbar',xlim=[],ylim=presRange,label='before',col='k')
            for d in y.keys():
                print('Temperature, after (red): ' + d)
                plot_prof(dataX=y[d]['temp'],dataY=y[d]['pres'],xlab='Temperature, degC',ylab='Pressure, dbar',xlim=[],ylim=presRange,label='after',col='r')
            plt.title('Temperature profiles',fontsize=24)
            plt.legend(fontsize=18)
            # salinity
            plt.subplot(122)
            for d in x.keys():
                print('Salinity, before (black): ' + d)
                try:
                    plot_prof(dataX=x[d]['psal'],dataY=x[d]['pres'],xlab='Salinity, psu',ylab='Pressure, dbar',xlim=[],ylim=presRange,label='before',col='k')
                except:
                    pass
            for d in y.keys():
                print('Salinity, after (red): ' + d)
                try:
                    plot_prof(dataX=y[d]['psal'],dataY=y[d]['pres'],xlab='Salinity, psu',ylab='Pressure, dbar',xlim=[],ylim=presRange,label='after',col='r')
                except:
                    pass
            plt.title('Salinity profiles',fontsize=24)
            plt.legend(fontsize=18)
            plt.show()

# **plot_SeaIce_argo_QC_temp_sal**
#After adding empty profiles to the float history as needed, we will generate a multi-panel figure showing: a) the time series of the float's position QC flag value and co-located fraction of sea-ice from SOSE, b) ocean temperature profiles in time as recorded by the float and, c) ocean salinity profiles in time as recorded by the float:
def plot_SeaIce_argo_QC_temp_sal(plevIntp,temp2d,psal2d,platformDf_plev,tdelta,seaice_val,platform_number):
    temp2d_new = np.empty((len(plevIntp),1))
    psal2d_new = np.empty((len(plevIntp),1))
    temp2d_new[:,0] = temp2d[:,0]
    psal2d_new[:,0] = psal2d[:,0]

    time_new = pd.Series(pd.to_datetime(platformDf_plev['date'][0]))
    for i in np.arange(1,len(platformDf_plev['date'])):
            nnans_days = (pd.to_datetime(platformDf_plev['date'][i])-pd.to_datetime(platformDf_plev['date'][i-1])).days
            # in the following we add 1 profile of nan in between profiles that are too far apart in time 
            # (if more profiles of nan are needed, then you need to edit np.empty and the number of iterations in j-loop)
            if nnans_days>tdelta:

                temp2d_new = np.append(temp2d_new,np.nan*np.empty((len(plevIntp),1)),axis=1)
                psal2d_new = np.append(psal2d_new,np.nan*np.empty((len(plevIntp),1)),axis=1)

                bfr_date = pd.to_datetime(platformDf_plev['date'][i-1])
                # this loop is now only looping once as we only need to add one profile of nans 
                for j in np.arange(1,2,1):
                    bfr_date += timedelta(days=1)
                    time_new = time_new.append(pd.Series(bfr_date))

            temp2d_new = np.append(temp2d_new,temp2d[:,np.newaxis,i],axis=1)
            psal2d_new = np.append(psal2d_new,psal2d[:,np.newaxis,i],axis=1)
            time_new = time_new.append(pd.Series(pd.to_datetime(platformDf_plev['date'][i])))

    plt.figure(figsize=(30, 20))
    for i in [1,2,3]:
        if i==1:
            ax = plt.subplot(311)
            color = 'tab:red'
            plt.plot(pd.to_datetime(platformDf_plev['date']),platformDf_plev['position_qc'],'sr',markersize=14)
            ax.set_ylabel('Position QC flag',size=24,labelpad=0)
            ax.yaxis.label.set_color('red')
            ax.tick_params(axis='y', colors='red')
            ax2 = ax.twinx()
            color = 'tab:blue'
            ax2.plot(pd.to_datetime(platformDf_plev['date']),np.stack( seaice_val ),'ob',markersize=12)
            ax2.set_ylabel('SOSE sea-ice fraction',size=24,labelpad=0)
            ax2.yaxis.label.set_color('blue')
            ax2.tick_params(axis='y', colors='blue')
            plt.title('Float #'+platform_number,size=24)
            plt.grid()

        if i==2:
            ax = plt.subplot(312)
            mt = np.ma.masked_where(np.isnan(temp2d_new),temp2d_new)
            plt.pcolor(pd.to_datetime(time_new),plevIntp,mt,cmap='plasma')
            #plt.pcolormesh(pd.to_datetime(platformDf_plev['date']),plevIntp,temp2d) # pd.Series(np.arange(0,79,1))
            #plt.title('Temperature, degC',size=24)
            cbar = plt.colorbar(orientation="horizontal", pad=0.2)
            cbar.ax.tick_params(labelsize=24)
            cbar.set_label('Temperature, degC',fontsize=24)

        if i==3:
            ax = plt.subplot(313)
            ms = np.ma.masked_where(np.isnan(psal2d_new),psal2d_new)
            plt.pcolor(pd.to_datetime(time_new),plevIntp,ms,cmap='viridis_r')
            # plt.pcolor(pd.to_datetime(platformDf_plev['date']),plevIntp,psal2d) # pd.Series(np.arange(0,79,1))
            #plt.title('Salinity, psu',size=24)
            cbar = plt.colorbar(orientation="horizontal", pad=0.2)
            cbar.ax.tick_params(labelsize=24)
            cbar.set_label('Salinity, PSU',fontsize=24)

        if i==2 or i==3:
            plt.gca().invert_yaxis()
            ax.set_ylabel('Pressure, dbar',size=24,labelpad=0)


        # change font
        for tick in ax.xaxis.get_majorticklabels():  # example for xaxis
            tick.set_fontsize(24) 
        for tick in ax.yaxis.get_majorticklabels():  # example for xaxis
            tick.set_fontsize(24)
        for tick in ax2.yaxis.get_majorticklabels():  # example for xaxis
            tick.set_fontsize(24)
        plt.xlim([min(pd.to_datetime(platformDf_plev['date'])),max(pd.to_datetime(platformDf_plev['date']))])

    plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4, 
                        hspace=0.4)
    plt.show()

#map sea-ice and argo float locations
#This function downloads SOSE sea-ice (see: get_SOSE_SeaIce), parses the data (see: parse_into_df_SeaIce), then downloads and parses Argo profiles (see: get_selection_profiles) and generates a map of sea-ice fraction and Argo float locations.
   # presRange = '[0,50]', defines pressure range of interest for Argo profiles.
  
   # plev_presRange = 30, is a value in the middle of the pressure range of interest.
  
   # color_map = plt.cm.get_cmap('Blues'), select colormap of preference.

def map_seaice_argo(date_ALL,yreg_ALL,xreg_ALL,delta_argo,presRange,plev_presRange,color_map):
    reversed_color_map = color_map.reversed()
    for l in np.arange(0,len(date_ALL)):
        fig = plt.figure(figsize=(15,15))
        ax = plt.axes(projection=ccrs.SouthPolarStereo())
        gl = ax.gridlines(draw_labels=True,color='black')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        ax.coastlines()
        ax.add_feature(cft.LAND)
        ax.add_feature(cft.OCEAN)
        ax.set_extent([-180, 180, -90, -35], crs=ccrs.PlateCarree())

        # dates for argo profiles
        bfr_startDate = datetime.strptime(date_ALL[0], '%Y-%m-%d')
        bfr_startDate -= timedelta(days=delta_argo)
        bfr_endDate = datetime.strptime(date_ALL[0], '%Y-%m-%d')
        bfr_endDate += timedelta(days=delta_argo)

        str_flag = True
        for i in tqdm(np.arange(0,len(xreg_ALL)-1)):
            time.sleep(0.01)
            for j in np.arange(0,len(yreg_ALL)-1):
                try:
                    selectionSeaIce  = get_SOSE_sea_ice(xreg=xreg_ALL[i:i+2],yreg=yreg_ALL[j:j+2],date=date_ALL[l],printUrl=False)
                    df_SeaIce       = parse_into_df_SeaIce(selectionSeaIce)
                    plt.scatter(df_SeaIce['lon'],df_SeaIce['lat'],transform=ccrs.Geodetic(),s=5,cmap=reversed_color_map,c=df_SeaIce['value'])

                except:
                    pass
                # add Argo profiles to the plot and print float number
                # region
                shape = [[[max(xreg_ALL[i:i+2]),min(yreg_ALL[j:j+2])],[max(xreg_ALL[i:i+2]),max(yreg_ALL[j:j+2])],
                          [min(xreg_ALL[i:i+2]),max(yreg_ALL[j:j+2])],[min(xreg_ALL[i:i+2]),min(yreg_ALL[j:j+2])],
                          [max(xreg_ALL[i:i+2]),min(yreg_ALL[j:j+2])]]]
                strShape = str(shape).replace(' ', '')

                bfr_selectionProfiles     = get_selection_profiles(bfr_startDate, bfr_endDate, strShape, presRange,printUrl=False)

                if len(bfr_selectionProfiles) > 0:
                    bfr_selectionDf = parse_into_df(bfr_selectionProfiles)
                    bfr_selectionDf = bfr_selectionDf.drop_duplicates('profile_id')
                    plt.plot(bfr_selectionDf['lon'],bfr_selectionDf['lat'],'+r',transform=ccrs.Geodetic())                
                    bfr2 = bfr_selectionDf.loc[bfr_selectionDf['position_qc'] == 8]
                    if len(bfr2) > 0:
                        plt.plot(bfr2['lon'],bfr2['lat'],'or',transform=ccrs.Geodetic())
                        if str_flag:
                            print('>>>> Profile ID for profiles under ice (i.e. platformNumber_profileNumber) <<<<')
                            str_flag = False
                        print(bfr2['profile_id'])
        plt.title('SOSE sea-ice fraction and Argo profile locations (+: QC =1, o: QC = 8)',fontsize=20)              
        colorbar =plt.colorbar(orientation='vertical',fraction=0.05,pad=0.02)
        colorbar.ax.tick_params(labelsize=15)
        colorbar.set_label('sea-ice fraction', fontsize=16)

### plot_argo_QC_location###
   # This function plots the location of a user specified Argo float colorcoded by the value of its QC flag at each position. The function uses parsed data(parse_into_df_plev) queried via get_profile_platform. 
   
   # lon: variable name for longitudes of each recorded position of the Argo float. E.g. platformDf_plev['lon']
   
   # lat: variable name for latitudes of each recorded position of the Argo float. E.g. platformDf_plev['lat']
   
   # qc: variable name for the QC flag value at each recorded position of the Argo float. E.g. platformDf_plev['position_qc']


def plot_argo_QC_location(platform_number,lon,lat,qc):
    fig, (ax) = plt.subplots(1,1, figsize=(15,10))
    plt.scatter(lon,lat,c=qc,s=80)
    plt.colorbar()
    # change font
    for tick in ax.xaxis.get_majorticklabels():  # example for xaxis
        tick.set_fontsize(24) 
    for tick in ax.yaxis.get_majorticklabels():  # example for xaxis
        tick.set_fontsize(24)
    ax.set_ylabel('Degree latitude',size=24,labelpad=0)
    ax.set_xlabel('Degree longitude',size=24,labelpad=0)
    plt.title('Float #'+platform_number,size=24)
    plt.grid()
    plt.show()
