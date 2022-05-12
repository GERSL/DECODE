# Interpolate the tidal information
import gdal
import pandas as pd
import os
import osgeo.gdalconst as gdalconst
pathIn = '/home/xiy19029/DECODE_github/prepareWL/'
pathOut = '/shared/zhulab/Yang/dailyWL/'

if not os.path.exists(pathOut):
    os.makedirs(pathOut)

dailyWLInfo = 'dailyWL_latlon_Alaska.csv'
cols = pd.read_csv(pathIn+dailyWLInfo, nrows=1).columns
WLInfo = pd.read_csv(pathIn+dailyWLInfo, usecols=cols[3:])
lat_lon = pd.read_csv(pathIn+dailyWLInfo, usecols=[1, 2])
flag = 0  # generate the vrt file one time
for index in WLInfo.columns:
    # print(wlInfo[index])
    # iterate the doy
    doySpecific = pd.DataFrame(WLInfo[index])
    doywl = doySpecific.rename(columns={index: 'wl'}, inplace=False)
    # ptDoyWl = pd.concat([lat_lon, doywl*1000], axis=1, sort=False)
    ptDoyWl = pd.concat([lat_lon, doywl*1000], axis=1)

    ptDoyWl.rename(index={0: "lat", 1: "lon", 2: "wl"})
    ptDoyWl.to_csv("wlDoyInfo.csv", index=False)  # export the temporal doy water level info to a csv file
    if flag == 0:
        with open("wlDoyInfo.vrt", 'w') as fn_vrt:
            fn_vrt.write('<OGRVRTDataSource>\n')
            fn_vrt.write('\t<OGRVRTLayer name="wlDoyInfo">\n')
            fn_vrt.write('\t\t<SrcDataSource>wlDoyInfo.csv</SrcDataSource>\n')
            fn_vrt.write('\t\t<GeometryType>wkbPoint</GeometryType>\n')
            fn_vrt.write('\t\t<LayerSRS>WGS84</LayerSRS>\n')
            fn_vrt.write('\t\t<GeometryField encoding="PointFromColumns" x="lon" y="lat" z="wl"/>\n')
            fn_vrt.write('\t</OGRVRTLayer>\n')
            fn_vrt.write('</OGRVRTDataSource>\n')
    else:
        flag = 1
        #    gdal.Grid(pathOut+'wlDoy' + index + '.tif', 'wlDoyInfo.vrt', outputType = gdalconst.GDT_Int32, algorithm = 'invdist:power=2',
        #             outputBounds=[-78, 36, -64, 46], width=1400, height=1000)
 #   gdal.Grid(pathOut + 'wlDoy' + index + '.tif', 'wlDoyInfo.vrt', outputType=gdalconst.GDT_Int32,algorithm='invdist:power=2',
  #            outputBounds=[-100, 22, -74, 38], width=2600, height=1600)
    gdal.Grid(pathOut + 'wlDoy' + index + '.tif', 'wlDoyInfo.vrt', outputType=gdalconst.GDT_Int32,algorithm='invdist:power=2',
              outputBounds=[-180, 50, -130, 74], width=5000, height=2400)
    print(pathOut + 'wlDoy' + index + '.tif')