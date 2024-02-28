import itertools
import os
# math \ plotting:
import numpy as np
from scipy import stats as st
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from cycler import cycler
# spectral iamges:
import spectral
from spectral import envi
# geospatial:
import rasterio
from rasterio import features as feats 
from rasterio.transform import from_origin
import geopandas as gpd
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import utm
import folium
from shapely.geometry import shape


class geo:
    def __init__(self, image_hdr_filename, verbose = True):
        
        # load the image
        im = envi.open(image_hdr_filename)
        coefIm = im.load()

        # build the band data
        for band_idx in range(coefIm.shape[2]):
            name = im.metadata['band names'][band_idx]
            output_tiff_fname = 'abund_'+name+'.tiff'
            abund = np.array(np.squeeze(coefIm[:,:,band_idx]))
            abund[np.where(abund<0)] = 0
            abund = (255*(abund/np.max(abund))).astype(np.uint8)
            
            # plot data with image
            if verbose:            
                plt.figure(figsize=(12,4))
                ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=2, rowspan=1)
                ax2 = plt.subplot2grid((2, 3), (1, 0), colspan=2, rowspan=1)
                ax3 = plt.subplot2grid((2, 3), (0, 2), colspan=1, rowspan=2)
                plt.tight_layout()
                ax1.imshow(np.rot90(abund>0));
                ax1.axes.get_xaxis().set_ticks([])
                ax1.axes.get_yaxis().set_ticks([])
                ax2.imshow(np.rot90(abund));
                ax2.axes.get_xaxis().set_ticks([])
                ax2.axes.get_yaxis().set_ticks([])
                ax3.hist(abund[np.where(abund>0)].flatten(), bins=256);
                ax3.set_title('Abundance Histogram')

            # save the geotiff of abundance
            im_gdal = gdal.Open(im.filename)
            self.write_geotiff(abund, im_gdal, output_tiff_fname, nbands=4, nodata=False, color_idx = band_idx, dtype=gdal.GDT_Byte)

    def write_geotiff(self, array, gdal_obj, outputpath, dtype=gdal.GDT_UInt16, options=0, color_idx = 0, nbands=4, nodata=False):
        """
        Description:
        Writes a geotiff from a Numpy array with appended georeferencing from parent geotiff.
        
        Parameters:
        array: numpy array to write as geotiff
        gdal_obj: object created by gdal.Open() using a tiff that has the SAME CRS, transformation, and resolution as the array you're writing
        outputpath: path including filename.tiff
        dtype (OPTIONAL): datatype to save as
        nodata (default: False): set to any value you want to use for nodata; if False, nodata is not set
        """

        gt = gdal_obj.GetGeoTransform()

        # Prepare destination file
        width = np.shape(array)[1]
        height = np.shape(array)[0]
        driver = gdal.GetDriverByName("GTiff")
        if options != 0:
            dest = driver.Create(outputpath, width, height, 4, dtype, options)
        else:
            dest = driver.Create(outputpath, width, height, 4, dtype)

        # Prepare color table
        mode = st.mode(array.flatten()[array.flatten()>20])
        val1 = int(np.min([2*mode[0],255]))
        val2 = int(np.min([4*mode[0],255]))
        r,g,b = self.cmap(color_idx)
        colors = gdal.ColorTable()
        colors.CreateColorRamp(0, (r, g, b, 50), val1, (r, g, b, 75))
        colors.CreateColorRamp(val1,  (r, g, b, 75), val2,  (r, g, b, 128)) 
        colors.CreateColorRamp(val2,  (r, g, b, 128), 255,  (r, g, b, 128))

        # Write output raster
        arrayR = np.zeros((array.shape[0],array.shape[1]), np.uint8)
        arrayG = np.zeros((array.shape[0],array.shape[1]), np.uint8)
        arrayB = np.zeros((array.shape[0],array.shape[1]), np.uint8)
        arrayA = np.zeros((array.shape[0],array.shape[1]), np.uint8)
        for i in range(1,256):
            RGBA = colors.GetColorEntry(i)
            idx = np.where(array==i)
            arrayR[idx] = RGBA[0]
            arrayG[idx] = RGBA[1]
            arrayB[idx] = RGBA[2]
            arrayA[idx] = RGBA[3]
        dest.GetRasterBand(1).WriteArray(arrayR)
        dest.GetRasterBand(2).WriteArray(arrayG)
        dest.GetRasterBand(3).WriteArray(arrayB)
        dest.GetRasterBand(4).WriteArray(arrayA)
        dest.GetRasterBand(1).SetRasterColorInterpretation(gdal.GCI_RedBand)
        dest.GetRasterBand(2).SetRasterColorInterpretation(gdal.GCI_GreenBand)
        dest.GetRasterBand(3).SetRasterColorInterpretation(gdal.GCI_BlueBand)
        dest.GetRasterBand(4).SetRasterColorInterpretation(gdal.GCI_AlphaBand)

        # Set transform and projection
        dest.SetGeoTransform(gt)
        wkt = gdal_obj.GetProjection()
        srs = osr.SpatialReference()
        srs.ImportFromWkt(wkt)
        dest.SetProjection(srs.ExportToWkt())
        
        # Close output raster dataset 
        dest = None
        
    def build_colormap(self):
        # build a list of colors to iterate through
        
        self.colors = []
        cmap = matplotlib.colormaps['Set1']  # type: matplotlib.colors.ListedColormap
        for c in cmap.colors:
            self.colors.append(c)
        cmap = matplotlib.colormaps['Dark2']  # type: matplotlib.colors.ListedColormap
        for c in cmap.colors:
            self.colors.append(c)
        cmap = matplotlib.colormaps['tab10']  # type: matplotlib.colors.ListedColormap
        for c in cmap.colors:
            self.colors.append(c)
    
    def cmap(self, i):
        # choose next color in colormap
        if 'self.colors' not in locals():
            self.build_colormap()
        RGB_float = self.colors[i%len(self.colors)]
        RGB_uint8 = np.zeros(3, dtype=np.uint8)
        for i in range(2):
            RGB_uint8[i] = int(255*RGB_float[i])
        return RGB_uint8
        