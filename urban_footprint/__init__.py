from pathlib import Path

import pyproj
import numpy as np
import pandas as pd
import geopandas as gpd
import pickle

import rasterio
import rasterio.mask

import shapely
from shapely.ops import transform
from shapely.geometry import Point, Polygon
import osmnx as ox

from scipy.ndimage import generic_filter

from .libs.bufferByPercentage import *

###

# list of raster data providers with tile schemas
RASTER_PROVIDERS = {
    'ghsl': {'tile_schema': Path('/home/users/lbarleta/leo_local/urban_footprint/datasets/ghsl_tile_schema.geojson')},
    'provider2': {'url': ''}
}  # TODO add relative path

STANDARD_BUFFER = 5000  # in meters
STANDARD_CRS = pyproj.CRS('epsg:3857')
NODATA = 255
VERBOSE = True

# TMP_FOLDER = Path('tmp')


METHODS = {
    'nyu': 'NYU Urban Extent (Angel at al, 20xx)',
    'nighttime': 'Nighttime Light',
}

RADIUS = 564  # meters

ghsl_dict = {
    1975: {'built': [6], 'not_built': [2, 3, 4, 5]},
    1990: {'built': [5, 6], 'not_built': [2, 3, 4]},
    2000: {'built': [4, 5, 6], 'not_built': [2, 3]},
    2014: {'built': [3, 4, 5, 6], 'not_built': [2]}
}

###

control = 0


class Utils:
    @staticmethod
    def load_raster(raster_path, ret='dataset'):
        rst = rasterio.open(raster_path)
        return rst if ret == 'dataset' else rst.read()

    @staticmethod
    def save_raster(raster_data, raster_metadata, file):
        # add a third dimensions in case its 2d
        raster_data = np.array([raster_data]) if raster_data.ndim == 2 else raster_data

        with rasterio.open(file, 'w', **raster_metadata) as dst:
            dst.write(raster_data)
        return True

    @staticmethod
    def save_geometry_to_vector(geometry, crs, file):
        gpd.GeoSeries.from_wkt(
            [geometry.wkt],
            name='geometry',
            crs=crs).to_frame().to_file(file)
        #
        # gpd.GeoDataFrame(
        #     [geometry],
        #     columns=['geometry'],
        #     crs=crs
        # ).to_file(file)

    @staticmethod
    def clip_raster(raster, mask):
        out_image, out_transform = rasterio.mask.mask(raster,
                                                      [shapely.geometry.mapping(mask)],
                                                      crop=True,
                                                      nodata=NODATA)

        out_meta = raster.meta
        out_meta.update({"driver": "GTiff",
                         "height": out_image.shape[1],
                         "width": out_image.shape[2],
                         "transform": out_transform,
                         "nodata": NODATA})

        return out_image, out_meta

    @staticmethod
    def reclassify_raster(raster, reclass_dict):
        """
        reclassify a raster according to dictionary
        """
        reclassified = np.full_like(raster, fill_value=NODATA)

        for new_val, old_vals in reclass_dict.items():
            for i in old_vals:
                reclassified = np.where(raster == i, new_val, reclassified)

        return reclassified

    @staticmethod
    def reclassify_raster_minmax(raster, minmax_dict, fill_val=NODATA):
        """
        reclassify a raster according to dictionary min and max values
        """
        reclassified = np.full_like(raster, fill_value=fill_val)

        for new_val, old_vals in minmax_dict.items():
            cond1 = raster >= old_vals['min'] if old_vals['min'] == 0 else raster > old_vals['min']
            cond2 = raster <= old_vals['max']
            reclassified = np.where(cond1 & cond2, new_val, reclassified)

        return reclassified.astype(int)

    @staticmethod
    def mean_func(input, nodata=NODATA):
        '''
        This function calculates the mean of a given
        np.array, excluding a given 'nodata' value
        '''
        global control

        new_input = np.delete(input, np.where(input == nodata))
        unique = np.unique(new_input)

        control += 1
        if control % 100000 == 0:
            print(f'      Reached {control:,}')

        return nodata if len(unique) == 0 else np.mean(new_input)


class BaseRaster:
    def __init__(self,
                 file_path=None,
                 dir_path=None,
                 data_provider=None):
        if VERBOSE: print("Loading raster dataset...")

        # load and clip files to the observation area
        try:
            if file_path:
                self.raw_raster = Utils.load_raster(file_path)

            # # This option is not implemented; merge large files before using the library
            # elif dir_path:
            #     # load raster files in the folder
            #     # TODO clip parts of each raster and THEN merge them all, to save processing
            #     base_raster = self.load_raster_from_dir(dir_path, observation_area)
            # elif data_provider:
            #     if data_provider in RASTER_PROVIDERS:
            #         base_raster = self.download_data(data_provider, observation_area)
            # TODO CRS check and, potentialluy, conversion

            if VERBOSE: print("Raster loaded successfully.")
        except Exception as e:
            print("Unable to load raster data.")
            print(e)

    #
    #
    #
    # def bbox_to_polygon(self, r):
    #     return Polygon([(r.left, r.top),
    #                     (r.right, r.top),
    #                     (r.right, r.bottom),
    #                     (r.left, r.bottom)])
    #
    # def load_raster_from_dir(self, dir_path, obs):
    #     dir_path = Path(dir_path) if isinstance(dir_path, Path) is False else dir_path
    #     if dir_path.is_dir():
    #         # loading rasters and getting their bounds
    #         raster_paths = [i for i in dir_path.glob('*.tif')]
    #         rasters_bboxes = {}
    #         for raster_path in raster_paths:
    #             with rasterio.open(raster_path) as rst:
    #                 # TODO validate CRS
    #                 rasters_bboxes[raster_path.name] = {
    #                     'bounds': self.bbox_to_polygon(rst.bounds),
    #                     'raster_file': raster_path
    #                 }
    #
    #         # checking overlap with observation area
    #         # gdf = gpd.GeoDataFrame.from_dict(rasters_bboxes, orient='index')
    #         # gdf.set_geometry('bounds', crs=STANDARD_CRS, inplace=True, drop=True)
    #         # gdf.geometry.intersects(obs_area.geometry)
    #
    #         # removing rasters that don't overlap
    #         overlap = []
    #         for index, row in rasters_bboxes.items():
    #             if row['bounds'].intersects(obs.geometry) is True:
    #                 overlap.append(row['raster_file'])
    #
    #         if len(overlap) > 1:
    #             # merge!
    #             pass
    #         elif len(overlap) == 1:
    #             # TODO test this
    #             with rasterio.open(overlap[0]) as rst:
    #                 self.base_raster = rst.read(0)
    #
    #     print(overlap)
    #     return overlap
    #     # find overlap
    #     # merge if necessary
    #     # self.base_raster
    #     return True
    #
    # # TODO refactor the download with load from directory method
    # def download_data(self, provider, obs):
    #     # downloading tile schema mapping tiled-files
    #     print(RASTER_PROVIDERS[provider])
    #     if 'tile_schema' in RASTER_PROVIDERS[provider]:
    #         # check intersection between observation area and tile schema
    #         tile_schema = gpd.read_file(RASTER_PROVIDERS[provider]['tile_schema'])
    #         to_download = tile_schema[tile_schema.geometry.intersects(obs.geometry)].ftp_path.tolist()
    #     elif 'url' in RASTER_PROVIDERS[provider]:
    #         to_download = [RASTER_PROVIDERS[provider]['url']]
    #
    #     # TODO Continue from here:
    #     # Find the actual tile schema for the GHS_BUILT
    #     # Download
    #     # Merge
    #     # Clip
    #     # Reclassify
    #
    #     print(to_download)
    #     return False


class ObservationArea:
    def __init__(self, polygon, name, base_dir, orig_crs=None):
        self.name = name
        self.raster_data = None
        self.raster_meta = None
        self.processed_rasters = None

        if VERBOSE: print("Loading observation area...")
        try:
            if polygon.is_valid:

                # check MultiPolygon
                if polygon.geom_type == 'MultiPolygon':
                    if len(polygon.geoms) > 1 and VERBOSE: print('   MultiPolygon detected! Using just the fist Polygon')
                    polygon = polygon.geoms[0]

                if orig_crs and orig_crs != STANDARD_CRS:
                    if VERBOSE: print("   Adjusting projection")
                    project = pyproj.Transformer.from_proj(orig_crs, STANDARD_CRS)

                    reprojected = transform(project.transform, polygon)  # apply projection
                    self.geometry = reprojected
                else:
                    self.geometry = polygon

                # TODO Check for MULTIPOLYGON

            if VERBOSE: print("Observation area loaded successfully!")
        except Exception as e:
            print("Invalid shapely Polygon object.")
            print(e)

    def process_raster(self,
                   base_raster,
                   buffer=STANDARD_BUFFER,
                   multi_year=None):
        """
        Clip base raster to the size of the observation area and reclassify according to dict
        """
        # clip raster to observation area
        self.raster_data, self.raster_meta = Utils.clip_raster(base_raster.raw_raster, self.geometry.buffer(buffer))

        # disaggregate multi-year data
        if multi_year:
            new_classes = {}
            self.processed_rasters = {}
            for year, classes in multi_year.items():
                new_classes[0] = classes['not_built']  # set not built to 0
                new_classes[1] = classes['built']  # set built to 1
                self.processed_rasters[year] = Utils.reclassify_raster(self.raster_data, new_classes)

        if VERBOSE: print("Raster pre-processed!")


class NyuFootprint:
    urban_classes = {
        1: {'min': 0, 'max': 0.25, 'label': 'rural'},
        2: {'min': 0.25, 'max': 0.5, 'label': 'suburban'},
        3: {'min': 0.5, 'max': 1, 'label': 'urban'},
    }

    def __init__(self, name, raster_data, raster_meta):
        self.name = name
        self.raster_data = raster_data
        self.raster_meta = raster_meta

        # Variables that will be calculated in the process
        self.urbanization_level = None
        self.urban_bands = None
        self.fringe = None
        self.captured = None
        self.footprint = None
        self.metrics = {}
        self.final = None

    def standard_method(self):
        self.calculate_urban_level()
        self.categorize_urban_level()
        self.rectify_open_spaces()
        self.rectify_components()

        #TODO check article to see what they do with the bands
        #TODO rename this
        self.final = pd.concat([
            self.urban_bands,
            self.captured,
            self.fringe
        ])

        if self.footprint:
            self.calculate_metrics()

    def calculate_urban_level(self, radius=RADIUS):
        global control

        # calculate the radius
        raster_resolution = self.raster_meta['transform'][0]
        actual_radius = round(radius / raster_resolution)

        if VERBOSE: print(
            f'   The actual radius is {int(actual_radius * raster_resolution):,} meters instead of {radius} meters, '
            f'given the resolution of the raster.')

        # filtering the pixels within the radius
        kernel = np.zeros((2 * actual_radius + 1, 2 * actual_radius + 1))
        y, x = np.ogrid[-actual_radius:actual_radius + 1, -actual_radius:actual_radius + 1]
        mask = x ** 2 + y ** 2 <= actual_radius ** 2
        kernel[mask] = 1

        if VERBOSE: print(f'   There are {self.raster_data.shape[1] * self.raster_data.shape[2]:,} calculations '
                          f'to be made. It might take a feel minutes.')

        control = 0
        self.urbanization_level = generic_filter(self.raster_data[0].astype(float),
                                                 Utils.mean_func,
                                                 footprint=kernel)

    def categorize_urban_level(self, minmax_dict=None):
        minmax_dict = self.urban_classes if minmax_dict is None else minmax_dict
        reclassified = Utils.reclassify_raster_minmax(self.urbanization_level, minmax_dict)

        # polygonizing raster
        features = [i for i in rasterio.features.shapes(reclassified,
                                                        mask=reclassified != NODATA,
                                                        transform=self.raster_meta['transform'])]
        collection = [(minmax_dict[f[1]]['label'], shapely.geometry.shape(f[0])) for f in features]

        self.urban_bands = gpd.GeoDataFrame(collection,
                                            columns=['class', 'geometry'],
                                            crs=self.raster_meta['crs'])
        self.urban_bands.set_index('class', inplace=True)

        # excluding 'rural' classification as it doesn't have a particular shape
        if 'rural' in self.urban_bands.index:
            self.urban_bands.drop('rural', inplace=True)

    def rectify_open_spaces(self, dist_fringe=100, area_captured=200):
        # Open Spaces
        # # Fringe open spaces: open spaces (rural) within 100 m of urban and suburban areas
        # # Captured open spaces: small (less than 200 ha) areas surrounded by urban, suburban, or fringe open spaces
        # # # TODO Question: what about large bodies of water within the city?
        # # Rural open spaces: all the rest, not included
        # Note: fringe and captured open spaces might overlap but they are merged at the end to become the footprint

        dissolved = self.urban_bands.dissolve()

        # include fringe open spaces
        fringe = (dissolved
                  .buffer(dist_fringe, single_sided=True)  # apply a buffer
                  .difference(dissolved))  # calculate the different as the fringe area

        fringe = fringe[fringe.is_empty == False].to_frame(name='geometry')
        fringe['class'] = 'fringe'
        fringe.set_index('class', inplace=True)

        self.fringe = fringe.explode(index_parts=False).geometry.to_frame()

        # include captured open spaces
        holes = dissolved.explode(index_parts=False).interiors.to_list()  # for each component, list the "holes"

        captured = gpd.GeoDataFrame(
            [('captured', Polygon(item)) for sublist in holes for item in sublist],
            columns=['class', 'geometry'],
            crs=self.raster_meta['crs']
        )

        captured = captured[captured.area / 10000 <= area_captured]
        # captured['class'] = 'captured'
        captured.set_index('class', inplace=True)

        self.captured = captured[captured.area / 10000 <= area_captured]

    def rectify_components(self, buffer_prop=1.25):
        # Urban Clusters
        # # Urban clusters: union of contiguous urban and suburban areas and fringe and captured urban spaces
        # # Main urban cluster: largest (in terms of area) clusters
        # # Inclusion rule: union of clusters which buffer -- equaling to one-quarter of the total area of a
        # # given cluster -- are connected to the main cluster

        components = (pd.concat([self.urban_bands, self.fringe, self.captured])
                      .dissolve()
                      .explode(index_parts=False)
                      .reset_index(drop=True)
                      )

        # create a buffer with one-quarter of the area of the clusters
        buffered = components.geometry.apply(bufferByPercentage, args=[buffer_prop, False])

        # finding main cluster and those which intercepts with it
        if len(components) > 0:
            #TODO throw a warning
            self.footprint = (
                components[
                    buffered.intersects(components.at[components.area.idxmax(), 'geometry'])
                ].unary_union
            )

    def calculate_metrics(self):
        overlap = self.urban_bands[self.urban_bands.intersects(self.footprint)].reset_index()
        self.overlap = overlap

        self.metrics['urban_area_ha'] = round(overlap[overlap['class'] == 'urban'].area.sum() / 10000, 2)
        self.metrics['suburban_area_ha'] = round(overlap[overlap['class'] == 'suburban'].area.sum() / 10000, 2)
        self.metrics['centroid'] = self.footprint.centroid


class StreetNetwork:
    def __init__(self, query=None):
        self.query = query

    # def download_road_network(self, clip_area):
