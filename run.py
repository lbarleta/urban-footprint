from urban_footprint import *
from slugify import slugify
from pathlib import Path
import pickle

'''
Under development

TODO:
- Refactor code below to include most of the steps inside the Python package
- Add option to save intermediary files
- Improve VERBOSE mode
- Create notebook with instructions
- Improve comments and documentation
- Include option to automatically download (and merge) raster files from the cloud
'''

if __name__ == '__main__':

    ############
    # Settings #
    ############

    # Path to base directory for inputs and outputs
    base_dir = Path('/')
    
    # Path to GHS_BUILT raster (make sure it's in the right projection)
    # Originally developed for GHS-BUILT-LAU2STAT R2022A: https://data.jrc.ec.europa.eu/dataset/94d62a61-25d0-42fd-9e1e-a41f877cf788
    # TODO: test and adapt to the newest release (GHS-BUILT-S - R2023A)
    ghs_file = base_dir / 'ghs_built_3857.tif'

    # Path to vector file with observation areas
    obs_file = base_dir / 'observations_3857.geojson'

    VERBOSE = True
    
    ####################
    # Batch Processing #
    ####################

    if VERBOSE: print(f'Start batch processing')

    # load raster and vector files
    ghs_raster = BaseRaster(ghs_file)
    observations = gpd.read_file(obs_file) 


    for i, obs in observations.iterrows():
        city_name = slugify(obs['NAME_MAIN'])
        city_dir = base_dir / 'output' / city_name
        Path.mkdir(city_dir, parents=True, exist_ok=True)

        if VERBOSE: print(f'Starting {obs["NAME_MAIN"]}')

        # load observation area and do preprocessing (clipping and reclassification)
        obs_area = ObservationArea(
            polygon=obs.geometry,
            name=city_name,
            base_dir=city_dir
        )

        # preprocessing raster to split multi-year data
        obs_area.process_raster(ghs_raster, multi_year=ghsl_dict)

        # create footprint
        nyu_footprint = {}
        for year, data in obs_area.processed_rasters.items():
            nyu_footprint[year] = NyuFootprint(name=obs_area.name,
                                               raster_data=data,
                                               raster_meta=obs_area.raster_meta)

            # apply the standard NYU Urban Extent method
            nyu_footprint[year].standard_method()

            if nyu_footprint[year].footprint:
                Utils.save_geometry_to_vector(
                    nyu_footprint[year].footprint,
                    nyu_footprint[year].raster_meta['crs'],
                    city_dir / f'{city_name}_footprint_{year}.geojson'
                )

        # Save as .pkl for later analysis
        with open(city_dir / f'{city_name}_footprint.pkl', 'wb') as output:
            pickle.dump(
                nyu_footprint,
                output,
                protocol=pickle.HIGHEST_PROTOCOL
            )

    print('>>END<<')