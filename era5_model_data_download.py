
''' 

https://cds.climate.copernicus.eu/datasets/reanalysis-era5-complete?tab=d_download
https://apps.ecmwf.int/codes/grib/param-db/ 

***https://apps.ecmwf.int/data-catalogues/era5/?class=ea

c = cdsapi.Client()
c.retrieve('reanalysis-era5-complete', { # Requests follow MARS syntax
                                        # Keywords 'expver' and 'class' can be dropped. They are obsolete
                                        # since their values are imposed by 'reanalysis-era5-complete'
    'date'    : '2013-01-01',            # The hyphens can be omitted
    'levelist': '1/10/100/137',          # 1 is top level, 137 the lowest model level in ERA5. Use '/' to separate values.
    'levtype' : 'ml',
    'param'   : '130',                   # Full information at https://apps.ecmwf.int/codes/grib/param-db/
                                        # The native representation for temperature is spherical harmonics
    'stream'  : 'oper',                  # Denotes ERA5. Ensemble members are selected by 'enda'
    'time'    : '00/to/23/by/6',         # You can drop :00:00 and use MARS short-hand notation, instead of '00/06/12/18'
    'type'    : 'an',
    'area'    : '80/-50/-25/0',          # North, West, South, East. Default: global
    'grid'    : '1.0/1.0',               # Latitude/longitude. Default: spherical harmonics or reduced Gaussian grid
    'format'  : 'netcdf',                # Output needs to be regular lat-lon, so only works in combination with 'grid'!
}, 'ERA5-ml-temperature-subarea.nc')     # Output file. Adapt as you wish.

'''

#!/usr/bin/env python
import cdsapi


# variable_sf_ml = {
#     "10m_u_component_of_wind": '165',
#     "10m_v_component_of_wind": '166',
#     "2m_dewpoint_temperature": '168',
#     "2m_temperature": '167',
#     "mean_sea_level_pressure": '151',
#     "sea_surface_temperature": '34',
#     "surface_pressure": '134',
#     "total_precipitation": '228',
#     "skin_temperature": '235',
#     "surface_latent_heat_flux": '137',
#     "top_net_solar_radiation_clear_sky": '176',
#     "snow_depth": '141',
#     "soil_temperature_level_1": '139',
#     "soil_temperature_level_2": '170',
#     "soil_temperature_level_3": '171',
#     "soil_temperature_level_4": '172',
#     "soil_type": '152',
#     "volumetric_soil_water_layer_1": '40',
#     "volumetric_soil_water_layer_2": '41',
#     "volumetric_soil_water_layer_3": '42',
#     "volumetric_soil_water_layer_4": '43',
#     "leaf_area_index_high_vegetation": '33',
#     "geopotential": '129',
#     "land_sea_mask": '172',
#     "sea_ice_cover": '31',
# }

variable_sf_ml = {
    "10m_u_component_of_wind": '165',
    "10m_v_component_of_wind": '166',
    "2m_dewpoint_temperature": '168',
    "2m_temperature": '167',
    "mean_sea_level_pressure": '151',
    "sea_surface_temperature": '34',
    "surface_pressure": '134',
    "total_precipitation": '228',
    "skin_temperature": '235',
    "surface_latent_heat_flux": '147',
    "top_net_solar_radiation_clear_sky": '208',
    "snow_depth": '141',
    "soil_temperature_level_1": '139',
    "soil_temperature_level_2": '170',
    "soil_temperature_level_3": '183',
    "soil_temperature_level_4": '236',
    "soil_type": '43',
    "volumetric_soil_water_layer_1": '39',
    "volumetric_soil_water_layer_2": '40',
    "volumetric_soil_water_layer_3": '41',
    "volumetric_soil_water_layer_4": '42',
    "leaf_area_index_high_vegetation": '67',
    #"geopotential": '129',
    "land_sea_mask": '172',
    "sea_ice_cover": '31',
    "Sea_ice_area_fraction": '262001',
}

Geopotential = '129'
log_sur_pressure = '152'
Specific_humidity = '133'
Temperature = '130'
U_component_of_wind = '131'
V_component_of_wind = '132'

# temperature = '130'
# u_wind = '131'
# v_wind = '132'
# rel_humidity = '157'
# geopotential = '129'
# temperature_2m = '167'
# rel_humidity_2m = '165'
# sea_surf_temp = '151'
# surf_skin_temp = '151'

variable_ml= (f'{Geopotential}/{Temperature}/{U_component_of_wind}/{V_component_of_wind}/{Specific_humidity}/{log_sur_pressure}')
#Relative humidity is not available on model levels in ERA5, but can be derived from specific humidity and temperature.
#Geopotential has to be downloaded separately, as it is not available on model levels in ERA5.
#grib_to_netcdf -o /mnt/hdd2/S_K_B/ERA5_model_data/ERA5-GS-ml-2022-01-05-2022-01-11.nc /mnt/hdd2/S_K_B/ERA5_model_data/ERA5-GS-ml-2022-01-05-2022-01-11.grib
#ds_zs = xr.open_dataset('/mnt/hdd2/S_K_B/ERA5_model_data/ERA5-GS-ml-2022-01-05-2022-01-11.nc')
#ds_zs.z.isel(level=[0])

v_sfc_l = '/'.join(variable_sf_ml.values())

all_var_ml = (f'{v_sfc_l}/{variable_ml}')

start_data = '2022-01-05'
end_data = '2022-01-11'
area = '42/67/22/88' # North, West, South, East. Default: global

model_folder = '/mnt/hdd4/S_K_B2/WRF/op/wrf_wps_29may25/model_levels/' #'/mnt/hdd2/S_K_B/ERA5_model_data/'
model_file = (f'{model_folder}ERA5-ztuvqp_ml-{start_data}-{end_data}.grib')

print(f"Downloading model data from {start_data} to {end_data} for area {area} and saving to {model_file}")

dataset = 'reanalysis-era5-complete'

request = {
    "class"   : "ea",
    'date'    : start_data + '/to/' + end_data,
    'time'    : '00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00',
   # '00/to/23/by/1',
    'area'    : area,
    'param'   : variable_ml,  # '129/130/131/132/133/152',
    'stream'  : 'oper',
    'type'    : 'an',
    'levelist': '/'.join(str(i) for i in range(1, 138)),
    'levtype' : 'ml',
    'grid'    : '0.25/0.25',
    # 'data_format'  : 'grib', # grib
    # 'download_format' : "unarchived",
}

c = cdsapi.Client()
c.retrieve(dataset, request, model_file)
print(f"Model data downloaded and saved to {model_file}")

#%%
# "31.128/34.128/39.128/40.128/41.128/42.128/43.128/67.128/134.128/139.128/141.128/151.128/"
# "165.128/166.128/167.128/168.128/170.128/172.128/183.128/235.128/236.128"

