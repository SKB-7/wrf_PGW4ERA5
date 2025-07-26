#%%
import xarray as xr
import eccodes as ec
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset, num2date, date2num
from netCDF4 import Variable
import numpy as np

#%%
# era_file = xr.open_dataset('/mnt/hdd2/S_K_B/ERA5_model_data/ERA5-qtuv-ml-2022-01-05-2022-01-11.grib', decode_cf=False)

# with open('/mnt/hdd2/S_K_B/ERA5_model_data/ERA5-qtuv-ml-2022-01-05-2022-01-11.grib', 'rb') as f:        
#     gid = ec.codes_grib_new_from_file(f)
#     pv = ec.codes_get_array(gid, "pv")
#     n_half_levels = len(pv) // 2
#     a_half = pv[:n_half_levels]
#     b_half = pv[n_half_levels:]
#     # a_half = a_half[:-1]
#     # b_half = b_half[:-1]
    
#     # Clean up
#     ec.codes_release(gid)

# with Dataset("/mnt/hdd2/S_K_B/ERA5_model_data/ak.nc", "w") as nc_ak:
#     nc_ak.createDimension("level1", len(a_half))
#     ak_var = nc_ak.createVariable("ak", "f8", ("level1",))
#     ak_var[:] = a_half
#     ak_var.units = "Pa"
#     ak_var.long_name = "A half-level coefficients"

# with Dataset("/mnt/hdd2/S_K_B/ERA5_model_data/bk.nc", "w") as nc_bk:
#     nc_bk.createDimension("level1", len(b_half))
#     bk_var = nc_bk.createVariable("bk", "f8", ("level1",))
#     bk_var[:] = b_half
#     bk_var.units = "1"
#     bk_var.long_name = "B half-level coefficients"

# #%%
# ds_qtuv= xr.open_dataset('/mnt/hdd2/S_K_B/ERA5_model_data/ERA5-qtuv_ml-2022-01-05-2022-01-11.nc')
# ds_z= xr.open_dataset('/mnt/hdd2/S_K_B/ERA5_model_data/ERA5-allvar_ml-2022-01-05-2022-01-11.nc')
# ds_sfc = xr.open_dataset('/mnt/hdd4/S_K_B2/WRF/op/wrf_wps_09apr25/surface_level/ERA5-20220105-20220111-sl.nc')
#%%

qtuv_path = '/mnt/hdd2/S_K_B/ERA5_model_data/ERA5-qtuv_ml-2022-01-05-2022-01-11.nc'
z_path = '/mnt/hdd2/S_K_B/ERA5_model_data/ERA5-allvar_ml-2022-01-05-2022-01-11.nc'
sfc_path = '/mnt/hdd4/S_K_B2/WRF/op/wrf_wps_09apr25/surface_level/ERA5-20220105-20220111-sl.nc'
merged_path = '/mnt/hdd2/S_K_B/ERA5_model_data/ERA5-merged_ml-2022-01-05-2022-01-11.nc'
ak_path = '/mnt/hdd2/S_K_B/ERA5_model_data/ak.nc'
bk_path = '/mnt/hdd2/S_K_B/ERA5_model_data/bk.nc'

# Overwrite output file if it exists
if os.path.exists(merged_path):
    os.remove(merged_path)

dim_rename_map = {
    'time': 'valid_time',
    'lon': 'longitude',
    'lat': 'latitude'
}

np2ncdtype = {
    'float32': 'f4',
    'float64': 'f8',
    'int32': 'i4',
    'int16': 'i2',
    'int8': 'i1',
    'uint8': 'u1'
}

with Dataset(merged_path, 'w') as nc_out:
    with Dataset(qtuv_path, 'r') as nc_qtuv, Dataset(z_path, 'r') as nc_z, Dataset(sfc_path, 'r') as nc_sfc:

        # Copy dimensions from qtuv
        for name, dim in nc_qtuv.dimensions.items():
            nc_out.createDimension(name, len(dim) if not dim.isunlimited() else None)

        # Copy variables from qtuv
        for name, var in nc_qtuv.variables.items():
            var_dims = var.dimensions
            print(f"Processing variable: {name} with dimensions: {var_dims}, 'dtype': {var.dtype}")
            var_data = var[:]

            if 'model_level' in var_dims and var.shape[var_dims.index('model_level')] == 137:
                print('var shape:', var.shape)
                # Insert a level 0 with zeros
                # axis = var_dims.index('model_level')
                # new_shape = list(var.shape)
                # new_shape[axis] = 1  # shape for level 0

                # level0 = np.zeros(new_shape, dtype=var.dtype)
                # var_data = np.concatenate([level0, var_data], axis=axis)
                # # Update dimension size if needed
                # if 'model_level' not in nc_out.dimensions:
                #     nc_out.createDimension('model_level', 138)
            if name == "model_level":
                var_data = np.arange(var_data.shape[0], dtype='i4')
                nc_dtype = 'i4'
                out_var = nc_out.createVariable(name, nc_dtype, var_dims)
                out_var[:] = var_data
            else:
                dtype_str = str(var_data.dtype)
                nc_dtype = np2ncdtype.get(dtype_str, 'f4')
                out_var = nc_out.createVariable(name, nc_dtype, var_dims)
                out_var[:] = var_data

            # out_var = nc_out.createVariable(name, nc_dtype, var.dimensions)
            # Set attributes safely
            for attr in var.ncattrs():
                if attr == "_FillValue":
                    continue
                out_var.setncattr(attr, var.getncattr(attr))

            if "_FillValue" in var.ncattrs():
                try:
                    fill_val = var.getncattr("_FillValue")
                    if isinstance(fill_val, np.ndarray):
                        fill_val = fill_val.item()
                    out_var.setncattr("_FillValue", type(var_data.flat[0])(fill_val))
                except Exception as e:
                    print(f"Skipping _FillValue for {name} due to: {e}")
            # out_var[:] = var[:]
        
        if "valid_time" in nc_qtuv.variables:
            time_var_in = nc_qtuv.variables["valid_time"]
            time_data = time_var_in[:]
            time_units = time_var_in.units
            calendar = getattr(time_var_in, 'calendar', 'standard')

            # Convert and re-encode
            time_datetimes = num2date(time_data, units=time_units, calendar=calendar)
            encoded_time = date2num(time_datetimes, units="hours since 1970-01-01 00:00:00", calendar="standard")

            # Modify the existing variable
            time_out = nc_out.variables["valid_time"]
            time_out[:] = encoded_time
            time_out.units = "hours since 1970-01-01 00:00:00"
            time_out.calendar = "standard"
            time_out.long_name = "Forecast valid time"


        # Handle and flatten z variables
        for name, var in nc_z.variables.items():
            if name in nc_out.variables:
                continue  # skip duplicates

            var_dims = list(var.dimensions)
            var_data = var[:]

            # If there's a singleton model_level dimension, drop it
            if 'model_level' in var_dims and var.shape[var_dims.index('model_level')] == 1:
                # Remove the dimension
                var_data = np.squeeze(var_data, axis=var_dims.index('model_level'))
                var_dims.remove('model_level')

            out_var = nc_out.createVariable(name, var.datatype, tuple(var_dims))
            out_var.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
            out_var[:] = var_data

        # # Copy surface variables from sfc
        for name, var in nc_sfc.variables.items():
            if name in nc_out.variables:
                continue

            # Skip coordinate variables (e.g., time, lat, lon themselves)
            if len(var.dimensions) == 1 and var.dimensions[0] == name:
                continue

            # Rename dimensions
            new_dims = tuple(dim_rename_map.get(d, d) for d in var.dimensions)

            # Create any missing renamed dimensions
            for new_dim, orig_dim in zip(new_dims, var.dimensions):
                if new_dim not in nc_out.dimensions:
                    nc_out.createDimension(new_dim, len(nc_sfc.dimensions[orig_dim]))

            # Create variable
            out_var = nc_out.createVariable(name, var.datatype, new_dims)
            out_var.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
            out_var[:] = var[:]

        # Optional: copy global attributes
        nc_out.setncatts(nc_qtuv.__dict__)

# %%
# Open all files
with Dataset(merged_path, "a") as nc_merged, \
     Dataset(ak_path, "r") as nc_ak, \
     Dataset(bk_path, "r") as nc_bk:

    # Copy ak
    if "ak" not in nc_merged.variables:
        # Ensure dimension exists
        if "level1" not in nc_merged.dimensions:
            nc_merged.createDimension("level1", len(nc_ak.dimensions["level1"]))
        var = nc_ak.variables["ak"]
        out_var = nc_merged.createVariable("ak", var.datatype, var.dimensions)
        out_var.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
        out_var[:] = var[:]

    # Copy bk
    if "bk" not in nc_merged.variables:
        if "level1" not in nc_merged.dimensions:
            nc_merged.createDimension("level1", len(nc_bk.dimensions["level1"]))
        var = nc_bk.variables["bk"]
        out_var = nc_merged.createVariable("bk", var.datatype, var.dimensions)
        out_var.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
        out_var[:] = var[:]
# %%
