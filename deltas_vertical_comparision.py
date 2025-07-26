import xarray as xr
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

c=xr.open_dataset('/mnt/hdd2/S_K_B/ERA5_model_data/cas1/cas2022010500.nc')
c1=xr.open_dataset('/mnt/hdd2/S_K_B/pre_processed_files/step03/cas2022010500.nc')
p=xr.open_dataset('/mnt/hdd4/S_K_B2/WRF/op/wrf_wps_09apr25/pressure_levels/ERA5-20220105-20220111-pl.nc')
p=p.sel(time='2022-01-05')
pressure = c.model_level  # or pgw.level
c_mean = c.mean(dim=["latitude", "longitude"])
c1_mean = c1.mean(dim=["latitude", "longitude"])
p_mean = p.mean(dim=["lat", "lon"])
# plt.figure(figsize=(12, 5))

# plt.plot(c_mean.t[0], pressure, label='CTRL')
# plt.plot(c1_mean.t[0], pressure, label='PGW')
# plt.gca().invert_yaxis()
# plt.title("Domain-Averaged Temperature")
# plt.xlabel("Temperature (K)")
# plt.ylabel("Pressure (hPa)")
# plt.legend()

pa_era = xr.open_dataset('/mnt/hdd2/S_K_B/pre_processed_files/step03/pa_era_cas2022010500.nc')
pa_hl_era = xr.open_dataset('/mnt/hdd2/S_K_B/pre_processed_files/step03/pa_hl_era_cas2022010500.nc')

# Temperature (shape: level, lat, lon)
T_ctrl = c['t'].isel(valid_time=0)
T_pgw = c1['t'].isel(time=0)

# Define standard pressure levels in hPa
plevs = np.array([1000, 975, 950, 925,
900, 875, 850, 825,
800, 775, 750,
700, 650,
600, 550,
500, 450,
400, 350,
300, 250, 225,
200, 175, 150, 125,
100, 70, 50, 30, 20,
10, 7, 5, 3, 2, 1])


# Interpolation function: interpolates along "level"
def interpolate_column(p_model_1d, T_1d, p_out):
    f = interp1d(p_model_1d, T_1d, bounds_error=False, fill_value=np.nan)
    return f(p_out)

# Apply using xarray
def interp_to_plevs(T, p_model, plevs):
    return xr.apply_ufunc(
        interpolate_column,
        p_model, T,
        input_core_dims=[["model_level"], ["model_level"]],
        output_core_dims=[["plev"]],
        vectorize=True,
        dask="parallelized",
        kwargs={"p_out": plevs},
        output_dtypes=[T.dtype]
    )

# Apply interpolation
T_ctrl_interp = interp_to_plevs(T_ctrl, pa_era.pa_era[0]/100, plevs)
T_pgw_interp = interp_to_plevs(T_pgw, pa_era.pa_era[0]/100, plevs)
T_ctrl_interp = T_ctrl_interp.transpose("plev", "latitude", "longitude")
T_pgw_interp = T_pgw_interp.transpose("plev", "latitude", "longitude")

T_ctrl_interp = xr.DataArray(
    T_ctrl_interp,
    dims=("plev", "latitude", "longitude"),
    coords={"plev": plevs, "latitude": T_ctrl.latitude, "longitude": T_ctrl.longitude}
)

T_pgw_interp = xr.DataArray(
    T_pgw_interp,
    dims=("plev", "latitude", "longitude"),
    coords={"plev": plevs, "latitude": T_pgw.latitude, "longitude": T_pgw.longitude}
)

T_ctrl_mean = T_ctrl_interp.mean(dim=["latitude", "longitude"])
T_pgw_mean = T_pgw_interp.mean(dim=["latitude", "longitude"])

plt.figure(figsize=(8, 6))
plt.plot(T_ctrl_mean, plevs, label="CTRL")
plt.plot(T_pgw_mean, plevs, label="PGW")
plt.plot(p_mean.t[0][::-1], plevs, label="ERA5")
plt.gca().invert_yaxis()
plt.xlabel("Temperature (K)")
plt.ylabel("Pressure (hPa)")
plt.title("Domain-Averaged Temperature Profile (Standard P Levels)")
plt.legend()
plt.grid(True)
