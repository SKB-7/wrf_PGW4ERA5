import cdsapi
import os
from concurrent.futures import ThreadPoolExecutor

# Global area: North/West/South/East
area = '90/-180/-90/180'
# start_data = '1979-01-01'
# end_data = '2024-12-31'

# Output folder
output_folder = '/mnt/hdd4/S_K_B2/era/mslp_global/'
os.makedirs(output_folder, exist_ok=True)

# Dataset and variable
dataset = 'derived-era5-single-levels-daily-statistics'
variable = 'mean_sea_level_pressure'

def download_mslp(year):
    """Download MSLP data for a given year (globally)."""
    c = cdsapi.Client()
    output_file = os.path.join(output_folder, f'ERA5_mslp_global_{year}.nc')
    print(f"📤 Submitting download for {year} → {output_file}")
    
    request = {
        'product_type': 'reanalysis',
        # 'format': 'grib',
        'variable': [variable],
        'year': str(year),
        'month': [f"{m:02d}" for m in range(1, 13)],
        'day': [f"{d:02d}" for d in range(1, 32)],
        # 'time': [f"{h:02d}:00" for h in range(24)],
        'area': area,
        'daily_statistic': 'daily_mean',
        'time_zone': 'utc+00:00',
        'frequency': '1_hourly'
    }

    try:
        c.retrieve(dataset, request, output_file)
        print(f"✅ Download complete: {output_file}")
    except Exception as e:
        print(f"❌ Failed for {year}: {e}")

# List of years to download
years = list(range(1979, 2025))
max_workers = 6  # You can adjust based on your system/network

# Asynchronous download
with ThreadPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(download_mslp, year) for year in years]

print("🚀 All global MSLP download jobs submitted.")

#     "class"   : "ea",
#     'date'    : start_data + '/to/' + end_data,
#     'time'    : '00/to/23/by/1',
# # '00/to/23/by/1',
#     'area'    : area,
#     'param'   : variable,  # '129/130/131/132/133/152',
#     'stream'  : 'oper',
#     'type'    : 'an',
#     'levtype' : 'sfc',
#     'grid'    : '0.25/0.25',
#     'format'  : 'netcdf', # grib