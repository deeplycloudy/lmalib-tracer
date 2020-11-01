import glob
import numpy as np
import datetime
import xarray as xr
import pyproj as proj4

from lmalibtracer.coords.sbu import get_sbu_proj

from pyxlma.lmalib.io import read as lma_read

from pyxlma.lmalib.flash.cluster import cluster_flashes
from pyxlma.lmalib.flash.properties import flash_stats, filter_flashes
from pyxlma.lmalib.grid import  create_regular_grid, assign_regular_bins, events_to_grid

import sys
filenames = sys.argv[1:]

duration_min = 10
resolution_m = 2500
chi2max = 1.0
stationsmin = 6
latlon_grid=False # False uses the SBU stereographic coordinate grid.

print("Reading LMA")
lma_data, starttime = lma_read.dataset(filenames)

good_events = (lma_data.event_stations >= 6) & (lma_data.event_chi2 <= chi2max)
lma_data = lma_data[{'number_of_events':good_events}]


print("Clustering flashes")
ds = cluster_flashes(lma_data)
print("Calculating flash stats")
ds = flash_stats(ds)
ds = filter_flashes(ds, flash_event_count=(5, None))
# ds0 = ds.copy()

print("Setting up grid spec")

dttuple = [starttime, starttime+datetime.timedelta(minutes=duration_min)]
# dttuple = lma_data.Datetime.min(), lma_data.Datetime.max()
tstring = 'LMA {}-{}'.format(dttuple[0].strftime('%H%M'),
                                      dttuple[1].strftime('%H%M UTC %d %B %Y '))
print(tstring, dttuple)

grid_dt = np.asarray(60, dtype='m8[s]')
grid_t0 = np.asarray(dttuple[0]).astype('datetime64[ns]')
grid_t1 = np.asarray(dttuple[1]).astype('datetime64[ns]')
time_range = (grid_t0, grid_t1+grid_dt, grid_dt)

# Change the dictionaries below to a consistent set of coordinates
# and adjust grid_spatial_coords in the call to events_to_grid to
# change what is gridded (various time series of 1D, 2D, 3D grids)

if latlon_grid:
    # Houston
    # center = 29.7600000, -95.3700000
    lat_range = (27.75, 31.75, 0.025)
    lon_range = (-97.37, -93.37, 0.025)
    alt_range = (0, 18e3, 1.0e3)


    grid_edge_ranges ={
        'grid_latitude_edge':lat_range,
        'grid_longitude_edge':lon_range,
    #     'grid_altitude_edge':alt_range,
        'grid_time_edge':time_range,
    }
    grid_center_names ={
        'grid_latitude_edge':'grid_latitude',
        'grid_longitude_edge':'grid_longitude',
    #     'grid_altitude_edge':'grid_altitude',
        'grid_time_edge':'grid_time',
    }

    event_coord_names = {
        'event_latitude':'grid_latitude_edge',
        'event_longitude':'grid_longitude_edge',
    #     'event_altitude':'grid_altitude_edge',
        'event_time':'grid_time_edge',
    }

    flash_ctr_names = {
        'flash_init_latitude':'grid_latitude_edge',
        'flash_init_longitude':'grid_longitude_edge',
    #     'flash_init_altitude':'grid_altitude_edge',
        'flash_time_start':'grid_time_edge',
    }
    flash_init_names = {
        'flash_center_latitude':'grid_latitude_edge',
        'flash_center_longitude':'grid_longitude_edge',
    #     'flash_center_altitude':'grid_altitude_edge',
        'flash_time_start':'grid_time_edge',
    }
else:
    # Project lon, lat to SBU map projection
    sbu_lla, sbu_map, x_edge, y_edge = get_sbu_proj()
    sbu_dx = x_edge[1] - x_edge[0]
    sbu_dy = y_edge[1] - y_edge[0]
    lma_sbu_xratio = resolution_m/sbu_dx
    lma_sbu_yratio = resolution_m/sbu_dy
    trnsf_to_map = proj4.Transformer.from_crs(sbu_lla, sbu_map)
    trnsf_from_map = proj4.Transformer.from_crs(sbu_map, sbu_lla)
    lmax, lmay = trnsf_to_map.transform(#sbu_lla, sbu_map,
                                 ds.event_longitude.data,
                                 ds.event_latitude.data)
    lma_initx, lma_inity = trnsf_to_map.transform(#sbu_lla, sbu_map,
                                 ds.flash_init_longitude.data,
                                 ds.flash_init_latitude.data)
    lma_ctrx, lma_ctry = trnsf_to_map.transform(#sbu_lla, sbu_map,
                                 ds.flash_center_longitude.data,
                                 ds.flash_center_latitude.data)
    ds['event_x'] = xr.DataArray(lmax, dims='number_of_events')
    ds['event_y'] = xr.DataArray(lmay, dims='number_of_events')
    ds['flash_init_x'] = xr.DataArray(lma_initx, dims='number_of_flashes')
    ds['flash_init_y'] = xr.DataArray(lma_inity, dims='number_of_flashes')
    ds['flash_ctr_x'] = xr.DataArray(lma_ctrx, dims='number_of_flashes')
    ds['flash_ctr_y'] = xr.DataArray(lma_ctry, dims='number_of_flashes')

    grid_edge_ranges ={
        'grid_x_edge':(x_edge[0],x_edge[-1]+.001,sbu_dx*lma_sbu_xratio),
        'grid_y_edge':(y_edge[0],y_edge[-1]+.001,sbu_dy*lma_sbu_yratio),
    #     'grid_altitude_edge':alt_range,
        'grid_time_edge':time_range,
    }
    grid_center_names ={
        'grid_x_edge':'grid_x',
        'grid_y_edge':'grid_y',
    #     'grid_altitude_edge':'grid_altitude',
        'grid_time_edge':'grid_time',
    }

    event_coord_names = {
        'event_x':'grid_x_edge',
        'event_y':'grid_y_edge',
    #     'event_altitude':'grid_altitude_edge',
        'event_time':'grid_time_edge',
    }

    flash_ctr_names = {
        'flash_init_x':'grid_x_edge',
        'flash_init_y':'grid_y_edge',
    #     'flash_init_altitude':'grid_altitude_edge',
        'flash_time_start':'grid_time_edge',
    }
    flash_init_names = {
        'flash_ctr_x':'grid_x_edge',
        'flash_ctr_y':'grid_y_edge',
    #     'flash_center_altitude':'grid_altitude_edge',
        'flash_time_start':'grid_time_edge',
    }


print("Creating regular grid")
grid_ds = create_regular_grid(grid_edge_ranges, grid_center_names)
if latlon_grid:
    pass
else:
    ctrx, ctry = np.meshgrid(grid_ds.grid_x, grid_ds.grid_y)
    hlon, hlat = trnsf_from_map.transform(ctrx, ctry)
    if False:
        houston_grid=xr.open_dataset('/data/Houston/realtime-tracer/lmalib-tracer/test/sbu-radar-grids/houston_grid.nc')
        # Confirm x,y as in the grid match lon lat
        try:
            assert np.allclose(hlon,houston_grid.lon.data)
            assert np.allclose(hlat,houston_grid.lat.data)
            print("Grid matches SBU grid")
        except AssertionError:
            max_delta = np.abs(hlon-houston_grid.lon.data).max()
            print(max_delta)
            max_delta = np.abs(hlat-houston_grid.lat.data).max()
            print(max_delta)
    else:
        print("Skipping check of grid match to SBU grid")
    # Add lon lat to the dataset, too.
    ds['lon'] = xr.DataArray(hlon, dims=['grid_y', 'grid_x'],
                    attrs={'standard_name':'longitude'})
    ds['lat'] = xr.DataArray(hlat, dims=['grid_y', 'grid_x'],
                    attrs={'standard_name':'latitude'})

print("Finding grid position for flashes")
pixel_id_var = 'event_pixel_id'
ds_ev = assign_regular_bins(grid_ds, ds, event_coord_names,
    pixel_id_var=pixel_id_var, append_indices=True)
# ds_flctr = assign_regular_bins(grid_ds, ds, flash_ctr_names,
#     pixel_id_var='flash_ctr_pixel_id', append_indices=True)
# flctr_gb = ds.groupby('flash_ctr_pixel_id')
# ds_flini = assign_regular_bins(grid_ds, ds, flash_init_names,
#     pixel_id_var='flash_init_pixel_id', append_indices=True)
# flini_gb = ds.groupby('flash_init_pixel_id')

# print('===== ev_gb')
# for event_pixel_id, dsegb in ev_gb:
#     print(dsegb)
#     break
# print('===== flctr_gb')
# for event_pixel_id, dsfgb in flctr_gb:
#     print(dsfgb)
#     break

print("Gridding data")
if latlon_grid:
    grid_spatial_coords=['grid_time', None, 'grid_latitude', 'grid_longitude']
    event_spatial_vars = ('event_altitude', 'event_latitude', 'event_longitude')
else:
    grid_spatial_coords=['grid_time', None, 'grid_y', 'grid_x']
    event_spatial_vars = ('event_altitude', 'event_y', 'event_x')

# print(ds_ev)
# print(grid_ds)
grid_ds = events_to_grid(ds_ev, grid_ds, min_points_per_flash=3,
                         pixel_id_var=pixel_id_var,
                         event_spatial_vars=event_spatial_vars,
                         grid_spatial_coords=grid_spatial_coords)


both_ds = xr.combine_by_coords((grid_ds, ds))

print("Writing data")
duration_sec = (dttuple[1]-dttuple[0]).total_seconds()
if latlon_grid:
    date_fmt = "LYLOUT_%y%m%d_%H%M%S_{0:04d}.nc".format(int(duration_sec))
else:
    date_fmt = "LYLOUT_%y%m%d_%H%M%S_{0:04d}_map{1:d}m.nc".format(
                    int(duration_sec), resolution_m)
outfile = dttuple[0].strftime(date_fmt)

comp = dict(zlib=True, complevel=5)
encoding = {var: comp for var in both_ds.data_vars}
both_ds.to_netcdf(outfile, encoding=encoding)