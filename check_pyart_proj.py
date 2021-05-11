"""
Sample pyart gridding call:
    grid = pyart.map.grid_from_radars(
        radar, grid_shape=(31, 501, 501),
        grid_limits=((0, 15000), (-250000,250000), (-250000, 250000)),
        fields=['reflectivity', 'differential_reflectivity', 'KDP_CSU', 'D0', 'NW', 'MU', 'MW', 'MI'],
        gridding_algo='map_gates_to_grid',
        h_factor=0., nb=0.6, bsp=1., min_radius=200.)
    return grid
We want to track down the coord transforms in pyart to make LMA match this.

Our use of map_gates_to_grid assumes the default projection
https://github.com/ARM-DOE/pyart/blob/master/pyart/map/gates_to_grid.py#L90
Which is aeqd in from _find_projparams
https://github.com/ARM-DOE/pyart/blob/master/pyart/map/gates_to_grid.py#L191
    # parse grid_projection
    if grid_projection is None:
        grid_projection = {
            'proj': 'pyart_aeqd', '_include_lon_0_lat_0': True}
    projparams = grid_projection.copy()
    if projparams.pop('_include_lon_0_lat_0', False):
        projparams['lon_0'] = grid_origin_lon
        projparams['lat_0'] = grid_origin_lat
    return projparams

Then it calls geographic_to_cartesian to find actual coords (gate_x, gate_y)

geographic_to_cartesian is
https://github.com/ARM-DOE/pyart/blob/ea9450a4dc58c070c7e51dd56106e05241276a47/pyart/core/transforms.py#L336
This uses built-in aeqd (geographic_to_cartesian_aeqd), not pyproj:
https://github.com/ARM-DOE/pyart/blob/ea9450a4dc58c070c7e51dd56106e05241276a47/pyart/core/transforms.py#L369

So we should be able to use aeqd in lmatools to get the same coordinates, with a spherical earth of 6370997.
"""


import numpy as np
import pyproj as proj4
from pyart.core.transforms import geographic_to_cartesian_aeqd

# Default earth radius from PyART's geographic_to_cartesian_aeqd
earth_R = 6370997.

khgx_lat, khgx_lon = 29.4719, -95.0792

lat = np.arange(28,31,.1)
lon = np.arange(-97,-93,.1)
lats, lons = np.meshgrid(lat, lon)

pyart_x, pyart_y = geographic_to_cartesian_aeqd(lons, lats,
							khgx_lon, khgx_lat, R=earth_R)



# Configure proj aeqd https://proj.org/operations/projections/aeqd.html
proj_map = proj4.crs.CRS(proj='aeqd', R=earth_R,
                     lat_0=khgx_lat, lon_0=khgx_lon)
proj_lla = proj4.crs.CRS(proj='latlong', R=earth_R)
trnsf_to_map = proj4.Transformer.from_crs(proj_lla, proj_map)
proj_x, proj_y = trnsf_to_map.transform(lons, lats)


dx = pyart_x-proj_x
dy = pyart_y-proj_y

print("Max x, y in pyart", pyart_x.max(), pyart_y.max())
print("Max x, y in proj", proj_x.max(), proj_y.max())

print("Should be effectively zero")
print(np.abs(dx).max())
print(np.abs(dy).max())
