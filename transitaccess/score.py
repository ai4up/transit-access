import sys
import logging

import numpy as np
import pandas as pd
import geopandas as gpd
import h3pandas

ID_SEP = '@@'

logging.basicConfig(format='%(asctime)s | %(levelname)s | %(message)s', level=logging.INFO, stream=sys.stdout)
logger = logging.getLogger(__name__)


def transit_access(stops, loc, decay_param=0.5):
    """
    Calculate the TransitAccess score for a set of locations.

    TransitAccess is a measure of the accessibility of a location using the public transportation network.
    It is calculated as the sum of the closeness centrality of all reachable transit stops weighted by the service frequency of the respective transit line.

    Parameters
    ----------
    stops : geopandas.GeoDataFrame
        GeoDataFrame of transit stops with precalculated closeness centrality and service frequency.
    loc : geopandas.GeoSeries
        GeoSeries of locations for which the TransitAccess score should be calculated.
    decay_param : float
        Decay parameter for the gaussian decay function used to weight the access to each stop by distance.


    Returns
    -------
    score : pandas.Series
        TransitAccess score for each location.


    Examples
    --------
    >>> G = network.transit_graph('some-dir/city-GTFS.zip', 'EPSG:26914')
    >>> stops = utils.nodes_to_gdf(G)
    >>> locations = geopandas.GeoSeries.from_xy([13.351798027529089], [52.49615200183667], crs='EPSG:4326')
    >>> locations = locations.to_crs(G.graph['crs'])
    >>> score.transit_access(stops, locations)
    0    0.00123
    """

    dm = _distance_matrix(loc, stops.geometry)
    score = _calculate_gravity_score(
        distance_matrix=dm,
        supply=stops['frequency'].values,
        decay_param=decay_param,
    )

    return score


def transit_access_for_grid(stops, area=None, h3_res=9):
    """
    Calculate the TransitAccess score for a hexagonal grid.

    TransitAccess is a measure of the accessibility of a location using the public transportation network.
    It is calculated as the sum of the closeness centrality of all reachable transit stops weighted by the service frequency of the respective transit line.

    Parameters
    ----------
    stops : geopandas.GeoDataFrame
        GeoDataFrame of transit stops with precalculated closeness centrality and service frequency.
    area : geopandas.GeoDataFrame
        GeoDataFrame of the area for which the TransitAccess score should be calculated.
        If None, the convex hull of the transit stops with a 2km buffer is used.
    h3_res : int
        Resolution of the hexagonal grid (see https://h3geo.org/docs/core-library/restable/).

    Returns
    -------
    score : geopandas.GeoDataFrame
        GeoDataFrame of the hexagonal grid with TransitAccess score for each hexagon.


    Examples
    --------
    >>> G = network.transit_graph('some-dir/city-GTFS.zip', 'EPSG:26914')
    >>> stops = utils.nodes_to_gdf(G)
    >>> score.transit_access_for_grid(stops, h3_res=8)
                               geometry                                                access_score
        h3_08
        8866e09107fffff	POLYGON ((987746.618 999962.825, 987300.330 99...	7.883864e-03
        8866e42f41fffff	POLYGON ((1007073.842 1015684.465, 1006628.554...	1.071192e-03
        8866e0911dfffff	POLYGON ((989646.084 1000076.113, 989199.868 9...	5.542570e-02
    """

    if area is None:
        logger.info('Calculating TransitAccess index for convex hull with 2km buffer around transit stops.')
        area = stops.dissolve().convex_hull.buffer(2000).to_frame('geometry')

    hex_grid = _create_hex_grid(h3_res, area)
    hex_grid['score_spatiotemporal'] = transit_access(stops, hex_grid.centroid)
    return hex_grid


def transit_access_for_neighborhood(stops, neighborhoods):
    """
    Calculate the TransitAccess score for a set of neighborhoods.

    TransitAccess is a measure of the accessibility of a location using the public transportation network.
    It is calculated as the sum of the closeness centrality of all reachable transit stops weighted by the service frequency of the respective transit line.
    The neighborhood score is calculated as the average of the TransitAccess scores of all h3 hexagon centroids within the neighborhood.

    Parameters
    ----------
    stops : geopandas.GeoDataFrame
        GeoDataFrame of transit stops with precalculated closeness centrality and service frequency.
    neighborhoods : geopandas.GeoDataFrame
        GeoDataFrame of the neighborhoods for which the TransitAccess score should be calculated.

    Returns
    -------
    score : geopandas.GeoDataFrame
        GeoDataFrame of the hexagonal grid with TransitAccess score for each hexagon.


    Examples
    --------
    >>> G = network.transit_graph('some-dir/city-GTFS.zip', 'EPSG:26914')
    >>> stops = utils.nodes_to_gdf(G)
    >>> score.transit_access_for_neighborhood(stops, gdf_zip_codes)
                 zip_code     geometry                                         	access_score
        0	10115	POLYGON ((389163.210 5821872.935, 389321.827 5...	0.379396
        1	10117	POLYGON ((389678.019 5820987.307, 389683.298 5...	0.450508
        2	10119	POLYGON ((391390.987 5820861.120, 391546.214 5...	0.364538
    """

    hex_grid = transit_access_for_grid(stops, neighborhoods)
    hex_grid['geometry'] = hex_grid.centroid
    neighborhoods = _mean_per_area(neighborhoods, hex_grid, 'score_spatiotemporal')
    return neighborhoods


def _distance_matrix(loc, stops_loc):
    if loc.crs != stops_loc.crs:
        raise Exception(f'Coordinate reference system (CRS) of provided locations ({loc.crs}) and transit stops ({stops_loc.crs}) is not the same.')

    if not loc.crs.is_projected:
        raise Exception(f'Please use projected coordinate reference system (CRS) to ensure accurate distance calculation.')

    dm = loc.apply(lambda g: stops_loc.distance(g)) / 1000
    dm.to_pickle(f'distance-matrix-{len(loc)}-{len(stops_loc)}.pkl')

    return dm


def _create_hex_grid(h3_res, area):
    ori_crs = area.crs
    hexbin = area.dissolve().to_crs('EPSG:4326')[['geometry']].h3.polyfill_resample(h3_res)
    hexbin = hexbin.rename_axis(f'h3_{h3_res:02d}').drop(columns=['index'])
    hexbin = hexbin.to_crs(ori_crs)
    return hexbin


def _mean_per_area(area, points, col):
    points = points[[col, 'geometry']].to_crs(area.crs)
    matched = gpd.sjoin(points, area, how='left', predicate='within')
    mean = matched.groupby('index_right')[col].mean().reset_index().set_index('index_right')
    area = pd.merge(area, mean, how='left', left_index=True, right_index=True)
    return area


def _calculate_gravity_score(distance_matrix, supply, decay_param):
    # weight access to each stop by distance using a gaussian decay
    access_to_each_stop = (supply * _gaussian_decay(distance_matrix, decay_param)).T

    # only consider a single stop per route and direction
    # intuition: closeness to two stops of same transit line is not better than to one stop with identical distance
    access_to_each_stop['route_id_w_direction'] = _route_id_w_direction(access_to_each_stop)
    access_to_each_route = access_to_each_stop.groupby('route_id_w_direction').transform(_sum_two_largest) # assuming direction of route is not encoded in index, but the stop very likely serves two directions
    # access_to_each_route = access_to_each_stop.groupby('route_id_w_direction').max() # assuming direction of route is encoded in index

    # summing up the access to all reachable stops
    access = access_to_each_route.sum()
    return access


def _route_id_w_direction(stops):
    return stops.index.str.split(ID_SEP, n=1).str[1]


def _gaussian_decay(distance_array, sigma):
    return np.exp(-(distance_array**2 / (2.0 * sigma**2)))


def _sum_n_largest(x, n):
    try:
        return np.partition(x, -n)[-n:].sum()
    except ValueError:
        return x.max()


def _sum_two_largest(x):
    return _sum_n_largest(x, n=2)
