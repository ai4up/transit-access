# transit-access

Calculate accessibility to transit for any location using a gravity-based model with a Gaussian decay function.


## Install
```bash
pip install git+https://github.com/ai4up/transit-access@main
```


## Usage
**Prerequisite**

Use [gtfs2nx](https://github.com/ai4up/gtfs2nx/) package to convert GTFS feed to GeoDataFrame of transit stops.
```Python
G = gtfs2nx.transit_graph('some-dir/city-GTFS.zip', 'EPSG:26914')
stops = utils.nodes_to_gdf(G)
```
**API**
```Python
# For specific lat lon location
locations = geopandas.GeoSeries.from_xy([13.351798027529089], [52.49615200183667], crs='EPSG:4326').to_crs(stops.crs)
score.transit_access(stops, locations)

# For H3 Uber hexagonal grid
score.transit_access_for_grid(stops, h3_res=8)

# For spatial units, e.g. zip codes
gdf_neighborhoods = gpd.read_file(<some-file>)
score.transit_access_for_neighborhood(stops, gdf_neighborhoods)
```



## Approach
1. Create routable [NetworkX](https://github.com/networkx/networkx) graph with realistic transfer times from [GTFS](https://developers.google.com/transit/gtfs/) feeds using [gtfs2nx](https://github.com/ai4up/gtfs2nx/) package.
1. Transit score per stop
    * Calculation depends on notion of transit accessibility (e.g. frequency weighted by closeness centrality)
1. Calculate accessibility score per location
    * Gravity model following a Gaussian decay with distance to stops
    * Consider only nearest stop per route
1. Mapping accessibility to neighborhoods
    * Calculate mean over multiple locations arranged in hexagonal grid (Uber H3 with 10 resolution)
