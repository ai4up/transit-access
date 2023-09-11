# transit-access

Calculate transit accessibility of any location using GTFS feeds and optionally [Overture POI data](https://docs.overturemaps.org/themes/places/place/).

Supported notions:
* Accessibility **to** transit
* Accessibility **to** transit system
* Accessibility to opportunities **by** transit


## Install
```bash
pip install git+https://github.com/ai4up/gtfs2nx@main
```


## Usage

**Accessibility to transit**
* Using service frequency (departures per hour)
```
transitaccess.access_to_transit(G)
```

**Accessibility to transit system**
* Using closeness centrality and average Euclidean travel speed
```
transitaccess.access_to_transit_system(G)
```

**Accessibility to opportunities by transit**
* Using [Overture POI data](https://docs.overturemaps.org/themes/places/place/) to calculate weighted average travel time
```
transitaccess.access_by_transit(G)
```


## Approach
1. Create routable [NetworkX](https://github.com/networkx/networkx) graph with realistic transfer times from [GTFS](https://developers.google.com/transit/gtfs/) feeds.
1. Transit score per stop
    * Calculation depends on notion of transit accessibility (e.g. frequency weighted by closeness centrality)
1. Calculate accessibility score per location
    * Gravity model following a Gaussian decay with distance to stops
    * Consider only nearest stop per route
1. Mapping accessibility to neighborhoods
    * Calculate mean over multiple locations arranged in hexagonal grid (Uber H3 with 10 resolution)
