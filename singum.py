import geopandas as gpd
import pandas as pd 
from load_data import load_lion_gdf
from scipy.spatial import Voronoi
from shapely.geometry import LineString, MultiLineString, Point, Polygon
from collections import defaultdict
import math
import itertools
from load_data import nyc_boundaries, nyc_raster, downsample_raster, manhattan_census
from rasterstats import zonal_stats
from fiona import listlayers
import numpy as np

#read data
gdf = load_lion_gdf()
# nyc_boundaries()
#manhattan_census()
#nyc_raster()
#downsample_raster(
    #input_path="nyc_pop_density.tif",
    #output_path="nyc_pop_density_downsampled.tif",
    #scale_factor=4  # 4x coarser in width & height → 16x fewer pixels
#)

# ----------------------------------- DRAWING SINGUMS --------------------------------------

#find coordinates of each node + whether it's start/end node 
def extract_coord(geom, from_end=True):
    if isinstance(geom, LineString):
        coords = list(geom.coords)
    elif isinstance(geom, MultiLineString):
        coords = list(geom.geoms[0].coords) if from_end else list(geom.geoms[-1].coords)
    else:
        return None
    return coords[0] if from_end else coords[-1]

#for each row segment, get beg. point of segment
gdf["From_Coord"] = gdf.apply(lambda row: extract_coord(row.geometry, True), axis=1)
#for each row segment, get end point of segment 
gdf["To_Coord"] = gdf.apply(lambda row: extract_coord(row.geometry, False), axis=1)

#map nodes to segments and coordinates
node_to_segments = defaultdict(list)
node_to_coord = {}

#when this for loop ends every node should be mapped to its start nodes and end nodes
for _, row in gdf.iterrows(): #each row in gdf is a road segment
    node_to_segments[row.NodeIDFrom].append((row.geometry, row.To_Coord)) #for each start node of segment add end and geo
    node_to_segments[row.NodeIDTo].append((row.geometry, row.From_Coord)) #for each end point add start point and geo
    node_to_coord[row.NodeIDFrom] = row.From_Coord #assign coordinate to node
    node_to_coord[row.NodeIDTo] = row.To_Coord 


#build array of coordinates and match index to NodeID
node_coords = np.array([coord for coord in node_to_coord.values()]) #need this to use voronoi package
node_ids = list(node_to_coord.keys()) #node_coords[i] now coresp. to node_ids[i]
id_to_index = {node_ids[i]: i for i in range(len(node_ids))} #each node id is now mapped to its respective index in node coords

#build Voronoi polygons, ignoring infinite regions 
unique_coords, unique_indices = np.unique(node_coords, axis=0, return_index=True)
unique_node_ids = [node_ids[i] for i in unique_indices]
vor = Voronoi(unique_coords)

poly_geoms = []
poly_ids = []
skipped = 0

for point_idx, region_idx in enumerate(vor.point_region):
    verts = vor.regions[region_idx]
    if -1 in verts or len(verts) == 0:
        skipped += 1
        continue
    try:
        polygon = Polygon([vor.vertices[i] for i in verts])
        if polygon.is_valid and not polygon.is_empty:
            poly_geoms.append(polygon)
            poly_ids.append(unique_node_ids[point_idx])
    except:
        skipped += 1

print("Skipped nodes:", skipped)
print("Total unique nodes:", len(unique_node_ids))
assert len(poly_geoms) + skipped == len(unique_node_ids)


#create GeoDataFrame
voronoi_gdf = gpd.GeoDataFrame({
    "NodeID": poly_ids,
    "geometry": poly_geoms
}, crs=gdf.crs)  # make sure CRS matches original
voronoi_gdf = voronoi_gdf.dissolve(by="NodeID", as_index=False)

manhattan_gdf = gpd.read_file("manhattan_boundaries.geojson")

filtered_voronoi = voronoi_gdf[voronoi_gdf.geometry.intersects(manhattan_gdf.iloc[0].geometry)]
filtered_voronoi.to_file("singum_voronoi.geojson", driver="GeoJSON")

print("Before clipping:", len(poly_geoms))
print("After clipping:", len(filtered_voronoi))
print("Unique NodeIDs after clipping:", filtered_voronoi['NodeID'].nunique())



'''
# ---------------------------------------ADDING POPULATION VALUES -------------------------------------------

SINGUM_PATH = "singum_nodes.geojson"
RASTER_PATH = "nyc_pop_density_downsampled.tif"
OUTPUT_PATH = "singum_with_population.geojson"
CHUNK_SIZE = 100  # number of polygons to process at once

singum_gdf = gpd.read_file(SINGUM_PATH)

#simplify singum geometry 
singum_gdf["geometry"] = singum_gdf["geometry"].simplify(0.0001, preserve_topology=True)
singum_gdf = singum_gdf[singum_gdf.is_valid & ~singum_gdf.is_empty]

#chunked zonal stats ===
def compute_population_in_chunks(gdf, raster_path, chunk_size=100):
    pop_values = []
    for i in range(0, len(gdf), chunk_size):
        print(f"processing chunk {i} to {i+chunk_size}")
        chunk = gdf.iloc[i:i+chunk_size]
        stats = zonal_stats(
            chunk,
            raster_path,
            stats=["sum"],
            geojson_out=True,
            all_touched=True
        )
        pop_chunk = [s["properties"]["sum"] if s["properties"]["sum"] is not None else 0 for s in stats]
        pop_values.extend(pop_chunk)
    return pop_values

#run zonal stats
pop_vals = compute_population_in_chunks(singum_gdf, RASTER_PATH, chunk_size=CHUNK_SIZE)
singum_gdf["population"] = pop_vals

singum_gdf.to_file(OUTPUT_PATH, driver="GeoJSON")

'''
#------------------------------------CLASSIFYING SINGUMS BASED ON TYPE ----------------------------------------
'''
#paths
gdb_path = "buildingtype/pluto.gdb"
pluto_layer = "MapPLUTO_25v1_1_clipped"
singum_path = "singum_nodes.geojson"
output_path = "singum_with_use_type.geojson"

#Read parcel polygons from GDB 
pluto_gdf = gpd.read_file(gdb_path, layer=pluto_layer)
pluto_gdf["BBL"] = pluto_gdf["BBL"].astype(str)

# === 2. Classify each lot as 'residential' or 'work' ===
def classify_use(row):
    res = row.get("ResArea", 0) or 0
    work = (
        (row.get("ComArea", 0) or 0) +
        (row.get("OfficeArea", 0) or 0) +
        (row.get("RetailArea", 0) or 0) +
        (row.get("GarageArea", 0) or 0) +
        (row.get("FactryArea", 0) or 0)
    )
    return "residential" if res >= work else "work"

pluto_gdf["use_type"] = pluto_gdf.apply(classify_use, axis=1)
pluto_gdf = pluto_gdf[~pluto_gdf["use_type"].isna()]

# === 3. Load Singum polygons ===
singum_gdf = gpd.read_file(singum_path)
singum_gdf = singum_gdf.to_crs(pluto_gdf.crs)

# === 4. Spatial join: Singums × intersecting parcels ===
intersections = gpd.overlay(singum_gdf, pluto_gdf, how="intersection")
dominant_use = intersections.groupby("NodeID")["use_type"].agg(lambda x: x.mode()[0])
singum_gdf = singum_gdf.merge(dominant_use, on="NodeID", how="left")

# === 6. Save the result ===

singum_gdf.to_file(output_path, driver="GeoJSON")
print(f"✅ Done. Result saved to {output_path}")

'''
