import geopandas as gpd
import pandas as pd 
from load_data import load_lion_gdf
from shapely.geometry import LineString, MultiLineString, Point, Polygon
from collections import defaultdict
import math
import itertools
from load_data import nyc_boundaries, nyc_raster, downsample_raster, manhattan_census
from rasterstats import zonal_stats

#read data
gdf = load_lion_gdf()
nyc_boundaries()
#manhattan_census()
nyc_raster()
downsample_raster(
    input_path="nyc_pop_density.tif",
    output_path="nyc_pop_density_downsampled.tif",
    scale_factor=4  # 4x coarser in width & height → 16x fewer pixels
)

# ----------------------------------- DRAWING SINGUMS --------------------------------------

#find coordinates of nodes 
def extract_coord(geom, from_end=True):
    if isinstance(geom, LineString):
        coords = list(geom.coords)
    elif isinstance(geom, MultiLineString):
        coords = list(geom.geoms[0].coords)
    else:
        return None
    return coords[0] if from_end else coords[-1]

gdf["From_Coord"] = gdf.apply(lambda row: extract_coord(row.geometry, True), axis=1)
gdf["To_Coord"] = gdf.apply(lambda row: extract_coord(row.geometry, False), axis=1)

#map nodes to segments and coordinates
node_to_segments = defaultdict(list)
node_to_coord = {}

for _, row in gdf.iterrows():
    node_to_segments[row.NodeIDFrom].append((row.geometry, row.To_Coord))
    node_to_segments[row.NodeIDTo].append((row.geometry, row.From_Coord))
    node_to_coord[row.NodeIDFrom] = row.From_Coord
    node_to_coord[row.NodeIDTo] = row.To_Coord

#find perp bisector
def perp_bisector(p1, p2, length=50):
    mx, my = (p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2
    dx, dy = p2[0] - p1[0], p2[1] - p1[1]
    if dx == dy == 0:
        return None
    nx, ny = -dy, dx
    norm = math.hypot(nx, ny)
    nx, ny = nx / norm, ny / norm
    return LineString([
        (mx - nx * length, my - ny * length),
        (mx + nx * length, my + ny * length)
    ])

#calculate intersections to form polygons
singum_node_ids = []
singum_polygons = []

for node_id, segments in node_to_segments.items():
    center = node_to_coord[node_id]
    bisectors = []

    for geom, other_coord in segments:
        bis = perp_bisector(center, other_coord)
        if bis:
            bisectors.append(bis)
    skipped_nodes = []
    if len(bisectors) < 3:
        skipped_nodes.append({
            "NodeID" : node_id, 
            "geometry" : Point(center)})
        continue

    #find all intersections
    intersections = []
    for b1, b2 in itertools.combinations(bisectors, 2):
        pt = b1.intersection(b2)
        if isinstance(pt, Point) and pt.is_valid:
            intersections.append(pt)

    if len(intersections) >= 3:
        sorted_pts = sorted(intersections, key=lambda p: math.atan2(p.y - center[1], p.x - center[0]))
        poly = Polygon([(p.x, p.y) for p in sorted_pts])
        if poly.is_valid:
            singum_node_ids.append(node_id)
            singum_polygons.append(poly)
    

#build geodataframe
singum_gdf = gpd.GeoDataFrame(columns=["NodeID", "geometry"], geometry="geometry", crs=gdf.crs)
for nid, poly in zip(singum_node_ids, singum_polygons):
    singum_gdf.loc[len(singum_gdf)] = [nid, poly]
skipped_gdf = gpd.GeoDataFrame(skipped_nodes, crs=gdf.crs) 

#save
singum_gdf.to_file("singum_nodes.geojson", driver="GeoJSON")
skipped_gdf.to_file("skipped_nodes.geojson", driver = "GeoJSON")



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
