import geopandas as gpd
import osmnx as ox

# for bounding box change it.
G = ox.graph_from_place("New Jersey, USA", network_type="drive")
print("done downloading")
ox.save_graphml(G, "new_jersey_drive.graphml")

# fsr i think data for speed limits don't work,
# so save both the graphml and gdfs
nodes, edges = ox.graph_to_gdfs(G)
nodes.to_file("nj_nodes.geojson", driver="GeoJSON")
edges.to_file("nj_edges.geojson", driver="GeoJSON")
print(nodes.columns)
print(edges.columns)