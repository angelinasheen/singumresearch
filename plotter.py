from shapely.geometry import Polygon, MultiPolygon
import geopandas as gpd
import matplotlib.pyplot as plt

gdf = gpd.read_file("manhattan_boundaries.geojson")

fig, ax = plt.subplots()

gdf.plot(ax=ax, facecolor=None, edgecolor=None, linewidth=1)

plt.title("visualize singums")
plt.show()