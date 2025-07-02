from shapely.geometry import Polygon, MultiPolygon
import geopandas as gpd
import matplotlib.pyplot as plt

manhattan_gdf = gpd.read_file("manhattan_boundaries.geojson")

fig, ax = plt.subplots()

manhattan_gdf.plot(ax=ax, facecolor=None, edgecolor=None, linewidth=1)

# just plot each singum i guess
singum_gdf = None
# singum_gdf.plot(ax=ax, facecolor="black", edgecolor="red", linewidth=1)

plt.title("visualize singums")
plt.show()