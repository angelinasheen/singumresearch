from shapely.geometry import Polygon, MultiPolygon
import geopandas as gpd
import matplotlib.pyplot as plt
import load_data

nyc_gdf = gpd.read_file("nyc.geojson")

fig, ax = plt.subplots()

nyc_gdf.plot(ax=ax, facecolor="Green", edgecolor=None, linewidth=1)


# just plot each singum i guess
singum_gdf = gpd.read_file("singum_voronoi.geojson")
singum_gdf.plot(ax=ax, facecolor=None, edgecolor="red", linewidth=.3, alpha=.5)

plt.title("visualize singums")
plt.show()