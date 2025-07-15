import geopandas as gpd
# import rasterio
# from rasterio.mask import mask
# import rasterio
# from rasterio.enums import Resampling


def load_lion_gdf():
    path = "lion/lion.gdb"
    gdf = gpd.read_file(path, layer="lion")
    columns = ["geometry", "POSTED_SPEED", "TrafDir", "NodeIDFrom", "NodeIDTo", "SegmentID"]
    return gdf[columns]

def nyc_boundaries():
    path = "nyc_boundary.geojson"
    gdf = gpd.read_file(path)
    manhattan_gdf = gdf[gdf["district"].str.contains("MN", case=False)]
    manhattan_gdf.to_file("manhattan_boundary.geojson", driver="GeoJSON")

def manhattan_census():
    path = "nyc_2020_census_blocks.geojson"
    gdf = gpd.read_file(path)
    manhattan_gdf = gdf[gdf["boroname"].str.contains("Manhattan", case=False)]
    manhattan_gdf.to_file("manhattan_census_blocks.geojson", driver="GeoJSON")

def nyc_raster():
    nyc_boundary = gpd.read_file("nyc_boundary.geojson")  
    nyc_boundary = nyc_boundary.to_crs("EPSG:4326") 

    #clip raster
    with rasterio.open("usa_density.tif") as src:
        out_image, out_transform = mask(src, nyc_boundary.geometry, crop=True)
        out_meta = src.meta.copy()

    out_meta.update({
        "height": out_image.shape[1],
        "width": out_image.shape[2],
        "transform": out_transform
    })

    #save clipped nyc pop raster 
    with rasterio.open("nyc_pop_density.tif", "w", **out_meta) as dest:
        dest.write(out_image)

def downsample_raster(input_path, output_path, scale_factor=2):
    with rasterio.open(input_path) as src:
        # Calculate new dimensions
        new_width = src.width // scale_factor
        new_height = src.height // scale_factor

        # Read and resample data
        data = src.read(
            out_shape=(src.count, new_height, new_width),
            resampling=Resampling.average  # or .bilinear / .nearest
        )

        # Scale transform
        new_transform = src.transform * src.transform.scale(
            (src.width / new_width),
            (src.height / new_height)
        )

        # Write to new file
        profile = src.profile
        profile.update({
            "height": new_height,
            "width": new_width,
            "transform": new_transform
        })

        with rasterio.open(output_path, "w", **profile) as dst:
            dst.write(data)

