import geopandas as gpd
from shapely.geometry import LineString, MultiLineString, Point, Polygon
import math
# summation of all the sides n_i n_j * len * K / A
# Q^I = K * (delta V) / len
# delta V = Q * len / speed limit
# CODE NOT TESTED AT ALL
def tensor(singum):
    singum = gpd.GeoSeries(singum)
    # requires geometry, flow through each link, flow through each side, center
    # center could be calculated here but its easier to just provide it.
    exterior = singum.geometry.exterior
    center = singum["center"]
    flows = singum["POSTED_SPEED"]
    res = [[0, 0], [0, 0]]
    A = exterior.area # this is probably not going to give a right result.
    for i in range(-1, len(exterior) - 1):
        a, b = exterior[i], exterior[i + 1]
        mid = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
        # the stupid thing is this is essentially re-finding the roads but whatever
        L = math.hypot(mid[0] - center[0], mid[1] - center[1])
        theta = math.atan2(mid[1] - center[1], mid[0] - center[0])
        K = flows[i]
        prod = L * K
        # get x comp
        res[0][0] += prod * math.cos(theta)
        res[1][1] += prod * math.sin(theta) # pls bro pls tell me im doing this right
        # [0][1] means as we go x how does y change and is that slope..
        res[0][1] += prod * math.tan(theta)
        res[1][0] += prod / math.tan(theta) # pls bro pls be right
        # hopes and dreams
    for i in range(2):
        for j in range(2):
            res[i][j] /= A
    return res