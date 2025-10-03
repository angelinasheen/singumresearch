import numpy as np
import math
import geopandas as gpd
from shapely.validation import make_valid
from scipy.interpolate import griddata
from collections import defaultdict
from scipy.sparse import csr_matrix
from shapely.geometry import LineString
from shapely import union_all
import matplotlib.pyplot as plt
from simulations import apply_dirichlet_row_replace, solve_with_component_pinning, greedy_path_on_nodes, dijkstra_path, path_time, nan_gaussian_blur


TARGET_CRS = "EPSG:2263"
OUT_GEOJSON    = "manhattan_routes.geojson"
NODES_PHI = "nodes_with_phi_2263.geojson"
ROADS     = "filtered_roads_manhattan_2263.geojson"
ROUTES    = "manhattan_routes.geojson" 
OUT_PNG   = "manhattan_routes_viz.png"
VORO_IN   = "singum_voronoi.geojson"            # your existing file (must have 'NodeID')
NODES_GDF = "nodes_points_2263.geojson"         # from Step C (columns: NodeID, index, geometry)
GRAPH_TXT = "graph_edges.txt"                   # from Step C
MAN_GDF   = "manhattan_boundary_2263.geojson"   # Manhattan polygon(s)

VORO_OUT      = "singum_voronoi_with_K.geojson"
K_CENTROIDS_CSV = "K_centroids_2263.csv"

MPH_TO_FTPS = 1.46667   # mph -> feet/second
SEC_TO_MIN  = 1.0 / 60.0
CONST_MPH   = 25.0  
START_NODEID = 21465   # your start NodeID
END_NODEID   = 78136   # your end NodeID   



roads = gpd.read_file("filtered_roads.geojson")
manhattan   = gpd.read_file(MAN_GDF)

#make sure roads and manhattan boundary are on the same coord system 
if manhattan.crs is None and roads.crs is not None:
    manhattan = manhattan.set_crs(roads.crs)

if roads.crs.to_string() != TARGET_CRS:
    roads = roads.to_crs(TARGET_CRS)

if manhattan.crs is None or manhattan.crs.to_string() != TARGET_CRS:
    manhattan = manhattan.to_crs(TARGET_CRS)

roads.to_file("filtered_roads_2263.geojson", driver="GeoJSON")
manhattan.to_file("manhattan_boundary_2263.geojson", driver="GeoJSON")


#clip roads to manhattan boundary 
roads = gpd.read_file("filtered_roads_2263.geojson")
manhattan   = gpd.read_file("manhattan_boundary_2263.geojson")

manhattan_union = make_valid(union_all(list(manhattan.geometry)))
before_total = len(roads)
roads = gpd.clip(roads, gpd.GeoDataFrame(geometry=[manhattan_union], crs=manhattan.crs))
roads = roads[~roads.geometry.is_empty].copy()
roads = roads.explode(index_parts=False, ignore_index=True)  # ensures LineString per row
roads = roads[roads.geom_type == "LineString"].copy()
roads.to_file("filtered_roads_manhattan_2263.geojson", driver="GeoJSON")

def load_graph(graph_txt):
    with open(graph_txt, "r") as f:
        V, E = map(int, f.readline().split())
        list_nodes = list(map(int, f.readline().split()))
        efrom, eto, etime = [], [], []
        for line in f:
            u, v, w = line.split()
            efrom.append(int(u)); eto.append(int(v)); etime.append(float(w))
    return V, np.array(list_nodes, int), np.array(efrom, int), np.array(eto, int), np.array(etime, float)



def build_node():
    #constants
    ft_per_sec = CONST_MPH * MPH_TO_FTPS
    sec_per_ft = 1.0 / ft_per_sec
    min_per_ft = sec_per_ft * SEC_TO_MIN

    #build node indexes 
    rep_node = {}     # {NodeID -> idx}
    list_nodes = []   # idx -> NodeID
    coord_map = {}    # idx -> (x,y)
    V = 0
    edges = []        # (u_idx, v_idx, time_min)

    for _, r in roads.iterrows():
        u_id = int(r["NodeIDFrom"])
        v_id = int(r["NodeIDTo"])

        #endpoints of road
        p0 = r.geometry.coords[0]
        p1 = r.geometry.coords[-1]

        #register nodes
        for nid, pt in ((u_id, p0), (v_id, p1)):
            if nid not in rep_node:
                rep_node[nid] = V
                list_nodes.append(nid)
                coord_map[V] = np.array(pt, dtype=float)
                V += 1

        u = rep_node[u_id]; v = rep_node[v_id]

        #length in feet
        length_ft = float(r.geometry.length)
        time_min  = length_ft * min_per_ft

        # Directionality: W=with (u->v), A=against (v->u), T=two-way
        traf = (r.get("TrafDir") or "T").strip().upper()
        if traf in ("W", "T"):
            edges.append((u, v, time_min))
        if traf in ("A", "T"):
            edges.append((v, u, time_min))

    # Write graph file (V E, then NodeIDs, then u v w per edge)
    with open(GRAPH_TXT, "w") as f:
        E = len(edges)
        f.write(f"{V} {E}\n")
        f.write(" ".join(map(str, list_nodes)) + "\n")
        for u, v, w in edges:
            f.write(f"{u} {v} {w}\n")

    #save node points 
    node_xy = np.vstack([coord_map[i] for i in range(V)])
    nodes_gdf = gpd.GeoDataFrame(
        {"NodeID": list_nodes, "index": list(range(V))},
        geometry=gpd.points_from_xy(node_xy[:,0], node_xy[:,1]),
        crs=roads.crs
    )
    nodes_gdf.to_file(NODES_GDF, driver="GeoJSON")


def compute_route():
    nodes = gpd.read_file(NODES_PHI).sort_values("index")
    coords = np.c_[nodes.geometry.x.values, nodes.geometry.y.values]  # (N,2)
    phi    = nodes["phi"].to_numpy(dtype=float)

    V, list_nodes, efrom, eto, etime = load_graph(GRAPH_TXT)
    assert V == coords.shape[0] == phi.shape[0], "Node count mismatch across inputs."

    node_index = {nid: i for i, nid in enumerate(list_nodes)}
    src = node_index[START_NODEID]
    dst = node_index[END_NODEID]

    sing_nodes = greedy_path_on_nodes(phi, efrom, eto, src, dst, coords=coords)
    dijk_nodes, dijk_time = dijkstra_path(coords, efrom, eto, etime, src, dst)
    time_map = {(u, v): t for u, v, t in zip(efrom, eto, etime)}
    sing_time = path_time(sing_nodes, time_map)


    def path_to_linestring(path):
        return LineString([coords[i] for i in path]) if path and len(path) >= 2 else None

    sing_line = path_to_linestring(sing_nodes)
    dijk_line = path_to_linestring(dijk_nodes)

    out = gpd.GeoDataFrame(
        {
            "name":      ["singum", "dijkstra"],
            "time_min":  [float(sing_time), float(dijk_time)],
            "start_id":  [int(START_NODEID)]*2,
            "end_id":    [int(END_NODEID)]*2,
        },
        geometry=[sing_line, dijk_line],
        crs=nodes.crs,
    )
    out = out[out.geometry.notnull()]  # drop None if any
    out.to_file(OUT_GEOJSON, driver="GeoJSON")

def calculate_voltages():
    nodes = gpd.read_file(NODES_GDF).sort_values("index")
    coords = np.c_[nodes.geometry.x.values, nodes.geometry.y.values]  

    with open("graph_edges.txt", "r") as f:
        V, E = map(int, f.readline().split())
        list_nodes = list(map(int, f.readline().split()))  # external NodeIDs, length V
        efrom, eto, etime = [], [], []
        for line in f:
            u, v, w = line.split()
            efrom.append(int(u)); eto.append(int(v)); etime.append(float(w))
    efrom = np.array(efrom, int); eto = np.array(eto, int); etime = np.array(etime, float)
    assert V == coords.shape[0], "Node count mismatch between nodes file and graph file."

    Gpair = defaultdict(float)  
    for u, v, t in zip(efrom, eto, etime):
        if t <= 0: continue
        a, b = (u, v) if u < v else (v, u)
        Gpair[(a, b)] += 1.0 / max(t, 1e-12)

    rows, cols, data = [], [], []
    for (a, b), G in Gpair.items():
        rows += [a, b, a, b]
        cols += [a, b, b, a]
        data += [ G,  G, -G, -G]
    L = csr_matrix((np.array(data, float), (np.array(rows, int), np.array(cols, int))), shape=(V, V))

    node_index = {nid: i for i, nid in enumerate(list_nodes)}
    start_idx = node_index[START_NODEID]
    end_idx   = node_index[END_NODEID]

    ML, b = apply_dirichlet_row_replace(L, [start_idx], 1.0)
    ML, b = apply_dirichlet_row_replace(ML, [end_idx],   0.0)
    phi = solve_with_component_pinning(ML, b, pinned_nodes=np.r_[start_idx, end_idx])

    nodes["phi"] = phi
    nodes.to_file(NODES_PHI, driver="GeoJSON")

def graph():
    # ----- load nodes with phi -----
    nodes = gpd.read_file(NODES_PHI).sort_values("index")
    coords = np.c_[nodes.geometry.x.values, nodes.geometry.y.values]
    phi    = nodes["phi"].to_numpy(float)

    # ----- load graph (inside with) -----
    with open(GRAPH_TXT, "r") as f:
        V, E = map(int, f.readline().split())
        list_nodes = list(map(int, f.readline().split()))
        efrom, eto, etime = [], [], []
        for line in f:
            u, v, w = line.split()
            efrom.append(int(u)); eto.append(int(v)); etime.append(float(w))
    efrom = np.array(efrom, int); eto = np.array(eto, int); etime = np.array(etime, float)

    assert V == coords.shape[0] == phi.shape[0], "Node count mismatch."

    # map NodeID -> compact index
    node_index = {nid: i for i, nid in enumerate(list_nodes)}
    src = node_index[START_NODEID]
    dst = node_index[END_NODEID]

    # routes + times
    time_map   = {(u, v): t for u, v, t in zip(efrom, eto, etime)}
    sing_nodes = greedy_path_on_nodes(phi, efrom, eto, src, dst, coords=coords)
    dijk_nodes, dijk_time = dijkstra_path(coords, efrom, eto, etime, src, dst)
    sing_time  = path_time(sing_nodes, time_map)

    # ----- plot basemap + routes -----
    fig, ax = plt.subplots(figsize=(8.2, 10.5))
    try:
        man = gpd.read_file(MAN_GDF)
        man_union = make_valid(union_all(list(man.geometry)))
        gpd.GeoDataFrame(geometry=[man_union], crs=man.crs).boundary.plot(ax=ax, color="black", linewidth=1.0)
    except Exception:
        man_union = None
    try:
        roads = gpd.read_file(ROADS)
        roads.plot(ax=ax, color="lightgray", linewidth=0.35)
    except Exception:
        pass

    # draw paths
    if sing_nodes and len(sing_nodes) >= 2:
        C = coords[np.array(sing_nodes)]
        ax.plot(C[:,0], C[:,1], '-', lw=3, color='cyan', label=f"Singum greedy ({sing_time:.1f} min)")
    if dijk_nodes and len(dijk_nodes) >= 2:
        C = coords[np.array(dijk_nodes)]
        ax.plot(C[:,0], C[:,1], '-', lw=3, color='magenta', label=f"Dijkstra ({dijk_time:.1f} min)")

    # start/end dots
    ax.plot(coords[src,0], coords[src,1], 'go', ms=8, label="Start")
    ax.plot(coords[dst,0], coords[dst,1], 'ro', ms=8, label="End")

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel("x (ft, EPSG:2263)"); ax.set_ylabel("y (ft, EPSG:2263)")
    ax.set_title("Manhattan: Singum (greedy) vs Dijkstra")
    ax.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=220)
    print(f"Singum time:   {sing_time:.3f} min")
    print(f"Dijkstra time: {dijk_time:.3f} min")
    print("Saved", OUT_PNG)

    # ====== Voronoi + K enrichment (once) ======
    vor = gpd.read_file(VORO_IN)[["NodeID","geometry"]]
    man = gpd.read_file(MAN_GDF)
    man_union = make_valid(union_all(list(man.geometry)))

    # clip to Manhattan and make one poly per NodeID
    vor = gpd.overlay(vor, gpd.GeoDataFrame(geometry=[man_union], crs=man.crs),
                      how="intersection", keep_geom_type=True)
    vor = vor.dissolve(by="NodeID", as_index=False)
    vor["area_ft2"] = vor.geometry.area

    # undirected conductance G = sum(1/t)
    from collections import defaultdict
    Gpair = defaultdict(float)
    for u, v, t in zip(efrom, eto, etime):
        if t > 0:
            a, b = (u, v) if u < v else (v, u)
            Gpair[(a, b)] += 1.0 / t

    # per-node tensors (pre-area)
    Kxx = np.zeros(V); Kxy = np.zeros(V); Kyy = np.zeros(V)
    for (a, b), G in Gpair.items():
        pa, pb = coords[a], coords[b]
        d = pb - pa
        L = np.linalg.norm(d)
        if L <= 0 or G <= 0:
            continue
        n = d / L
        factor = 0.5 * (L * L) * G   # (L^2/2) * G
        nx, ny = float(n[0]), float(n[1])
        for u in (a, b):
            Kxx[u] += factor * nx * nx
            Kxy[u] += factor * nx * ny
            Kyy[u] += factor * ny * ny

    # normalize by Voronoi area aligned by NodeID
    nodeid_to_idx = {nid: i for i, nid in enumerate(list_nodes)}
    idx_to_nodeid = np.array(list_nodes, int)

    area = np.zeros(V)
    vor = vor[vor["NodeID"].isin(nodeid_to_idx)].copy()
    area_idx = vor["NodeID"].map(nodeid_to_idx).to_numpy(int)
    area[area_idx] = vor["area_ft2"].to_numpy(float)
    safe_area = np.where(area > 0, area, np.nan)

    Kxx /= safe_area; Kxy /= safe_area; Kyy /= safe_area

    # attach to polygons + export
    Ktab = gpd.GeoDataFrame({"NodeID": idx_to_nodeid, "Kxx": Kxx, "Kxy": Kxy, "Kyy": Kyy}, geometry=None)
    vor = vor.merge(Ktab, on="NodeID", how="left")

    vor.to_file(VORO_OUT, driver="GeoJSON")
    cent = vor.geometry.centroid
    cent_xy = np.c_[cent.x.values, cent.y.values]
    np.savetxt(K_CENTROIDS_CSV,
               np.c_[cent_xy[:,0], cent_xy[:,1],
                     vor["Kxx"].to_numpy(float),
                     vor["Kxy"].to_numpy(float),
                     vor["Kyy"].to_numpy(float)],
               delimiter=",",
               header="x_ft,y_ft,Kxx,Kxy,Kyy",
               comments="")
    print("Wrote:", VORO_OUT, "and", K_CENTROIDS_CSV)

if __name__ == "__main__":
    build_node()
    calculate_voltages()
    compute_route()
    graph()


  

