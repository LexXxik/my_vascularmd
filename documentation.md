# vascularmd — Codebase Documentation

> **Scope of this document (Stage 1).** This is a read-only survey of the entire
> `vascularmd` codebase. It describes the overall capability of the code and then
> gives a concise, per-function description of every file. It performs no code
> changes. It is the reference used by the two Stage 2 planning documents
> ([`plan_for_nas_tagging.md`](plan_for_nas_tagging.md) and
> [`plan_for_extending_volume_mesh.md`](plan_for_extending_volume_mesh.md)).

---

## 1. High-level overview

`vascularmd` turns **vascular centerlines** (points `(x, y, z, r)` + directed
connectivity) into **high-quality structured hexahedral meshes** of the vessel
lumen, ready for CFD (and, via the two Stage 2 plans, FSI). It is backed by two
publications:

- **Modeling & meshing engine** — Decroocq et al., *"Modeling and hexahedral
  meshing of cerebral arterial networks from centerlines"* (Medical Image
  Analysis, 2023). The parametric spline/bifurcation model and O-grid hexahedral
  mesher (`ArterialTree`, `Nfurcation`, `Spline`, `Model`).
- **Interactive software** — Decroocq et al., *"A Software to Visualize, Edit,
  Model and Mesh Vascular Networks"* (IEEE EMBC 2022). The GUI (`Editor`) that
  wraps the engine, letting non-expert users (e.g. clinicians) correct centerlines
  against medical images and mesh them without programming.

### What it can currently do

| Capability | Entry point |
|---|---|
| Import centerlines (swc / edg+nds / txt ideal-model / vtk / tre) | `ArterialTree(...)` constructor |
| Parametric **spline model** of every branch (penalized B-splines, automatic control-point count, tangent continuity across branches) | `ArterialTree.model_network()` |
| **Bifurcation / n-furcation modeling** (VMTK-style key points, rounded apex) | `Nfurcation` class |
| **Surface mesh** — structured quads over lumen | `ArterialTree.mesh_surface()` |
| **Volume mesh** — structured hexahedra with a radial **O-grid** cross-section pattern | `ArterialTree.mesh_volume()` |
| Mesh quality check (scaled Jacobian) | `ArterialTree.check_mesh()` |
| Post-processing: **flow extensions**, **inlet/outlet capping**, **virtual stenosis** from templates, deform-to-target | `add_extensions`, `close_surface`, `deform_surface_to_template`, `deform_surface_to_mesh` |
| Export: **VTK** (surface `.vtk`, volume `.vtk`/`.vtu`), **OpenFOAM** polyMesh + boundary conditions, **Nastran `.nas`** (via meshio/pyvista) | `write_surface_mesh`, `write_volume_mesh`, `Simulation`, `pv.save_meshio` |
| **Tagged Nastran `.nas` for COMSOL FSI** — 2 cell domains (fluid / wall) + 5 boundary selections (interface, inlet, outlet, outer wall, wall-ends), emitted natively at write time | `notebooks/generate_nas_mesh.ipynb` |
| Interactive **GUI** for editing (move points, add/remove branches, change angle, add pathology, mesh, check) across four representation modes | `Editor` class |
| **Medical-image overlay** — confront the network against a NIFTI image (local 2D slices, ray-traced along the current viewpoint) to correct centerlines to the real anatomy | `Editor.load_image / compute_slice` |
| Preprocessing from image segmentation to directed centerlines | `Others/segmentation_to_centerlines.py` |

### Known limitations (from README + code)

- Input centerlines must be **directed** in the flow direction.
- No complex/non-planar multifurcations (coplanar n-furcations only).
- Meshing can fail on hard geometry; failures are fixed by hand-editing input in the GUI.
- The `mesh_volume` mesher fills the **lumen only** (single fluid domain); it grows no
  dedicated vessel-wall band. For FSI, [`notebooks/generate_nas_mesh.ipynb`](notebooks/generate_nas_mesh.ipynb)
  works around this by **reinterpreting the outermost O-grid layer `a` as the solid wall**
  (fluid = `b`+`c`), yielding a two-domain mesh with no extra cells; a real wall band of
  controllable thickness remains the fidelity upgrade in
  [`plan_for_extending_volume_mesh.md`](plan_for_extending_volume_mesh.md).
- The library's built-in `.nas` writers are **untagged** (one implicit domain). Tagging is
  now done at generation time by [`notebooks/generate_nas_mesh.ipynb`](notebooks/generate_nas_mesh.ipynb)
  (see §2.10), following [`plan_for_nas_boundary_tagging.md`](plan_for_nas_boundary_tagging.md).

### The four-graph data model (central design idea)

All geometry/topology lives in the attributes of **NetworkX directed graphs**, one
per processing stage. Each is derived from the previous and can be converted back:

1. **`full_graph`** — one node per raw centerline data point `(x,y,z,r)`; edges = raw connectivity.
2. **`topo_graph`** — collapses chains of degree-2 points into single edges; nodes are typed
   `end` / `bif` / `reg` / `sink`; the collapsed points are stored on the edge as `coords`.
3. **`model_graph`** — same topology, but each edge now carries a `Spline` object and each
   bifurcation node carries an `Nfurcation` object; adds `sep` nodes (bifurcation ends),
   reference vectors `ref`, and rotation angles `alpha`.
4. **`crsec_graph`** — each node/edge carries the **cross-section node coordinates** (`crsec`)
   plus connectivity indices (`connect`) and an integer `id` (base offset into the global
   vertex array) used when the flat mesh arrays are assembled.

The surface and volume meshes themselves are stored on the `ArterialTree` as flat NumPy
arrays (`_surface_mesh`, `_volume_mesh`) plus an optional **`link_graph`** that maps each
face/cell back to `[edge_start, edge_end, cross_section_index]` — the natural substrate for
tagging (see the NAS tagging plan).

These back-end graphs surface directly in the GUI as its **four representation modes**
(*data → topology → model → mesh*, per the EMBC 2022 paper): the user edits at whichever
stage is most convenient — dragging raw data points, trimming branches or rotating angles
in topology mode, moving spline control points / adjusting smoothing in model mode, or
checking cell quality in mesh mode — and each mode can be overlaid on the medical image
and on the other modes simultaneously.

### End-to-end pipeline

```
centerline file
  └─ __load_file / __set_topo_graph        full_graph → topo_graph
       └─ model_network()                  topo_graph → model_graph  (splines + Nfurcations)
            └─ compute_cross_sections(N,d)  model_graph → crsec_graph (surface node rings)
                 ├─ mesh_surface()          crsec_graph → _surface_mesh  (quads)
                 └─ mesh_volume(ratio,a,b)  crsec_graph → _volume_mesh   (O-grid hexahedra)
                      └─ Simulation / write_volume_mesh / pv.save_meshio  → OpenFOAM / VTK / NAS
```

### Module map

| File | Lines | Role |
|---|---|---|
| `ArterialTree.py` | ~5200 | Central class. Owns the 4 graphs, modeling, cross-sections, surface & volume meshing, O-grid pattern, post-processing, I/O. |
| `Nfurcation.py` | ~1480 | Bifurcation/n-furcation object: key points, shape & trajectory splines, separation section, apex smoothing, its own surface mesh. |
| `Spline.py` | ~1470 | `Spline` wrapper around a geomdl B-spline: evaluation, arc-length, projection, intersection, approximation entry point. |
| `Model.py` | ~710 | Penalized least-squares B-spline fitting (the math behind `Spline.approximation`) + model-selection criteria. |
| `Simulation.py` | ~760 | Export volume+boundary mesh to OpenFOAM polyMesh and write BC files; boundary-patch classification (wall/inlet/outlet). |
| `utils.py` | ~1100 | Geometry helpers, 2D O-grid smoothing intersections, mesh neighbor queries, quality metrics, multiprocessing meshing helpers. |
| `Editor.py` | ~3150 | vpython GUI wrapping every `ArterialTree` capability. |
| `Others/segmentation_to_centerlines.py` | ~1610 | Image/skeleton → directed centerline graph → swc/vtk (preprocessing). |
| `main.py` | 93 | Example driver script. |
| `Install/module_check.py` | 34 | Dependency installer. |

---

## 2. Per-file, per-function descriptions

Legend: methods prefixed `__` are private (name-mangled); a leading `_attr` denotes a stored graph/mesh attribute.

### 2.1 `ArterialTree.py` — the central class

**Constructor / getters / setters**

- `__init__(patient_name, database_name, filename=None, automatic_resampling=True)` — Initialise; if a filename is given, load it into `full_graph`, derive `topo_graph`, and optionally resample data density.
- `get_full_graph / get_topo_graph / get_model_graph / get_crsec_graph()` — Return the corresponding graph (warn if not computed).
- `get_surface_mesh()` — Return the surface mesh as a `pv.PolyData` (`vertices`, `faces`).
- `get_volume_mesh(formt="pyvista")` — Return the volume mesh as a `pv.UnstructuredGrid` (default) or a `meshio.Mesh` (`"meshio"`).
- `get_surface_link / get_volume_link()` — Return the `link_graph` array mapping each face/cell to `[edge_start, edge_end, crsec_index]`.
- `get_number_of_faces / get_number_of_cells()` — Sizes of the surface/volume meshes.
- `get_bifurcations()` — List of all `Nfurcation` objects in the model graph.
- `get_inlet_outlet()` — For CFD: list of `[center, tangent, "inlet"/"outlet"]` at every terminal node (in-degree 0 = inlet).
- `check_full_graph()` — Validate that every node's (in,out) degree is an allowed vascular pattern.
- `reset_model_mesh()` — Null out model, crsec, surface and volume (used after any topology edit).
- `set_full_graph / set_topo_graph / set_model_graph / set_crsec_graph(G, ...)` — Replace a graph and cascade-recompute the downstream graphs.
- `automatic_resampling(mind=0.4, maxd=0.6)` — Resample each branch's data points to a target density (points per unit length) for cleaner modeling/editing.
- `__load_file(filename)` — Dispatch on extension (`.swc/.vtk/.vtp/.txt/.tre` or `[edg, nds]` list) to the matching loader; validate.
- `__set_topo_graph()` — Collapse degree-2 chains into edges, label nodes `end/bif/reg/sink`, store `full_id` mappings, relabel to integers.

**Modeling (`topo_graph` → `model_graph`)**

- `model_network(radius_model=True, criterion="AIC", ...)` — Top-level: DFS the tree, build an `Nfurcation` at every `bif` node, fit a `Spline` to every plain vessel, then compute reference-vector rotations. Produces `model_graph`.
- `__model_furcation(n, ...)` — Extract one bifurcation's parameters from data (approximate daughter branches, find the **apex** by spline intersection, cut mother/daughter splines at fixed radius-multiples, build the `Nfurcation`, splice `sep`/`bif` nodes and edges into the model graph). Handles the merge/combine cases when branches are too short (→ trifurcation).
- `__reorder_branches(n, splines, ...)` — For >2 daughters, order branches by inter-branch angle so a coplanar furcation is meshed consistently.
- `__model_vessel(e, ...)` — Fit a single vessel spline between two nodes with end-point/tangent constraints (handles `sink` merges).
- `__compute_rotations(n)` — Propagate a reference vector along a branch chain and compute the per-edge twist angle `alpha` + `connect` index so that cross-section node #0 is continuous across bifurcations (this is what keeps the hex mesh watertight and untwisted).
- `__combine_nfurcation(n)` — Merge two adjacent bifurcations that share a very short segment into a single combined furcation mesh.

**Cross-sections (`model_graph` → `crsec_graph`)**

- `compute_cross_sections(N, d, parallel=True)` — Build the `crsec_graph`. For every bifurcation call `Nfurcation.compute_cross_sections`; for every plain edge call `__segment_crsec`. Stores `crsec` (node ring coordinates), `center`, and `connect` (node-index rotation). `N` = nodes per ring (multiple of 4/8), `d` = longitudinal density. `parallel=True` uses a process `Pool`.
- `recompute_cross_sections(n)` — Locally recompute one branch's cross-sections after an edit (walks to bounding `sep`/`end` nodes).
- `furcation_cross_sections(n)` — Serial version of the bifurcation cross-section step for node `n`.
- `vessel_cross_sections(e)` — Serial version of the vessel cross-section step for edge `e`.
- `__add_node_id_surface(nds, edg)` — Assign each crsec node/edge its base **`id`** (running offset) into the global surface-vertex array; set `_nb_nodes`.
- `__add_node_id_volume(num_a, num_b, nds, edg)` — Same, but with the **O-grid node count per cross-section** (`nb_nds_ogrid`, and a larger count at bifurcations). **This is one of the count formulas that must change when adding an outer wall layer.**
- `__segment_crsec(spl, num, N, v0, alpha)` — Compute `num+2` rings of `N` nodes along a vessel spline, transporting/rotating the reference vector `v0` by `alpha`.
- `__single_crsec(spl, t, v, N)` — One ring: rotate `v` around the tangent in `N` steps and project each direction to the vessel surface (`spl.project_time_to_surface`).

**Surface meshing (`crsec_graph` → `_surface_mesh`)**

- `mesh_surface(edg=[], link=True)` — Assemble the structured quad surface mesh: for every edge, write the first ring, the edge rings, the last ring (respecting `connect` rotation and normal flip), emitting quad faces `[4, i, j, k, l]`. Fills `link_graph` `[e0, e1, crsec_index]` per face. Stores/returns the mesh.

**Volume meshing (`crsec_graph` → `_volume_mesh`) — the O-grid**

- `mesh_volume(layer_ratio=[0.2,0.4,0.4], num_a=4, num_b=4, edg=[], link=True)` — Assemble the structured **hexahedral** volume mesh. For each edge it computes the O-grid **vertices** and **faces** of every cross-section (`ogrid_pattern_vertices/faces`, or `bif_ogrid_pattern_*` at bifurcations), then forms a hexahedron per O-grid quad by joining cross-section *i* to cross-section *i+1* (`cells = [8, faceA(4), faceB(4)]`). Handles connect-rotation + normal flip. Fills `link_graph` per cell. Stores `[cells, cell_types, vertices, link_graph]`. **This is the primary function targeted by both Stage 2 plans.**
- `write_pyvista_mesh_from_vtk(vertices, cells)` — Build a VTK unstructured grid of hexahedra explicitly and round-trip it through a temp `.vtu` (legacy helper).
- `ogrid_pattern(center, crsec)` — Convenience: return `(vertices, faces)` of the plain O-grid for a ring using the stored `_layer_ratio/_num_a/_num_b`.
- `ogrid_pattern_vertices(center, crsec, layer_ratio, num_a, num_b)` — **Core geometry.** For each of the four quarters of a ring, build the central square (side = `layer_ratio[2]·radius`) then cast a "ray" from every surface node inward: `layer a` (surface→`pb`, `num_a` cells, thickness `layer_ratio[0]·r`), `layer b` (`pb`→square, `num_b` cells, thickness `layer_ratio[1]·r`), `layer c` (inside square). Returns `nb_nds_ogrid` node coordinates. **The surface node `crsec[i]` is the first vertex of each ray — i.e. the outermost point. An outer wall layer is added by prepending vertices _before_ it.**
- `ogrid_pattern_faces(N, num_a, num_b)` — Quad connectivity of one O-grid cross-section (ray-by-ray, closing the circle and sharing central-square edges). Returns `(N/8)·(2·(num_a+num_b+2)+N/8)·4` quads. **The count/ordering here directly determines the hex count and per-cell radial position.**
- `bif_ogrid_pattern(center, crsec, N, nbif, ...)` — Convenience wrapper returning bifurcation-section `(vertices, faces)`.
- `bif_ogrid_pattern_faces(N, nbif, num_a, num_b)` — O-grid quad connectivity for the **half-sections** of a bifurcation separation plane (shared edges between the `nbif` half-sections; more complex than the plain case).
- `bif_ogrid_pattern_vertices(center, crsec, N, nbif, layer_ratio, num_a, num_b)` — Vertex coordinates for the bifurcation O-grid half-sections (same a/b/c layering logic applied per half-section).
- `__reorder_faces(h, f, N, num_a, num_b)` — Reorder the bifurcation half-section faces so the O-grid of the incoming vessel connects to the correct rotated half-section of the bifurcation.

**Graph conversions (going back up the chain)**

- `topo_to_full(replace=True)` — Rebuild `full_graph` from `topo_graph` (re-expanding collapsed edge points). A second, commented alternative version is kept in the source.
- `model_to_full(replace=True)` — Sample the spline of each model edge into a dense `full_graph`, and rederive `topo_graph`.

**Post-processing**

- `close_surface(edg=[], layer_ratio=[0.2,0.3,0.5], num_a=10, num_b=10)` — Cap open ends by appending an O-grid disk (`ogrid_pattern_vertices/faces`) at each `end` node (used to make the surface watertight; see the STL notebook).
- `open_surface(edg=[])` — Re-mesh the surface without caps.
- `add_extensions(edg=[], size=6)` — Add straight cylindrical **flow extensions** at inlets/outlets (new spline + node/edge + recompute cross-sections). Size is a multiple of the radius.
- `remove_extensions(edg=[], prop=3)` — Remove previously added extensions.
- `load_pathology_template(path)` — Load a stenosis/aneurysm template (cross-section text files + `info.txt`) as `(template, radius, center)`.
- `deform_surface_to_template(edg, t0, t1, template, ..., method="bicubic", rotate=0)` — Deform the surface between times `t0..t1` on a branch to match a loaded pathology template (adds a virtual stenosis).
- `deform_surface_to_mesh(mesh, edges=[], search_dist=40)` — Move every cross-section surface node onto a **target mesh** by ray-casting outward from the section center (used in the STL notebook to fit the inner model onto the outer wall; `search_dist` must be limited to avoid rays hitting the wrong far surface).
- `__intersection(mesh, center, coord, search_dist=40)` — Ray-trace helper for the two `deform_*` methods.
- `extract_vessel(e, volume=False) / extract_furcation(n, volume=False)` — Return a standalone mesh of a single vessel/furcation.
- `subgraph(nodes)` — Return a sub-`ArterialTree` for the given nodes.
- `make_inlet(n)` — Re-root the directed tree at node `n`.
- `crop_branch_degree(max_deg)` — Delete branches beyond a branching depth.
- `merge_branch(n, mode="topo")` — Merge a branch across a node (used when bifurcations combine).
- `remove_branch(e, preserve_shape=True, from_node=False)` — Remove a branch, optionally re-modeling the affected bifurcation (`remove_bifurcation_branch` nested helper).
- `translate / scale(val)` — Rigid transforms of all data points.
- `set_minimim_radius(val) / transform_radius() / modify_branch_radius(edg, eps)` — Bulk radius edits.
- `resample(p) / low_sample(p)` — Change data-point density.
- `delete_data_point / delete_data_edge / add_data_edge / add_data_point(...)` — Fine-grained centerline edits.
- `modify_control_point_coords/radius`, `modify_data_point_coords/radius` — Move a single control/data point.
- `rotate_branch(normal, edg, alpha)` — Rotate a branch about an axis (change bifurcation angle).
- `smooth_spline(edg, l, radius=False)` — Increase spline smoothing `lambda` on an edge.
- `convert_edge_mode(edg, mode1, mode2)` — Convert an edge id between the graph levels.
- `angle(distance=None, mode="topo")` — Compute branch angles.
- `count_nodes()` — Node counts per graph.
- `add_noise_centerline(std, normal=False) / add_noise_radius(std)` — Add synthetic noise (for robustness studies).
- `check_mesh(thres=0, edg=[], mode="crsec")` — **Quality check.** Re-mesh each vessel/bifurcation individually, compute pyvista `scaled_jacobian`, and flag any cell ≤ `thres`; return a per-surface-face boolean field plus lists of failed edges/bifs. Uses the `link_graph` to map bad cells back to surface faces.
- `show(points=False, centerline=False, control_points=False)` — Matplotlib 3D preview.

**I/O**

- `__swc_to_graph / __edg_nds_to_graph / __txt_to_graph / __vtk_to_graph / __tre_to_graph(filename)` — Loaders producing a `full_graph`.
- `write_volume_mesh(output="volume_mesh.vtk")` — Write the volume mesh via **meshio** (`("hexahedron", cells[:,1:])`). **The natural place to add tagged/`.nas` export.**
- `write_surface_mesh(output=...)` — Save the surface `pv.PolyData`.
- `write_vtk(type, filename)` — Write centerlines (`full/topo/spline`) as a VTK polyline.
- `write_edg_nds(filename) / write_swc(filename)` — Export centerlines back to edg+nds / swc.
- `write_data_as_circles()` — Debug view of data points as oriented circles.

---

### 2.2 `Nfurcation.py` — bifurcation / n-furcation object

Represents one coplanar furcation defined either from splines or from cross-sections
(`model="spline"|"crsec"`). Builds a smooth, watertight junction mesh.

- `__init__(model, args)` — Store end sections, apex sections and apex points; derive shape splines, key points, geometric center `X`, separation points `SP`, common points `CT`, and trajectory splines `tspl`. Set default `N=24`, `d=0.2`, smoothing radius `R`.
- `check_input_parameters()` — Verify (and translate if needed) that each apical cross-section is consistent with its apex point.
- Getters: `get_spl / get_tspl / get_endsec / get_apexsec / get_X / get_CT / get_SP / get_crsec / get_AP / get_tAP / get_N / get_d / get_R / get_n`.
- `get_reference_vectors()` — Reference vector for each end section (aligns cross-section node #0 with the separation plane).
- `get_angles()` — Furcation angle(s) between daughter branches (degrees) + the vectors.
- `get_crsec_normals()` — Surface normals of the bifurcation nodes in crsec structure (used for apex smoothing).
- `get_curves()` — All longitudinal curves crossing the separation plane, plus a per-curve local frame (for apex smoothing).
- `curves_to_crsec(curve_set)` — Inverse of `get_curves` — write smoothed curves back into the crsec structure.
- `set_crsec / set_R` — Setters.
- `__set_apexsec / set_apexsec_radius(radius, branch_id)` — Set/adjust apical cross-section radius while keeping the apex fixed.
- `rotate_apex_section(ind, angle)` — Rotate an apical cross-section about the apex point.
- `correct_apex_section(ind)` — (Commented out) dichotomy correction of an apex section that dips below another branch's surface.
- `__set_spl()` — Fit the **shape splines** (one per branch) from end+apex sections.
- `__set_AP / __set_tAP()` — Find apex point(s) via spline intersection and their parametric times.
- `__set_key_pts()` — VMTK-style key points on each branch.
- `__set_X()` — Geometric center of the furcation.
- `__set_SP()` — The two separation points on the surface.
- `__set_CT()` — The two "common" points (top/bottom of the separation plane); recenters `X`.
- `__set_tspl()` — The **trajectory splines** from each end section to `X`.
- `mesh_to_crsec(mesh)` — Read node coordinates back out of a meshed surface into the crsec structure (after smoothing/deformation).
- `show(nodes=False)` — Matplotlib preview of splines, key points, and optionally the mesh nodes.
- `mesh_surface()` — Build the bifurcation's own quad surface mesh from its crsec (end sections, connecting nodes, separation plane).
- `compute_cross_sections(N, d, end_ref=None)` — **Main mesher.** Build the separation-plane nodes and end-section rings, the connectivity indices, and the connecting node rings; then relax, apex-smooth, relax again. Stores `_crsec = [end_crsec, bif_crsec, nds, connect_index]`.
- `__bifurcation_connect(tind, ind, P0, P1, n)` — The `n` connecting nodes between an end node and its separation node, following a fitted trajectory spline and projected to the surface.
- `__separation_section(N)` — Nodes of the separation plane (the `CT`/`SP`/`AP` skeleton subdivided into `N/2`-node half-arcs).
- `__separation_segment(P1, P2, n)` — `n` surface nodes on the arc between two separation-plane key points.
- `__end_sections(N, end_ref=None)` — The `N`-node rings at each branch end.
- `relaxation(n_iter=5)` — Alternate Laplacian smoothing and re-projection to improve cell shape.
- `smooth(n_iter)` — One Laplacian smoothing pass of the surface mesh.
- `optimal_smooth_radius(degree)` — Apex smoothing radius as a function of the furcation angle.
- `smooth_apex(radius)` — Round the sharp apex using the inscribed-circle polyline smoothing in a local frame.
- `deform(mesh)` — Project the bifurcation nodes onto an external target surface (used inside `relaxation`).
- `__intersection(mesh, center, coord)` — Ray-trace helper for `deform`.
- `send_to_surface(O, n, ind)` — Project a point in direction `n` onto the multi-spline surface (checks all shape splines).
- `__projection(O, n, c0, c1, ind)` — Bisection projection onto one shape spline's surface.

---

### 2.3 `Spline.py` — spline wrapper

Wraps a `geomdl` B-spline (dimension 3 = spatial, 4 = spatial + radius) and stores the
underlying penalized `Model`(s).

- `__init__(control_points=None, knot=None, order=None)` — Build a cubic B-spline (uniform knot by default).
- Getters: `get_spl / get_knot / get_control_points / get_points / get_times / get_length / get_nb_control_points / get_model / get_data / get_lbd`.
- `set_spl(spl) / set_control_points(P)` — Replace the underlying curve.
- `__set_length_tab()` — Precompute an arc-length table + KD-tree of sample points + a `delta` step (used by all length/projection queries).
- `_set_lambda_model(lbd)` — Change the smoothing `lambda` of the spatial/radius models and rebuild control points.
- `first_derivative / second_derivative / tangent(t, radius=False)` — Evaluate derivatives / unit tangent at time(s) `t`.
- `point(t, radius=False)` — Evaluate the point (optionally including radius) at `t` (scalar or list).
- `radius(t)` — Radius component at `t`.
- `interpolation(D)` — Interpolate data points exactly (centripetal).
- `approximation(D, end_constraint, end_values, derivatives, radius_model=True, ..., criterion="None")` — **Main fitting entry point.** Choose control-point count (fixed / RMSE / Akaike), fit a spatial `Model` and a separate radius `Model`, optionally enforce minimum tangent magnitude and a curvature-vs-radius constraint; store both models + control points.
- `__nb_control_points_accuracy(data, thres, criterion) / __nb_control_points_akaike(data)` — Choose the number of control points.
- `__optimize_model(model, criterion, max_distance)` — Golden-section search for the smoothing `lambda` under a fidelity criterion and a max-distance cap.
- `__constraint_curvature(spatial_model, radius_model)` — Increase `lambda` until the curvature radius stays above the vessel radius (avoids self-intersecting tubes).
- `__uniform_knot(p, n)` — Uniform knot vector.
- `split_time(t) / split_length(l)` — Split into two splines at a time/length (re-fitting sub-models).
- `reverse / copy_reverse()` — Reverse orientation.
- `length / mean_radius()` — Arc length / mean radius.
- `curvature(t) / curvature_radius(T)` — Curvature and its reciprocal.
- `length_to_time(L) / time_to_length(T)` — Convert between arc length and parameter using the precomputed table.
- `resample_time(n, t0, t1, include_start)` — `n` equally-arc-spaced times (used to place cross-sections).
- `transport_vector(v, t0, t1)` — Parallel-transport a normal vector along the spline (keeps cross-sections un-twisted).
- `project_time_to_surface(v, t)` — Point on the tube surface at time `t` in direction `v` (`center + r·v̂`).
- `project_point_to_surface(pt) / project_point_to_centerline(pt)` — Nearest surface / nearest centerline parameter for an arbitrary 3D point (KD-tree + segment projection).
- `distance(D)` — Per-point `(time, distance)` to the spline.
- `estimated_point / estimated_first_derivative / estimated_tangent / estimated_curvature(D)` — Spline estimates at the projections of `D`.
- `RMSE / RMSEder / RMSEcurv(...)` — Fit-error metrics.
- `first_intersection / first_intersectionv2(spl,...)` — Find the furthest surface-intersection (apex) between two tubes.
- `intersection(spl, v0, t0, t1)` — Bisection surface intersection along a transported vector.
- `first_intersection_centerline / first_intersection_centerlinev2(spl,...)` — Centerline-based intersection (radius vs distance sign change).
- `show(knot=False, control_points=True, data=[])` — Matplotlib preview.

---

### 2.4 `Model.py` — penalized B-spline fitting

The linear-algebra engine behind `Spline.approximation`.

- `__init__(D, n, p, end_constraint, end_values, derivatives, lbd, knot=None, t=None)` — Set up a penalized least-squares B-spline of `n` control points, degree `p`, smoothing weight `lbd`, with optional fixed end points/tangents; build matrices and solve.
- `get_data / get_length / get_n / get_t / get_order / get_knot / get_lbd` — Getters.
- `set_lambda(lbd)` — Change smoothing weight and re-solve.
- `get_magnitude()` — Magnitudes of the end tangents (with sign).
- `quality(criterion="CV")` — Model-selection / fidelity score: `AIC / AICC / SBC / CV / GCV / SSE / RMSE / max_dist / RMSEder` (uses the hat matrix `H`).
- `__solve_system()` — Solve `(NᵀN + λΔ) P = NᵀD − λQ` (pseudo-inverse), re-insert the fixed end control points; fall back to resampling on numerical failure.
- `__compute_matrices()` — Build the basis-function matrix `N`, difference-operator penalty `Δ = UᵀU`, and the constraint vectors `Q1/Q2/Pt`, for both the "fixed-derivative" and "fixed-tangent" formulations.
- `uniform_knot / uniform_averaging_knot / averaging_knot(t)` — Knot-vector constructions.
- `chord_length_parametrization()` — Chord-length parameter values for the data.
- `plot_basis_functions()` — Debug plot.
- `__basis_functions(t)` — B-spline basis values at `t` (Cox–de Boor recurrence).
- `__basis_functions_derivative(t)` — First-derivative basis values (via geomdl helpers).

---

### 2.5 `Simulation.py` — OpenFOAM export & boundary conditions

- `__init__(arterial_tree, output_dir)` — Require a volume mesh; classify inlet/outlet terminal nodes; build the boundary-patch MultiBlock and cache the pv volume mesh.
- `__boundary_patches()` — Build a `pv.MultiBlock` of boundary surfaces: `wall` = full surface mesh; one `inlet_i`/`outlet_j` O-grid disk per terminal (via `ogrid_pattern`). **This is the existing wall/inlet/outlet face classification — directly reusable for face tagging.**
- `file_description / file_header / top_separator / bottom_separator()` — OpenFOAM file boilerplate.
- `write_mesh_files()` — **The polyMesh writer.** For every hex cell, extract its 6 faces via VTK; a face with no neighbor is a boundary face (matched to the nearest boundary-patch midpoint via a KD-tree locator to get its patch id), otherwise it is internal (owner < neighbour). Write `points`, `faces`, `owner`, `neighbour`, `boundary`. **The face-owner/boundary logic here is the template for extracting and tagging boundary faces in the NAS plan.**
- `write_mesh_files_numpy()` — Same output, pre-counting faces and filling NumPy arrays (faster/lower-memory variant).
- `write_pressure_boundary_condition_file()` — Write `0/p` (zeroGradient on wall/inlet, fixedValue at outlets).
- `write_velocity_boundary_condition_file(vel_mag)` — Write `0/U` (noSlip wall, fixedValue inlet from spline tangent × magnitude, zeroGradient outlet).
- `import_field(filename, fieldname, outfieldname)` — Load a simulation result field back onto the mesh.
- `display_crsec(step, fieldname)` — (Stub) cross-section field viewer.

---

### 2.6 `utils.py` — helpers

**Geometry:** `cart2pol / pol2cart`, `rotate_vector(v, axis, theta)` (Rodrigues), `angle(v1, v2, axis, signed)`, `directed_angle / directed_angle_negative`, `linear_interpolation / lin_interp(p0, p1, num)` (the workhorse used all through the O-grid), `length_polyline(D)`, `resample(D, num)`, `averaging_knot / optimal_knot`, `order_points` (shortest-path reordering via KD-tree + Dijkstra).

**2D O-grid smoothing** (inscribed-circle apex rounding): `intersection_segment_segment / intersection_cercle_cercle / intersection_segment_cercle(...)` (find where offset polyline pieces intersect) and `smooth_polyline(data, radius, show=False)` (iteratively remove intersections to round a corner).

**Mesh neighbor queries:** `neighbor_faces(vertices, faces)`, `neighbor_faces_id`, `neighbor_faces_normals`, `neighbor_vertices_id`, `neighbor_vertices_coords` — adjacency used by smoothing/deformation.

**Validation:** `distance(mesh1, mesh2)` (implicit distance), `quality(mesh, metric='scaled_jacobian')` (mean/min/max cell quality).

**Misc:** `split_tubes(...)` (VMTK centerline → per-tube swc), `evaluate_model(spl1, data, spl2)` (model comparison metrics), `chord_length_parametrization(D)`.

**Multiprocessing meshing helpers** (module-level so they pickle): `parallel_bif(bif, N, d, end_ref)`, `parallel_apex(spl)`, `segment_crsec(spl, num, N, v0, alpha)` and `single_crsec(spl, t, v, N)` (standalone copies of the `ArterialTree` cross-section methods), and a standalone `intersection(...)`.

---

### 2.7 `Editor.py` — vpython GUI

A large interactive editor — the subject of the EMBC 2022 software paper. Built on the
`vpython` package (chosen over Blender/Paraview for interaction-oriented simplicity), it
holds an `ArterialTree` and mirrors every capability as buttons/sliders. The layout is a
3D viewer + message bar + edit menu + import/export field + four independently toggled
representation panes (data / topology / model / mesh). Editing is guided by an optional
medical-image overlay so the network can be matched to the true patient anatomy at every
step. Typical use cases from the paper: repairing noisy open-database centerlines (e.g.
BraVa) into CFD-valid meshes, and rapidly generating mesh databases by varying branch
angle or topology from a single input network. Only the meshing-relevant methods are
highlighted here.

- `__init__(tree, width, height)` — Build the whole vpython scene, sliders and widgets.
- `reset_scene / do_nothing / node_color / barycenter_finder` — Scene helpers.
- `crop_network / get_closest_data_point / update_visibility_*` — Selection/visibility.
- `load_image / compute_slice / slice_opacity` — Overlay a medical image slice.
- `create_template_pathology / show_hide_template / save_pathology_template / next_template / locate_pathology / add_pathology` — Virtual-pathology workflow.
- `update_mesh_representation / show / hide / lock / unlock / disable` — Display modes.
- `update_save_directory / update_save_filename / load_pathology_template / save` — I/O widgets.
- `apply_function / create_elements / update_graph / converted_selection / refresh_display` — Bridge between GUI selection and `ArterialTree` calls; rebuild vpython objects from graphs/meshes.
- `output_message / keyboard_control / update_edition_mode / resample_nodes / move / move_edges / update_spline / drop / select / unselect` — Interaction/edit handlers.
- `select_smooth_parameter` — Smoothing slider.
- `update_mesh_parameters(b)` — **Read the meshing parameters from the GUI:** `N`, `d`, `a/b/c` (normalised to sum 1 = `layer_ratio`), `num_a` ("layers in the boundary"), `num_b` ("layers in the intermediate part"), plus max-distance modeling caps. (Confirms the user-facing meaning of the O-grid parameters.)
- `deform_mesh` — Call `tree.deform_surface_to_mesh(target)`.
- `check_mesh` — Call `tree.check_mesh` and colour failed faces red.
- `close_mesh / disable_close_check` — Call `tree.close_surface / open_surface`.
- `manage_extensions` — Call `tree.add_extensions / remove_extensions`.
- `mesh_volume` — Call `tree.mesh_volume(layer_ratio=[a,b,c], num_a, num_b, edg=selection)`.
- `mesh_surface` — Call `tree.compute_cross_sections(N,d)` then `tree.mesh_surface`.
- `update_edge_size / update_node_size / update_opacity_state` — Rendering tweaks.

---

### 2.8 `Others/segmentation_to_centerlines.py` — preprocessing

Standalone tools to get **directed** centerlines from images/skeletons/third-party tools.

- **I/O / skeletonisation:** `swc_to_networkx`, `load_nii`, `write_centerline_nii`, `extract_radius_from_image` (distance transform → radius), `sknw_to_networkx` (skeletonize + sknw), `voreen_to_networkx`, `vesselvio_to_networkx`, `register_coords`, `networkx_to_vtk`, `voreen_to_vtk`, `vesselvio_to_vtk`, `sknw_to_vtk`, `networkx_to_swc`.
- **Graph conversion between representations:** `sknw_topo_to_undirected_graph`, `topo_to_undirected_graph`, `topo_to_directed_graph`, `undirected_graph_to_topo`, `orient_tree`, `directed_graph_to_topo`, `get_ordered_nds`, `add_clique_multifurcation`, `max_node_id`.
- **Attributes / cleanup:** `add_length`, `add_branch_id`, `add_angle`, `angle`, `add_weight`, `find_branch_degree`, `add_orientation`, `undirected_graph_to_oriented_tree` (the main pipeline: orient + prune + build cliques for multifurcations), `find_inlet`, `crop_branch_degree`, `remove_multifurcations`, `remove_small_end_branches`, `remove_small_branches`, `keep_largest_cc`.
- **Analysis:** `graph_summary`, `find_cycles`, `degree_distribution`, `degree_distribution_directed`, `number_loop`, `find_bulges`, `find_connected_components`.
- **Top-level driver:** `write_swc(outfile, path_img, algo="sknw", ...)` — image → oriented tree → swc.

---

### 2.9 `main.py` & `Install/module_check.py`

- `main.py` — Example driver: load a centerline, `model_network`, `compute_cross_sections`, `mesh_surface`, `mesh_volume`, save VTK, write OpenFOAM (`Simulation`), or launch the `Editor`; plus a single-bifurcation `Nfurcation` example.
- `Install/module_check.py` — Try-import each dependency (`numpy, pyvista, scipy, networkx, nibabel, geomdl, vpython`) and `pip install` the missing ones.

---

### 2.10 `notebooks/generate_nas_mesh.ipynb` — native tagged Nastran (`.nas`) export for COMSOL FSI

Implements [`plan_for_nas_boundary_tagging.md`](plan_for_nas_boundary_tagging.md) end-to-end:
it emits, **at write time and with only native `ArterialTree` functions**, the domain and
boundary property-ids (PIDs) a COMSOL fluid–structure-interaction model needs — so the setup
is scriptable via LiveLink with **zero manual face-picking**. Under the proof-of-concept
paradigm the existing O-grid **layer `a`** (outermost band) is reinterpreted as the **solid
wall** and parts `b`+`c` as the **fluid**, giving a full two-domain FSI mesh from today's
lumen-only mesher with no extra cells. It builds the mesh once (`model_network` →
`compute_cross_sections` → `mesh_volume`) then writes three files of increasing richness:

- **One-domain** — every `CHEXA` gets PID 1 via `mesh.cell_data["nastran:ref"]` (meshio's
  Nastran writer emits it into the element PID field; no `pv.save_meshio` hop needed).
- **Two-domain** — per-cell fluid/wall split. A faithful mirror of `ogrid_pattern_faces`
  (`ogrid_faces_with_region`, asserted equal to the library's own faces) flags each O-grid
  quad's radial band (`is_layer_a`, `k ≤ num_a`); tiled across every slab (validated against
  `get_volume_link()` and by centroid radius) it gives PID 2 = wall (layer a), PID 1 = fluid.
- **Full FSI** — the two `CHEXA` domains **plus five `CQUAD4` boundary selections**, each a
  *coincident shell* that reuses the hex faces' corner node IDs and carries its own PID:

  | PID | Group | Native derivation |
  |---|---|---|
  | 101 | fluid↔wall interface (a\|b) | faces shared by the two domains — split the mesh per PID (`threshold`), `extract_surface` each, intersect. Interior → FSI coupling |
  | 111 | inlet fluid cap | inlet terminal's O-grid disk (`ogrid_pattern_faces` + node `id`), fluid quads only (`~is_layer_a`) → velocity BC |
  | 131 | outlet fluid cap(s) | same at each outlet terminal → pressure BC |
  | 201 | outer wall surface | all boundary faces **minus** every terminal cap (= layer-a outer surface, the lumen surface) → free / external pressure |
  | 221 | wall-end annulus rings | the wall (`is_layer_a`) quads of the terminal disks → Fixed Constraint (removes rigid-body motion) |

  Terminals are identified by the same criterion as `get_inlet_outlet()` (inlet = `end` node
  with in-degree 0). Every cap is reconstructed **exactly** from `ogrid_pattern_faces` offset
  by the section's base vertex id (`crsec_graph.nodes[n]['id']`) — the same bookkeeping
  `mesh_volume` uses — so no geometric tolerance is involved. Each cell self-verifies:
  boundary faces partition exactly into `{outer ∪ caps}`, the interface is disjoint from the
  boundary, per-PID counts match the file, tags survive a meshio round-trip, and cap radii
  order correctly (fluid inside annulus). meshio shares one global element-ID sequence across
  the `CHEXA` and `CQUAD4` blocks, so EIDs never collide. On import COMSOL turns each PID into
  a named domain/boundary selection (**with "Allow partitioning of shells" off**).

---

## 3. Context relevant to the Stage 2 plans

Discovered while reading the repository (grounds the two plans in the user's real workflow):

- **Target application is COMSOL FSI.** The dataset (`Data/BG0014.CNG/`) stores each BraVa
  bifurcation as an **inner/outer SWC pair** (`bifurcation_X.swc` = lumen,
  `bifurcation_X_outer.swc` = wall) that share identical `(x,y,z)` topology and differ only
  in radius; the wall thickness is a uniform **0.5** (e.g. r 0.62 → 1.12).
- **`.nas` export already works** but *was* **untagged** at the time the plans were written:
  `notebooks/vtk_to_nas.ipynb` called `pv.save_meshio(nas_path, volume_mesh)`, producing a
  single implicit domain with **no property/component IDs**. This is the gap now closed by
  [`notebooks/generate_nas_mesh.ipynb`](notebooks/generate_nas_mesh.ipynb) (§2.10), which tags
  domains and boundaries at write time.
- **The wall is currently produced as two separate surfaces**, not as hex cells.
  `notebooks/bifurcations_to_stl.ipynb` models inner and outer separately, deforms the inner
  onto the outer via `deform_surface_to_mesh`, caps both with `close_surface`, and relies on
  **COMSOL to subtract them** (`wall domain = outer solid − inner solid`). There is no
  structured hex wall layer today.
- **A per-cell/per-face back-reference already exists:** the `link_graph`
  (`get_volume_link` / `get_surface_link`) tags each cell/face with
  `[edge_start, edge_end, crsec_index]`. It does **not** yet encode radial position
  (fluid-vs-wall) — an extension both plans build on.
- **Boundary (inlet/outlet/wall) classification already exists** in
  `Simulation.__boundary_patches()` and `get_inlet_outlet()`.

These five facts are the starting point for
[`plan_for_nas_tagging.md`](plan_for_nas_tagging.md) and
[`plan_for_extending_volume_mesh.md`](plan_for_extending_volume_mesh.md).
