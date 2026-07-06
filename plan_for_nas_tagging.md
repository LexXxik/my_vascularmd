# Plan — Tagging cells and faces in the Nastran (`.nas`) export

> **Stage 2 deliverable (planning only — no code is changed here).**
> Companion to [`documentation.md`](documentation.md) and
> [`plan_for_extending_volume_mesh.md`](plan_for_extending_volume_mesh.md).
> Goal: make the exported `.nas` carry **per-cell tags** (element IDs / property
> IDs) and **per-face tags** so that COMSOL (and other Nastran readers) can select
> individual regions and boundaries of the mesh.

---

## 1. Where things stand today

- The volume mesh is a flat array on `ArterialTree._volume_mesh = [cells, cell_types, vertices, link_graph]`
  (built in [`mesh_volume`](ArterialTree.py) at `ArterialTree.py:2003`). `cells` rows are
  `[8, g0..g7]` VTK hexahedra.
- Export today is **untagged**. `notebooks/vtk_to_nas.ipynb` does
  `pv.save_meshio(nas_path, volume_mesh)`, and `ArterialTree.write_volume_mesh`
  (`ArterialTree.py:5024`) does the same through `meshio`. The produced
  `Output/bifurcation_116_volume_mesh.nas` (meshio v4.4.6) contains exactly:

  | Card | Count | Meaning |
  |---|---|---|
  | `GRID*` | 93 764 | nodes |
  | `CHEXA` | 89 232 | hex cells — **all with the same (default) property id** |
  | (none) | — | no `PSOLID`, no `CQUAD4`, no components/sets |

  So there is currently **one implicit domain and no boundaries** — nothing to select in COMSOL.

- **A per-cell back-reference already exists but is unused by the exporter.** The
  `link_graph` (returned by `get_volume_link`, `ArterialTree.py:157`; built inside
  `mesh_volume`) tags each cell with `[edge_start, edge_end, crsec_index]` — i.e. which
  branch and which longitudinal cross-section the cell came from. It does **not** yet encode
  radial position (which O-grid layer) or a fluid/wall flag.

- **Boundary classification already exists.** `Simulation.__boundary_patches()`
  (`Simulation.py:61`) builds `wall` / `inlet_i` / `outlet_j` surfaces, and
  `Simulation.write_mesh_files()` (`Simulation.py:166`) already walks every hex face,
  detects boundary faces (a face with no neighbouring cell), and matches each to its patch.
  That is exactly the machinery needed to emit tagged boundary faces.

**Conclusion:** the data needed for tagging is essentially already computed; what is missing
is (a) a richer per-cell tag array, and (b) a writer that puts tags into the `.nas`.

---

## 2. Can individual **faces** be tagged? — Yes

This was called out as something to investigate. The answer is **yes**, and there are two
standard Nastran mechanisms; COMSOL reads both:

1. **Element (cell) tags → `PID` on the `CHEXA` card.**
   `CHEXA EID PID G1 G2 … G8`. The `PID` field groups solid elements. Each distinct `PID`
   (declared by a `PSOLID PID MID` card) becomes a **selectable domain** in COMSOL. This is
   how you tag *cells*.

2. **Face tags → coincident shell elements (`CQUAD4`) with their own `PID`.**
   Nastran does not attach a label to "face 3 of solid element 57" in a way COMSOL imports as
   a boundary. The portable idiom is to emit a **`CQUAD4 EID PID g1 g2 g3 g4`** shell element
   lying exactly on each boundary quad (reusing the hexahedron's own corner node IDs), grouped
   by `PID`. On import COMSOL turns each shell `PID` into a **boundary selection**
   (inlet_0, outlet_1, wall, fluid–wall interface, …). This is the accepted way to carry named
   boundaries through a Nastran file, and it is directly analogous to what
   `Simulation.write_mesh_files` already does to build OpenFOAM patches.

   *(A second option — Nastran `SET`/`Component` grouping of element+face pairs — is more
   solver-specific and less reliably imported by COMSOL, so the CQUAD4 route is recommended.)*

So: **cells are taggable via `CHEXA` PID; faces are taggable via coincident `CQUAD4` PID.**
Both live in the same bulk-data file and can be emitted together.

---

## 3. What we can tag (the tag sources)

Even on the current lumen-only mesh, these per-cell tags are available or cheaply derivable:

| Tag | Source | Status |
|---|---|---|
| **Branch id** (which vessel/bifurcation) | `link_graph[:, 0:2]` (`edge_start,edge_end`) | already computed |
| **Longitudinal index** (cross-section number along the branch) | `link_graph[:, 2]` | already computed |
| **Radial region** (boundary layer *a* / intermediate *b* / core *c*) | the radial index `k` in `ogrid_pattern_faces` (`ArterialTree.py:2640`) — `k < num_a` → a, `< num_a+num_b` → b, else c | **needs a small addition** (emit `k` as a parallel array) |
| **Fluid vs wall domain** | radial region once the wall layer exists | depends on [`plan_for_extending_volume_mesh.md`](plan_for_extending_volume_mesh.md) |

Per-face (boundary) tags:

| Tag | Source | Status |
|---|---|---|
| **wall / inlet_i / outlet_j** | `Simulation.__boundary_patches()` + the no-neighbour-face test in `write_mesh_files` | already computed |
| **fluid–wall interface** | the old lumen boundary faces after a wall is added | depends on the wall plan |

This plan delivers the **mechanism** (rich per-cell tag + tagged writer). It works
immediately with branch/region tags on the current mesh, and automatically gains the
fluid/wall/interface tags once the wall layer from the other plan exists.

---

## 4. Design

Two pieces: **(4.1)** enrich the per-cell tag record, and **(4.2)** a tagged `.nas` writer.

### 4.1 Enrich the per-cell tag (extend `link_graph`)

The cell loop in `mesh_volume` already writes `link_graph[nb_cells:…]`. The cheapest robust
change is to make the O-grid face builder also return a **radial-region code per quad**, then
carry it into a widened `link_graph` (or a parallel `cell_region` array).

- In `ogrid_pattern_faces` (`ArterialTree.py:2591`) and `bif_ogrid_pattern_faces`
  (`ArterialTree.py:2268`): alongside each `faces[nb_faces] = [...]`, record the radial index
  `k` (already the loop variable at `ArterialTree.py:2640`). Return `(faces, region)` where
  `region[i] ∈ {0=a, 1=b, 2=c}` (bucket `k` by `num_a`, `num_b`).
- In `mesh_volume`: widen `link_graph` from 3 to 4+ columns to
  `[edge_start, edge_end, crsec_index, radial_region]` (and later a `domain` column for
  fluid/wall). Because a hexahedron is built from one O-grid quad extruded between two
  cross-sections, the quad's `region` is exactly the cell's radial region — no extra geometry.

This keeps all tag logic inside the mesher (single source of truth) and leaves `link_graph` as
the one place every downstream tool (this exporter, `check_mesh`, the GUI) reads cell metadata.

### 4.2 A tagged `.nas` writer

**Recommendation: hand-write the `.nas`**, in the exact style the repo already uses for
OpenFOAM (`Simulation.write_mesh_files`), swc, and vtk. Reasons:

- Full, version-independent control over `CHEXA` `PID` (cell tags) and `CQUAD4` `PID`
  (face tags) — `pv.save_meshio` gave us no control (hence the untagged file today), and
  meshio's Nastran tag mapping varies by version.
- The hex node order the mesh already uses (VTK hexahedron: bottom `0-1-2-3`, top `4-5-6-7`)
  is the **same** convention as Nastran `CHEXA G1..G8`; the existing meshio file already
  imports correctly, so we can reuse `cells[:,1:]` verbatim.
- It fits the codebase (no new dependency; siblings write their own ASCII formats).

Proposed new method, e.g. `ArterialTree.write_volume_mesh_nas(output, cell_tag="domain", face_tags=True)`:

**Bulk-data layout to emit**

```
$ vascularmd tagged Nastran export
BEGIN BULK
$ --- properties (one PSOLID per cell tag / domain) ---
PSOLID, 1, 1            $ e.g. fluid
PSOLID, 2, 1            $ e.g. wall
$ --- shell properties (one PSHELL per boundary tag) ---
PSHELL, 101, 1, 0.0     $ wall
PSHELL, 102, 1, 0.0     $ inlet_0
PSHELL, 103, 1, 0.0     $ outlet_0  ...
$ --- nodes ---
GRID*, <id>, , x, y, *  / z          (large-field, same as meshio output)
$ --- solid cells, PID = chosen cell tag ---
CHEXA, <eid>, <pid>, g1, g2, g3, g4, g5, g6, +
+, g7, g8
$ --- boundary faces as coincident shells, PID = boundary tag ---
CQUAD4, <eid>, <pid>, n1, n2, n3, n4
ENDDATA
```

**How each part is produced**

1. **`GRID`** — from `_volume_mesh[2]` (vertices). Reuse meshio's `GRID*` large-field
   formatting (the current `.nas` shows the exact column layout) so numeric precision matches
   what already imports cleanly.
2. **`CHEXA` + PID** — from `_volume_mesh[0]` (cells) and the chosen tag column of the widened
   `link_graph`. Map the tag value → a small integer `PID`; keep a `{pid: name}` legend. This
   is the **per-cell tag**. Selectable tag choices: `domain` (fluid/wall), `radial_region`
   (a/b/c), or `branch` (per-vessel), chosen by the `cell_tag` argument.
3. **`PSOLID`** — one per distinct `PID` used by the `CHEXA` cards.
4. **Boundary `CQUAD4` + PID** — obtain the boundary faces **in volume-node IDs** by reusing
   the face-walk already written in `Simulation.write_mesh_files` (`Simulation.py:247–276`):
   for each hex, each face with `GetCellNeighbors == 0` is a boundary face; classify it to a
   patch via the existing nearest-`boundary_patches`-midpoint locator. Emit one `CQUAD4` per
   boundary face with `PID = patch id`. This is the **per-face tag**.
   - This reuse matters: the surface mesh from `mesh_surface` uses a *different* node-ID scheme
     (`__add_node_id_surface`) than the volume mesh (`__add_node_id_volume`), so surface-mesh
     face indices cannot be dropped straight into the `.nas`. The topological no-neighbour-face
     extraction gives boundary quads already expressed in the **volume** node IDs the `CHEXA`
     cards use — the only consistent choice.
5. **Interface faces (FSI):** once the wall layer exists, the fluid–wall interface is the set
   of faces shared by a fluid cell and a wall cell (both have neighbours, but the neighbour's
   `domain` tag differs). Detect them in the same face-walk (`nn == 1` and
   `domain[cell] != domain[neighbour]`) and give them their own `PID`. COMSOL needs this
   interface as an explicit boundary for FSI coupling.

**Legend / provenance.** Also write a small sidecar (e.g. `<output>.tags.txt` or `$`-comment
header) listing `pid → name` for solids and shells, so the COMSOL operator knows which
selection is which.

---

## 5. Two-stage delivery

**Stage A — tag what already exists (independent of the wall plan).**
1. Add the `region` return to `ogrid_pattern_faces` / `bif_ogrid_pattern_faces` and widen
   `link_graph` in `mesh_volume`.
2. Implement `write_volume_mesh_nas(...)` with `CHEXA` PID from a chosen cell tag
   (`branch` or `radial_region`) and boundary `CQUAD4` PID from the reused boundary-face walk.
3. Result: a `.nas` where COMSOL can already select each vessel/bifurcation, each O-grid
   radial band, and each inlet/outlet/wall boundary — on the current lumen-only mesh.

**Stage B — add the FSI tags (after the wall layer exists).**
4. Add the `domain` (fluid/wall) column to `link_graph` from the wall plan's radial region.
5. Emit `PSOLID` fluid vs wall, and the fluid–wall **interface** `CQUAD4` group.
6. Result: a `.nas` ready for COMSOL FSI (fluid domain, wall domain, interface, inlets,
   outlets, wall exterior) — replacing today's "subtract two STL/VTK surfaces in COMSOL" step.

---

## 6. Risks / things to verify during implementation

- **Node-ID base.** Nastran IDs are 1-based; meshio's current file starts `GRID* 1`. Emit
  `GRID`/`CHEXA`/`CQUAD4` with 1-based IDs and offset the `cells` node indices (`+1`)
  accordingly. Verify against the existing importable file.
- **Large-field vs small-field.** meshio used `GRID*` (large field, 16-col) because
  coordinates need the precision. Match that for `GRID`; `CHEXA`/`CQUAD4` integer cards fit
  small field but must respect the 8-char field width / continuation rules (`+` continuation
  for `CHEXA`'s 8 nodes, exactly as meshio already emits).
- **`CQUAD4` winding.** Give interface/wall shells outward-consistent winding. The repo already
  had to solve surface orientation (`fix_orientation` in `bifurcations_to_stl.ipynb`,
  `compute_normals(auto_orient_normals=True)`); reuse that convention so COMSOL normals are sane.
- **Face de-duplication.** Each internal face is shared by two cells; only emit a `CQUAD4`
  once per boundary/interface face (the `cell_id < neighbour_cell_id` guard already used in
  `Simulation.write_mesh_files` handles this).
- **`check_mesh` side effect.** `check_mesh` re-runs `mesh_volume` per segment and overwrites
  `_volume_mesh` (noted in `vtk_to_nas.ipynb`). Tag/export the mesh object captured right after
  `mesh_volume`, or re-mesh before exporting.
- **Bifurcation region codes.** `bif_ogrid_pattern_faces` (`ArterialTree.py:2268`) has more
  complex per-ray bookkeeping than the plain case; deriving the radial `k` there needs care.
  Branch/boundary tags work regardless; the radial-region tag at bifurcations should be
  validated visually (colour cells by tag in ParaView) before trusting it.
- **COMSOL round-trip.** Validate by importing the tagged `.nas` into COMSOL (or re-reading
  with meshio and checking `cell_sets`) and confirming each `PID` shows up as a distinct
  domain/boundary selection with the expected cell/face counts.

---

## 7. Summary

- Faces **can** be tagged (coincident `CQUAD4` `PID`); cells **can** be tagged (`CHEXA` `PID`).
- The tag data is mostly already computed (`link_graph`, `Simulation.__boundary_patches`); the
  only new per-cell info is the radial region, a one-line capture of the existing `k` index.
- Recommended writer is a **hand-written tagged `.nas`** method, matching the repo's existing
  OpenFOAM/vtk/swc writers, giving full control over cell and face `PID`s that
  `pv.save_meshio` does not.
- Deliver in two stages: branch/region + boundary tags now (lumen-only mesh); fluid/wall/
  interface tags after [`plan_for_extending_volume_mesh.md`](plan_for_extending_volume_mesh.md)
  adds the wall layer.
