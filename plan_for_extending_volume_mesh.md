# Plan — Extending the volume mesh with an outer wall layer

> **Stage 2 deliverable (planning only — no code is changed here).**
> Companion to [`documentation.md`](documentation.md) and
> [`plan_for_nas_tagging.md`](plan_for_nas_tagging.md).
> Goal: add a **new outermost layer of hexahedral cells** to the volume mesh —
> a structured vessel-wall shell outside the current lumen surface — with a
> **user-controllable thickness** (and controllable number of sub-layers).

---

## 1. What "add the 4th / outermost layer" means here

### The current radial structure of the volume mesh

`mesh_volume` (`ArterialTree.py:2003`) builds each cross-section with the **O-grid** pattern
(`ogrid_pattern_vertices`, `ArterialTree.py:2515`). Going **from the lumen wall inward to the
centre**, a cross-section has three radial *parts* (the `layer_ratio = [a, b, c]`,
`a+b+c = 1`, radial fractions of the vessel radius):

| Part | What it is | Radial thickness | # cell layers |
|---|---|---|---|
| **a** — boundary layer | thin cells hugging the lumen wall | `layer_ratio[0]·r` | `num_a` |
| **b** — intermediate | transition to the core | `layer_ratio[1]·r` | `num_b` |
| **c** — core (O-grid square) | the central square that removes the axis singularity | `layer_ratio[2]·r` | (square block) |

Each surface node `crsec[i]` is the **outermost** point of its radial "ray"
(`ArterialTree.py:2561–2577`: the ray is built `lin_interp(crsec[i], pb, num_a+2) + … inward`).
So today **the outermost surface of the volume mesh *is* the lumen wall** — there is nothing
outside it. The mesh fills the blood domain only.

### The requested extension

Add a **fourth radial part — call it `d` (the wall)** — *outside* part `a`, i.e. outside the
current lumen surface, subdivided into `num_d` cell layers of total thickness `wall_thickness`.
After the change, the radial structure from the outside in becomes:

```
   NEW  wall exterior surface
    │   part d  (wall)         num_d cell layers, thickness = wall_thickness   ← NEW outermost
   ─┼── lumen surface  (crsec[i])  ← was the outermost surface, becomes the fluid–wall interface
    │   part a  (boundary)     num_a layers
    │   part b  (intermediate) num_b layers
    │   part c  (core square)
   centre
```

The single volume mesh then contains **both** the fluid (lumen: parts a+b+c) and the solid wall
(part d) as conformal, node-matched hexahedra — the arrangement COMSOL FSI wants, replacing the
current "model inner and outer surfaces separately and subtract them in COMSOL"
(`notebooks/bifurcations_to_stl.ipynb`).

> **On the "4th layer / thickness of the 5th layer" wording.** The task phrases this as adding
> the *4th* layer and controlling the *5th*. Interpreting the three O-grid *parts* (a, b, c) as
> layers 1–3, the **new wall is part `d` — the "4th layer", added as the outermost band**, and
> its thickness is the controllable quantity. The plan exposes **two** wall parameters —
> `wall_thickness` and `num_d` (number of wall sub-layers) — which covers both "add the layer"
> and "control its thickness" regardless of how the parts are counted. This interpretation is
> stated up front so the numbering ambiguity does not affect the implementation.

---

## 2. Why this is tractable — how `mesh_volume` forms cells

The key architectural fact (from [`documentation.md`](documentation.md) §2.1): **a hexahedron
is one O-grid quad of cross-section *i* joined to the same quad of cross-section *i+1*.**
Concretely, `mesh_volume` does (`ArterialTree.py:2101`):

```
cells = [8,  id_sectionA + f_ogrid[:,1:],  id_sectionB + f_ogrid[:,1:]]
```

where `f_ogrid = ogrid_pattern_faces(...)` is the **2-D cross-section quad pattern**, and the
per-section vertices come from `ogrid_pattern_vertices(...)`.

**Therefore: if the 2-D O-grid pattern (vertices + faces) is extended to include the outer wall
band, the 3-D hex extrusion produces the wall cells automatically** — `mesh_volume`'s
section-to-section join is unchanged. The work concentrates in the per-cross-section pattern
functions plus the node-count bookkeeping, not in the extrusion logic.

The radial index `k` in `ogrid_pattern_faces` (`ArterialTree.py:2640`,
`faces=[4, ray_prec[k], ray_ind[k], ray_ind[k+1], ray_prec[k+1]]`) is exactly the radial layer
number (k=0 at the wall surface, increasing inward). Wall cells are simply the new quads added
at the outer end of each ray — trivially identifiable for tagging (see
[`plan_for_nas_tagging.md`](plan_for_nas_tagging.md)).

---

## 3. Two implementation options

### Option 1 — Native O-grid extension (integrated, structured, harder)

Extend the O-grid pattern itself so every cross-section carries the wall band.

**3.1 Where the outer wall nodes come from.** For each surface node `crsec[i]`, create `num_d`
new nodes going *outward* to `outer_i = crsec[i] + wall_thickness · n̂_i`, where `n̂_i` is the
outward surface normal. On a plain vessel cross-section the ring lies in the plane normal to the
spline tangent and `crsec[i] − center` already points radially outward, so
`n̂_i = (crsec[i] − center)/‖·‖` gives a **uniform-thickness** wall — directly satisfying
"control the thickness". (Anatomically-faithful alternative: ray-cast onto the existing
`_outer.swc` surface — see Option 1 variant below.)

**3.2 `ogrid_pattern_vertices` (`ArterialTree.py:2515`).** For each ray, **prepend** the wall
nodes before the current surface-first ray:
```
ray = lin_interp(outer_i, crsec[i], num_d+1)[:-1]   # NEW: num_d nodes outside the surface
      + lin_interp(crsec[i], pb, num_a+2)[:-1]        # existing part a
      + lin_interp(pb, square, num_b+2)[:-1]          # existing part b
      + …                                             # existing core
```
Node count per ray grows by `num_d`, so the per-section node count becomes
`nb_nds_ogrid = N·(num_a + num_b + num_d + 3) + ((N−4)/4)²` — the same shape as today with
`num_d` slotting in exactly like `num_a`/`num_b` (strong evidence the original design
anticipated extra radial bands).

**3.3 `ogrid_pattern_faces` (`ArterialTree.py:2591`).** Because each ray gains `num_d` nodes at
the outer end, the existing `k`-loop emits `num_d` extra quads per ray-pair automatically; total
per-section faces grows by `N·num_d`. Capture the new radial region for tagging here.

**3.4 Propagate the new count everywhere `nb_nds_ogrid` is hard-coded.** The formula
`N·(num_a+num_b+3)+((N−4)/4)²` (and the bifurcation form
`nb_nds_bif_ogrid`, `ArterialTree.py:1791`, `2114`, `2374`) appears in **all** of:
`__add_node_id_volume` (`:1772`), `mesh_volume` (`:2039`, `:2114`),
`bif_ogrid_pattern_vertices` (`:2373`), `bif_ogrid_pattern_faces` (`:2285`),
`close_surface` (`:2810`). Every one must take `num_d` and use the extended formula. **This
fan-out is the main error surface** — the counts must stay consistent or the flat-array packing
in `mesh_volume` corrupts silently.

**3.5 Bifurcations (the hard part).** `bif_ogrid_pattern_vertices/faces` (`ArterialTree.py:2350`,
`2268`) build **half-sections with shared edges** and are reordered by `__reorder_faces`
(`:2450`). The wall band must be prepended there too, the shared-edge bookkeeping and
`nb_nds_bif_ogrid` updated, and — critically — the **outward normal at the concave apex** can
make a uniform normal offset self-intersect for a thick wall. Recommended: at bifurcations,
define the outer nodes by **ray-casting onto the existing `_outer` surface**
(reuse `deform_surface_to_mesh`, `ArterialTree.py:3250`, which the STL notebook already uses with
`search_dist ≈ 3`) rather than pure normal offset. Stage this **after** vessels work.

**Signature change:** `mesh_volume(layer_ratio=[a,b,c], num_a, num_b, num_d=0, wall_thickness=0.0, …)`
with `num_d=0` reproducing today's mesh exactly (safe default / backward compatible). Surface the
two new parameters in `Editor.update_mesh_parameters` (`Editor.py:2877`) next to `num_a`/`num_b`.

**Pros:** fully structured/conformal by construction; wall inherits the O-grid quality; one code
path; radial tag falls out for free. **Cons:** touches the most delicate, count-sensitive code;
bifurcation half-section logic is intricate; highest risk.

### Option 2 — Boundary-face extrusion post-process (simpler, lower-risk, recommended first)

Leave the O-grid untouched. After `mesh_volume`, **extrude the lumen boundary faces outward** to
build the wall as stacked hexahedra, then append them to `_volume_mesh`.

**Steps:**
1. Build the lumen volume mesh as today (`mesh_volume`).
2. Extract its **lateral** boundary faces in volume-node IDs (the no-neighbour-face walk already
   in `Simulation.write_mesh_files`, `Simulation.py:247`), excluding inlet/outlet caps
   (classified via `Simulation.__boundary_patches` / `get_inlet_outlet`). These faces are exactly
   the lumen surface, and their nodes are shared with the outermost lumen hexes.
3. Compute a consistent outward normal per boundary node (reuse the notebook's
   `compute_normals(auto_orient_normals=True, consistent_normals=True)` + signed-volume flip from
   `bifurcations_to_stl.ipynb`).
4. For each boundary node, create `num_d` offset copies outward at
   `node + (j/num_d)·wall_thickness·n̂`, `j = 1..num_d`.
5. For each boundary quad, emit `num_d` stacked hexahedra connecting layer `j` to `j+1` (inner
   layer 0 = the existing lumen surface nodes → **conformal / node-matched with the lumen mesh,
   no gaps**). Append vertices + cells to `_volume_mesh`; tag them `domain = wall`.

**Pros:** does not touch the O-grid or bifurcation code at all; works **uniformly at
bifurcations** because it operates on the final surface (the hub is just more boundary quads);
reuses code that already exists and is tested (boundary extraction, normal orientation);
thickness and `num_d` are natural parameters. **Cons:** wall nodes are new points welded onto the
lumen surface (not part of the O-grid's own radial grading); at sharp concavities a large
`wall_thickness` can pinch — mitigate by capping thickness or using the `_outer` surface as the
extrusion target (via `deform_surface_to_mesh`) instead of a fixed normal offset.

**Recommendation:** implement **Option 2 first** — it delivers a tagged fluid+wall FSI mesh
quickly and robustly for the user's bifurcation dataset, and it dovetails with the existing
inner/outer workflow. Treat **Option 1** as the longer-term, fully-integrated solution once
Option 2 validates the downstream COMSOL FSI pipeline.

---

## 4. Controlling the thickness (the explicit requirement)

Expose two parameters wherever the wall is built (both options):

- **`wall_thickness`** — absolute distance (e.g. `0.5`, matching this dataset's uniform inner→
  outer radius gap), *or* a multiple of local radius `r` for a proportional wall. Absolute is
  recommended since the BraVa data uses a constant 0.5.
- **`num_d`** — number of radial cell layers inside the wall (analogue of `num_a`/`num_b`),
  controlling through-thickness resolution for the solid-mechanics solve.

Optional refinement: a radial grading ratio for `num_d` layers (finer near the interface), same
idea as `layer_ratio` for the fluid side.

An alternative that reuses the dataset directly: instead of a numeric thickness, take the wall
from each bifurcation's **`_outer.swc`** (the outer surface already models the true wall) and set
each wall node by ray-casting the lumen boundary onto that outer surface — giving an
anatomically-varying thickness. Provide both modes (`wall_thickness=<float>` for uniform, or
`wall_target=<outer surface>` for anatomical).

---

## 5. Tagging hook (ties to the other plan)

Once the wall cells exist, add a **`domain`** code (fluid=0 / wall=1) to the per-cell record
(the widened `link_graph` from [`plan_for_nas_tagging.md`](plan_for_nas_tagging.md) §4.1):

- Option 1: `domain = wall` when the radial region is the new part `d` (from the `k` index).
- Option 2: the appended extrusion cells are all `domain = wall`; the original cells stay `fluid`.

The **fluid–wall interface** (needed by COMSOL FSI) is the set of faces shared by a fluid cell
and a wall cell — i.e. the *old* lumen boundary faces. They are found in the same boundary-face
walk (a face whose two cells have different `domain`). The NAS writer then emits them as a
dedicated `CQUAD4` boundary group. This is precisely the interface COMSOL couples across.

---

## 6. Staged delivery

1. **Vessels, Option 2.** Boundary-face extrusion on straight/curved vessels; `wall_thickness`
   + `num_d`; append + tag `wall`. Validate node-matching (no duplicate/cracked interface).
2. **Bifurcations, Option 2.** Same post-process over the whole network including hubs; validate
   watertightness and thickness at the apex (reuse the notebook's validation:
   `select_enclosed_points`, signed volume, feature-edge/manifold checks).
3. **Tag + export.** Add the `domain` column; emit fluid/wall/interface via the tagged `.nas`
   writer; import into COMSOL and confirm two domains + interface + boundaries.
4. **(Optional) Option 1.** Native O-grid wall band for a fully-structured, single-code-path
   mesh: vessels first (extend `ogrid_pattern_*` + all `nb_nds_ogrid` sites), then the
   bifurcation half-sections.

---

## 7. Risks / what to check

- **Count-formula fan-out (Option 1).** `nb_nds_ogrid` / `nb_nds_bif_ogrid` are recomputed
  independently in ≥6 places (`ArterialTree.py:1783, 1791, 2039, 2114, 2285, 2373, 2810`). A
  single missed site silently mis-packs the flat vertex/cell arrays. Centralise the formula in
  one helper and call it everywhere.
- **Conformity at the interface.** The wall's inner layer must reuse the **exact** lumen boundary
  node IDs (not new coincident points) so fluid and wall share faces. Option 2 guarantees this by
  extruding the actual boundary faces; verify no duplicated nodes at the interface (a `.clean()`
  or an explicit ID check).
- **Self-intersection of the offset.** At concave regions (bifurcation apex, tight curves) a
  fixed normal offset can fold for large `wall_thickness`. Cap thickness to a small multiple of
  local radius, or extrude onto the `_outer` surface via `deform_surface_to_mesh`
  (`search_dist` limited as the STL notebook found necessary — default 40 lets rays hit the wrong
  surface; ~3 is safe for a 0.5 wall).
- **Normal orientation.** vascularmd's surface comes out inward-facing by convention (documented
  in `bifurcations_to_stl.ipynb`); reuse `fix_orientation` / `auto_orient_normals` so the wall
  extrudes *outward*, not into the lumen.
- **Inlet/outlet handling.** Extrude only lateral wall faces; the wall's own inlet/outlet rims
  become new boundary groups. Don't extrude the O-grid end caps from `close_surface`.
- **Mesh quality.** Check the new wall hexes with the existing `check_mesh` / `utils.quality`
  (`scaled_jacobian`) — thin near-wall cells and apex cells are the likely low-quality spots.
- **`check_mesh` overwrites `_volume_mesh`.** As in the NAS notebook, capture/re-mesh before
  export.
- **Backward compatibility.** Default `num_d=0` / `wall_thickness=0` must reproduce today's mesh
  byte-for-byte so existing CFD-only workflows and `main.py` are unaffected.

---

## 8. Summary

- The requested change is to add a **new outermost radial band of hexahedra (the vessel wall,
  part `d`)** outside the current lumen surface, with **controllable `wall_thickness` and
  `num_d` sub-layers** — turning the single-domain lumen mesh into a conformal fluid+wall FSI
  mesh.
- It is tractable because `mesh_volume` builds hexes by extruding the 2-D O-grid pattern between
  cross-sections: extend the pattern (Option 1) *or* extrude the finished boundary faces
  (Option 2) and the wall cells follow.
- **Option 2 (boundary-face extrusion) is recommended first**: low risk, reuses tested code
  (boundary extraction, normal orientation, deform-to-outer), and works at bifurcations
  uniformly. **Option 1 (native O-grid band)** is the fully-integrated long-term form.
- The wall's radial region feeds the `domain` (fluid/wall) tag and the fluid–wall interface,
  which [`plan_for_nas_tagging.md`](plan_for_nas_tagging.md) then writes into the `.nas` for
  COMSOL — together replacing the current subtract-two-surfaces approach with a single tagged,
  structured hex mesh.
