# Summaries & quick-implementation recommendations

Brief summaries of the Stage 2 plans, each with the **fastest path to a working result**.
Full detail lives in [`plan_for_nas_tagging.md`](plan_for_nas_tagging.md),
[`plan_for_nas_tagging_layer_a.md`](plan_for_nas_tagging_layer_a.md) and
[`plan_for_extending_volume_mesh.md`](plan_for_extending_volume_mesh.md); the codebase
reference is [`documentation.md`](documentation.md).

The two plans interlock: tagging provides the *mechanism* to label cells/faces in the
`.nas`; the wall extension provides the *fluid/wall/interface* regions worth labelling.
Tagging can ship first, on its own.

---

## 1. NAS tagging — [`plan_for_nas_tagging.md`](plan_for_nas_tagging.md)

**Summary**

- Today's `.nas` (via `pv.save_meshio`) is **untagged** — 89 232 `CHEXA` with one implicit
  property, no boundaries — so COMSOL can't select any sub-region.
- **Cells** are taggable via `CHEXA` `PID`; **faces** are taggable via coincident `CQUAD4`
  shell elements with their own `PID` (COMSOL reads both as domains / boundary selections).
  So yes, individual faces *can* be tagged.
- The tag data mostly already exists: `link_graph` gives per-cell `[branch, crsec]`;
  `Simulation.__boundary_patches` + the boundary-face walk in `Simulation.write_mesh_files`
  already classify wall / inlet / outlet. Only the **radial region** (a one-line capture of
  the existing `k` index in `ogrid_pattern_faces:2640`) is missing.
- Recommends a **hand-written tagged `.nas`** method (matching the repo's OpenFOAM/vtk/swc
  writers) for full `PID` control that `pv.save_meshio` doesn't give.

**Fastest path to a result**

1. **~½ day — per-branch domains, zero mesher changes.** Write a small `.nas` function reusing
   `_volume_mesh` (cells + vertices) and the **existing** `link_graph` branch columns as the
   `CHEXA` `PID`. COMSOL immediately gets one selectable domain per vessel/bifurcation.
2. **+~1 day — radial + boundary tags.** Capture `k` in `ogrid_pattern_faces` → widen
   `link_graph` with a radial-region column (enables a/b/c and, later, fluid/wall). Lift the
   no-neighbour-face walk from `Simulation.write_mesh_files` → emit boundary `CQUAD4` with
   inlet/outlet/wall `PID`.
3. Validate by re-reading with meshio / importing to COMSOL and checking each `PID`'s cell/face
   counts. Watch the gotchas already listed: 1-based IDs, `GRID*` large-field format, `CHEXA`
   `+` continuation, and `check_mesh` overwriting `_volume_mesh` before export.

---

## 2. Boundary-layer-only (`num_a`) tagging — [`plan_for_nas_tagging_layer_a.md`](plan_for_nas_tagging_layer_a.md)

**Summary**

- A narrow specialization of the general tagging plan: give **only the `num_a` boundary-layer
  cells** (the near-wall band) a distinct `.nas` `PID`; every other cell stays in one default
  domain. Useful for near-wall WSS / boundary-layer zones, a distinct material, or isolating the
  wall-adjacent fluid for FSI.
- The `num_a` band is the **outermost radial run** of the O-grid (the `k = 0…` cells before the
  ray reaches `pb` in `ogrid_pattern_faces:2640`), so it separates cleanly with a single binary
  `is_layer_a` flag captured in the mesher and carried on `link_graph`.

**Fastest path to a result**

1. **~½ day, no core change (fallback).** Geometric per-cell test: tag part-`a` when a cell's
   centroid radial distance from its local centerline exceeds `(1 − layer_ratio[0]) · R`. Fine for
   straight/curved vessels; weaker at bifurcation hubs.
2. **~1 day, clean.** Add an `is_layer_a` flag in `ogrid_pattern_faces` → `link_graph`, then set a
   distinct `CHEXA` `PID` for those cells in the hand-written tagged `.nas` writer. Exact, and it
   unlocks full a/b/c tagging later.
3. Mind the off-by-one (surface→`pb` is `num_a+1` cells — label at generation, don't hard-code
   `k < num_a`) and validate the tagged band visually + by cell count.

---

## 3. Extending the volume mesh — [`plan_for_extending_volume_mesh.md`](plan_for_extending_volume_mesh.md)

**Summary**

- Today the hex mesh fills the **lumen only** (O-grid parts a + b + core). The goal is a new
  **outermost layer of hexahedra** (the vessel wall) with controllable thickness and sub-layer
  count — turning the single-domain mesh into a conformal fluid + wall mesh for COMSOL FSI,
  replacing the current "model inner & outer surfaces separately and subtract in COMSOL".
- Tractable because a hex = one O-grid quad of section *i* joined to section *i+1*, so extending
  the 2-D cross-section pattern makes the 3-D wall cells follow automatically.
- Two options: **Option 1** — native O-grid extension (prepend an outer band to every ray;
  elegant and fully structured, but must thread a new layer count through **8 count-formula
  sites** and the intricate bifurcation half-section code). **Option 2** — post-process that
  extrudes the finished lumen boundary faces outward into stacked hexes (leaves the O-grid
  untouched, works uniformly at bifurcations, reuses tested boundary-extraction + normal-
  orientation code).

**Fastest path to a result**

1. **~2–3 days — Option 2 (boundary-face extrusion), recommended first.** After `mesh_volume`:
   extract lateral boundary faces (reuse the `Simulation.write_mesh_files` walk), offset each
   boundary node outward by `wall_thickness` along an oriented normal (reuse the notebook's
   `fix_orientation` / `auto_orient_normals`), stack `num_d` hexes per face (inner layer reuses
   the existing surface node IDs → conformal, no gaps), append + tag `wall`.
2. **Validate** with the notebook's existing watertight / signed-volume / manifold checks plus
   `check_mesh` (`scaled_jacobian`). Cap thickness at concave apexes, or extrude onto the
   existing `_outer.swc` surface, to avoid self-intersection.
3. **Option 1 (native band) later** if a fully-structured single-code-path mesh is wanted — do
   vessels first (extend `ogrid_pattern_*` + centralise the count formula), bifurcations after.

---

## 4. Fastest end-to-end (tagged fluid + wall `.nas` for COMSOL FSI)

Combine the two quick paths:

1. **Wall** via Option 2 extrusion (~2–3 days) → cells carry a `domain` = fluid/wall flag.
2. **Tag** via the hand-written `.nas` writer (~1–2 days): `CHEXA` `PID` = fluid vs wall;
   `CQUAD4` groups for inlets, outlets, wall exterior, and the **fluid–wall interface** (faces
   whose two cells differ in `domain`).
3. Result: a single tagged, structured hex mesh COMSOL imports directly as fluid domain + wall
   domain + interface + boundaries — no more subtracting two surfaces. Total ≈ 1 working week,
   reusing existing/tested code throughout.

**Sequencing note.** Ship NAS tagging Stage A (branch/region/boundary) first — it's independent
and immediately useful on the current lumen mesh. The wall + FSI interface tags slot in once the
extrusion lands.
