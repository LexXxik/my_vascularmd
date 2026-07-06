# Plan — Tagging *only* the `num_a` boundary layer in the `.nas` export

> **Stage 2 deliverable (planning only — no code is changed here).**
> Focused specialization of [`plan_for_nas_tagging.md`](plan_for_nas_tagging.md); shares its
> tagged-`.nas` writer and region-capture mechanism. Codebase reference:
> [`documentation.md`](documentation.md).
> Goal: give the cells of the **`num_a` boundary layer only** a distinct property id (`PID`)
> in the exported `.nas`, so COMSOL can select the near-wall band as its own domain while every
> other cell stays in a single default domain.

---

## 1. Why tag only the `num_a` layer

`num_a` is the **boundary layer** of the O-grid — the thin band of cells hugging the lumen
surface (`layer_ratio[0]`, `num_a` cells; see [`documentation.md`](documentation.md) and
`ogrid_pattern_vertices`, `ArterialTree.py:2515`). Selecting *only* this band is useful for:

- near-wall CFD post-processing (wall-shear-stress / boundary-layer analysis on a defined zone),
- assigning a distinct material / porous / wall-function treatment to the wall-adjacent cells,
- (for FSI) isolating the fluid layer that couples to the wall interface.

This plan is deliberately narrow: **one distinguished region (`a`) vs. everything else.** It is
a strict subset of the general tagging plan — if that plan's full radial-region tagging is built,
"num_a-only" is just choosing which region gets the non-default `PID`.

---

## 2. Where the `num_a` cells are in the mesh

- A hexahedron = one O-grid quad of cross-section *i* joined to the same quad of *i+1*
  (`mesh_volume`, `ArterialTree.py:2101`), so a cell's radial position **is** its base quad's
  radial position.
- In `ogrid_pattern_faces` (`ArterialTree.py:2591`) each quad is built at a radial index `k`
  (`faces=[4, ray_prec[k], ray_ind[k], ray_ind[k+1], ray_prec[k+1]]`, `ArterialTree.py:2640`),
  where **k = 0 is the outermost cell (at the lumen surface)** and increasing `k` goes inward.
- The ray in `ogrid_pattern_vertices` is built as
  `lin_interp(crsec[i], pb, num_a+2)[:-1]  +  lin_interp(pb, square, num_b+2)[:-1]  +  core`
  (`ArterialTree.py:2571`). The **first segment (surface → `pb`) is exactly part `a`** — so the
  `num_a` band is the **outermost contiguous run of `k` values**, i.e. the cells generated before
  the ray reaches `pb`.

> **Off-by-one caution.** `lin_interp(crsec, pb, num_a+2)[:-1]` yields `num_a+1` nodes before
> `pb`, so the surface→`pb` span is `num_a+1` cells, not `num_a`. Rather than hard-code a
> `k < num_a` threshold (which would be wrong by one), **label part `a` at generation time**,
> where the code knows it is still emitting the first `lin_interp` segment. Confirm the exact
> count directly from `ogrid_pattern_vertices` during implementation.

---

## 3. Design

### 3.1 Mark the part-`a` cells (small mesher change)

Reuse §4.1 of [`plan_for_nas_tagging.md`](plan_for_nas_tagging.md), but only a binary flag is
needed:

- In `ogrid_pattern_faces` / `bif_ogrid_pattern_faces` (`ArterialTree.py:2591`, `2268`): emit a
  parallel `is_layer_a` array (True while the ray is still in the surface→`pb` segment, False
  after). Return `(faces, is_layer_a)`.
- In `mesh_volume` (`ArterialTree.py:2003`): carry the flag into a new `link_graph` column (or a
  parallel `cell_region` array), so every hex inherits its base quad's flag. No extra geometry —
  it is a by-product of the loop that already exists.

### 3.2 Selective `.nas` export

Use the hand-written tagged `.nas` writer from [`plan_for_nas_tagging.md`](plan_for_nas_tagging.md)
§4.2, with a trivial `PID` rule:

- `CHEXA … PID` = **1** for cells with `is_layer_a == True`, **2** for all other cells.
- Emit `PSOLID 1 1` (boundary layer) and `PSOLID 2 1` (rest). COMSOL then shows the boundary
  layer as one selectable domain and everything else as another.
- Everything else about the writer (1-based IDs, `GRID*` large-field format, `CHEXA` `+`
  continuation) is unchanged from the general plan.

### 3.3 Optional — tag the faces bounding the `num_a` band

If the near-wall *surfaces* are wanted (not just the cells), emit coincident `CQUAD4` shells
(same mechanism as the general plan §4.2 step 4) for:

- the **lumen surface** (outer face of the `num_a` band — the `k=0` faces), and
- the **a/b interface** (inner face of the band — the faces at the `part a → part b` transition).

Both are obtainable from the boundary/region-aware face walk; give each its own `PID`.

---

## 4. Fastest path / fallback

- **Quickest, no core change (fallback).** Tag geometrically after meshing: for each cell, take
  its centroid, find the local centerline point via its `link_graph` cross-section index
  (`crsec` center / spline), and mark it part-`a` if its radial distance from the axis exceeds
  `(1 − layer_ratio[0]) · R_local`. Works for a first result on straight/curved vessels; less
  reliable at bifurcation hubs and near the core-square corners.
- **Clean approach (recommended).** The §3.1 generation-time flag — exact and robust, and it
  also unlocks full a/b/c tagging later for free.

---

## 5. Risks / what to verify

- **Off-by-one on the part-`a` range** — see §2 caution; label at generation, don't guess `k`.
- **Bifurcation region codes** — `bif_ogrid_pattern_faces` (`ArterialTree.py:2268`) has more
  intricate per-ray bookkeeping than the plain case; validate the flag visually (colour cells by
  `PID` in ParaView) at a bifurcation before trusting it.
- **`check_mesh` overwrites `_volume_mesh`** (noted in `notebooks/vtk_to_nas.ipynb`) — tag/export
  the mesh captured right after `mesh_volume`, or re-mesh before exporting.
- **Validation** — re-read the `.nas` with meshio / import to COMSOL and confirm `PID 1` contains
  exactly the boundary-layer cell count (≈ `num_a`-per-ray × rays × cross-section-pairs) and forms
  a single connected shell around the lumen.

---

## 6. Summary

- The `num_a` boundary layer is the **outermost radial band** of the O-grid (the `k=0…` cells
  before the ray reaches `pb`), so it is cleanly separable.
- Add a binary `is_layer_a` flag in the O-grid face builder → `link_graph`, then set a distinct
  `CHEXA` `PID` for those cells in the hand-written tagged `.nas` writer — one selectable
  near-wall domain in COMSOL, everything else default.
- Reuses the mechanism and writer from [`plan_for_nas_tagging.md`](plan_for_nas_tagging.md); a
  geometric centroid-radius test is the no-core-change fallback for a first result.
