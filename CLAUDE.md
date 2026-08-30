# CLAUDE.md

Context dump for picking this project back up after a context reset. Read this
first, before touching any script.

## What this project is

Single-cell (fixed RNA profiling / 10x Flex, probe-based) RNA-seq analysis of
ALS across multiple tissue compartments, heading toward a manuscript's
figures. Will eventually be merged with sibling repos for Visium and HMIF
analyses of the same cohort.

**Cohort**: ~30 human donors (Control, sALS, C9orf72-ALS/"C9ORF72-ALS"),
collected across 6 processing batches/pools, from 3 sites (WashU = "AU-",
Barrow = "GWF", Georgetown = "GBB" sample-ID prefixes). 3 tissues per donor
(mostly): motor cortex ("_b" suffix), cervical spinal cord ("_s"), skeletal
muscle ("_m") — ~90 samples total. Sample IDs look like `AU-072_b`,
`GWF19-47_s`, `GBB-23-11_m`. Cell barcodes throughout the pipeline are
`<sample_id>_<10x barcode>`, unique across the whole cohort.

CellBender was run upstream (ambient RNA removal) before this repo's pipeline
starts; `00_cellbender_plotting.R` summarizes those metrics.

**HPC**: Northwestern Quest cluster. Project root:
`/projects/b1169/boles/als_cns_scrnaseq` (every script does
`setwd()` there first). SLURM allocations have moved around between
`b1169`/`b1169` and `b1042`/`genomics` over the course of this project —
check the most recent `jobs/*.sh` for whichever is current, don't assume.
`module load R/4.4.0` + `module load hdf5/1.14.1-2-gcc-12.3.0` in every job
script.

## Rules for working on this repo

- **Never push to `main` directly.** Always: create a feature branch, commit,
  push, open a PR. The user reviews and merges on GitHub. This has held for
  every single change across ~50 PRs — no exceptions, even for "small" fixes.
- **One task = one branch = one PR.** Don't bundle unrelated changes.
- Match the established code style exactly: `message2()` header comments
  before major sections, `setwd()` + relative paths (never hardcoded old
  project paths), `suppressMessages({ library(...) })` blocks, no
  `Sys.time()`/timestamp tracking, comments explain *why* not *what*.
- **Convert per-sample/per-tissue/per-cell-type loops into SLURM job
  arrays** whenever the units of work are independent — this has been the
  standing default since script 04. Prefer a flat `jobs/*.txt` params file
  (comma-separated, no header) over a hardcoded in-script list when the
  targets came from the user as pushed data; use a hardcoded R list when the
  targets are heterogeneous (e.g. some single cell types, some multi-type
  groups — see `15_subclustering2.R`) or when the user described them inline
  in a message rather than as a file.
  - Every array script reads `SLURM_ARRAY_TASK_ID`, `stop()`s with a clear
    message if it's unset (not meant to run outside an array) or out of
    range, and cross-references it against a small lookup table/list. This
    fail-fast pattern is load-bearing — copy it exactly for new scripts.
- **When uncertain, ask or flag — don't guess silently**, especially for:
  scope/column-name ambiguities, whether to redo expensive compute, and
  anything that would be costly to get wrong across many parallel tasks.
  This has come up repeatedly and the user consistently prefers being
  asked. Conversely, don't ask about things confidently resolvable from
  established precedent in this exact codebase — over-asking is also a
  failure mode.
- **Diagnose before guessing a fix**, especially for anything that smells
  like a package/environment/version issue. Add diagnostics and ask the user
  to rerun rather than blind-guessing a third fix. This project has burned
  real time on confident-but-wrong fixes (see Known gotchas below) — when a
  fix is genuinely unverifiable (no R environment available in this session,
  see below), say so explicitly rather than implying more confidence than
  warranted.
- **No R/Seurat/BPCells runtime is available in this session** — it's a
  generic cloud dev container for git/code editing, not the HPC cluster. All
  work here is static code reading/reasoning against docs, stack traces the
  user pastes in, and this codebase's own established patterns. Every fix is
  an informed hypothesis until the user runs it. Say this plainly when
  proposing anything non-trivial, especially version-dependent Seurat/BPCells
  internals.
- Don't add scope beyond what's asked (no speculative abstractions, no extra
  QC steps, no unrequested cleanup) — but do surface genuine bugs found along
  the way (the user has consistently welcomed this, e.g. the `factor()`
  level-swap bug, undefined variables, stale path references).

## Directory conventions

- `data/<script_name>/<tissue>/[<cell_type>/]` — BPCells matrix dirs +
  `.rds` metadata/reductions. Almost always `write_matrix_dir()` for a
  matrix (named `bpcells_data`, or `bpcells` for the two `_obj_reassembly`
  scripts) + separate `saveRDS()` calls for metadata/reductions — essentially
  never a whole-object `saveRDS()`.
- `plots/<script_name>/<tissue>/[<cell_type>/]` — PNGs, `dpi` varies by plot
  density, `bg = "white"` on anything that can render transparent.
- `tab_data/<script_name>/<tissue>/[<cell_type>/]` — CSVs (markers,
  annotation templates, modularity tables).
- `jobs/<script_name>.sh` — SLURM submission scripts, one per R script.
  `jobs/<N>_params.txt` for flat per-task parameter files.
- `.gitignore` only covers `data/`, `plots/`, `tab_data/`, `logs/` at the top
  level — never create a new top-level output folder outside these four
  without adding a gitignore rule, or its contents will get committed.
- A couple of scripts deviate slightly from `<datatype>/<script>/<tissue>/`
  (e.g. `15_subclustering2.R`'s round-1/round-2 label-comparison plots stay
  at the flat `plots/15_annotation2/` from before that script was renamed) —
  check the actual script before assuming the pattern is perfectly uniform.

## Core BPCells/Seurat gotcha: what's actually in the "counts" layer

`CreateSeuratObject()` requires a `counts` argument, but from
`09_*_integration.R` onward, most scripts only have *normalized* data on
disk (the previous step's `$data` layer) — so the established pattern is:

```r
obj <- CreateSeuratObject(counts = data_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat   # same matrix, both slots
```

**The "counts" layer in these objects is not real counts.** This is
harmless for anything that only reads `$data` (ScaleData, Harmony,
FindAllMarkers, DimPlot, etc.), but `FindVariableFeatures()`'s default
`"vst"` selection method needs genuine raw counts to fit correctly. Whenever
a script needs to run `FindVariableFeatures()` fresh (i.e. anything doing its
own re-normalization from scratch, not just reusing an already-selected
variable-feature set), it pulls real raw counts from
`data/06_obj_reassembly/bpcells` (uint32_t, true CellRanger/probe counts,
whole cohort, subset to whatever cells are needed by barcode) instead of
reusing whatever's already loaded. This pattern first appeared in
`13_subclustering1.R` and has been repeated in `15_subclustering2.R`,
`17_obj_reassembly.R`, and `18_full_integration.R` — treat it as the
standard whenever a new script needs fresh variable features/scaling.

## Established naming conventions (reductions, graphs, resolutions)

- PCA reduction: `"pca"`. Harmony-integrated reduction: `"harmony"`. UMAP:
  `"harmony_umap"` (reduction.key auto-derived). Neighbor object (from
  `FindNeighbors(..., return.neighbor = T)`): `"RNA.nn"`. Graphs (from
  `FindNeighbors(..., compute.SNN = T)`): `"RNA_nn"`/`"RNA_snn"`.
- The neighbor graph/UMAP block is called **twice** — once with
  `return.neighbor = T` (Neighbor object, needed by `neighbor_purity()` and
  as `RunUMAP()`'s `nn.name`), once with `compute.SNN = T` (Graph object,
  needed by clustering and `graph_modularity()`). This exact 2-call block
  (see `10_clustering.R`) is copy-pasted verbatim into every script that
  needs a fresh neighbor graph/UMAP (`13`, `15`, `17`'s UMAP-only variant,
  `18`).
- Leiden clustering: `FindClusters(algorithm = 4, method = "igraph")` —
  needs Python `leidenalg` via `reticulate` (`py_module_available()` check,
  currently commented out in `10_clustering.R` after the user confirmed it
  wasn't needed at runtime).
- Cluster-quality resolution sweeps use one of two established lists:
  `10_clustering.R`'s (19 values, `0.2` to `12`) for full-tissue clustering,
  or `13_subclustering1.R`'s (13 values, `0.2` to `5`) for cell-type/group
  subclustering. When a new script says "the same resolutions as before"
  without specifying, infer from which script it's most directly modeled on
  (this has come up and been resolved by inference + explicit flagging, not
  by asking, when only one precedent plausibly matched).
- `RunPCA(npcs = ...)` scales with population size, not fixed project-wide:
  100 for full-tissue objects (`07_norm_pca.R`, `17_obj_reassembly.R`,
  `18_full_integration.R`), 50 for cell-type/group subsets (`13`, `15`).
  `dims = 1:20` for Harmony/neighbors/UMAP is constant everywhere regardless
  of npcs.
- `message2()` is redefined locally in every script (not sourced from a
  shared file) — this project has no shared utils file; small helper
  functions (`graph_modularity()`, `approx_silhouette()`, `neighbor_purity()`)
  are likewise copy-pasted into each script that needs them, not centralized.

## Multi-round annotation system

Cell-type labels accumulate across 3 rounds, each a new metadata column,
never overwriting the previous:

- **`cell_type1`** — round 1, coarse, one label per Leiden cluster at a
  fixed per-tissue resolution (brain=1, sc=1, muscle=1.2 — this exact
  mapping recurs everywhere as a `tissues` data.frame with a `resolution`
  column). Built from `tab_data/12_annotation1/<tissue>_annotations.csv`
  via scCustomize's `Pull_Cluster_Annotation()` → `Rename_Clusters()`.
- **`cell_type2`** — round 2, refines specific `cell_type1` labels by
  subclustering each one individually (`13_subclustering1.R`, one SLURM task
  per (tissue, cell type) pair from `jobs/13_params.txt`) and re-annotating
  (`14_findmarkers2.R`). Folded onto the full object cell-by-cell, matched
  directly by cluster *value* (not scCustomize's `Rename_Clusters()`, which
  relies on positional alignment to sorted factor levels — matching by value
  avoids a whole class of off-by-one risk).
- **`cell_type3`** — round 3, re-integrates/re-clusters selected
  `cell_type2` groups — some single cell types, some deliberately-combined
  multi-type groups to double-check separation between related types (see
  the hardcoded `subclustering_targets` list, reused verbatim across `15`,
  `16`, `17`) — via `15_subclustering2.R`, then folds the result back with a
  **blank-row fallback**: an unfilled `annotations.csv` cluster row (the
  user's way of saying "the existing label already looked right") keeps that
  cell's `cell_type2` value instead of going empty. Cell types never chosen
  as a `subclustering_targets` entry also just keep their `cell_type2` value
  untouched. `cell_type3 == "Remove"` cells are dropped entirely in
  `17_obj_reassembly.R`, which is also where the final per-tissue,
  fully-annotated, cleaned object gets produced (fresh
  Normalize/FindVariableFeatures/ScaleData/RunPCA/Harmony from real counts,
  no reclustering — `cell_type3` is the terminal annotation).

The annotation CSV schema is always `cluster,cell_type` (scCustomize's
`Create_Cluster_Annotation_File()` default column names — never overridden).

## Pipeline stage-by-stage

| Script | Purpose | Array | Notes |
|---|---|---|---|
| `00_cellbender_plotting.R` | Summarize CellBender ambient-RNA metrics | — | Reads raw CellBender output dirs directly, not part of the main chain |
| `01_obj_creation.R` | Per-sample Seurat objects from CellRanger/probe output, saved as per-sample BPCells | array (per sample) | See BPCells segfault saga below |
| `02_qc1.R` | Pre-filter QC stats, builds `tissue`/`batch`/`group`/`id` metadata | — | See factor-level-swap and orig.ident bugs below |
| `03_qc2.R` | Applies QC filters, per-sample thresholds | — | Establishes the `tissues` data.frame convention (title + file) reused everywhere after |
| `04_doubletfinder.R` | DoubletFinder per sample (adjusted + unadjusted) | array (per sample, 1-90) | Materializes BPCells→dgCMatrix (DoubletFinder samples with replacement, BPCells' `[` forbids duplicates). Produces `DF.unadj`/`DF.adj`/`pANN*` columns |
| `05_qc3.R` | Pre/post-filter QC summary plots | — | Plotting only |
| `06_obj_reassembly.R` | Merges all 90 per-sample objects into one whole-cohort object | — | Real raw counts (`uint32_t`) live here: `data/06_obj_reassembly/bpcells` + `metadata.rds`. Referenced by every later script that needs true counts |
| `07_norm_pca.R` | Per-tissue Normalize/FindVariableFeatures/ScaleData/RunPCA(npcs=100)/JackStraw | array (3, per tissue) | Merged/replaced the old sketch-based 07+08 |
| `08_pca_evaluation.R`, `08_sketch_pca.R` | PCA diagnostics; sketch-based PCA (abandoned) | — | Sketch approach abandoned entirely, see below |
| `09_brain/muscle/sc_integration.R` | Per-tissue Harmony integration | 3 separate scripts, not an array | Went CCA→Harmony→(briefly CCA again)→Harmony for good — see CCA saga below |
| `10_clustering.R` | Per-tissue neighbor graph/UMAP/Leiden sweep + cluster-quality scoring | array (3, per tissue) | `graph_modularity()`/`approx_silhouette()`/`neighbor_purity()` defined here first — see silhouette OOM saga below |
| `11_findmarkers.R` | `FindAllMarkers()` at the chosen per-tissue resolution | array (3, per tissue) | `transpose_storage_order()` fix originates here — see below |
| `12_annotation1.R` | Diagnostic plots (cluster UMAP, top-5-marker dot plot, per-gene feature/violin plots) for manual round-1 annotation | array (3, per tissue) | Output: `tab_data/12_annotation1/<tissue>_annotations.csv` |
| `13_subclustering1.R` | Re-clusters each `cell_type1` individually | array (20, `jobs/13_params.txt`) | First script to pull real raw counts from 06 — see counts gotcha above |
| `14_findmarkers2.R` | `FindAllMarkers()` + top-5 dot plot + diagnostic plots per `13` subcluster | array (20, same params file) | Resolution picked by max graph modularity (established as the reliable criterion — see below) |
| `15_subclustering2.R` (renamed from `15_annotation2.R`) | Folds `cell_type1`→`cell_type2` onto full object; re-clusters selected `cell_type2` groups (some multi-type) into `cell_type3`'s source data | array (13, hardcoded `subclustering_targets`) | The file was renamed mid-development as its scope grew from "just fold annotations" to "also re-subcluster" |
| `16_subclustering2_inspection.R` | User's interactive/manual annotation session for `15`'s output | — | Not automated; not part of the array chain, but its `subclustering_targets` list and `annotations.csv` schema are load-bearing for `17` |
| `17_obj_reassembly.R` | Final per-tissue object: folds `cell_type3`, drops `"Remove"` cells, cleans metadata, fresh Normalize/PCA(npcs=100)/Harmony/UMAP (no reclustering) | array (3, per tissue) | Terminal per-tissue output |
| `18_full_integration.R` | Cross-tissue integration on top of `17`'s output | array (2: all 3 tissues, and brain+sc only) | Ported from an earlier version of the project; ~cleaned up for BPCells/Harmony/cell_type3 |
| `demographics_figure.R` | Cohort demographics heatmap for the manuscript | — | Standalone, not part of the processing chain |

**Not yet built** (as of this writing): anything past `18_full_integration.R`
— hdWGCNA, MiloR, or other downstream analyses the user has mentioned as
eventual uses for `18`'s output, but hasn't asked for scripts for yet.

## Known gotchas / troubleshooting sagas (don't re-diagnose these from scratch)

- **BPCells `write_matrix_dir()` segfault** (`01_obj_creation.R` era): took
  several wrong hypotheses (non-integer values, a `gcc` module conflict)
  before diagnosing via package-version/build-date diagnostics that
  BPCells/SeuratObject were compiled under a stale R 4.2.3 while Matrix was
  under R 4.4.0 — an ABI mismatch, fixed by package reinstall (user-side,
  not a code fix). Lesson institutionalized above: diagnose before guessing.
- **`factor(x, labels = ...)` without explicit `levels =`** (`02_qc1.R`):
  silently sorted alphabetically, swapping the Motor cortex/Cervical spinal
  cord tissue labels. Always pass explicit `levels =` with `labels =`.
- **`orig.ident` overwritten, losing the tissue suffix** (`02_qc1.R`): fixed
  by adding a separate `id` column instead of mutating `orig.ident` in place.
- **DoubletFinder + BPCells**: `paramSweep()` samples cell pairs with
  replacement; BPCells' `[` rejects duplicate column selection. Fixed with
  `as(mat, "dgCMatrix")` before handing off to DoubletFinder.
- **Sketch-based workflow abandoned**: `SketchData()` prefixes sketch-assay
  cell names with a leading `_` to disambiguate them from the full-data
  assay's cells coexisting in the same object — this caused a fatal barcode
  mismatch in the old `10_*_integration.R` (CCA-on-sketch) scripts that was
  never fully fixed, because the user separately found (from earlier-project
  notes) that sketched integration gave visibly worse results than full-data
  integration and abandoned sketching entirely. `07`/`08` were rewritten to
  drop sketching before this was resolved — don't resurrect the sketch path
  without reason.
- **CCA vs. Harmony, round 1**: switched `09_*_integration.R` to CCA per
  user request; hit three sequential bugs before it ran at all —
  (1) a `reference=` computation using `Layers(obj, search = "data")`
  turned out to be a red herring (the crash persisted after "fixing" it);
  (2) the real cause was a Seurat/BPCells bug where `IntegrateLayers()`
  can't cleanly densify a BPCells-backed *split* assay internally (matches
  `satijalab/seurat#8004`/`#7113`) — fixed by explicit
  `as(obj[["RNA"]]$data, "dgCMatrix")` before `split()`; (3) CCA's internal
  legacy-Assay step then needed real `scale.data` (dense, not the sparse
  matrix it tried to synthesize), fixed by an explicit `ScaleData()` call.
  Only after all three fixes did CCA actually run and get fairly compared
  against Harmony.
- **CCA vs. Harmony, round 2**: once comparable, CCA and Harmony gave
  similar integration quality, and Harmony was meaningfully faster (PCA-
  embedding-only vs. CCA's dense expression access) — `09_*_integration.R`
  reverted to Harmony for good. If asked to reconsider CCA again, this
  history is why Harmony is the default.
- **BPCells column-major vs. row-major storage**: `09_integration`'s
  `bpcells_data` is stored column-major (cell-major, right for PCA/
  clustering), but `FindAllMarkers()` needs row-major (per-gene) access.
  Without a pre-transpose, Seurat re-transposes small chunks live during the
  test loop (very slow, with an explicit warning naming the fix). Fixed by
  `BPCells::transpose_storage_order()` once, up front, writing to a scratch
  temp dir (not back to the shared `bpcells_data`) — first done in
  `11_findmarkers.R`, then `14_findmarkers2.R`, `15_subclustering2.R`.
  Do this in any new script that runs `FindAllMarkers()` on loaded BPCells
  data.
- **`cluster::silhouette()` OOM**: needs a full cell × cell distance matrix,
  O(cells²) — OOM-killed `10_clustering.R` on the larger tissues. Replaced
  with 3 approximations that never build that matrix:
  `approx_silhouette()` (centroid-based, O(cells × clusters)),
  `neighbor_purity()` (reuses the already-computed Neighbor object),
  `graph_modularity()` (reuses the already-computed sparse SNN graph via
  `igraph`). In practice, silhouette and neighbor purity turned out to be
  **inversely correlated with resolution** (they reward fewer, coarser
  clusters almost by construction — real biology is usually a continuum,
  and coarse partitions look "cleaner" by pure geometric separation even
  when finer real structure exists). **Graph modularity is the metric
  actually used to pick resolutions** (`14`, `15` both pick
  `resolution = res_tests[which.max(modularity)]`) — it measures whether a
  partition explains the graph's actual edges, not just embedding-space
  compactness, and empirically produces a genuine peak rather than a
  monotonic curve.
- **scCustomize's `Pull_Cluster_Annotation()`/`Rename_Clusters()` pairing**:
  the standard way to apply an `annotations.csv` (default columns `cluster`,
  `cell_type`) to a Seurat object's `Idents()`. `Rename_Clusters()` matches
  by *position* against `levels(Idents(object))` — safe when used exactly
  as documented (set `Idents()` to the right resolution column first), but
  risky to reimplement by hand. When folding annotations onto a *different*
  object than the one the CSV was built from (e.g. round-2/round-3 folding
  onto the full tissue object), this project matches by cluster *value*
  directly instead (`match(cluster_col, annot_csv$cluster)`), not by
  position — deliberately avoiding `Rename_Clusters()`'s positional
  assumption in that specific cross-object case.

## Where things actually stand right now

`18_full_integration.R` is the most recently added script (cross-tissue
integration on top of the fully-annotated per-tissue objects from `17`).
Nothing downstream of it has been built yet. The user has mentioned hdWGCNA
and MiloR as intended downstream uses of `18`'s output but hasn't requested
scripts for either. If picking this project back up, ask what's next rather
than assuming — the pipeline's shape has changed direction more than once
already (sketch→full-data, CCA→Harmony→(briefly CCA)→Harmony), and it's
cheap to check before building on an assumption.
