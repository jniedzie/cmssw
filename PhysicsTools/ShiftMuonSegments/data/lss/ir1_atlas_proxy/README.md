# IR1/ATLAS proxy magnetic-field maps

These five files are exact copies of the active field maps in the frozen
engineering fixture at:

```text
shift_cmssw_workflow/models/lss5_ir1_atlas_proxy/source
```

The source model is the legacy IR1/ATLAS Beam-2 LineBuilder deck, not an IR5
or Run-3 CMS model. The files are packaged here only so CMSSW configurations
can resolve them through `edm::FileInPath` without EOS or checkout-specific
absolute paths. `FIELD_MAP_SHA256SUMS` must continue to match the workflow
source bundle.

The map file alone does not define a magnet placement. Every CMSSW field
element must additionally provide its geometry-derived local bounds, origin,
proper rotation, and FLUKA field scale. Do not infer these from the map grid:
the grid is a transverse interpolation domain, not the physical magnet volume.

`shiftLssIr1AtlasProxy_cff.py` records the 17 active map-backed lattice-cell
placements from the same deck. Its factory requires an explicit model-to-CMS
translation and proper rotation; there is intentionally no default IR1-to-IR5
alignment. Cylindrical and annular cell bounds are retained rather than
replacing them with their enclosing boxes. The factory also records the 12
active constant-field cells, including the two excluded beam holes in the
MCBCH outer cell. Map-backed and constant-field families can be disabled
independently for controlled validation, but both are enabled by default.
