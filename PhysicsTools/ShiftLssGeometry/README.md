# SHIFT external LSS geometry

This package adds a bounded LSS material model to the standard CMSSW DD4hep
geometry. It never substitutes a FLUKA detector description for the CMS
detector. The source reconstructs the normal `Extended` conditions payload,
retains every pre-existing CMS world placement, and then attaches one
explicitly transformed external GDML child.

The customization is deliberately off by default. It requires a packaged
GDML file, a proper model-to-CMS rotation, an origin, and a positive protected
CMS boundary. Unsupported ROOT 6.36 `multiUnion` input, a volume crossing that
boundary, a changed standard world placement, and newly introduced geometry
overlaps fail closed.

`shift_cmssw_workflow/run_step1_generation.sh` enables the extension only when
`SHIFT_LSS_GEOMETRY_MODE=external`. The required values are:

```bash
SHIFT_LSS_GDML_FILE=PhysicsTools/ShiftLssGeometry/data/validated_lss.gdml
SHIFT_LSS_MODEL_ORIGIN_CM=x,y,z
SHIFT_LSS_MODEL_TO_CMS=r00,r01,r02,r10,r11,r12,r20,r21,r22
SHIFT_LSS_MINIMUM_ABS_Z_CM=2750
```

The GDML must be installed below the CMSSW `src` directory. Production mode
always enables overlap checks. Leaving `SHIFT_LSS_GEOMETRY_MODE` unset keeps
the unmodified `DB:Extended` geometry path.

Build and run the bounded positive and negative runtime checks with:

```bash
scram b -j 4 PhysicsTools/ShiftLssGeometry
edmPluginRefresh "$CMSSW_BASE/lib/$SCRAM_ARCH"
cmsRun PhysicsTools/ShiftLssGeometry/test/run_bounded_external_geometry_cfg.py
cmsRun PhysicsTools/ShiftLssGeometry/test/run_protected_region_rejection_cfg.py
cmsRun PhysicsTools/ShiftLssGeometry/test/run_multi_union_rejection_cfg.py
```

The last two commands must fail with `UnsupportedGeometry`.

The frozen IR1/ATLAS proxy is not a Run-3 IR5 model. The audited converter can
now produce a self-contained, ROOT-compatible bounded GDML that passes its
standalone ROOT overlap, extrusion, gap, and source-containment gates. Full
CMSSW/DD4hep attachment, protected-region, and overlap checks pass at an
artificial +1 km placement. The proxy still has no authoritative IR1-to-IR5/CMS
transform and does not implement the FLUKA magnetic fields. It must therefore
not be enabled for detector simulation or reconstruction. The bounded fixture
is for software validation only, not physics.
