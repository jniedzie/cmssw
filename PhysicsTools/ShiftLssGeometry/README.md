# SHIFT external LSS geometry

This package adds a bounded LSS material model to the standard CMSSW DD4hep
geometry. It never substitutes a FLUKA detector description for the CMS
detector. The source reconstructs the normal `Extended` conditions payload,
retains every pre-existing CMS world placement, and then attaches one
explicitly transformed external GDML child.

The customization is deliberately off by default. It requires a packaged
GDML file, its converter-recorded origin in the source model, a proper
model-to-CMS rotation, a model origin, and a positive protected CMS boundary.
The explicit artifact origin is required because the bounded converter
recentres its GDML; the runtime placement is
`modelOriginInCms + modelToCms * artifactOriginInModel`. Unsupported ROOT 6.36
`multiUnion` input, a volume crossing that
boundary, a changed standard world placement, and newly introduced geometry
overlaps fail closed.

`shift_cmssw_workflow` enables the extension only when
`SHIFT_LSS_MATERIAL_MODE=external`. Step 1 and Step 4 use the same values:

```bash
SHIFT_LSS_GDML_FILE=PhysicsTools/ShiftLssGeometry/data/validated_lss.gdml
SHIFT_LSS_GDML_SHA256=recorded_64_character_digest
SHIFT_LSS_ARTIFACT_ORIGIN_IN_MODEL_CM=converter_x,converter_y,converter_z
SHIFT_LSS_MODEL_ORIGIN_CM=x,y,z
SHIFT_LSS_MODEL_TO_CMS=r00,r01,r02,r10,r11,r12,r20,r21,r22
SHIFT_LSS_MINIMUM_ABS_Z_CM=2750
SHIFT_LSS_MATERIAL_BOUNDARY_ABS_Z_CM=14800
SHIFT_LSS_GEANT4E_MAXIMUM_PATH_LENGTH_CM=20000
```

The GDML must be installed below the CMSSW `src` directory. Production mode
always enables overlap checks. Leaving `SHIFT_LSS_MATERIAL_MODE` unset keeps
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

The frozen IR1/ATLAS proxy is not a Run-3 IR5 model. Its installed bounded GDML
and manifest pass standalone ROOT overlap, extrusion, gap, and source-
containment gates; CMSSW/DD4hep attachment, protected-region, overlap, and
combined field-integration tests also pass. The separately implemented proxy
field and all four material/field switch combinations complete a fixed-seed
one-event chain through NanoAODSIM. This establishes software wiring only. The
proxy has no authoritative IR1-to-IR5/CMS transform, source side, machine
state, or reviewed CMS polarity, so it must not be used for a CMS physics
claim.
