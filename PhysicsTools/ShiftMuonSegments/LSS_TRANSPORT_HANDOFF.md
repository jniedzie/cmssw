# LSS transport handoff

Last updated: 2026-08-28

## Repository state

- CMSSW checkout: `/afs/cern.ch/work/j/jniedzie/private/shift_cmssw/CMSSW_17_0_0_pre4/src`
- Development branch: `shift-muon-segments-pre4`
- Remote tracking branch: `my-cmssw/shift-muon-segments-pre4`
- Current commit: `3061ca590ae` (`Add FLUKA 2D field map support and IR1 Atlas proxy configuration`)
- Preceding LSS commit: `b26eb705a28` (`Add ShiftLssMagneticFieldESProducer and related configurations for LSS transport`)
- The checkout was clean and synchronized with the remote when this note was written.

The workflow repository remains on `main`. Its relevant committed IR1 proxy
fixture and extraction/conversion work is under:

```text
shift_cmssw_workflow/models/lss5_ir1_atlas_proxy
```

## Implemented CMSSW pieces

1. `ShiftLssMagneticFieldESProducer` wraps a separately labelled standard CMS
   field and publishes a configurable composite field. Simulation and SHIFT
   reconstruction can therefore consume the same empty-label EventSetup field.
2. Field elements support uniform fields, analytic quadrupoles, and the FLUKA
   `QUAD`, `QUADINT`, `INTER2D`, and `KICKINT` 2D map formats. The implementation
   reproduces the FLUKA symmetry, analytic-core, row ordering, interpolation,
   scaling, and local/global rotation behavior.
3. Spatial bounds support boxes, cylinders, annuli, and excluded cylindrical
   holes. This represents all 29 active IR1 deck assignments: 17 map-backed
   cells and 12 constant-field cells.
4. `shiftLssIr1AtlasProxyFieldElements` constructs those assignments from the
   frozen deck. It requires an explicit `modelOriginCm` and proper
   `modelToCms` rotation. There is intentionally no default IR1-to-IR5 mapping.
   Mapped and constant field families can be toggled independently.
5. The five referenced maps are packaged below
   `PhysicsTools/ShiftMuonSegments/data/lss/ir1_atlas_proxy` and verified against
   `FIELD_MAP_SHA256SUMS`.
6. `customiseShiftLssTransport` configures the shared field label, material
   boundary, and Geant4e momentum, step, and path limits. The remaining active
   `1100 cm` target-constraint boundary was replaced by the configured
   `lssTransport.materialBoundaryAbsZCm` value.

## Important limitations

- This is a legacy 6.5 TeV Beam-2 IR1/ATLAS LineBuilder model. It is not an
  IR5/CMS model and not an authoritative Run-3 model. Use it only as a
  provisional sensitivity proxy.
- The magnetic field is implemented, but the converted tunnel/magnet/material
  geometry is not yet installed in CMSSW. Field-only propagation does not model
  momentum loss or multiple scattering.
- A FLUKA-to-GDML conversion preserves geometry and materials only. It does not
  carry field assignments, FLUKA physics, source, biasing, or scoring. The
  separately implemented field must be validated against the converted
  placements.
- No IR1-proxy-to-CMS coordinate transform has been approved. Do not enable the
  factory in production with a guessed identity or reflection transform.
- Full `cmsRun` validation was not possible on 2026-08-28 because even an
  `EmptySource` baseline configuration segfaulted at startup in both the shared
  and disposable CMSSW environments, independently of the LSS plugin.
- `condor_q` also failed because `wtfis.cern.ch` could not be resolved. Recheck
  the scheduler before editing live campaign configuration or rebuilding the
  shared CMSSW release.

## Validation completed

The focused Python configuration tests passed:

```bash
cd /tmp/jniedzie/shift_lss_build_20260828/src
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval "$(scram runtime -sh)"
python3 PhysicsTools/ShiftMuonSegments/test/test_shift_lss_configuration.py
```

All five tests passed. A disposable build also succeeded and registered
`shiftLssMagneticFieldESProducer`:

```bash
scram b -j 4 PhysicsTools/ShiftMuonSegments
```

`python3 -m py_compile` and `git diff --check` passed as well.

Useful reference values independently calculated from the frozen maps and the
FLUKA interpolation rules are:

```text
MQXA  local (10.25, -3.75), scale +1.8928521203:
      Bx = -2.448081013108509 T, By = -0.515688631654532 T
MBXW  local (10, -10), scale -1.1999831604:
      Bx = -0.1496211003376344 T, By = -1.6437825324760151 T
MQTL  local (-12, 3), scale -0.37950611892:
      Bx = -0.323903695368136 T, By = +0.8455491124744595 T
MQXB  local (8.25, 7.75), scale -1.8928521203:
      Bx = -1.3095434607059622 T, By = -0.9075209544852374 T
```

These should become direct C++ field-sampling regression checks once a stable
EventSetup test can run.

## Restart order

1. Confirm the branch and cleanliness with `git status --short --branch` and
   `git log -3 --oneline --decorate`. Continue on `shift-muon-segments-pre4`.
2. Check `condor_q` before touching the shared release or campaign files. Do
   not rebuild shared AFS CMSSW while jobs use it.
3. Diagnose the baseline `cmsRun` startup segfault with a pristine minimal
   CMSSW configuration. Keep that investigation separate from LSS plugin bugs.
4. Run the workflow GDML conversion using the recorded lattice-AABB workaround.
   Inspect its `conversion_report.json`, material inventory, lattice count,
   world bounds, and overlaps. Do not drop lattice cards to make conversion
   succeed.
5. Add the converted geometry to CMSSW behind an explicit configuration switch.
   Ensure both `g4SimHits` and Geant4e reconstruction see the same Geant4
   material geometry; field configuration alone is insufficient.
6. Establish and document an explicit IR1-model-to-CMS placement transform and
   source side. Validate magnet centers and axes against geometry before using
   the proxy factory.
7. Add direct field sampling at the reference points above and at element
   boundaries, annuli, excluded holes, overlaps, and outside-field locations.
8. Dump the resolved CMSSW configuration and verify that simulation and every
   SHIFT reconstruction propagator consume the intended field label, geometry,
   material boundary, and Geant4e limits.
9. Run a bounded 1-20 event test. Check logs, readable ROOT output, forward
   generator transport, backward reconstruction transport, charge signs,
   momentum closure, material path, and failure counters before considering a
   larger campaign.

Keep electronics integration windows, BX assignment, readout, and trigger
rules unchanged. This work concerns transport geometry, material, and magnetic
fields; it must not be used to loosen Run-3 detector or trigger constraints.
