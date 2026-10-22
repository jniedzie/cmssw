# Code style: AI-code. AI-generated paired nominal/particle-only controls.
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

here = Path(__file__).resolve().parent
configuration = Path(sys.argv[1]).resolve() if len(sys.argv) == 2 else here/"process_inventory_cfg.py"
binary = shutil.which("testShiftDarkPhotonProcessInventory")
if binary is None:
    raise RuntimeError("Build the private package tests before running this control")

# This standalone binary links the Geant4 shared libraries. The bundled
# Simulation plugin also embeds Geant4; loading both creates separate thread
# state. Use the ordinary individual CMS plugins for these geometry-free tests.
runtime = dict(os.environ)
runtime["LD_LIBRARY_PATH"] = ":".join(
    path for path in runtime.get("LD_LIBRARY_PATH", "").split(":") if "/biglib/" not in path
)

with tempfile.TemporaryDirectory(prefix="shift-dp-processes-") as temporary:
    results = {}
    for mode in ("nominal", "noop", "dark-photon"):
        output = Path(temporary)/(mode+".json")
        subprocess.run([binary,str(configuration),mode,str(output)],
                       check=True,timeout=180,env=runtime)
        results[mode] = json.loads(output.read_text())
    nominal, noop, dark = (results[mode] for mode in ("nominal","noop","dark-photon"))
    assert len(nominal["sm_particles"]) > 100
    assert nominal["sm_particles"] == noop["sm_particles"] == dark["sm_particles"]
    assert nominal["signal_decay_process_pdgs"] == noop["signal_decay_process_pdgs"] == []
    assert sorted(dark["signal_decay_process_pdgs"]) == [-321,130,310,321]
    assert all("/biglib/" not in path for result in results.values() for path in result["runtime_library_paths"])
    if nominal["nominal_factory"] == "SimG4Core/Physics/FTFP_BERT_EMH":
        for result in results.values():
            for pdg in (22,11,-11):
                particle = next(row for row in result["sm_particles"] if row["pdg_id"] == pdg)
                assert particle["g4hepem_tracking_manager"]
                assert "G4HepEmTrackingManager" in particle["tracking_manager_class"]
    assert nominal["constructors"] == noop["constructors"] == dark["constructors"][:-1]
    assert dark["constructors"][-1] == {"name":"ShiftDarkPhoton","type":0}
    assert nominal["dark_photon"] is None and noop["dark_photon"] is None
    point = dark["dark_photon"]
    assert point["pdg_id"] == 32 and point["mass_gev"] == 0.3
    assert point["twice_spin"] == 2 and point["charge_e"] == 0.
    assert not point["stable"] and not point["short_lived"]
    assert len(point["processes"]) == 2
    assert sorted(row["name"] for row in point["processes"]) == ["Decay","Transportation"]
    for field in ("massGeV", "properCTauMm"):
        invalid = Path(temporary)/(field+"_zero.py")
        invalid.write_text(configuration.read_text()
                           + f"\nphysics.DarkPhoton.{field} = cms.double(0.0)\n")
        failure = subprocess.run([binary,str(invalid),"dark-photon",str(Path(temporary)/"invalid.json")],
                                 capture_output=True,text=True,timeout=180,env=runtime)
        assert failure.returncode != 0 and "DarkPhoton requires" in failure.stderr
    conflict = subprocess.run([binary,str(configuration),"conflict",str(Path(temporary)/"conflict.json")],
                              capture_output=True,text=True,timeout=180,env=runtime)
    assert conflict.returncode != 0 and "Conflicting Geant4 PDG32" in conflict.stderr
    # Preserve full compact inventories when a qualification caller requests it.
    destination = os.environ.get("SHIFT_DP_PROCESS_TEST_OUTPUT")
    if destination:
        target = Path(destination)
        target.mkdir(exist_ok=False)
        for mode,result in results.items():
            (target/(mode+".json")).write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
