<!-- Code style: human-code. -->
# Optional dark-photon particle in the nominal CMS physics list

`SimG4Core/Physics/FTFP_BERT_EMM_DarkPhoton` calls the existing CMS
`FTFP_BERT_EMM` factory and returns its actual physics list, with one added
neutral spin-one PDG32 particle.
The existing nominal physics list supplies transportation and standard decay;
the extension checks that each occurs exactly once. For the four physical kaon
species, a subclass of standard decay also updates the nested daughter frames
when a preassigned dark photon is present. Ordinary kaon decays use the inherited
behavior. Process names, types, ordering and all other particle processes stay
unchanged.

The selected signal point must provide an explicit tracked `DarkPhoton` PSet
with finite positive `massGeV` and `properCTauMm`. Use the physical proper
lifetime from the point's pinned width metadata. A mean LAB flight distance is
not proper c-tau. The configuration rejects conflicting particle definitions.
Existing campaigns keep their original physics-list selection.

The geometry-free process-inventory test compares nominal, a no-op factory wrapper,
and the optional extension in separate processes. It checks every defined
ordinary particle's processes and ordering. The standalone test uses individual
CMS plugins and Geant4 shared libraries, avoiding a second embedded Geant4 copy
from CMSSW's Simulation bundle. This is not a detector transport, material,
trigger, kaon-coherence or source-admission qualification.

Build and test in a private release. Do not rebuild a shared CMSSW release
while production uses it. Full simulation and production require a separately
reviewed frozen configuration and the remaining source and transport checks.
