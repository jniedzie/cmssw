import FWCore.ParameterSet.Config as cms


shiftLssMagneticField = cms.ESProducer(
    "ShiftLssMagneticFieldESProducer",
    # The base and output labels must differ to avoid an EventSetup cycle.
    baseFieldLabel=cms.string("shiftLssBaseMagneticField"),
    outputLabel=cms.string(""),
    # False replaces the base field inside an element and delegates to it
    # everywhere else. True adds each configured element to the base field.
    addBaseField=cms.bool(False),
    # False gives the first matching element priority; True sums overlaps.
    sumOverlaps=cms.bool(False),
    elements=cms.VPSet(),
)


def shiftLssUniformFieldElement(
    name,
    minimumCm,
    maximumCm,
    fieldTesla,
    originCm=(0.0, 0.0, 0.0),
    localToGlobal=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
):
    return cms.PSet(
        name=cms.string(name),
        type=cms.string("uniform"),
        minimumCm=cms.vdouble(*minimumCm),
        maximumCm=cms.vdouble(*maximumCm),
        originCm=cms.vdouble(*originCm),
        localToGlobal=cms.vdouble(*localToGlobal),
        fieldTesla=cms.vdouble(*fieldTesla),
        gradientTeslaPerCm=cms.double(0.0),
    )


def shiftLssQuadrupoleFieldElement(
    name,
    minimumCm,
    maximumCm,
    gradientTeslaPerCm,
    originCm=(0.0, 0.0, 0.0),
    localToGlobal=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
):
    return cms.PSet(
        name=cms.string(name),
        type=cms.string("quadrupole"),
        minimumCm=cms.vdouble(*minimumCm),
        maximumCm=cms.vdouble(*maximumCm),
        originCm=cms.vdouble(*originCm),
        localToGlobal=cms.vdouble(*localToGlobal),
        fieldTesla=cms.vdouble(0.0, 0.0, 0.0),
        gradientTeslaPerCm=cms.double(gradientTeslaPerCm),
    )
