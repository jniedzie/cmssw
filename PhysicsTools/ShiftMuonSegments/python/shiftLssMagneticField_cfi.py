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
    boundsShape="box",
    boundsCenterCm=(0.0, 0.0, 0.0),
    innerRadiusCm=0.0,
    outerRadiusCm=0.0,
    excludedCylindersCm=(),
):
    return cms.PSet(
        name=cms.string(name),
        type=cms.string("uniform"),
        minimumCm=cms.vdouble(*minimumCm),
        maximumCm=cms.vdouble(*maximumCm),
        originCm=cms.vdouble(*originCm),
        boundsShape=cms.string(boundsShape),
        boundsCenterCm=cms.vdouble(*boundsCenterCm),
        innerRadiusCm=cms.double(innerRadiusCm),
        outerRadiusCm=cms.double(outerRadiusCm),
        excludedCylindersCm=cms.vdouble(*excludedCylindersCm),
        localToGlobal=cms.vdouble(*localToGlobal),
        fieldTesla=cms.vdouble(*fieldTesla),
        gradientTeslaPerCm=cms.double(0.0),
        mapFile=cms.string(""),
        fieldScaleTesla=cms.double(1.0),
    )


def shiftLssQuadrupoleFieldElement(
    name,
    minimumCm,
    maximumCm,
    gradientTeslaPerCm,
    originCm=(0.0, 0.0, 0.0),
    localToGlobal=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    boundsShape="box",
    boundsCenterCm=(0.0, 0.0, 0.0),
    innerRadiusCm=0.0,
    outerRadiusCm=0.0,
    excludedCylindersCm=(),
):
    return cms.PSet(
        name=cms.string(name),
        type=cms.string("quadrupole"),
        minimumCm=cms.vdouble(*minimumCm),
        maximumCm=cms.vdouble(*maximumCm),
        originCm=cms.vdouble(*originCm),
        boundsShape=cms.string(boundsShape),
        boundsCenterCm=cms.vdouble(*boundsCenterCm),
        innerRadiusCm=cms.double(innerRadiusCm),
        outerRadiusCm=cms.double(outerRadiusCm),
        excludedCylindersCm=cms.vdouble(*excludedCylindersCm),
        localToGlobal=cms.vdouble(*localToGlobal),
        fieldTesla=cms.vdouble(0.0, 0.0, 0.0),
        gradientTeslaPerCm=cms.double(gradientTeslaPerCm),
        mapFile=cms.string(""),
        fieldScaleTesla=cms.double(1.0),
    )


def shiftLssFlukaMap2DFieldElement(
    name,
    minimumCm,
    maximumCm,
    mapFile,
    fieldScaleTesla,
    originCm=(0.0, 0.0, 0.0),
    localToGlobal=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    boundsShape="box",
    boundsCenterCm=(0.0, 0.0, 0.0),
    innerRadiusCm=0.0,
    outerRadiusCm=0.0,
    excludedCylindersCm=(),
):
    return cms.PSet(
        name=cms.string(name),
        type=cms.string("flukaMap2D"),
        minimumCm=cms.vdouble(*minimumCm),
        maximumCm=cms.vdouble(*maximumCm),
        originCm=cms.vdouble(*originCm),
        boundsShape=cms.string(boundsShape),
        boundsCenterCm=cms.vdouble(*boundsCenterCm),
        innerRadiusCm=cms.double(innerRadiusCm),
        outerRadiusCm=cms.double(outerRadiusCm),
        excludedCylindersCm=cms.vdouble(*excludedCylindersCm),
        localToGlobal=cms.vdouble(*localToGlobal),
        fieldTesla=cms.vdouble(0.0, 0.0, 0.0),
        gradientTeslaPerCm=cms.double(0.0),
        mapFile=cms.string(mapFile),
        fieldScaleTesla=cms.double(fieldScaleTesla),
    )
