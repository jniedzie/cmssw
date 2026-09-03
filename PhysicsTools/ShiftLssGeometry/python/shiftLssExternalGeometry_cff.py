import FWCore.ParameterSet.Config as cms


def customiseShiftLssExternalGeometry(
    process,
    *,
    gdmlFile,
    artifactOriginInModelCm,
    modelOriginCm,
    modelToCms,
    minimumAbsZCm,
    geometryLabel="Extended",
    detectorElementName="shiftLssExternal",
    externalMotherVolumeName="cms:CMSE",
    overlapToleranceCm=0.001,
    checkOverlaps=True,
):
    """Append bounded LSS material to, and never replace, standard CMS geometry."""
    if not hasattr(process, "DDDetectorESProducerFromDB"):
        raise RuntimeError(
            "Cannot extend LSS geometry: standard DDDetectorESProducerFromDB is absent"
        )
    if hasattr(process, "shiftLssGeometryESSource"):
        raise RuntimeError("SHIFT LSS external geometry is already configured")
    if len(artifactOriginInModelCm) != 3 or len(modelOriginCm) != 3 or len(modelToCms) != 9:
        raise ValueError(
            "artifactOriginInModelCm, modelOriginCm, and modelToCms require 3, 3, and 9 values"
        )
    if minimumAbsZCm <= 0.0 or overlapToleranceCm <= 0.0:
        raise ValueError("minimumAbsZCm and overlapToleranceCm must be positive")

    source = process.DDDetectorESProducerFromDB
    if not source.fromDB.value():
        raise RuntimeError("SHIFT LSS extension requires the standard geometry conditions payload")
    if source.label.value() != geometryLabel:
        raise ValueError(
            "geometryLabel must match DDDetectorESProducerFromDB.label so the same CMS payload is used"
        )

    process.shiftLssGeometryESSource = cms.ESSource(
        "ShiftLssGeometryESSource",
        gdmlFile=cms.FileInPath(gdmlFile),
        geometryLabel=cms.string(geometryLabel),
        detectorElementName=cms.string(detectorElementName),
        externalMotherVolumeName=cms.string(externalMotherVolumeName),
        artifactOriginInModelCm=cms.vdouble(*artifactOriginInModelCm),
        modelOriginCm=cms.vdouble(*modelOriginCm),
        modelToCms=cms.vdouble(*modelToCms),
        minimumAbsZCm=cms.double(minimumAbsZCm),
        overlapToleranceCm=cms.double(overlapToleranceCm),
        checkOverlaps=cms.bool(checkOverlaps),
    )
    del process.DDDetectorESProducerFromDB
    process.shiftLssGeometryContract = cms.PSet(
        standardCmsGeometryPreserved=cms.bool(True),
        externalExtensionOnly=cms.bool(True),
        geometryConditionsLabel=cms.string(geometryLabel),
        externalMotherVolumeName=cms.string(externalMotherVolumeName),
        gdmlFile=cms.string(gdmlFile),
        artifactOriginInModelCm=cms.vdouble(*artifactOriginInModelCm),
        modelOriginCm=cms.vdouble(*modelOriginCm),
        modelToCms=cms.vdouble(*modelToCms),
        minimumAbsZCm=cms.double(minimumAbsZCm),
    )
    return process
