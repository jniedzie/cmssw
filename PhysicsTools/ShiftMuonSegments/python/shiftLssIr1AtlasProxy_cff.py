"""Provisional magnetic-field placements from the frozen ATLAS IR1 model.

The caller must supply the coordinate transform into CMS. This prevents the
IR1 proxy from being mistaken for an aligned IR5 model.
"""

import math

from PhysicsTools.ShiftMuonSegments.shiftLssMagneticField_cfi import (
    shiftLssFlukaMap2DFieldElement,
    shiftLssUniformFieldElement,
)


_DATA_DIRECTORY = "PhysicsTools/ShiftMuonSegments/data/lss/ir1_atlas_proxy"

# name, map, FLUKA field scale [T or T/cm], origin [cm], local AABB [cm]
_MAP_PLACEMENTS = (
    ("MQXA.1R1.outer", "MQXA.dat", 1.8928521203, (0.0, 0.0, 2615.0), (-24.5, -24.5, -332.0), (24.5, 24.5, 332.0)),
    ("MQXA.1R1.aperture", "MQXA.dat", 1.8928521203, (0.0, 0.0, 2615.0), (-3.3, -3.3, -332.0), (3.3, 3.3, 332.0)),
    ("MQXB.A2R1", "MQXB.dat", -1.8928521203, (0.0, 0.0, 3480.0), (-20.1, -20.1, -287.0), (20.1, 20.1, 287.0)),
    ("MQXB.B2R1", "MQXB.dat", -1.8928521203, (0.0, 0.0, 4130.0), (-20.1, -20.1, -287.0), (20.1, 20.1, 287.0)),
    ("MQXA.3R1.outer", "MQXA.dat", 1.8928521203, (0.0, 0.0, 5015.0), (-24.5, -24.5, -332.0), (24.5, 24.5, 332.0)),
    ("MQXA.3R1.aperture", "MQXA.dat", 1.8928521203, (0.0, 0.0, 5015.0), (-3.325, -3.325, -332.0), (3.325, 3.325, 332.0)),
    ("MBXW.A4R1", "MBXW.dat", -1.1999831604, (0.0, 0.0, 6132.2), (-40.25, -30.1, -194.2), (40.25, 30.1, 192.6)),
    ("MBXW.B4R1", "MBXW.dat", -1.1999831604, (0.0, 0.0, 6558.8), (-40.25, -30.1, -194.2), (40.25, 30.1, 192.6)),
    ("MBXW.C4R1", "MBXW.dat", -1.1999831604, (0.0, 0.0, 6985.4), (-40.25, -30.1, -194.2), (40.25, 30.1, 192.6)),
    ("MBXW.D4R1", "MBXW.dat", -1.1999831604, (0.0, 0.0, 7412.0), (-40.25, -30.1, -194.2), (40.25, 30.1, 192.6)),
    ("MBXW.E4R1", "MBXW.dat", -1.1999831604, (0.0, 0.0, 7838.6), (-40.25, -30.1, -194.2), (40.25, 30.1, 192.6)),
    ("MBXW.F4R1", "MBXW.dat", -1.1999831604, (0.0, 0.0, 8265.2), (-40.25, -30.1, -194.2), (40.25, 30.1, 192.6)),
    ("MQY.4R1", "MQYana.dat", -0.72198822788, (0.0, 0.0, 16955.3), (-24.75, -24.75, -179.6), (24.75, 24.75, 179.6)),
    ("MQML.5R1", "MQTL.dat", -0.37950611892, (0.0, 0.0, 19649.0), (-45.7, -53.7, -250.4), (45.7, 37.7, 247.4)),
    ("MQML.6R1", "MQTL.dat", 0.11863777627, (0.0, 0.0, 22839.0), (-45.7, -53.7, -250.4), (45.7, 37.7, 247.4)),
    ("MQM.A7R1", "MQTL.dat", -1.6090171308, (0.0, 0.0, 26170.4), (-45.7, -53.7, -170.0), (45.7, 37.7, 170.0)),
    ("MQM.B7R1", "MQTL.dat", -1.6090171308, (0.0, 0.0, 26547.1), (-45.7, -53.7, -170.0), (45.7, 37.7, 170.0)),
)

_CYLINDRICAL_BOUNDS = {
    "MQXA.1R1.outer": ((0.0, 0.0, 0.0), 3.45, 24.5),
    "MQXA.1R1.aperture": ((0.0, 0.0, 0.0), 0.0, 3.3),
    "MQXB.A2R1": ((0.0, 0.0, 0.0), 0.0, 20.1),
    "MQXB.B2R1": ((0.0, 0.0, 0.0), 0.0, 20.1),
    "MQXA.3R1.outer": ((0.0, 0.0, 0.0), 3.45, 24.5),
    "MQXA.3R1.aperture": ((0.0, 0.0, 0.0), 0.0, 3.325),
    "MQY.4R1": ((0.0, 0.0, 0.0), 0.0, 24.75),
    "MQML.5R1": ((0.0, -8.0, 0.0), 0.0, 45.7),
    "MQML.6R1": ((0.0, -8.0, 0.0), 0.0, 45.7),
    "MQM.A7R1": ((0.0, -8.0, 0.0), 0.0, 45.7),
    "MQM.B7R1": ((0.0, -8.0, 0.0), 0.0, 45.7),
}

# name, field [T], origin [cm], z bounds [cm], center, inner/outer radii,
# and excluded (x, y, radius) cylinders in the local model frame.
_CONSTANT_PLACEMENTS = (
    ("MCBXV.1R1.outer", (0.36866803721, 0.0, 0.0), (0.0, 0.0, 2984.2), -33.9, 36.1, (0.0, 0.0, 0.0), 4.5, 17.5, ()),
    ("MCBXV.1R1.aperture", (0.36866803721, 0.0, 0.0), (0.0, 0.0, 2984.2), -33.9, 36.1, (0.0, 0.0, 0.0), 0.0, 3.3, ()),
    ("MCBXV.2R1.outer", (0.36866803721, 0.0, 0.0), (0.0, 0.0, 3801.9), -33.9, 36.1, (0.0, 0.0, 0.0), 4.5, 17.5, ()),
    ("MCBXV.2R1.aperture", (0.36866803721, 0.0, 0.0), (0.0, 0.0, 3801.9), -33.9, 36.1, (0.0, 0.0, 0.0), 0.0, 3.325, ()),
    ("MCBXV.3R1.outer", (0.36866803721, 0.0, 0.0), (0.0, 0.0, 5381.4), -33.9, 36.1, (0.0, 0.0, 0.0), 4.5, 17.5, ()),
    ("MCBXV.3R1.middle", (0.36866803721, 0.0, 0.0), (0.0, 0.0, 5381.4), -33.925, 36.125, (0.0, 0.0, 0.0), 3.5, 4.45, ()),
    ("MCBXV.3R1.aperture", (0.36866803721, 0.0, 0.0), (0.0, 0.0, 5381.4), -33.9, 36.1, (0.0, 0.0, 0.0), 0.0, 3.325, ()),
    ("MBRC.4R1", (0.0, 2.5904398384, 0.0), (0.0, 0.0, 15790.0), -490.7, 490.7, (0.0, 0.0, 0.0), 0.0, 35.0, ()),
    ("MCBYH.4R1", (-1.3248807929, 0.0, 0.0), (0.0, 0.0, 16573.5), -54.8, 54.8, (0.0, 0.0, 0.0), 0.0, 24.75, ()),
    ("MCBCH.5R1.outer", (-0.28389660314, 0.0, 0.0), (0.0, 0.0, 19344.8), -45.0, 45.0, (0.0, -8.0, 0.0), 0.0, 45.7, (-9.7, 0.0, 2.8, 9.7, 0.0, 2.8)),
    ("MCBCH.5R1.negativeX", (-0.28389660314, 0.0, 0.0), (0.0, 0.0, 19344.8), -45.0, 45.0, (-9.7, 0.0, 0.0), 0.0, 2.69, ()),
    ("MCBCH.5R1.positiveX", (-0.28389660314, 0.0, 0.0), (0.0, 0.0, 19344.8), -45.0, 45.0, (9.7, 0.0, 0.0), 0.0, 2.69, ()),
)


def _validated_transform(modelOriginCm, modelToCms):
    if len(modelOriginCm) != 3 or len(modelToCms) != 9:
        raise ValueError("modelOriginCm and modelToCms must contain 3 and 9 values")
    values = tuple(float(value) for value in (*modelOriginCm, *modelToCms))
    if not all(math.isfinite(value) for value in values):
        raise ValueError("IR1 proxy coordinate transform must contain only finite values")
    origin = values[:3]
    rotation = values[3:]
    for row in range(3):
        for other in range(3):
            dot = sum(rotation[3 * row + column] * rotation[3 * other + column] for column in range(3))
            if abs(dot - (1.0 if row == other else 0.0)) > 1.0e-9:
                raise ValueError("modelToCms must be an orthonormal rotation")
    determinant = (
        rotation[0] * (rotation[4] * rotation[8] - rotation[5] * rotation[7])
        - rotation[1] * (rotation[3] * rotation[8] - rotation[5] * rotation[6])
        + rotation[2] * (rotation[3] * rotation[7] - rotation[4] * rotation[6])
    )
    if abs(determinant - 1.0) > 1.0e-9:
        raise ValueError("modelToCms must have determinant +1")
    return origin, rotation


def shiftLssIr1AtlasProxyFieldElements(
    *,
    modelOriginCm,
    modelToCms,
    fieldScale=1.0,
    includeMappedFields=True,
    includeConstantFields=True,
):
    """Build mapped IR1 field elements in an explicitly supplied CMS frame."""
    model_origin, rotation = _validated_transform(modelOriginCm, modelToCms)
    field_scale = float(fieldScale)
    if not math.isfinite(field_scale) or field_scale == 0.0:
        raise ValueError("fieldScale must be finite and nonzero")
    if not includeMappedFields and not includeConstantFields:
        raise ValueError("at least one IR1 proxy field family must be enabled")

    def cms_origin(origin):
        return tuple(
            model_origin[axis] + sum(rotation[3 * axis + local] * origin[local] for local in range(3))
            for axis in range(3)
        )

    elements = []
    if includeMappedFields:
        for name, map_name, scale, origin, minimum, maximum in _MAP_PLACEMENTS:
            element_parameters = {}
            if name in _CYLINDRICAL_BOUNDS:
                center, inner_radius, outer_radius = _CYLINDRICAL_BOUNDS[name]
                element_parameters = dict(
                    boundsShape="cylinderZ",
                    boundsCenterCm=center,
                    innerRadiusCm=inner_radius,
                    outerRadiusCm=outer_radius,
                )
            elements.append(
                shiftLssFlukaMap2DFieldElement(
                    name,
                    minimum,
                    maximum,
                    f"{_DATA_DIRECTORY}/{map_name}",
                    field_scale * scale,
                    originCm=cms_origin(origin),
                    localToGlobal=rotation,
                    **element_parameters,
                )
            )

    if includeConstantFields:
        for name, field, origin, minimum_z, maximum_z, center, inner_radius, outer_radius, exclusions in _CONSTANT_PLACEMENTS:
            minimum = (center[0] - outer_radius, center[1] - outer_radius, minimum_z)
            maximum = (center[0] + outer_radius, center[1] + outer_radius, maximum_z)
            elements.append(
                shiftLssUniformFieldElement(
                    name,
                    minimum,
                    maximum,
                    tuple(field_scale * component for component in field),
                    originCm=cms_origin(origin),
                    localToGlobal=rotation,
                    boundsShape="cylinderZ",
                    boundsCenterCm=center,
                    innerRadiusCm=inner_radius,
                    outerRadiusCm=outer_radius,
                    excludedCylindersCm=exclusions,
                )
            )
    return elements
