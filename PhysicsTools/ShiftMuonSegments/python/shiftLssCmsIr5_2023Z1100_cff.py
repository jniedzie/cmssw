# Generated, checksummed field placement for the CMS IR5 LSS model.
# Do not edit by hand. The caller must supply the reviewed model-to-CMS transform.

import math
import os

from PhysicsTools.ShiftMuonSegments.shiftLssMagneticField_cfi import (
    shiftLssFlukaMap2DFieldElement,
    shiftLssUniformFieldElement,
)


_DATA_DIRECTORY = 'PhysicsTools/ShiftMuonSegments/data/lss/cms_ir5_2023_z1100'
_FIELD_DOMAINS_SHA256 = 'c51f96336d9824c303adbd29a3aae4d43436f31dfa9628ac244be58ebeac64d9'
_ELEMENTS = (
    dict(name='MQXA.QXA1R__a', type='flukaMap2D', minimum=(-24.5, -24.5, -332), maximum=(24.5, 24.5, 332), origin=(-0, -0, 2615), bounds_shape='cylinderZ', bounds_center=(0, 0, 0), inner_radius=3.4500000000000002, outer_radius=24.5, excluded_cylinders=(), map_file='MQXA.dat', field_scale=1.9741973648239945),
    dict(name='MQXA.QXA1R__c', type='flukaMap2D', minimum=(-3.2999999999999998, -3.2999999999999998, -332), maximum=(3.2999999999999998, 3.2999999999999998, 332), origin=(-0, -0, 2615), bounds_shape='cylinderZ', bounds_center=(0, 0, 0), inner_radius=0, outer_radius=3.2999999999999998, excluded_cylinders=(), map_file='MQXA.dat', field_scale=1.9741973648239945),
    dict(name='DIPOLE.CBXH1R_a', type='uniform', minimum=(-17.5, -17.5, -33.899999999999999), maximum=(17.5, 17.5, 36.100000000000001), origin=(-0, -0, 2984.1999999999998), bounds_shape='cylinderZ', bounds_center=(0, 0, 0), inner_radius=4.5, outer_radius=17.5, excluded_cylinders=(), field=(0, 0.31232785923831791, 0)),
    dict(name='DIPOLE.CBXH1R_c', type='uniform', minimum=(-3.2999999999999998, -3.2999999999999998, -33.899999999999999), maximum=(3.2999999999999998, 3.2999999999999998, 36.100000000000001), origin=(-0, -0, 2984.1999999999998), bounds_shape='cylinderZ', bounds_center=(0, 0, 0), inner_radius=0, outer_radius=3.2999999999999998, excluded_cylinders=(), field=(0, 0.31232785923831791, 0)),
    dict(name='MQXB.QXBA2R', type='flukaMap2D', minimum=(-20.100000000000001, -20.100000000000001, -287), maximum=(20.100000000000001, 20.100000000000001, 287), origin=(-0, -0, 3480), bounds_shape='cylinderZ', bounds_center=(0, 0, 0), inner_radius=0, outer_radius=20.100000000000001, excluded_cylinders=(), map_file='MQXB.dat', field_scale=-1.9916423070862641),
    dict(name='DIPOLE.CBXH2R_a', type='uniform', minimum=(-17.5, -17.5, -33.899999999999999), maximum=(17.5, 17.5, 36.100000000000001), origin=(-0, -0, 3801.9000000000001), bounds_shape='cylinderZ', bounds_center=(0, 0, 0), inner_radius=4.5, outer_radius=17.5, excluded_cylinders=(), field=(0, 0.31232785923831791, 0)),
    dict(name='DIPOLE.CBXH2R_c', type='uniform', minimum=(-3.3250000000000002, -3.3250000000000002, -33.899999999999999), maximum=(3.3250000000000002, 3.3250000000000002, 36.100000000000001), origin=(-0, -0, 3801.9000000000001), bounds_shape='cylinderZ', bounds_center=(0, 0, 0), inner_radius=0, outer_radius=3.3250000000000002, excluded_cylinders=(), field=(0, 0.31232785923831791, 0)),
    dict(name='MQXB.QXBB2R', type='flukaMap2D', minimum=(-20.100000000000001, -20.100000000000001, -287), maximum=(20.100000000000001, 20.100000000000001, 287), origin=(-0, -0, 4130), bounds_shape='cylinderZ', bounds_center=(0, 0, 0), inner_radius=0, outer_radius=20.100000000000001, excluded_cylinders=(), map_file='MQXB.dat', field_scale=-1.9916423070862641),
    dict(name='MQXA.QXA3R__a', type='flukaMap2D', minimum=(-24.5, -24.5, -332), maximum=(24.5, 24.5, 332), origin=(-0, -0, 5015), bounds_shape='cylinderZ', bounds_center=(0, 0, 0), inner_radius=3.4500000000000002, outer_radius=24.5, excluded_cylinders=(), map_file='MQXA.dat', field_scale=1.983922303527639),
    dict(name='MQXA.QXA3R__c', type='flukaMap2D', minimum=(-3.3250000000000002, -3.3250000000000002, -332), maximum=(3.3250000000000002, 3.3250000000000002, 332), origin=(-0, -0, 5015), bounds_shape='cylinderZ', bounds_center=(0, 0, 0), inner_radius=0, outer_radius=3.3250000000000002, excluded_cylinders=(), map_file='MQXA.dat', field_scale=1.983922303527639),
    dict(name='DIPOLE.CBXH3R_a', type='uniform', minimum=(-17.5, -17.5, -33.899999999999999), maximum=(17.5, 17.5, 36.100000000000001), origin=(-0, -0, 5381.3999999999996), bounds_shape='cylinderZ', bounds_center=(0, 0, 0), inner_radius=4.5, outer_radius=17.5, excluded_cylinders=(), field=(0, 0.31232785923831791, 0)),
    dict(name='DIPOLE.CBXH3R_c', type='uniform', minimum=(-4.4500000000000002, -4.4500000000000002, -33.899999999999999), maximum=(4.4500000000000002, 4.4500000000000002, 36.100000000000001), origin=(-0, -0, 5381.3999999999996), bounds_shape='cylinderZ', bounds_center=(0, 0, 0), inner_radius=3.5, outer_radius=4.4500000000000002, excluded_cylinders=(), field=(0, 0.31232785923831791, 0)),
    dict(name='DIPOLE.CBXH3R_e', type='uniform', minimum=(-3.3250000000000002, -3.3250000000000002, -33.899999999999999), maximum=(3.3250000000000002, 3.3250000000000002, 36.100000000000001), origin=(-0, -0, 5381.3999999999996), bounds_shape='cylinderZ', bounds_center=(0, 0, 0), inner_radius=0, outer_radius=3.3250000000000002, excluded_cylinders=(), field=(0, 0.31232785923831791, 0)),
    dict(name='MBXW.BXWA4R', type='flukaMap2D', minimum=(-40.25, -30.100000000000001, -194.19999999999999), maximum=(40.25, 30.100000000000001, 192.59999999999999), origin=(-0, -0, 6132.1999999999998), bounds_shape='box', bounds_center=(0, 0, 0), inner_radius=0, outer_radius=0, excluded_cylinders=(), map_file='MBXW.dat', field_scale=-1.2553669967398966),
    dict(name='MBXW.BXWB4R', type='flukaMap2D', minimum=(-40.25, -30.100000000000001, -194.19999999999999), maximum=(40.25, 30.100000000000001, 192.59999999999999), origin=(-0, -0, 6558.8000000000002), bounds_shape='box', bounds_center=(0, 0, 0), inner_radius=0, outer_radius=0, excluded_cylinders=(), map_file='MBXW.dat', field_scale=-1.2553669967398966),
    dict(name='MBXW.BXWC4R', type='flukaMap2D', minimum=(-40.25, -30.100000000000001, -194.19999999999999), maximum=(40.25, 30.100000000000001, 192.59999999999999), origin=(-0, -0, 6985.3999999999996), bounds_shape='box', bounds_center=(0, 0, 0), inner_radius=0, outer_radius=0, excluded_cylinders=(), map_file='MBXW.dat', field_scale=-1.2553669967398966),
    dict(name='MBXW.BXWD4R', type='flukaMap2D', minimum=(-40.25, -30.100000000000001, -194.19999999999999), maximum=(40.25, 30.100000000000001, 192.59999999999999), origin=(-0, -0, 7412), bounds_shape='box', bounds_center=(0, 0, 0), inner_radius=0, outer_radius=0, excluded_cylinders=(), map_file='MBXW.dat', field_scale=-1.2553669967398966),
    dict(name='MBXW.BXWE4R', type='flukaMap2D', minimum=(-40.25, -30.100000000000001, -194.19999999999999), maximum=(40.25, 30.100000000000001, 192.59999999999999), origin=(-0, -0, 7838.6000000000004), bounds_shape='box', bounds_center=(0, 0, 0), inner_radius=0, outer_radius=0, excluded_cylinders=(), map_file='MBXW.dat', field_scale=-1.2553669967398966),
    dict(name='MBXW.BXWF4R', type='flukaMap2D', minimum=(-40.25, -30.100000000000001, -194.19999999999999), maximum=(40.25, 30.100000000000001, 192.59999999999999), origin=(-0, -0, 8265.2000000000007), bounds_shape='box', bounds_center=(0, 0, 0), inner_radius=0, outer_radius=0, excluded_cylinders=(), map_file='MBXW.dat', field_scale=-1.2553669967398966),
    dict(name='DIPOLE.BRC4R', type='uniform', minimum=(-20, -20, -490.69999999999999), maximum=(20, 20, 490.69999999999999), origin=(-0, -0, 15790), bounds_shape='cylinderZ', bounds_center=(0, 0, 0), inner_radius=0, outer_radius=20, excluded_cylinders=(), field=(0, 2.7099984571548212, 0)),
    dict(name='DIPOLE.CBYH4R', type='uniform', minimum=(-20, -20, -54.799999999999997), maximum=(20, 20, 54.799999999999997), origin=(-0, -0, 16573.5), bounds_shape='cylinderZ', bounds_center=(0, 0, 0), inner_radius=0, outer_radius=20, excluded_cylinders=(), field=(0, -0.87986527590314367, 0)),
    dict(name='MQYana.QY4R', type='flukaMap2D', minimum=(-24.75, -24.75, -179.59999999999999), maximum=(24.75, 24.75, 179.59999999999999), origin=(-0, -0, 16955.299999999999), bounds_shape='cylinderZ', bounds_center=(0, 0, 0), inner_radius=0, outer_radius=24.75, excluded_cylinders=(), map_file='MQYana.dat', field_scale=-0.64397635566936118),
)


def _validated_transform(modelOriginCm, modelToCms):
    if len(modelOriginCm) != 3 or len(modelToCms) != 9:
        raise ValueError("modelOriginCm and modelToCms must contain 3 and 9 values")
    values = tuple(float(value) for value in (*modelOriginCm, *modelToCms))
    if not all(math.isfinite(value) for value in values):
        raise ValueError("coordinate transform must contain only finite values")
    origin, rotation = values[:3], values[3:]
    for row in range(3):
        for other in range(3):
            dot = sum(rotation[3 * row + column] * rotation[3 * other + column]
                      for column in range(3))
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


def shiftLssCmsIr5_2023Z1100FieldElements(
    *, modelOriginCm, modelToCms, fieldScale=1.0, dataDirectory=None
):
    model_origin, rotation = _validated_transform(modelOriginCm, modelToCms)
    field_scale = float(fieldScale)
    if not math.isfinite(field_scale) or field_scale == 0.0:
        raise ValueError("fieldScale must be finite and nonzero")
    map_directory = _DATA_DIRECTORY if dataDirectory is None else dataDirectory
    if dataDirectory is not None and not os.path.isabs(map_directory):
        raise ValueError("dataDirectory must be an absolute path")

    def cms_origin(origin):
        return tuple(
            model_origin[axis]
            + sum(rotation[3 * axis + local] * origin[local] for local in range(3))
            for axis in range(3)
        )

    result = []
    for element in _ELEMENTS:
        common = dict(
            originCm=cms_origin(element["origin"]),
            localToGlobal=rotation,
            boundsShape=element["bounds_shape"],
            boundsCenterCm=element["bounds_center"],
            innerRadiusCm=element["inner_radius"],
            outerRadiusCm=element["outer_radius"],
            excludedCylindersCm=element["excluded_cylinders"],
        )
        if element["type"] == "uniform":
            result.append(shiftLssUniformFieldElement(
                element["name"], element["minimum"], element["maximum"],
                tuple(field_scale * value for value in element["field"]), **common,
            ))
        else:
            result.append(shiftLssFlukaMap2DFieldElement(
                element["name"], element["minimum"], element["maximum"],
                f"{map_directory}/{element['map_file']}",
                field_scale * element["field_scale"], **common,
            ))
    return result
