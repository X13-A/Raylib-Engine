#pragma once

#include <raymath.h>
#include "MyStruct.h"

Polar CartesianToPolar(Vector2 cart, bool keepThetaPositive = true);

Vector2 PolarToCartesian(Polar polar);

Cylindrical CartesianToCylindrical(Vector3 cart, bool keepThetaPositive = true);

Vector3 CylindricalToCartesian(Cylindrical cylindrical);

Spherical CartesianToSpherical(Vector3 cart, bool keepThetaPositive = true);

Vector3 SphericalToCartesian(Spherical spherical);