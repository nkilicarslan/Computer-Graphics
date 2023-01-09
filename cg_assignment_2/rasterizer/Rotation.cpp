#include "Rotation.h"
#include "Helpers.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

Rotation::Rotation() {}

Rotation::Rotation(int rotationId, double angle, double x, double y, double z)
{
    this->rotationId = rotationId;
    this->angle = angle;
    this->ux = x;
    this->uy = y;
    this->uz = z;
}

ostream &operator<<(ostream &os, const Rotation &r)
{
    os << fixed << setprecision(3) << "Rotation " << r.rotationId << " => [angle=" << r.angle << ", " << r.ux << ", " << r.uy << ", " << r.uz << "]";

    return os;
}

// create Matrix4 for rotation with angle, ux, uy, uz
/*Value specifies α, ux, uy, and uz , which are the rotation parameters, i.e., the object is
 * rotated α degrees around the rotation axis which pass through points (ux, uy, uz ) and (0, 0, 0).
 * The positive angle of rotation is given as the counter-clockwise rotation along the direction (ux, uy, uz ).
 */
void Rotation::applyRotation(Matrix4 &matrix4) {
    // yanlissa slayttakini yap
    Matrix4 rotationMatrix = getIdentityMatrix();
    double c = cos(angle * M_PI / 180);
    double s = sin(angle * M_PI / 180);
    double t = 1 - c;
    rotationMatrix.val[0][0] = t * ux * ux + c;
    rotationMatrix.val[0][1] = t * ux * uy - s * uz;
    rotationMatrix.val[0][2] = t * ux * uz + s * uy;
    rotationMatrix.val[1][0] = t * ux * uy + s * uz;
    rotationMatrix.val[1][1] = t * uy * uy + c;
    rotationMatrix.val[1][2] = t * uy * uz - s * ux;
    rotationMatrix.val[2][0] = t * ux * uz - s * uy;
    rotationMatrix.val[2][1] = t * uy * uz + s * ux;
    rotationMatrix.val[2][2] = t * uz * uz + c;

    // multiply pMatrix4 with rotationMatrix and store the result in pMatrix4
    matrix4 = multiplyMatrixWithMatrix(rotationMatrix, matrix4);

}


