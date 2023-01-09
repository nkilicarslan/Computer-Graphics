#include "Scaling.h"
#include "Helpers.h"
#include <iostream>
#include <iomanip>

using namespace std;

Scaling::Scaling() {}

Scaling::Scaling(int scalingId, double sx, double sy, double sz)
{
    this->scalingId = scalingId;
    this->sx = sx;
    this->sy = sy;
    this->sz = sz;
}

ostream &operator<<(ostream &os, const Scaling &s)
{
    os << fixed << setprecision(3) << "Scaling " << s.scalingId << " => [" << s.sx << ", " << s.sy << ", " << s.sz << "]";

    return os;
}

void Scaling::applyScaling(Matrix4 &pMatrix4) const
{
    // create Matrix4 for scaling with sx, sy, sz
    Matrix4 scalingMatrix = getIdentityMatrix();
    scalingMatrix.val[0][0] = this->sx;
    scalingMatrix.val[1][1] = this->sy;
    scalingMatrix.val[2][2] = this->sz;

    // multiply pMatrix4 with scalingMatrix and store the result in pMatrix4
    pMatrix4 = multiplyMatrixWithMatrix(scalingMatrix, pMatrix4);
}
