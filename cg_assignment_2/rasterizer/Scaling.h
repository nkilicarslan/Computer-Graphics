#ifndef __SCALING_H__
#define __SCALING_H__

#include <iostream>
#include "Triangle.h"
#include "Matrix4.h"

using namespace std;

class Scaling
{
public:
    int scalingId;
    double sx, sy, sz;

    Scaling();
    Scaling(int scalingId, double sx, double sy, double sz);
    friend ostream &operator<<(ostream &os, const Scaling &s);

    void applyScaling(Matrix4 &pMatrix4) const;
};

#endif