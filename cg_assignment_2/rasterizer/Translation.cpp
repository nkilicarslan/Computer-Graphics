#include "Translation.h"
#include <iostream>
#include <iomanip>
#include "Helpers.h"

using namespace std;

Translation::Translation()
{
    this->translationId = -1;
    this->tx = 0.0;
    this->ty = 0.0;
    this->tz = 0.0;
}

Translation::Translation(int translationId, double tx, double ty, double tz)
{
    this->translationId = translationId;
    this->tx = tx;
    this->ty = ty;
    this->tz = tz;
}

ostream &operator<<(ostream &os, const Translation &t)
{
    os << fixed << setprecision(3) << "Translation " << t.translationId << " => [" << t.tx << ", " << t.ty << ", " << t.tz << "]";

    return os;
}

void Translation::applyTranslation(Matrix4 &matrix4) const {
    // create Matrix4 for translation with tx, ty, tz
    Matrix4 translationMatrix = getIdentityMatrix();
    translationMatrix.val[0][3] = this->tx;
    translationMatrix.val[1][3] = this->ty;
    translationMatrix.val[2][3] = this->tz;

    // multiply matrix4 with translationMatrix and store the result in matrix4
    matrix4 = multiplyMatrixWithMatrix(translationMatrix, matrix4);

}
