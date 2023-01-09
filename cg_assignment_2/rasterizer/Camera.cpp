#include "Camera.h"
#include "Matrix4.h"
#include <string>
#include <iostream>
#include <iomanip>
#include "Helpers.h"

using namespace std;

Camera::Camera() {}

Camera::Camera(int cameraId,
               int projectionType,
               Vec3 pos, Vec3 gaze,
               Vec3 u, Vec3 v, Vec3 w,
               double left, double right, double bottom, double top,
               double near, double far,
               int horRes, int verRes,
               string outputFileName)
{

    this->cameraId = cameraId;
    this->projectionType = projectionType;
    this->pos = pos;
    this->gaze = gaze;
    this->u = u;
    this->v = v;
    this->w = w;
    this->left = left;
    this->right = right;
    this->bottom = bottom;
    this->top = top;
    this->near = near;
    this->far = far;
    this->horRes = horRes;
    this->verRes = verRes;
    this->outputFileName = outputFileName;
}

Camera::Camera(const Camera &other)
{
    this->cameraId = other.cameraId;
    this->projectionType = other.projectionType;
    this->pos = other.pos;
    this->gaze = other.gaze;
    this->u = other.u;
    this->v = other.v;
    this->w = other.w;
    this->left = other.left;
    this->right = other.right;
    this->bottom = other.bottom;
    this->top = other.top;
    this->near = other.near;
    this->far = other.far;
    this->horRes = other.horRes;
    this->verRes = other.verRes;
    this->outputFileName = other.outputFileName;
}

// calculate the camera matrix return the result matrix
void Camera::getCameraMatrix(Matrix4 &matrix) const
{
    // move camera to origin
    // create Identity matrix
    Matrix4 cameraMatrix = getIdentityMatrix();
    // translate camera to origin
    cameraMatrix.val[0][3] = -pos.x;
    cameraMatrix.val[1][3] = -pos.y;
    cameraMatrix.val[2][3] = -pos.z;

    // create the rotation matrix
    Matrix4 rotationMatrix = getIdentityMatrix();
    rotationMatrix.val[0][0] = u.x;
    rotationMatrix.val[0][1] = u.y;
    rotationMatrix.val[0][2] = u.z;
    rotationMatrix.val[1][0] = v.x;
    rotationMatrix.val[1][1] = v.y;
    rotationMatrix.val[1][2] = v.z;
    rotationMatrix.val[2][0] = w.x;
    rotationMatrix.val[2][1] = w.y;
    rotationMatrix.val[2][2] = w.z;

    // multiply cameraMatrix with rotationMatrix and store the result in cameraMatrix
    cameraMatrix = multiplyMatrixWithMatrix(rotationMatrix, cameraMatrix);

    matrix = cameraMatrix;

}

void Camera::getOrthographicMatrix(Matrix4 &matrix4) const {
    // create the orthographic matrix
    Matrix4 orthographicMatrix = getIdentityMatrix();
    orthographicMatrix.val[0][0] = 2 / (right - left);
    orthographicMatrix.val[1][1] = 2 / (top - bottom);
    orthographicMatrix.val[2][2] = -2 / (far - near);
    orthographicMatrix.val[0][3] = -(right + left) / (right - left);
    orthographicMatrix.val[1][3] = -(top + bottom) / (top - bottom);
    orthographicMatrix.val[2][3] = -(far + near) / (far - near);
    matrix4 = orthographicMatrix;
}



ostream &operator<<(ostream &os, const Camera &c)
{
    const char *camType = c.projectionType ? "perspective" : "orthographic";

    os << fixed << setprecision(6) << "Camera " << c.cameraId << " (" << camType << ") => pos: " << c.pos << " gaze: " << c.gaze << endl
       << "\tu: " << c.u << " v: " << c.v << " w: " << c.w << endl
       << fixed << setprecision(3) << "\tleft: " << c.left << " right: " << c.right << " bottom: " << c.bottom << " top: " << c.top << endl
       << "\tnear: " << c.near << " far: " << c.far << " resolutions: " << c.horRes << "x" << c.verRes << " fileName: " << c.outputFileName;

    return os;
}

void Camera::getPerspectiveMatrix(Matrix4 &matrix4) const {
    // create the perspective matrix
    Matrix4 perspectiveMatrix = getIdentityMatrix();
    perspectiveMatrix.val[0][0] = 2 * near / (right - left);
    perspectiveMatrix.val[1][1] = 2 * near / (top - bottom);
    perspectiveMatrix.val[0][2] = (right + left) / (right - left);
    perspectiveMatrix.val[1][2] = (top + bottom) / (top - bottom);
    perspectiveMatrix.val[2][2] = -(far + near) / (far - near);
    perspectiveMatrix.val[2][3] = -2 * far * near / (far - near);
    perspectiveMatrix.val[3][2] = -1;
    perspectiveMatrix.val[3][3] = 0;
    matrix4 = perspectiveMatrix;
}

void Camera::getViewportMatrix(Matrix4 &matrix4) const {
    // create the viewport matrix

    Matrix4 viewportMatrix = getIdentityMatrix();

    viewportMatrix.val[0][0] = horRes / 2.0;
    viewportMatrix.val[1][1] = verRes / 2.0;
    viewportMatrix.val[0][3] = (horRes-1) / 2.0;
    viewportMatrix.val[1][3] = (verRes-1) / 2.0;
    viewportMatrix.val[2][2] = 0.5;
    viewportMatrix.val[2][3] = 0.5;

    matrix4 = viewportMatrix;

}


