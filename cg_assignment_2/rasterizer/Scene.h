#ifndef _SCENE_H_
#define _SCENE_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "Vec4.h"

using namespace std;

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	vector< vector<Color> > image;
	vector< Camera* > cameras;
	vector< Vec3* > vertices;
	vector< Color* > colorsOfVertices;
	vector< Scaling* > scalings;
	vector< Rotation* > rotations;
	vector< Translation* > translations;
	vector< Mesh* > meshes;

	Scene(const char *xmlPath);

	void initializeImage(Camera* camera);
	void forwardRenderingPipeline(Camera* camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera* camera);
	void convertPPMToPNG(string ppmFileName, int osType);

    void getTransformationMatrix(Matrix4 &transformation_matrix, Mesh *mesh);

    vector<Vec4> getTransformedPoints(Triangle &triangle, Matrix4 &matrix4);

    vector<Color> getColorsOfTriangle(Triangle triangle);

    void rasterization(vector<Vec4> pVector, vector<Color> pVector1, Matrix4 pMatrix4, Camera* camera);

    static void multiplyVec4WithScalar(Vec4 &vec4, double d);

    bool clipping(Vec4 &point1, Vec4 &point2, Color &color1, Color &color2);

    bool visible(double den, double num, double &tE, double &tL);

    void line_rasterization(Vec4 vec4, Vec4 vec41, Color color, Color color1, Matrix4 matrix4, Camera *pCamera);
};

#endif
