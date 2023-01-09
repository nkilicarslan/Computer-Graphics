#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <vector>
#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"


using namespace tinyxml2;
using namespace std;

/*
	Transformations, clipping, culling, rasterization are done here.
*/

void Scene::getTransformationMatrix(Matrix4 &transformation_matrix, Mesh *mesh) {
    for (int i = 0; i < mesh->numberOfTransformations; i++) {
        int transformationId = mesh->transformationIds[i];
        char transformationType = mesh->transformationTypes[i];
        switch (transformationType) {
            case 's':
                this->scalings[transformationId - 1]->applyScaling(transformation_matrix);
                break;
            case 'r':
                this->rotations[transformationId - 1]->applyRotation(transformation_matrix);
                break;
            case 't':
                this->translations[transformationId - 1]->applyTranslation(transformation_matrix);
                break;
            default:
                break;
        }
    }
}

vector<Vec4> Scene::getTransformedPoints(Triangle &triangle, Matrix4 &matrix4) {
    vector<Vec4> transformedPoints;
    Vec3 first_vec = *vertices[triangle.getFirstVertexId()-1];
    Vec3 second_vec = *vertices[triangle.getSecondVertexId()-1];
    Vec3 third_vec = *vertices[triangle.getThirdVertexId()-1];


    Vec4 first_vec4 = Vec4(first_vec.x, first_vec.y, first_vec.z, 1, first_vec.colorId);
    Vec4 second_vec4 = Vec4(second_vec.x, second_vec.y, second_vec.z, 1, second_vec.colorId);
    Vec4 third_vec4 = Vec4(third_vec.x, third_vec.y, third_vec.z, 1, third_vec.colorId);

    transformedPoints.push_back(multiplyMatrixWithVec4(matrix4, first_vec4));
    transformedPoints.push_back(multiplyMatrixWithVec4(matrix4, second_vec4));
    transformedPoints.push_back(multiplyMatrixWithVec4(matrix4, third_vec4));

    return transformedPoints;
}

vector<Color> Scene::getColorsOfTriangle(Triangle triangle) {
    vector<Color> colors;
    Vec3 first_vec = *vertices[triangle.getFirstVertexId()-1];
    Vec3 second_vec = *vertices[triangle.getSecondVertexId()-1];
    Vec3 third_vec = *vertices[triangle.getThirdVertexId()-1];

    colors.push_back(*colorsOfVertices[first_vec.colorId-1]);
    colors.push_back(*colorsOfVertices[second_vec.colorId-1]);
    colors.push_back(*colorsOfVertices[third_vec.colorId-1]);
    return colors;
}

bool Scene::visible(double den, double num, double &tE, double &tL) {
    if(den > 0){
        double t = num/den;
        if(t > tL) return false;
        else if(t > tE) tE = t;
    }
    else if(den < 0){
        double t = num/den;
        if(t < tE) return false;
        else if(t < tL) tL = t;
    }
    else if(num > 0) return false;
    return true;
}

bool Scene::clipping(Vec4 &point1, Vec4 &point2, Color &color1, Color &color2) {
    double tE = 0, tL = 1;
    bool visible_var = false;

    point1.x /= point1.t;
    point1.y /= point1.t;
    point1.z /= point1.t;
    point1.t /= point1.t;
    point2.x /= point2.t;
    point2.y /= point2.t;
    point2.z /= point2.t;
    point2.t /= point2.t;

    double dx = point2.x - point1.x;
    double dy = point2.y - point1.y;
    double dz = point2.z - point1.z;

    Color diff_color;
    diff_color.r = color2.r - color1.r;
    diff_color.g = color2.g - color1.g;
    diff_color.b = color2.b - color1.b;


    if(visible(dx,-1-point1.x,tE,tL)) {
        if(visible(-dx,point1.x-1,tE,tL)) {
            if(visible(dy,-1-point1.y,tE,tL)) {
                if(visible(-dy,point1.y-1,tE,tL)) {
                    if(visible(dz,-1-point1.z,tE,tL)) {
                        if(visible(-dz,point1.z-1,tE,tL)) {
                            visible_var = true;
                            if(tL < 1) {
                                point2.x = point1.x + tL*dx;
                                point2.y = point1.y + tL*dy;
                                point2.z = point1.z + tL*dz;
                                color2.b = color1.b + diff_color.b * tL;
                                color2.g = color1.g + diff_color.g * tL;
                                color2.r = color1.r + diff_color.r * tL;

                            }
                            if(tE > 0) {
                                point1.x = point1.x + tE*dx;
                                point1.y = point1.y + tE*dy;
                                point1.z = point1.z + tE*dz;
                                color1.b = color1.b + diff_color.b * tE;
                                color1.g = color1.g + diff_color.g * tE;
                                color1.r = color1.r + diff_color.r * tE;
                            }
                        }
                    }
                }
            }
        }
    }
    return visible_var;

}



void Scene::rasterization(vector<Vec4> pVector, vector<Color> color_Vector, Matrix4 pMatrix4, Camera* camera) {
    for (int i = 0; i < pVector.size(); ++i) {
        multiplyVec4WithScalar(pVector[i], 1/pVector[i].t);
        pVector[i] = multiplyMatrixWithVec4(pMatrix4, pVector[i]);
    }
    //take the minimum and maximum x and y values
    int minX = pVector[0].x;
    int maxX = pVector[0].x;
    int minY = pVector[0].y;
    int maxY = pVector[0].y;
    for (int i = 1; i < pVector.size(); ++i) {
        if (pVector[i].x < minX) {
            minX = pVector[i].x;
        }
        if (pVector[i].x > maxX) {
            maxX = pVector[i].x;
        }
        if (pVector[i].y < minY) {
            minY = pVector[i].y;
        }
        if (pVector[i].y > maxY) {
            maxY = pVector[i].y;
        }
    }
    if(minX < 0) minX = 0;
    if(minX > camera->horRes -1 )   minX = camera->horRes -1;
    if(minY < 0) minY = 0;
    if(minY > camera->verRes -1 )   minY = camera->verRes -1;
    if(maxX < 0) maxX = 0;
    if(maxX > camera->horRes -1 )   maxX = camera->horRes -1;
    if(maxY < 0) maxY = 0;
    if(maxY > camera->verRes -1 )   maxY = camera->verRes -1;
    // create Vec3 vector using the pVector by discarding t value
    vector<Vec3> vec3_Vector;
    for (int i = 0; i < pVector.size(); ++i) {
        Vec3 temp = Vec3(pVector[i].x, pVector[i].y, pVector[i].z, pVector[i].colorId);
        vec3_Vector.push_back(temp);
    }
    // according to rasterization algorithm on slides 8-29, midpoint algorithm with baricentric coordinates
    for(int i=minX; i<=maxX; i++){
        for(int j=minY; j<=maxY; j++){
            // get the baricentric coordinates of vec3_Vector[0], vec3_Vector[1], vec3_Vector[2]
            double alpha = ((vec3_Vector[1].y - vec3_Vector[2].y)*i + (vec3_Vector[2].x - vec3_Vector[1].x)*j + vec3_Vector[1].x*vec3_Vector[2].y - vec3_Vector[2].x*vec3_Vector[1].y)
                           / ((vec3_Vector[1].y - vec3_Vector[2].y)*vec3_Vector[0].x + (vec3_Vector[2].x - vec3_Vector[1].x)*vec3_Vector[0].y + vec3_Vector[1].x*vec3_Vector[2].y - vec3_Vector[2].x*vec3_Vector[1].y);
            double beta = ((vec3_Vector[2].y - vec3_Vector[0].y)*i + (vec3_Vector[0].x - vec3_Vector[2].x)*j + vec3_Vector[2].x*vec3_Vector[0].y - vec3_Vector[0].x*vec3_Vector[2].y)
                          / ((vec3_Vector[1].y - vec3_Vector[2].y)*vec3_Vector[0].x + (vec3_Vector[2].x - vec3_Vector[1].x)*vec3_Vector[0].y + vec3_Vector[1].x*vec3_Vector[2].y - vec3_Vector[2].x*vec3_Vector[1].y);
            double gamma = 1 - alpha - beta;
            // if the point is inside the triangle
            if(alpha >= 0 && beta >= 0 && gamma >= 0){
                // create the c0 c1 c2
                Color c0 = color_Vector[0];
                Color c1 = color_Vector[1];
                Color c2 = color_Vector[2];
                c0.r = c0.r * alpha;
                c0.g = c0.g * alpha;
                c0.b = c0.b * alpha;
                c1.r = c1.r * beta;
                c1.g = c1.g * beta;
                c1.b = c1.b * beta;
                c2.r = c2.r * gamma;
                c2.g = c2.g * gamma;
                c2.b = c2.b * gamma;
                // add them together
                image[i][j].r = c0.r + c1.r + c2.r;
                image[i][j].g = c0.g + c1.g + c2.g;
                image[i][j].b = c0.b + c1.b + c2.b;
            }
        }
    }
}



void Scene::forwardRenderingPipeline(Camera *camera)
{
    Matrix4 camera_transformation_matrix = getIdentityMatrix();
    Matrix4 projection_transformation_matrix = getIdentityMatrix();
    Matrix4 viewport_transformation_matrix = getIdentityMatrix();

    camera->getCameraMatrix(camera_transformation_matrix);
    camera->getViewportMatrix(viewport_transformation_matrix);

    if(camera->projectionType == 1) camera->getPerspectiveMatrix(projection_transformation_matrix);
    else camera->getOrthographicMatrix(projection_transformation_matrix);


    // for each mesh in the scene, apply transformations
    for (auto mesh : this->meshes) {
        // for each triangle in the mesh apply transformations
        Matrix4 modeling_matrix = getIdentityMatrix();
        getTransformationMatrix(modeling_matrix, mesh);
        // multiply camera matrix with transformation matrix
        Matrix4 viewing_without_vp = multiplyMatrixWithMatrix(projection_transformation_matrix,
                                                              multiplyMatrixWithMatrix(camera_transformation_matrix,
                                                                                       modeling_matrix));
        for (auto triangle: mesh->triangles) {
            vector<Vec4> transformedPoints = getTransformedPoints(triangle, viewing_without_vp);
            vector<Color> colors = getColorsOfTriangle(triangle);
            if(this->cullingEnabled){
                // iterate over the transformed points and create a Vec3 from Vec4
                Vec3 first_edge = Vec3(transformedPoints[1].x - transformedPoints[0].x, transformedPoints[1].y - transformedPoints[0].y, transformedPoints[1].z - transformedPoints[0].z, 0);
                Vec3 second_edge = Vec3(transformedPoints[2].x - transformedPoints[0].x, transformedPoints[2].y - transformedPoints[0].y, transformedPoints[2].z - transformedPoints[0].z, 0);
                Vec3 crossProduct = crossProductVec3(first_edge, second_edge);
                Vec3 normal = normalizeVec3(crossProduct);
                if(dotProductVec3(normal, Vec3(transformedPoints[0].x, transformedPoints[0].y, transformedPoints[0].z,0)) < 0){
                    // if the dot product is negative, then the triangle is facing away from the camera
                    continue;
                }


            }
            if(mesh->type == 1){
                // this is part for solid
                rasterization(transformedPoints, colors, viewport_transformation_matrix, camera);
            }
            else if(mesh->type == 0){
                // this is part for wireframe
                //take the all points and colors
                Vec4 first_point = transformedPoints[0];
                Vec4 first_point2 = transformedPoints[0];
                Vec4 second_point = transformedPoints[1];
                Vec4 second_point2 = transformedPoints[1];
                Vec4 third_point = transformedPoints[2];
                Vec4 third_point2 = transformedPoints[2];
                Color first_color = colors[0];
                Color first_color2 = colors[0];
                Color second_color = colors[1];
                Color second_color2 = colors[1];
                Color third_color = colors[2];
                Color third_color2 = colors[2];

                bool first_second = clipping(first_point,second_point,first_color, second_color);
                bool second_third = clipping(second_point2,third_point,second_color2, third_color);
                bool third_first  = clipping(third_point2,first_point2,third_color2, first_color2);

                if (first_second) {
                    line_rasterization(first_point, second_point, first_color, second_color, viewport_transformation_matrix, camera);
                }
                if (second_third) {
                    line_rasterization(second_point2, third_point, second_color2, third_color, viewport_transformation_matrix, camera);
                }
                if (third_first) {
                    line_rasterization(third_point2, first_point2, third_color2, first_color2, viewport_transformation_matrix, camera);
                }

            }
        }


    }

}

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
    const char *str;
    XMLDocument xmlDoc;
    XMLElement *pElement;

    xmlDoc.LoadFile(xmlPath);

    XMLNode *pRoot = xmlDoc.FirstChild();

    // read background color
    pElement = pRoot->FirstChildElement("BackgroundColor");
    str = pElement->GetText();
    sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

    // read culling
    pElement = pRoot->FirstChildElement("Culling");
    if (pElement != NULL) {
        str = pElement->GetText();

        if (strcmp(str, "enabled") == 0) {
            cullingEnabled = true;
        }
        else {
            cullingEnabled = false;
        }
    }

    // read cameras
    pElement = pRoot->FirstChildElement("Cameras");
    XMLElement *pCamera = pElement->FirstChildElement("Camera");
    XMLElement *camElement;
    while (pCamera != NULL)
    {
        Camera *cam = new Camera();

        pCamera->QueryIntAttribute("id", &cam->cameraId);

        // read projection type
        str = pCamera->Attribute("type");

        if (strcmp(str, "orthographic") == 0) {
            cam->projectionType = 0;
        }
        else {
            cam->projectionType = 1;
        }

        camElement = pCamera->FirstChildElement("Position");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

        camElement = pCamera->FirstChildElement("Gaze");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

        camElement = pCamera->FirstChildElement("Up");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

        cam->gaze = normalizeVec3(cam->gaze);
        cam->u = crossProductVec3(cam->gaze, cam->v);
        cam->u = normalizeVec3(cam->u);

        cam->w = inverseVec3(cam->gaze);
        cam->v = crossProductVec3(cam->u, cam->gaze);
        cam->v = normalizeVec3(cam->v);

        camElement = pCamera->FirstChildElement("ImagePlane");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
               &cam->left, &cam->right, &cam->bottom, &cam->top,
               &cam->near, &cam->far, &cam->horRes, &cam->verRes);

        camElement = pCamera->FirstChildElement("OutputName");
        str = camElement->GetText();
        cam->outputFileName = string(str);

        cameras.push_back(cam);

        pCamera = pCamera->NextSiblingElement("Camera");
    }

    // read vertices
    pElement = pRoot->FirstChildElement("Vertices");
    XMLElement *pVertex = pElement->FirstChildElement("Vertex");
    int vertexId = 1;

    while (pVertex != NULL)
    {
        Vec3 *vertex = new Vec3();
        Color *color = new Color();

        vertex->colorId = vertexId;

        str = pVertex->Attribute("position");
        sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

        str = pVertex->Attribute("color");
        sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

        vertices.push_back(vertex);
        colorsOfVertices.push_back(color);

        pVertex = pVertex->NextSiblingElement("Vertex");

        vertexId++;
    }

    // read translations
    pElement = pRoot->FirstChildElement("Translations");
    XMLElement *pTranslation = pElement->FirstChildElement("Translation");
    while (pTranslation != NULL)
    {
        Translation *translation = new Translation();

        pTranslation->QueryIntAttribute("id", &translation->translationId);

        str = pTranslation->Attribute("value");
        sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

        translations.push_back(translation);

        pTranslation = pTranslation->NextSiblingElement("Translation");
    }

    // read scalings
    pElement = pRoot->FirstChildElement("Scalings");
    XMLElement *pScaling = pElement->FirstChildElement("Scaling");
    while (pScaling != NULL)
    {
        Scaling *scaling = new Scaling();

        pScaling->QueryIntAttribute("id", &scaling->scalingId);
        str = pScaling->Attribute("value");
        sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

        scalings.push_back(scaling);

        pScaling = pScaling->NextSiblingElement("Scaling");
    }

    // read rotations
    pElement = pRoot->FirstChildElement("Rotations");
    XMLElement *pRotation = pElement->FirstChildElement("Rotation");
    while (pRotation != NULL)
    {
        Rotation *rotation = new Rotation();

        pRotation->QueryIntAttribute("id", &rotation->rotationId);
        str = pRotation->Attribute("value");
        sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

        rotations.push_back(rotation);

        pRotation = pRotation->NextSiblingElement("Rotation");
    }

    // read meshes
    pElement = pRoot->FirstChildElement("Meshes");

    XMLElement *pMesh = pElement->FirstChildElement("Mesh");
    XMLElement *meshElement;
    while (pMesh != NULL)
    {
        Mesh *mesh = new Mesh();

        pMesh->QueryIntAttribute("id", &mesh->meshId);

        // read projection type
        str = pMesh->Attribute("type");

        if (strcmp(str, "wireframe") == 0) {
            mesh->type = 0;
        }
        else {
            mesh->type = 1;
        }

        // read mesh transformations
        XMLElement *pTransformations = pMesh->FirstChildElement("Transformations");
        XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

        while (pTransformation != NULL)
        {
            char transformationType;
            int transformationId;

            str = pTransformation->GetText();
            sscanf(str, "%c %d", &transformationType, &transformationId);

            mesh->transformationTypes.push_back(transformationType);
            mesh->transformationIds.push_back(transformationId);

            pTransformation = pTransformation->NextSiblingElement("Transformation");
        }

        mesh->numberOfTransformations = mesh->transformationIds.size();

        // read mesh faces
        char *row;
        char *clone_str;
        int v1, v2, v3;
        XMLElement *pFaces = pMesh->FirstChildElement("Faces");
        str = pFaces->GetText();
        clone_str = strdup(str);

        row = strtok(clone_str, "\n");
        while (row != NULL)
        {
            int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

            if (result != EOF) {
                mesh->triangles.push_back(Triangle(v1, v2, v3));
            }
            row = strtok(NULL, "\n");
        }
        mesh->numberOfTriangles = mesh->triangles.size();
        meshes.push_back(mesh);

        pMesh = pMesh->NextSiblingElement("Mesh");
    }
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
    if (this->image.empty())
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            vector<Color> rowOfColors;

            for (int j = 0; j < camera->verRes; j++)
            {
                rowOfColors.push_back(this->backgroundColor);
            }

            this->image.push_back(rowOfColors);
        }
    }
    else
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            for (int j = 0; j < camera->verRes; j++)
            {
                this->image[i][j].r = this->backgroundColor.r;
                this->image[i][j].g = this->backgroundColor.g;
                this->image[i][j].b = this->backgroundColor.b;
            }
        }
    }
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
    if (value >= 255.0)
        return 255;
    if (value <= 0.0)
        return 0;
    return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
    ofstream fout;

    fout.open(camera->outputFileName.c_str());

    fout << "P3" << endl;
    fout << "# " << camera->outputFileName << endl;
    fout << camera->horRes << " " << camera->verRes << endl;
    fout << "255" << endl;

    for (int j = camera->verRes - 1; j >= 0; j--)
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
                 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
                 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
        }
        fout << endl;
    }
    fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
    string command;

    // call command on Ubuntu
    if (osType == 1)
    {
        command = "convert " + ppmFileName + " " + ppmFileName + ".png";
        system(command.c_str());
    }

        // call command on Windows
    else if (osType == 2)
    {
        command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
        system(command.c_str());
    }

        // default action - don't do conversion
    else
    {
    }
}

void Scene::multiplyVec4WithScalar(Vec4 &vec4, double d) {
    vec4.x *= d;
    vec4.y *= d;
    vec4.z *= d;
    vec4.t *= d;
}

void Scene::line_rasterization(Vec4 first_point, Vec4 second_point, Color first_color, Color second_color, Matrix4 viewport_matrix4, Camera *pCamera) {
    Vec4 first_point_viewport = multiplyMatrixWithVec4(viewport_matrix4, first_point);
    Vec4 second_point_viewport = multiplyMatrixWithVec4(viewport_matrix4, second_point);

    // calculate the slope of the line
    double slope = (second_point_viewport.y - first_point_viewport.y) / (second_point_viewport.x - first_point_viewport.x);
    int x;
    int y;
    double x1_x0 = abs(second_point_viewport.x - first_point_viewport.x);
    double y1_y0 = abs(second_point_viewport.y - first_point_viewport.y);
    if(slope >= 0){
        x = min(first_point_viewport.x, second_point_viewport.x);
        y = min(first_point_viewport.y, second_point_viewport.y);
    }
    else if(slope < -1){
        x = max(first_point_viewport.x, second_point_viewport.x);
        y = min(first_point_viewport.y, second_point_viewport.y);
    }
    else {
        x = min(first_point_viewport.x, second_point_viewport.x);
        y = max(first_point_viewport.y, second_point_viewport.y);
    }
    if(slope < 1 && slope >= 0){
        double d = 2 * y1_y0 - x1_x0;
        double dE = 2 * y1_y0;
        double dNE = 2 * (y1_y0 - x1_x0);
        while(x <= max(first_point_viewport.x, second_point_viewport.x)){
            image[x][y].r = round((first_color.r * abs(x - second_point_viewport.x) + second_color.r * abs(x - first_point_viewport.x)) / x1_x0);
            image[x][y].g = round((first_color.g * abs(x - second_point_viewport.x) + second_color.g * abs(x - first_point_viewport.x)) / x1_x0);
            image[x][y].b = round((first_color.b * abs(x - second_point_viewport.x) + second_color.b * abs(x - first_point_viewport.x)) / x1_x0);
            if(d <= 0){
                d += dE;
                x++;
            }
            else{
                d += dNE;
                x++;
                y++;
            }
        }
    }
    else if(slope >= 1){
        double d = 2 * x1_x0 - y1_y0;
        double dE = 2 * x1_x0;
        double dNE = 2 * (x1_x0 - y1_y0);
        while(y <= max(first_point_viewport.y, second_point_viewport.y)){
            image[x][y].r = round((first_color.r * abs(y - second_point_viewport.y) + second_color.r * abs(y - first_point_viewport.y)) / y1_y0);
            image[x][y].g = round((first_color.g * abs(y - second_point_viewport.y) + second_color.g * abs(y - first_point_viewport.y)) / y1_y0);
            image[x][y].b = round((first_color.b * abs(y - second_point_viewport.y) + second_color.b * abs(y - first_point_viewport.y)) / y1_y0);
            if(d <= 0){
                d += dE;
                y++;
            }
            else{
                d += dNE;
                x++;
                y++;
            }
        }
    }
    else if(slope < -1){
        double d = 2 * x1_x0 - y1_y0;
        double dE = 2 * x1_x0;
        double dNE = 2 * (x1_x0 - y1_y0);
        while(y < max(first_point_viewport.y, second_point_viewport.y)){
            image[x][y].r = round((first_color.r * abs(y - second_point_viewport.y) + second_color.r * abs(y - first_point_viewport.y)) / y1_y0);
            image[x][y].g = round((first_color.g * abs(y - second_point_viewport.y) + second_color.g * abs(y - first_point_viewport.y)) / y1_y0);
            image[x][y].b = round((first_color.b * abs(y - second_point_viewport.y) + second_color.b * abs(y - first_point_viewport.y)) / y1_y0);
            if(d <= 0){
                d += dE;
                y++;
            }
            else{
                d += dNE;
                x--;
                y++;
            }
        }
    }
    else{
        double d = 2 * y1_y0 - x1_x0;
        double dE = 2 * y1_y0;
        double dNE = 2 * (y1_y0 - x1_x0);
        while(x < max(first_point_viewport.x, second_point_viewport.x)){
            image[x][y].r = round((first_color.r * abs(x - second_point_viewport.x) + second_color.r * abs(x - first_point_viewport.x)) / x1_x0);
            image[x][y].g = round((first_color.g * abs(x - second_point_viewport.x) + second_color.g * abs(x - first_point_viewport.x)) / x1_x0);
            image[x][y].b = round((first_color.b * abs(x - second_point_viewport.x) + second_color.b * abs(x - first_point_viewport.x)) / x1_x0);
            if(d <= 0){
                d += dE;
                x++;
            }
            else{
                d += dNE;
                x++;
                y--;
            }
        }
    }
}
