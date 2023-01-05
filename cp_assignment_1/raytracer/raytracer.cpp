#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include "vector"
#include <chrono>
using namespace std;
using namespace parser;

struct Ray{
    Vec3f origin;
    Vec3f direction;
};


struct Hit {
    Material material;
    Vec3f pointofIntersect;
    Vec3i pixel;
    Vec3f normal;
};


Vec3f crossProduct(Vec3f const &firstVec, Vec3f const &secondVec){
    Vec3f tmpVec;

    tmpVec.x = firstVec.y*secondVec.z-secondVec.y*firstVec.z;
    tmpVec.y = secondVec.x*firstVec.z-firstVec.x*secondVec.z;
    tmpVec.z = firstVec.x*secondVec.y-secondVec.x*firstVec.y;

    return tmpVec;

}


double dotProduct(Vec3f const &firstVec, Vec3f const secondVec){

    Vec3f tmpVec;
    tmpVec.x = firstVec.x* secondVec.x;
    tmpVec.y = firstVec.y* secondVec.y;
    tmpVec.z = firstVec.z* secondVec.z;

    return tmpVec.x+tmpVec.y+tmpVec.z;
}


void vectorNormalizer(Vec3f &vec){


    double length = sqrt(pow(vec.x,2)+pow(vec.y,2)+pow(vec.z,2));
    vec.x = vec.x/length;
    vec.y = vec.y/length;
    vec.z = vec.z/length;

}


Ray rayCreater(const Camera &camera, int i, int j){

    double left = camera.near_plane.x;
    float right = camera.near_plane.y;
    float top = camera.near_plane.w;
    float bottom = camera.near_plane.z;

    Vec3f m = camera.position + camera.gaze*camera.near_distance;
    Vec3f u = crossProduct(camera.up,camera.gaze*(-1));
    vectorNormalizer(u);
    Vec3f q = m + u*left + camera.up*top;

    double su = (i+0.5)*(right-left)/camera.image_width;
    double sv = (j+0.5)*(top-bottom)/camera.image_height;

    Vec3f s = q + u*su - camera.up*sv;

    Ray ray;
    ray.origin = camera.position;
    ray.direction = s - camera.position;
    vectorNormalizer(ray.direction);

    return  ray;
}


float intersectionSphere(Scene const &scene, Ray const &ray, Sphere const &sphere) {
    double radius = sphere.radius;
    Vec3f center = scene.vertex_data[sphere.center_vertex_id - 1];

    float A = dotProduct(ray.direction, ray.direction);
    float B = 2 * dotProduct(ray.direction, ray.origin - center);
    float C = dotProduct(ray.origin - center, ray.origin - center) - pow(radius, 2);

    float discriminant = pow(B, 2) - 4 * A * C;

    if (discriminant < 0) return -1;
    else {
        float t1 = (-1*B + sqrtf(discriminant)) / (2 * A);
        float t2 = (-1*B - sqrtf(discriminant)) / (2 * A);
        float min = t1 < t2 ? t1 : t2;
        return min;
    }
}


float intersectionPointLight(Scene const &scene, Ray const &ray, Vec3f center){
    double radius = pow(10,-1);

    float A = dotProduct(ray.direction, ray.direction);
    float B = 2 * dotProduct(ray.direction, ray.origin - center);
    float C = dotProduct(ray.origin - center, ray.origin - center) - pow(radius, 2);

    float discriminant = pow(B, 2) - 4 * A * C;

    if (discriminant < 0) return -1;
    else {
        float t1 = (-1*B + sqrtf(discriminant)) / (2 * A);
        float t2 = (-1*B - sqrtf(discriminant)) / (2 * A);
        float min = t1 < t2 ? t1 : t2;
        return min;
    }
}


float intersectionTriangle(Scene const &scene, Ray const &ray, Triangle const &triangle){

    Vec3f triangleA = scene.vertex_data[triangle.indices.v0_id-1];
    Vec3f triangleB = scene.vertex_data[triangle.indices.v1_id-1];
    Vec3f triangleC = scene.vertex_data[triangle.indices.v2_id-1];

    float ax_ox = triangleA.x - ray.origin.x;
    float ax_cx = triangleA.x - triangleC.x;
    float dx = ray.direction.x;

    float ay_oy = triangleA.y - ray.origin.y;
    float ay_cy = triangleA.y - triangleC.y;
    float dy = ray.direction.y;

    float az_oz = triangleA.z - ray.origin.z;
    float az_cz = triangleA.z - triangleC.z;
    float dz = ray.direction.z;

    float ax_bx = triangleA.x - triangleB.x;
    float ay_by = triangleA.y - triangleB.y;
    float az_bz = triangleA.z - triangleB.z;

    float beta_nominator = ax_ox * (ay_cy * dz - dy * az_cz) - ax_cx * (ay_oy * dz - dy * az_oz) + dx * (ay_oy * az_cz - ay_cy * az_oz);
    float gama_nominator = ax_bx * (ay_oy * dz - dy * az_oz) - ax_ox * (ay_by * dz - dy * az_bz) + dx * (ay_by * az_oz - ay_oy * az_bz);
    float t_nominator = ax_bx * (ay_cy * az_oz - ay_oy * az_cz) - ax_cx * (ay_by * az_oz - ay_oy * az_bz) + ax_ox * (ay_by * az_cz - ay_cy * az_bz);
    float A = ax_bx * (ay_cy * dz - dy * az_cz) - ax_cx * (ay_by * dz - dy * az_bz) + dx * (ay_by * az_cz - ay_cy * az_bz);

    float beta = beta_nominator / A;
    float gama = gama_nominator / A;
    float t = t_nominator / A;

    double epsilon = pow(10, -6);
    if (beta >= 0 - epsilon && gama >= 0 - epsilon && (beta + gama) <= 1 + epsilon) return t;
    else return -1;

}


Vec3f calculateIrradience(Hit hit, PointLight pointLight){

    Vec3f irradience;
    Vec3f value_before_length = pointLight.position - hit.pointofIntersect;
    Vec3f tmp = value_before_length * value_before_length;
    float r = tmp.x + tmp.y + tmp.z;

    if(r !=0){
        irradience.x = pointLight.intensity.x/r;
        irradience.y = pointLight.intensity.y/r;
        irradience.z = pointLight.intensity.z/r;
    }

    return irradience;
}


Vec3f calculateDiffuse(Hit hit, PointLight pointLight, Vec3f irradiance){

    Vec3f W_i = pointLight.position - hit.pointofIntersect;
    vectorNormalizer(W_i);
    double cos_teta = dotProduct(hit.normal,W_i);
    double epsilon = pow(10, -6);
    if(cos_teta < 0 + epsilon) cos_teta = 0;

    return hit.material.diffuse * (irradiance * cos_teta);
}


Vec3f calculateSpecular(Hit hit, Ray ray, Vec3f irradiance, PointLight pointLight){

    Vec3f W_i = pointLight.position - hit.pointofIntersect;
    vectorNormalizer(W_i);
    Vec3f W_o = ray.origin - hit.pointofIntersect;
    vectorNormalizer(W_o);
    Vec3f h = W_i + W_o;
    vectorNormalizer(h);

    double cos_alpha_teta = dotProduct(hit.normal,h);
    double epsilon = pow(10, -6);
    if(cos_alpha_teta < 0 + epsilon) cos_alpha_teta = 0;
    double tmp = pow(cos_alpha_teta,hit.material.phong_exponent);

    return hit.material.specular * (irradiance*tmp);
}


int detectShadow(Scene const &scene, PointLight const &pointLight, Vec3f const &intersectionPoint, Hit const &hit) {

    Ray shadow;
    shadow.direction = pointLight.position - intersectionPoint;
    vectorNormalizer(shadow.direction);
    shadow.origin = intersectionPoint + hit.normal * scene.shadow_ray_epsilon;

    vector<Sphere> allSpheres = scene.spheres;
    vector<Triangle> allTriangles = scene.triangles;
    vector<Mesh> allMeshes = scene.meshes;

    float tLight = intersectionPointLight(scene, shadow, pointLight.position);

    for (int i = 0; i < allSpheres.size(); ++i) {
        double t = intersectionSphere(scene, shadow, allSpheres[i]);
        if(t >= 0 && t<= tLight){
            return 1;
        }
    }

    for (int i = 0; i < allTriangles.size(); ++i) {
        double t = intersectionTriangle(scene, shadow, allTriangles[i]);
        if(t>= 0 && t <=tLight){
            return 1;
        }
    }

    for (int i = 0; i < allMeshes.size(); i++) {
        for (int j = 0; j < allMeshes[i].faces.size(); ++j) {
            Triangle tmp_triangle;
            tmp_triangle.indices = allMeshes[i].faces[j];
            double t = intersectionTriangle(scene, shadow, tmp_triangle);
            if(t>= 0 && t <= tLight){
                return 1;
            }
        }
    }
    return -1;
}


Ray detectMirror(Scene const &scene, Ray const &ray, Hit const &hit){

    Ray result;
    Vec3f w0_direction = ray.direction * -1;
    float dot = dotProduct(hit.normal, w0_direction);
    Vec3f w_r_direction = ray.direction + (hit.normal * 2) * dot;
    vectorNormalizer(w_r_direction);
    result.origin = hit.pointofIntersect + hit.normal * scene.shadow_ray_epsilon;
    result.direction = w_r_direction;
    return result;
}


Hit sendRayToObjects(int recursion_number, Scene const &scene, Ray const &ray){

    Hit hit;
    hit.pixel = scene.background_color;
    float t = -1;
    vector<Sphere> allSpheres = scene.spheres;
    vector<Triangle> allTriangles = scene.triangles;
    vector<Mesh> allMeshes = scene.meshes;
    int indexSphere = 0, indexTriangle = 0, indexMesh = 0;
    int selector;
    int faceIndex = 0;
    for (int i = 0; i < allSpheres.size(); ++i) {
        float tmp_val = intersectionSphere(scene,ray,allSpheres[i]);
        if(t < 0 && tmp_val > 0){
            t = tmp_val;
            selector = 0;
            indexSphere = i;
        }
        else if(tmp_val > 0 && t>0 && tmp_val < t){
            t = min(t,tmp_val);
            selector = 0;
            indexSphere = i;
        }
    }
    for (int i = 0; i < allTriangles.size(); ++i) {
        float tmp_val = intersectionTriangle(scene,ray,allTriangles[i]);
        if(t < 0 && tmp_val > 0){
            t = tmp_val;
            selector = 1;
            indexTriangle = i;
        }
        else if(tmp_val > 0 && t>0 && tmp_val < t){
            t = min(t,tmp_val);
            selector = 1;
            indexTriangle = i;
        }
    }
    for (int i = 0; i < allMeshes.size(); ++i) {
        for (int j = 0; j < allMeshes[i].faces.size(); ++j) {
            Triangle tmp_triangle;
            tmp_triangle.indices.v0_id = allMeshes[i].faces[j].v0_id;
            tmp_triangle.indices.v1_id = allMeshes[i].faces[j].v1_id;
            tmp_triangle.indices.v2_id = allMeshes[i].faces[j].v2_id;
            float tmp_val = intersectionTriangle(scene, ray, tmp_triangle);
            if (t < 0 && tmp_val > 0) {
                t = tmp_val;
                selector = 2;
                indexMesh = i;
                faceIndex = j;
            } else if (tmp_val > 0 && t > 0 && tmp_val < t) {
                t = min(t, tmp_val);
                selector = 2;
                indexMesh = i;
                faceIndex = j;
            }
        }
    }
    if(recursion_number > 0 && t < 0) {
        hit.pixel.x = 0;
        hit.pixel.y = 0;
        hit.pixel.z = 0;
        return hit;
    }
    if(t < 0) return hit;   // base case

    // sphere intersection
    if(selector == 0){  // selector:0 sphere

        hit.material = scene.materials[allSpheres[indexSphere].material_id-1];
        hit.pointofIntersect = ray.origin + ray.direction*t;
        hit.normal = hit.pointofIntersect - scene.vertex_data[allSpheres[indexSphere].center_vertex_id - 1];
        vectorNormalizer(hit.normal);
    }

    // triangle intersection
    if(selector == 1){  // selector:1 triangle

        hit.material = scene.materials[allTriangles[indexTriangle].material_id-1];
        hit.pointofIntersect = ray.origin + ray.direction*t;
        Vec3f triangle_a = scene.vertex_data[allTriangles[indexTriangle].indices.v0_id-1];
        Vec3f triangle_b = scene.vertex_data[allTriangles[indexTriangle].indices.v1_id-1];
        Vec3f triangle_c = scene.vertex_data[allTriangles[indexTriangle].indices.v2_id-1];
        hit.normal = crossProduct(triangle_b-triangle_a,triangle_c-triangle_a);
        vectorNormalizer(hit.normal);
    }

    // mesh intersection
    if(selector == 2){  // selector:2 mesh

        Triangle tmp_triangle;
        tmp_triangle.indices = allMeshes[indexMesh].faces[faceIndex];
        hit.material = scene.materials[allMeshes[indexMesh].material_id-1];
        hit.pointofIntersect = ray.origin + ray.direction*t;
        Vec3f triangle_a = scene.vertex_data[allMeshes[indexMesh].faces[faceIndex].v0_id-1];
        Vec3f triangle_b = scene.vertex_data[allMeshes[indexMesh].faces[faceIndex].v1_id-1];
        Vec3f triangle_c = scene.vertex_data[allMeshes[indexMesh].faces[faceIndex].v2_id-1];
        hit.normal = crossProduct(triangle_b-triangle_a,triangle_c-triangle_a);
        vectorNormalizer(hit.normal);
    }

    Vec3f color = {0,0,0};
    for (int j = 0; j < scene.point_lights.size(); ++j) {
        int check_shadow = detectShadow(scene,scene.point_lights[j],hit.pointofIntersect,hit);
        if(check_shadow == -1) {
            Vec3f irradiance = calculateIrradience(hit, scene.point_lights[j]);
            Vec3f diffuse = calculateDiffuse(hit, scene.point_lights[j], irradiance);
            Vec3f specular = calculateSpecular(hit, ray, irradiance, scene.point_lights[j]);
            color = color + diffuse + specular;
        }
    }

    color = color + (hit.material.ambient * scene.ambient_light);
    if (hit.material.is_mirror && recursion_number < scene.max_recursion_depth) {
        Ray mirrorRay = detectMirror(scene, ray, hit);
        Hit mirrored_hit = sendRayToObjects((recursion_number + 1), scene, mirrorRay);
        color = color + (mirrored_hit.pixel * hit.material.mirror);
    }

    hit.pixel.x = int(round(color.x));
    hit.pixel.y = int(round(color.y));
    hit.pixel.z = int(round(color.z));

    return hit;
}


int main(int argc, char* argv[])
{
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    for (int i = 0; i < scene.cameras.size(); ++i) {

        int width = scene.cameras[i].image_width;
        int height = scene.cameras[i].image_height;

        unsigned char* image = new unsigned char [width * height * 3];

        int j = 0;
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                Ray ray = rayCreater(scene.cameras[i], x, y);
                Hit hit = sendRayToObjects(0, scene,ray);
                image[j++] = min(max(hit.pixel.x, 0), 255);
                image[j++] = min(max(hit.pixel.y, 0), 255);
                image[j++] = min(max(hit.pixel.z, 0), 255);
            }
        }
        write_ppm(scene.cameras[i].image_name.c_str(), image, width, height);
    }

    return 0;
}