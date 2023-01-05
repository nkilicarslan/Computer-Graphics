#ifndef __HW1__PARSER__
#define __HW1__PARSER__

#include <string>
#include <vector>

namespace parser
{
    //Notice that all the structures are as simple as possible
    //so that you are not enforced to adopt any style or design.
    struct Vec3f
    {
        float x, y, z;


    };

    inline Vec3f operator +(const Vec3f& vec1,const Vec3f& vec2){
        Vec3f tmpVec;
        tmpVec.x = vec1.x+vec2.x;
        tmpVec.y = vec1.y+vec2.y;
        tmpVec.z = vec1.z+vec2.z;

        return tmpVec;
    }

    inline Vec3f operator -(const Vec3f& vec1,const Vec3f& vec2){
        Vec3f tmpVec;
        tmpVec.x = vec1.x-vec2.x;
        tmpVec.y = vec1.y-vec2.y;
        tmpVec.z = vec1.z-vec2.z;

        return tmpVec;
    }



    inline Vec3f operator *(const Vec3f& vec1,const float &number){
        Vec3f tmpVec;
        tmpVec.x = vec1.x*number;
        tmpVec.y = vec1.y*number;
        tmpVec.z = vec1.z*number;

        return tmpVec;
    }

    inline Vec3f operator *(const Vec3f& vec1,const Vec3f &vec2){
        Vec3f tmpVec;
        tmpVec.x = vec1.x*vec2.x;
        tmpVec.y = vec1.y*vec2.y;
        tmpVec.z = vec1.z*vec2.z;

        return tmpVec;
    }

    inline Vec3f operator /(const Vec3f& vec1,const float &number){
        Vec3f tmpVec;
        tmpVec.x = vec1.x/number;
        tmpVec.y = vec1.y/number;
        tmpVec.z = vec1.z/number;

        return tmpVec;
    }


    struct Vec3i
    {
        int x, y, z;
    };

    inline Vec3f operator *(const Vec3i& vec1,const Vec3f &vec2){
        Vec3f tmpVec;
        tmpVec.x = vec1.x*vec2.x;
        tmpVec.y = vec1.y*vec2.y;
        tmpVec.z = vec1.z*vec2.z;

        return tmpVec;
    }

    struct Vec4f
    {
        float x, y, z, w;
    };

    struct Camera
    {
        Vec3f position;
        Vec3f gaze;
        Vec3f up;
        Vec4f near_plane;
        float near_distance;
        int image_width, image_height;
        std::string image_name;
    };

    struct PointLight
    {
        Vec3f position;
        Vec3f intensity;
    };

    struct Material
    {
        bool is_mirror;
        Vec3f ambient;
        Vec3f diffuse;
        Vec3f specular;
        Vec3f mirror;
        float phong_exponent;
    };

    struct Face
    {
        int v0_id;
        int v1_id;
        int v2_id;
    };

    struct Mesh
    {
        int material_id;
        std::vector<Face> faces;
    };

    struct Triangle
    {
        int material_id;
        Face indices;
    };

    struct Sphere
    {
        int material_id;
        int center_vertex_id;
        float radius;
    };

    struct Scene
    {
        //Data
        Vec3i background_color;
        float shadow_ray_epsilon;
        int max_recursion_depth;
        std::vector<Camera> cameras;
        Vec3f ambient_light;
        std::vector<PointLight> point_lights;
        std::vector<Material> materials;
        std::vector<Vec3f> vertex_data;
        std::vector<Mesh> meshes;
        std::vector<Triangle> triangles;
        std::vector<Sphere> spheres;

        //Functions
        void loadFromXml(const std::string &filepath);
    };
}

#endif