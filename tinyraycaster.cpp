#include <limits>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "geometry.h"

struct Light {
    Light(const Vec3f &p, const float &i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

struct Material {
    Material(const float &r, const Vec4f &a, const Vec3f &color, const float &spec) : refractive_index(r), albedo(a), diffuse(color), specular(spec) {}
    Material() : refractive_index(), albedo(1,0,0,0), diffuse(), specular() { }
    float refractive_index;
    Vec4f albedo;
    Vec3f diffuse;
    float specular;
};

struct Sphere {
    Sphere(const Vec3f &c, const float &r, const Material &m) : center(c), radius(r), material(m) {}

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const {
        Vec3f L = center - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0       = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }

    Material material;
    Vec3f center;
    float radius;
};

Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N*2.f*(I*N);
}

Vec3f refract(const Vec3f &I, const Vec3f &N, const float &refractive_index) { // Snell's law
    float cosi = - std::max(-1.f, std::min(1.f, I*N));
    float etai = 1, etat = refractive_index;
    Vec3f n = N;
    if (cosi < 0) {                      // We need to handle with care the two possible situations:
        cosi = -cosi;                    //    - When the ray is inside the object
        std::swap(etai, etat); n = -N;   //    - When the ray is outside.
    }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? Vec3f(0,0,0) : I*eta + n*(eta * cosi - sqrtf(k));
}

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, Vec3f &hit, Vec3f &N, Material &material) {
    float spheres_dist = std::numeric_limits<float>::max();
    for (size_t i=0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].intersect(orig, dir, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = orig + dir*dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }

    float checkerboard_dist = std::numeric_limits<float>::max();
    if (fabs(dir.y)>1e-3)  {
        float d = -(orig.y+4)/dir.y; // the checkerboard plane has equation y = -4
        Vec3f pt = orig + dir*d;
        if (d>0 && fabs(pt.x)<10 && pt.z<-10 && pt.z>-30 && d<spheres_dist)  {
            checkerboard_dist = d;
            hit = pt;
            N = Vec3f(0,1,0);
            material.diffuse = (int(.3*hit.x+1000) + int(.3*hit.z)) & 1 ? Vec3f(1,1,1) : Vec3f(1, .7, .3);
            material.diffuse = material.diffuse*.3;
        }
    }

    return std::min(spheres_dist, checkerboard_dist)<1000;
}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<Light> &lights, size_t depth) {
    Vec3f point, N;
    Material material;

    if (depth > 4 || !scene_intersect(orig, dir, spheres, point, N, material))
        return Vec3f(0.2, 0.7, 0.8); // background color

    Vec3f reflectionDirection = reflect(dir, N).normalize();
    Vec3f refractionDirection = refract(dir, N, material.refractive_index).normalize();
//  Vec3f reflectionRayOrig = reflectionDirection*N < 0 ? point - N*1e-3 : point + N*1e-3;
//  Vec3f refractionRayOrig = refractionDirection*N < 0 ? point - N*1e-3 : point + N*1e-3;
//  Vec3f reflectionColor = cast_ray(reflectionRayOrig, reflectionDirection, spheres, lights, depth + 1);
//  Vec3f refractionColor = cast_ray(refractionRayOrig, refractionDirection, spheres, lights, depth + 1);

    Vec2f reflectionColor = cast_ray(point, reflectionDirection, spheres, lights, depth + 1);
    Vec3f refractionColor = cast_ray(point, refractionDirection, spheres, lights, depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i=0; i<lights.size(); i++) {
        Vec3f light_dir      = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();

        Vec3f test_point = dir*N < 0 ? point + N*1e-3 : point - N*1e-3; // offset the original point to avoid occlusion by the object itself
        Vec3f tmppoint, tmpN;
        Material tmpmaterial;
        if (scene_intersect(test_point, light_dir, spheres, tmppoint, tmpN, tmpmaterial) && (tmppoint-test_point).norm() < light_distance)
            continue;

        diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir*N);
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*dir), material.specular)*lights[i].intensity;
    }
    return material.diffuse * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1] + reflectionColor*material.albedo[2] + refractionColor*material.albedo[3];
}


void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
    const int width  = 1024;
    const int height = 768;
    const int fov    = M_PI/2.;

    std::vector<Vec3f> framebuffer(width*height);

    float scale = tan(fov/2.);
    float aspect_ratio = width/(float)height;

    Vec3f orig(0,0,0);
    for (size_t j = 0; j < height; ++j) {
        for (size_t i = 0; i < width; ++i) {
            float x = (2 * (i + 0.5) / (float)width - 1) * aspect_ratio * scale;
            float y = (1 - 2 * (j + 0.5) / (float)height) * scale;
            Vec3f dir = Vec3f(x, y, -1).normalize();
            framebuffer[i+j*width] = cast_ray(orig, dir, spheres, lights,  0);
        }
    }

    std::ofstream ofs; // save framebuffer to file
    ofs.open("./out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c/max;
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

int main(int argc, char **argv) {
    std::vector<Sphere> spheres;
    std::vector<Light>  lights;

    Material      ivory(1.0, Vec4f(0.7,  0.3, 0.10, 0.0), Vec3f(0.4, 0.4, 0.3),   50.);
    Material      glass(1.5, Vec4f(0.0,  0.5, 0.15, 0.7), Vec3f(0.6, 0.7, 0.8),  125.);
    Material red_rubber(1.0, Vec4f(0.8,  0.2, 0.00, 0.0), Vec3f(0.9, 0.3, 0.3),   25.);
    Material     mirror(1.0, Vec4f(0.0, 10.0, 0.85, 0.0), Vec3f(1.0, 1.0, 1.0), 1425.);

    spheres.push_back(Sphere(Vec3f(-3  ,  0  , -15), 2,   ivory));
    spheres.push_back(Sphere(Vec3f( 2  , -1  , -16), 2.5, red_rubber));
    spheres.push_back(Sphere(Vec3f(-0.5, -1.5, -11), 2,   glass));
    spheres.push_back(Sphere(Vec3f( 6,    4,   -15), 3,    mirror));

    lights.push_back(Light(Vec3f(-20, 20,  20), 1.5));
    lights.push_back(Light(Vec3f( 30, 50, -12), 1.8));
    lights.push_back(Light(Vec3f( 30, 20,  30), 1.7));

    render(spheres, lights);

    return 0;
}

