#include <limits>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>

#include "geometry.h"

const int width  = 1024;
const int height = 768;
const int fov    = M_PI/2.;
const Vec3f background_color = Vec3f(0.235294, 0.67451, 0.843137);
const int max_depth = 4;

struct Light {
    Light(const Vec3f &p, const float &i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

struct Material {
    Material(const float &r, const float &kd, const float &ks, const float &kr, const float &kt, const Vec3f &color, const float &spec) : ior(r), Kd(kd), Ks(ks), Kr(kr), Kt(kt), diffuse(color), specular(spec) {}
    Material() : ior(0), Kd(1), Ks(0), Kr(0), Kt(0), diffuse(Vec3f(0.6, 0.7, 0.8)), specular(25.) {}
    float ior;  // refractive index
    float Kd;   // diffuse albedo
    float Ks;   // specular albedo
    float Kr;   // reflective albedo
    float Kt;   // refractive albedo
    Vec3f diffuse;
    float specular;
};

struct Sphere {
    Sphere(const Vec3f &c, const float &r, const Material &m) : center(c), radius(r), material(m) {}

    bool solve_quadratic_equation(const float &a, const float &b, const float &c, float &x0, float &x1) const {
        float discr = b*b - 4*a*c;
        if (discr < 0) return false;
        else if (fabs(discr) < 1e-3) x0 = x1 = -0.5*b/a;
        else {
            float q = (b > 0) ?            // The well known solution (-b +- sqrt(b^2 - 4ac)) / 2a is known to be non-robust
                -0.5*(b + sqrtf(discr)) :  // in computation when ac is very small compered to b^2,
                -0.5*(b - sqrtf(discr));   // because one is subtracting two very similar values.
            x0 = q / a;                    // It is better to use the lesser known solution 2c / (-b -+ sqrt(b^2 -4ac)) for the other root.
            x1 = c / q;                    // here the computation of q ensures that we are not subtracting two similar values.
        }
        if (x0 > x1) std::swap(x0, x1);
        return true;
    }

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &tnear) const {
        Vec3f L = orig - center;
        float a = dir*dir;
        float b = 2*(dir*L);
        float c = L*L - powf(radius, 2);
        float t0, t1;
        if (!solve_quadratic_equation(a, b, c, t0, t1)) return false;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        tnear = t0;
        return true;
    }

    Material material;
    Vec3f center;
    float radius;
};

template <typename t> t clamp(const t& low, const t& high, const t& value) {
    return value < low ? low : (value > high ? high : value);
}

// Compute reflection direction
Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N*2.f*(I*N);
}

// Compute refraction direction using Snell's law
Vec3f refract(const Vec3f &I, const Vec3f &N, const float &ior) {
    float cosi = -clamp(-1.f, 1.f, I*N);
    float etai = 1, etat = ior;
    Vec3f n = N;
    if (cosi < 0) {                      // We need to handle with care the two possible situations:
        cosi = -cosi;                    //    - When the ray is inside the object
        std::swap(etai, etat); n = -N;   //    - When the ray is outside.
    }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? Vec3f(0,0,0) : I*eta + n*(eta * cosi - sqrtf(k));
}

// Compute Fresnel equation (the amount of reflected light)
float fresnel(const Vec3f &I, const Vec3f &N, const float &ior) {
    float cosi = clamp(-1.f, 1.f, I*N);
    float etai = 1, etat = ior;
    if (cosi > 0) std::swap(etai, etat);
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi)); // Compute sini using Snell's law
    if (sint >= 1) {
        return 1;  // Total internal reflection
    } else {
        float cost = sqrtf(std::max(0.f, 1 - sint * sint));
        cosi = fabsf(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        return (Rs * Rs + Rp * Rp) / 2;
    }
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
        float d = -(orig.y+5)/dir.y;
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

    if (depth > max_depth || !scene_intersect(orig, dir, spheres, point, N, material)) {
        return background_color;
    }

    Vec3f reflectionDirection = reflect(dir, N).normalize();
    Vec3f refractionDirection = refract(dir, N, material.ior).normalize();
    Vec3f reflectionRayOrig = reflectionDirection*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vec3f refractionRayOrig = refractionDirection*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vec3f reflectionColor = cast_ray(reflectionRayOrig, reflectionDirection, spheres, lights, depth + 1);
    Vec3f refractionColor = cast_ray(refractionRayOrig, refractionDirection, spheres, lights, depth + 1);
    float kr = fresnel(dir, N, material.ior);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i=0; i<lights.size(); i++) {
        Vec3f light_dir      = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();

        {
            Vec3f test_point = dir*N < 0 ? point + N*1e-3 : point - N*1e-3; // offset the original point to avoid occlusion by the object itself
            Vec3f tmppoint, tmpN;
            Material tmpmaterial;
            if (scene_intersect(test_point, light_dir, spheres, tmppoint, tmpN, tmpmaterial) && (tmppoint-test_point).norm() < light_distance) {
                continue;
            }
        }

        diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir*N);
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*dir), material.specular)*lights[i].intensity;
    }
    return material.diffuse * diffuse_light_intensity * material.Kd + Vec3f(1., 1., 1.)*specular_light_intensity * material.Ks + reflectionColor*material.Kr*kr + refractionColor*material.Kt*(1-kr);
}


void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
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

    // save framebuffer to file
    std::ofstream ofs;
    ofs.open("./out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c/max;
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * clamp(0.f, 1.f, framebuffer[i][j]));
        }
    }
    ofs.close();
}

int main(int argc, char **argv) {
    std::vector<Sphere> spheres;
    std::vector<Light> lights;

    Material ivory(0, .7, 0.3, 0.03, 0., Vec3f(1., 1., .8)*.4, 50.);
    Material glass(1.5, 0, .5, .75, .8, Vec3f(0.6, 0.7, 0.8), 125.);
    Material red_rubber(0, 0.8, 0.2, 0., 0, Vec3f(0.9, 0.3, 0.3), 25.);
    Material mirror(0, 0.0, 10., .85, 0, Vec3f(1,1,1), 1425.);

    spheres.push_back(Sphere(Vec3f(-2  ,  0  , -12), 2  , ivory));
    spheres.push_back(Sphere(Vec3f(2  ,  -2  , -12), 1.5  , red_rubber));
    spheres.push_back(Sphere(Vec3f(-0.5, -1, -8 ), 1.5, glass));
    spheres.push_back(Sphere(Vec3f(5, 2, -15 ), 3, mirror));

    lights.push_back(Light(Vec3f(-20, 20,  20), 1.5));
    lights.push_back(Light(Vec3f( 30, 50, -12), 1.8));
    lights.push_back(Light(Vec3f( 30, 20, 30),   1.7));

    render(spheres, lights);

    return 0;
}

