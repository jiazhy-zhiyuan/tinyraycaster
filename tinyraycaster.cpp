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
    enum MaterialType { SLIGHTLY_GLOSSY, REFLECTION_REFRACTION, REFLECTION };
    Material(const MaterialType &mt, const float &r, const float &kd, const float &ks, const Vec3f &color, const float &spec) : type(mt), ior(r), Kd(kd), Ks(ks), diffuse(color), specular(spec) {}
    Material() : type(SLIGHTLY_GLOSSY), ior(1.5), Kd(0.8), Ks(0.2), diffuse(Vec3f(0.6, 0.7, 0.8)), specular(25.) {}
    MaterialType type;       // material properties
    float ior;               // refractive index
    float Kd;                // diffuse albedo
    float Ks;                // specular albedo
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
    float cosi = clamp(-1.f, 1.f, I*N);
    float etai = 1, etat = ior;
    Vec3f n = N;
    if (cosi < 0) {                      // We need to handle with care the two possible situations:
        cosi = -cosi;                    //    - When the ray is inside the object
    } else {                             //    - When the ray is outside.
        std::swap(etai, etat); n = -N;   // If the ray is outside, you need to make cosi positive cosi = -N.I
    }                                    // If the ray is inside, you need to invert the refractive indices and negate the normal N
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? Vec3f(0,0,0) : I*eta + n*(eta * cosi - sqrtf(k));
}

// Compute Fresnel equation (the amount of reflected light)
float fresnel(const Vec3f &I, const Vec3f &N, const float &ior) {
    float cosi = clamp(-1.f, 1.f, I*N);
    float etai = 1, etat = ior;
    if (cosi > 0) std::swap(etai, etat);
    // Compute sini using Snell's law
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
    // Total internal reflection
    if (sint >= 1) {
        return 1;
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
        if (d>0) checkerboard_dist = d;
        if (checkerboard_dist<spheres_dist)  {
            hit = orig + dir*checkerboard_dist;
            N = Vec3f(0,1,0);
            material.diffuse = (int(.3*hit.x+1000) + int(.3*hit.z)) & 1 ? Vec3f(1,1,1) : Vec3f(.5,.3,.3);
        }
    }

    return std::min(spheres_dist, checkerboard_dist)<1000;
}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<Light> &lights, size_t depth) {
    if (depth > max_depth) {
        return background_color;
    }

    Vec3f point, N;
    Material material;

    if (!scene_intersect(orig, dir, spheres, point, N, material)) {
        return background_color;
    }

    switch (material.type) {
        case Material::REFLECTION:
            {
                Vec3f reflectionDirection = reflect(dir, N).normalize();
                Vec3f reflectionRayOrig = reflectionDirection*N < 0 ? point - N*1e-3 : point + N*1e-3;
                return cast_ray(reflectionRayOrig, reflectionDirection, spheres, lights, depth + 1) * .9;
            }

        case Material::REFLECTION_REFRACTION:
            {
                Vec3f reflectionDirection = reflect(dir, N).normalize();
                Vec3f refractionDirection = refract(dir, N, material.ior).normalize();
                Vec3f reflectionRayOrig = reflectionDirection*N < 0 ? point - N*1e-3 : point + N*1e-3;
                Vec3f refractionRayOrig = refractionDirection*N < 0 ? point - N*1e-3 : point + N*1e-3;
                Vec3f reflectionColor = cast_ray(reflectionRayOrig, reflectionDirection, spheres, lights, depth + 1);
                Vec3f refractionColor = cast_ray(refractionRayOrig, refractionDirection, spheres, lights, depth + 1);
                float kr = fresnel(dir, N, material.ior);
                return reflectionColor*kr + refractionColor*(1 - kr);  // the conservation of energy means that the transmittance can be obtained as 1-kr
            }
        case Material::SLIGHTLY_GLOSSY:
            {
                float diffuse_light_intensity = 0, specular_light_intensity = 0;

                for (size_t i=0; i<lights.size(); i++) {
                    Vec3f light_dir = lights[i].position - point;
                    float light_distance = light_dir.norm();
                    light_dir.normalize();

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
                return material.diffuse * diffuse_light_intensity * material.Kd + Vec3f(1., 1., 1.)*specular_light_intensity * material.Ks;
            }
        default:
            assert(0);
    };
    return background_color;
}


void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
    Vec3f *framebuffer = new Vec3f[width*height];

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

    delete [] framebuffer;
}

int main(int argc, char **argv) {
    std::vector<Sphere> spheres;
    std::vector<Light> lights;

    Material rubber(Material::SLIGHTLY_GLOSSY,      1.3, 0.8, 0.05, Vec3f(0.9, 0.3, 0.3), 25.);
    Material metal(Material::SLIGHTLY_GLOSSY,       1.3, 0.8, 0.8, Vec3f(.8, .8, .9), 55.);
    Material glass(Material::REFLECTION_REFRACTION, 1.5, 0.8, 0.2, Vec3f(0.6, 0.7, 0.8), 25.);
    Material mirror(Material::REFLECTION, 1.5, 0.8, 0.2, Vec3f(0.6, 0.7, 0.8), 25.);

    spheres.push_back(Sphere(Vec3f(-1  ,  0  , -12), 2  , metal));
    spheres.push_back(Sphere(Vec3f(3  ,  -2  , -12), 1.5  , rubber));
    spheres.push_back(Sphere(Vec3f( 0.5, -1, -8 ), 1.5, glass));
    spheres.push_back(Sphere(Vec3f(3, 3, -13 ), 2.5, mirror));

    lights.push_back(Light(Vec3f(-20, 70,  20), 0.5));
    lights.push_back(Light(Vec3f( 30, 50, -12), 0.3));
    lights.push_back(Light(Vec3f( 30, 20, 30),   0.7));

    render(spheres, lights);

    return 0;
}

