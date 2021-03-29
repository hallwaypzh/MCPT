#pragma once

#ifndef __SCENE_H
#define __SCENE_H

#include <vector>

#include "Triangle.h"
#include "KDTree.h"
#include "TriangleMeshes.h"
#include "Ray.h"
#include "Hit.h"
#include "SkyBox.h"
#include "EnvMap.h"

struct LightSample {
    glm::vec3 lp;
    glm::vec3 ln;
    glm::vec3 le;
    float pdf;

    glm::vec3 getRadiance(glm::vec3 p) {
        glm::vec3 D = p - lp;
        float r = glm::length(D);
        return le * glm::dot(ln, D) / (r * r * r) / pdf;
    }
};

class Scene {
    std::vector<Triangle> triangles;
    std::vector<Triangle> luminaires;
    KdTreeAccel* kdTree;
    std::random_device rd;
    std::mt19937 e2;
    std::uniform_real_distribution<> u01;
    std::uniform_int_distribution<> ud;
    SkyBox *skyBox;
    EnvMap *envMap;
    glm::mat4x4 inv_view_mat;
public:
    Scene(TriangleMeshes &meshes) {
        meshes.initializeTriangles(triangles);
        luminaires = meshes.getLuminaires();
        kdTree = new KdTreeAccel(triangles);
        e2 = std::mt19937(rd());
        ud = std::uniform_int_distribution<> (0, luminaires.size() - 1);
        skyBox = NULL;
        envMap = NULL;
    }

    bool intersect(const Ray &ray, float &t, Hit &hit_record) const {
        return kdTree->Intersect(ray, t, hit_record);
    }

    bool intersect(const Ray &ray) const {
        return kdTree->IntersectP(ray);
    }

    bool intersect(const Ray &ray, float &t)  const {
        return kdTree->Intersect(ray, t);
    }

    void sampleLight(glm::vec3 &p, glm::vec3 &n, glm::vec3 &le, float &pdf) {
        pdf = 1.f / luminaires.size();
        int id = ud(e2);
        float u1 = u01(e2);
        float u2 = u01(e2);
        float pdf1;
        le = luminaires[id].getLe();
        luminaires[id].uniformSample(u1, u2, pdf1, p, n);
        pdf *= pdf1;
    }

    LightSample sampleLight(Sampler &sampler) {
        LightSample sample;
        int id = ud(e2);
        sample.pdf = 1.f / luminaires.size();
        float pdf1;
        glm::vec2 u = sampler.uniformSample2D();       
        sample.le = luminaires[id].getLe();
        luminaires[id].uniformSample(u.x, u.y, pdf1, sample.lp, sample.ln);
        sample.pdf *= pdf1;
        return sample;
    }

    glm::vec3 sampleSkybox(const Ray &ray) const {
        glm::vec3 d = glm::vec3(inv_view_mat * glm::vec4(ray.d, 0.f));
        //std::cout << d.x << " " << d.y << " " << d.z << std::endl;
        if (skyBox == NULL && envMap == NULL)
            return glm::vec3(0.7f);
        if (envMap != NULL) {
            return envMap->sampleLi1(d);
        }
        return skyBox->sampleEnvLight(d);
    }

    void createSkyBox(std::string dir, std::string name, int type) {
        envMap = new EnvMap(dir, name);
    }

    int getLightNums() { return luminaires.size(); }

    void setCam2World(glm::mat4x4 &mat) { this->inv_view_mat = mat; }
};

#endif