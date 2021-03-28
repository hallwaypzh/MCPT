#pragma once
#ifndef __TRIANGLE_H
#define __TRIANGLE_H

#include <glm/glm.hpp>

#include "Material.h"
#include "AABB.h"
#include "Ray.h"
#include "Hit.h"

class Triangle {
    glm::vec3 p[3];
    glm::vec3 n[3];
    glm::vec2 t[3];
    Material* material;
    glm::vec3 p2p0;
    glm::vec3 p2p1;
    glm::vec3 gn;
    float inv_area;
    float epsilon;
    
public:
    AABB bbox;
    Triangle(glm::vec3 &p0, glm::vec3 &p1, glm::vec3 &p2,
             glm::vec3 &n0, glm::vec3 &n1, glm::vec3 &n2,
             glm::vec2 &t0, glm::vec2 &t1, glm::vec2 &t2, Material *material) {
        p[0] = p0; p[1] = p1; p[2] = p2;
        n[0] = glm::normalize(n0); n[1] = glm::normalize(n1); n[2] = glm::normalize(n2);
        t[0] = t0; t[1] = t1; t[2] = t2;
        this->material = material;
        p2p0 = p0 - p2;
        p2p1 = p1 - p2;
        gn = glm::normalize(glm::cross(p2p0, p2p1));
        bbox.pmin = glm::min(glm::min(p0, p1), p2);
        bbox.pmax = glm::max(glm::max(p0, p1), p2);
        inv_area = 2.f / glm::length(glm::cross(p2p0, p2p1));
        epsilon = 1e-7;
    }

    bool intersect(const Ray &ray) const {
        // naive ray triangle intersection
        // will probably implement a sphosicate one later
        if (!bbox.intersect(ray))
            return false;
        glm::mat3x3 M(p2p0, p2p1, -ray.d);
        float d = glm::determinant(M);
        if (fabs(d) < epsilon) {
            return false;
        }

        glm::vec3 v = glm::inverse(M) * (ray.o - p[2]);
        if (v[2] < epsilon)
            return false;
        if (v[2] > ray.tmax - epsilon)
            return false;
        float a = v[0], b = v[1], c = v[0] + v[1];
        return (a > -epsilon && a < 1+epsilon) && (b > -epsilon && b < 1+epsilon) && (c > -epsilon && c < 1+epsilon);
    }

    bool intersect(const Ray &ray, float &t) const {
        if (!bbox.intersect(ray))
            return false;
        glm::mat3x3 M(p2p0, p2p1, -ray.d);
        float d = glm::determinant(M);
        if (fabs(d) < 1e-6) {
            return false;
        }

        glm::vec3 v = glm::inverse(M) * (ray.o - p[2]);
        if (v[2] < 1e-6 || v[2] > t - 1e-6)
            return false;
        float a = v[0], b = v[1], c = v[0] + v[1];
        if ((a > -1e-6 && a < 1+1e-6) && (b > -1e-6 && b < 1+1e-6) && (c > -1e-6 && c < 1+1e-6)) {
            t = v[2];
            return true;
        }
        return false;
    }

    bool intersect(const Ray &ray, float &t, Hit &hit) const {
        if (!bbox.intersect(ray))
            return false;
        glm::mat3x3 M(p2p0, p2p1, -ray.d);
        float d = glm::determinant(M);        
        if (fabs(d) < 1e-6) {
            return false;
        }
        glm::vec3 v = glm::inverse(M) * (ray.o - p[2]);
        if (v[2] < 1e-6 || v[2] > t + 1e-6)
            return false;
        float a = v[0], b = v[1], c = 1.f - (a + b);
        if ((a > -1e-6 && a < 1+1e-6) && (b > -1e-6 && b < 1+1e-6) && (c > -1e-6 && c < 1+1e-6)) {
            glm::vec3 pos = a * p[0] + b * p[1] + c * p[2];            
            glm::vec3 norm = glm::normalize(a * n[0] + b * n[1] + c * n[2]);
            hit.texcoord = a * this->t[0] + b * this->t[1] + c * this->t[2];
            hit.hitpoint = pos + 1e-6f * norm;
            hit.norm = norm;
            hit.gn = this->gn;
            hit.material = material;
            hit.constructFrame();
            t = v[2];
            return true;
        }
        return false;
    }

    void uniformSample(const float u, const float v, float &pdf, glm::vec3 &sp, glm::vec3 &sn) {
        pdf = inv_area;
        float s = std::sqrt(u);
        float t = v;
        sp = (1 - s) * p[0] + s * (1 - t) * p[1] + s * t * p[2];
        sn = (1 - s) * n[0] + s * (1 - t) * n[1] + s * t * n[2];
        sn = glm::normalize(sn);
    }

    glm::vec3 getLe() { return material->le; }
};

#endif