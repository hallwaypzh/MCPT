#pragma once
#ifndef __AABB_H
#define __AABB_H

#include <algorithm>

#include <glm/glm.hpp>
#include <glm/gtx/component_wise.hpp>

#include "Ray.h"

struct AABB {
    glm::vec3 pmax;
    glm::vec3 pmin;
    AABB() {}
    AABB(glm::vec3 pmin, glm::vec3 pmax) { this->pmax = pmax; this->pmin = pmin; }
    ~AABB() {}

    void Union(AABB &bbox) {
        pmin = glm::min(pmin, bbox.pmin);
        pmax = glm::max(pmax, bbox.pmax);
    }

    float SurfaceArea() const {
        glm::vec3 d = glm::abs(pmax - pmin);
        return 2 * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0]);
    }

    int MaximumExtent() const {
        glm::vec3 d = glm::abs(pmax - pmin);
        int axis = 0;
        if (d[axis] < d[1])
            axis = 1;
        if (d[axis] < d[2])
            axis = 2;
        return axis;
    }

    bool intersect(const Ray &ray, float &tmin, float &tmax) const {
        //  reference
        //  https://gamedev.stackexchange.com/questions/18436/most-efficient-aabb-vs-ray-collision-algorithms
        glm::vec3 v1 = pmin - ray.o;
        glm::vec3 v2 = pmax - ray.o;
        glm::vec3 tmp1 = v1 * ray.invd;
        glm::vec3 tmp2 = v2 * ray.invd;
        glm::vec3 tmp3 = glm::min(tmp1, tmp2);
        glm::vec3 tmp4 = glm::max(tmp1, tmp2);
        float t0 = glm::compMax(tmp3);
        float t1 = glm::compMin(tmp4);

        if (t1 < 0)
            return false;
        if (t0 > t1) 
            return false;
        if (t0 > ray.tmax)
            return false;
        tmin = t0;
        tmax = t1;
        return true;
    }
    bool intersect(const Ray &ray) const {
        //  reference
        //  https://gamedev.stackexchange.com/questions/18436/most-efficient-aabb-vs-ray-collision-algorithms
        glm::vec3 v1 = pmin - ray.o;
        glm::vec3 v2 = pmax - ray.o;
        
        glm::vec3 tmp1 = v1 * ray.invd;
        glm::vec3 tmp2 = v2 * ray.invd;
        glm::vec3 tmp3 = glm::min(tmp1, tmp2);
        glm::vec3 tmp4 = glm::max(tmp1, tmp2);
        float t0 = glm::compMax(tmp3);
        float t1 = glm::compMin(tmp4);
        
        if (t1 < 0)
            return false;
        if (t0 > t1) 
            return false;
        if (t0 > ray.tmax)
            return false;
        return true;
    }
};

#endif