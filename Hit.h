#pragma once

#ifndef __HIT_H
#define __HIT_H

#include <glm/glm.hpp>
#include "Ray.h"

struct Hit {
    glm::vec3 hitpoint;
    glm::vec3 norm;
    glm::vec3 gn; //geometry norm
    glm::vec2 texcoord;
    Material *material;
    glm::mat3x3 M;

    Hit() { }
    Hit(glm::vec3 &p, glm::vec3 &n, glm::vec2 &t, Material *m) : hitpoint(p), norm(n), texcoord(t), material(m) { 
        //  [T B N] right-handed coord frame
        //  https://www.scratchapixel.com/lessons/3d-basic-rendering/global-illumination-path-tracing/global-illumination-path-tracing-practical-implementation
        glm::vec3 N = glm::normalize(n);
        glm::vec3 T, B;
        if (std::fabs(N.x) > std::fabs(N.y))
            T = glm::normalize(glm::vec3(N.z, 0.f, -N.x));
        else
            T = glm::normalize(glm::vec3(0.f, N.z, -N.y));
        B = glm::cross(N, T);
        M = glm::mat3x3(T, B, N);
    }

    void constructFrame() {
        //  [T B N] right-handed coord frame
        //  https://www.scratchapixel.com/lessons/3d-basic-rendering/global-illumination-path-tracing/global-illumination-path-tracing-practical-implementation
        glm::vec3 N = norm;
        glm::vec3 T, B;
        if (std::fabs(N.x) > std::fabs(N.y))
            T = glm::normalize(glm::vec3(N.z, 0.f, -N.x));
        else
            T = glm::normalize(glm::vec3(0.f, N.z, -N.y));
        B = glm::cross(N, T);
        M = glm::mat3x3(T, B, N);
    }

    Ray castRay(glm::vec3 look_at) {
        glm::vec3 dir = glm::normalize(look_at - hitpoint);
        float tmax = std::fabs(look_at.x - hitpoint.x);
        return Ray(hitpoint, dir, tmax);
    }

    glm::vec3 world2local(glm::vec3 v) {
        // M in ortho, transpose of M = inverse of M
        return glm::transpose(M) * v;
    }

    glm::vec3 local2world(glm::vec3 v) {
        return M * v;
    }

    void sampleBRDF(Sampler &sampler, glm::vec3 wo, glm::vec3 &wi, glm::vec3 &f, float &pdf) {
        wo = glm::transpose(M) * wo;    
        material->sampleF(sampler, wo, wi, f, pdf, texcoord);
        wi = M * wi;
    }


    glm::vec3 F(const glm::vec3 &wi, const glm::vec3 &wo) {
        // both wi, wo in world space
        // both wi, wo should be unit vector
        //return material->F(world2local(wi), world2local(wo));
        return material->F(world2local(wi), world2local(wo), texcoord);
    }
};

#endif