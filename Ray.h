#pragma once

#ifndef __RAY_H
#define __RAY_H

#include <glm/glm.hpp>

struct Ray {
    glm::vec3 o;
    glm::vec3 d;
    glm::vec3 invd;
    float tmax;
public:
    Ray(glm::vec3 &origin, glm::vec3 &dir, float tmax = FLT_MAX) : o(origin) {
        d = glm::normalize(dir);
        invd = 1.f / d; 
        this->tmax = tmax;
    }

    Ray(const Ray &r) : o(r.o), d(r.d), invd(1.f / r.d), tmax(r.tmax) { }
    ~Ray() { }
};

#endif