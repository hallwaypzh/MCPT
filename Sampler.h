#pragma once

#ifndef __SAMPLER_H
#define __SAMPLER_H

#include <random>
#include <glm/glm.hpp>

struct Sampler {
    std::random_device rd;
    std::mt19937 e2;
    std::uniform_real_distribution<> u01;

    Sampler() {
        e2 = std::mt19937(rd());
        u01 = std::uniform_real_distribution<> (0, 1);
    }

    float uniformSample1D() {
        return u01(e2);
    }

    glm::vec2 uniformSample2D() {
        return glm::vec2(u01(e2), u01(e2));
    }
    
    int uniformSampleInt(int a, int b) {
        std::uniform_int_distribution<> dist(a, b);
        return dist(e2);
    }
};

#endif