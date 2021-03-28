#pragma once

#ifndef __ENVMAP_H
#define __ENVMAP_H

#include <glm/glm.hpp>
#include <vector>
#include <string>
#include <iostream>

struct EnvMap {
    //[An efficient representation for irradiance environment maps]
    std::string name;
    int W, H, C;
    std::vector<glm::vec3> texels;
    std::vector<float> Ls; 
    std::vector<std::vector<float> > distribution;
    EnvMap(std::string dir, std::string name);

    std::vector<float> sh[3];

    std::vector<float> Y(glm::vec3 d) {
        //[An efficient representation for irradiance environment maps]
        float C[9] = {0.282095, 0.488603, 0.488603,
                      0.488603, 1.092548, 1.092548,
                      0.315392, 1.092548, 0.546274};
        
        std::vector<float> v(9, 0.f);
        v[0] = C[0]; 
        v[1] = C[1] * d.y; 
        v[2] = C[2] * d.z;
        v[3] = C[3] * d.x;
        v[4] = C[4] * d.x * d.y; 
        v[5] = C[5] * d.y * d.z; 
        v[6] = C[6] * (3 * d.z * d.z - 1);
        v[7] = C[7] * d.x * d.z;
        v[8] = C[8] * (d.x + d.y) * (d.x - d.y);
        return v;
    }

    glm::vec3 sampleLi(glm::vec3 d);

    glm::vec3 sampleLi1(glm::vec3 d);
};

#endif