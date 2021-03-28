#pragma once

#ifndef __SKYBOX_H
#define __SKYBOX_H

#include <string>
#include <vector>
#include <glm/glm.hpp>
#include "Camera.h"
class SkyBox {
    //https://cglearn.codelight.eu/pub/computer-graphics/environment-mapping
    /*
    Layout
        [..][y+][..][..]
        [x-][z+][x+][z-]
        [..][y-][..][..]
    */
    std::string name;        
    std::vector<glm::vec3> env_map[6]; //[x-, x+, y-, y+, z-, z+]
    Camera* cameras[6];
    int size;
public:
    SkyBox(std::string dir, std::string name);
    glm::vec3 sampleEnvLight(glm::vec3 d);
    
};

#endif