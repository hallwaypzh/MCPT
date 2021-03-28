#pragma once

#ifndef __TEXTURE_H
#define __TEXTRUE_H

#include <string>
#include <vector>
#include <iostream>

#include <glm/glm.hpp>

struct Texture {
    //std::vector<float> data;
    std::vector<glm::vec3> data;
    std::string name;
    int w, h, c; //width, height, channels
    Texture() : w(0), h(0), c(0), name("Dummy") {}
    Texture(std::string dir, std::string name);
    glm::vec3 sampleColor(float u, float v) const;
    glm::vec3 getColor(int i, int j) const;

};

#endif