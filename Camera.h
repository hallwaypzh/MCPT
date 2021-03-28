#pragma once

#ifndef __CAMERA_H
#define __CAMERA_H

#include <glm/glm.hpp>
#include <vector>
#include "Ray.h"


struct Camera {
    glm::mat4x4 view_mat;
    float ar;
    float fovy;
    float inv_H;

    glm::mat3x3 M, MT;
    Camera(glm::vec3 &eye, glm::vec3 &center, glm::vec3 &up, float fovy, float aspect_ratio, int H);
    glm::mat4x4 getViewMat() { return this->view_mat; }
    Ray castRay(int i, int j);
    std::vector<Ray> castRay(int i, int j, int spp, const glm::vec2 &u);
    glm::vec2 getScreenCoord(glm::vec3 d);
};

#endif