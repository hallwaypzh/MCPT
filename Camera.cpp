#include "Camera.h"

#include <cmath>
#include <iostream>

Camera::Camera(glm::vec3 &eye, glm::vec3 &center, glm::vec3 &up, float fovy, float aspect_ratio, int H) {

    glm::vec3 w = glm::normalize(eye - center);
    glm::vec3 u = glm::normalize(glm::cross(up, w));
    glm::vec3 v = glm::cross(w, u);

    //M is orthographic, M's transpose is also M's inverse
    this->M = glm::mat3x3(u, v, w);
    this->MT = glm::transpose(M);
    view_mat = glm::mat4x4(MT);
    view_mat[3] = glm::vec4(- MT * eye, 1.f);

    this->fovy = tanf(fovy / 2.f);
    this->inv_H = 2.f / H;
    this->ar = aspect_ratio;
}

Ray Camera::castRay(int i, int j) {
    //sample one ray per pixel so far
    //choose top left corner of the pixel 
    float dx = (i * inv_H - ar) * fovy;
    float dy = (1 - j * inv_H) * fovy;
    glm::vec3 origin(0.f, 0.f, 0.f);
    glm::vec3 dir = glm::normalize(glm::vec3(dx, dy, -1.f));
    return Ray(origin, dir);
}

std::vector<Ray> Camera::castRay(int i, int j, int spp, const glm::vec2 &u) {
    //jitters sampling for each pixel
    std::vector<Ray> rays;
    spp = (int) (sqrt(spp) + 0.5);
    
    //size of each sub pixel
    float delta = fovy * inv_H / spp;
    float dx = delta * u[0];
    float dy = delta * u[1];
    float x0 = i * fovy * inv_H;
    float y0 = j * fovy * inv_H;
    glm::vec3 origin(0.f, 0.f, 0.f);
    for (int v = 0; v < spp; v++) {
        for (int u = 0; u < spp; u++) {
            float x1 = x0 + delta * u + dx;
            float y1 = y0 + delta * v + dy;
            x1 = x1 - ar * fovy;
            y1 = -y1 + fovy;
            glm::vec3 dir = glm::normalize(glm::vec3(x1, y1, -1.f));
            rays.push_back(Ray(origin, dir));
        }
    }
    return rays;
}

glm::vec2 Camera::getScreenCoord(glm::vec3 d) {
    d = glm::normalize(d);
    d = MT * d;
    if (d.z > 0) {        
        return glm::vec2(0.f, 0.f);
    }
    return glm::vec2((-d.x / d.z + 1) * .5f, (d.y / d.z + 1) * .5f);
}