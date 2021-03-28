#pragma once

#ifndef __TRIANGLEMESHES_H
#define __TRIANGLEMESHES_H

#include <vector>
#include <utility>
#include <tuple>
#include <string>
#include <map>
#include <glm/glm.hpp>

#include "Ray.h"
#include "Triangle.h"
#include "Material.h"

class TriangleMeshes {
    //first load geometry

    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec2> texcoord;

    std::vector<std::tuple<int, int, int> > vids;
    std::vector<std::tuple<int, int, int> > nids;
    std::vector<std::tuple<int, int, int> > tids;
    std::vector<int> mids;
    std::vector<Material> materials;
    std::vector<Triangle> luminaires;

    //TO-DO load textures & materials
    //TO-DO construct acceleration structure
public:
    TriangleMeshes(std::string &dir, std::string &file_name);
    ~TriangleMeshes() {}
    bool intersect(Ray &r);
    bool intersect(Ray &r, float &t, glm::vec4 &hit_point);
    void applyTransform(glm::mat4x4 &modelview);
    void initializeTriangles(std::vector<Triangle> &triangles);
    std::vector<Triangle> getLuminaires() { return luminaires; }
    std::vector<Material> getMaterials() { return materials; }
};

#endif