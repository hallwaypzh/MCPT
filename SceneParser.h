#pragma once

#ifndef __SCENEPARSER_H
#define __SCENEPARSER_H

#include <glm/glm.hpp>
#include <string>
#include "Camera.h"
#include "Scene.h"
#include "TriangleMeshes.h"

struct SceneParser {
    int height;
    int width;
    int supersampling;
    int rayperpixel;
    glm::vec3 position;
    glm::vec3 center;
    glm::vec3 up;
    float fovy;
    std::string scene_name;
    std::string output_name;
    std::string mesh_dir;
    std::string mesh_name;
    std::string skybox_dir;
    std::string skybox_name;
    bool has_skybox;
    int env_map_type;
    Scene* scene;
    SceneParser(std::string name);


    void parseScene(std::ifstream &fin);
    void parseCamera(std::ifstream &fin);
    void parseMeshes(std::ifstream &fin);
    void parseSkyBox(std::ifstream &fin);
    void parseOutput(std::ifstream &fin);

};

#endif