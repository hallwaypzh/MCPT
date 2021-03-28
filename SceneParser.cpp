#include "SceneParser.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <glm/glm.hpp>

SceneParser::SceneParser(std::string name) {
    std::ifstream fin;
    fin.open(name);
    std::string line;
    scene_name = name;
    has_skybox = false;
    env_map_type = -1;
    while (getline(fin, line)) {
        if (line.size() == 0) {
            continue;
        }
        std::stringstream ss(line);
        std::string discription;
        ss >> discription;
        if (discription == "Output") {
            parseOutput(fin);
        } else if (discription == "Camera") {
            parseCamera(fin);
        } else if (discription == "Meshes") {
            parseMeshes(fin);
        } else if (discription == "SkyBox") {
            parseSkyBox(fin);
        } else {
            std::cerr << "Bad line: " << line << std::endl;
        }
    }
    fin.close();

}

void SceneParser::parseCamera(std::ifstream &fin) {
    std::string line;
    while (getline(fin, line)) {
        if (line.size() == 0) {
            continue;
        }
        std::stringstream ss(line);
        std::string discription;
        ss >> discription;
        if (discription == "position") {
            ss >> position[0] >> position[1] >> position[2];
        } else if (discription == "center") {
            ss >> center[0] >> center[1] >> center[2];
        } else if (discription == "up") {
            ss >> up[0] >> up[1] >> up[2];
        } else if (discription == "fov") {
            ss >> fovy;
            fovy = fovy * M_PI / 180.f;
        }else if (discription == "}") {
            break;
        } else {
            std::cerr << "Bad line: " << line << std::endl;
        }
    }
    float aspect_ratio = width * 1.f / height;
}

void SceneParser::parseOutput(std::ifstream &fin) {
    std::string line;
    while (getline(fin, line)) {
        if (line.size() == 0) {
            continue;
        }
        std::stringstream ss(line);
        std::string discription;
        ss >> discription;
        if (discription == "ss") {
            ss >> supersampling;
        } else if (discription == "rpp") {
            ss >> rayperpixel;
        } else if (discription == "width") {
            ss >> width;
        } else if (discription == "height") {
            ss >> height;
        } else if (discription == "name") {
            ss >> output_name;
        } else if (discription == "}") {
            break;
        } else {
            std::cerr << "Bad line: " << line << std::endl;
        }
    }
}

void SceneParser::parseMeshes(std::ifstream &fin) {
    std::string line;
    while (getline(fin, line)) {
        if (line.size() == 0) {
            continue;
        }
        std::stringstream ss(line);
        std::string discription;
        ss >> discription;
        if (discription == "directory") {
            ss >> mesh_dir;
        } else if (discription == "name") {
            ss >> mesh_name;
        } else if (discription == "}") {
            break;
        } else {
            std::cerr << "Bad line: " << line << std::endl;
        }
    }
}

void SceneParser::parseSkyBox(std::ifstream &fin) {
    std::string line;
    has_skybox = true;
    env_map_type = 0;
    while (getline(fin, line)) {
        if (line.size() == 0) {
            continue;
        }
        std::stringstream ss(line);
        std::string discription;
        ss >> discription;
        if (discription == "directory") {
            ss >> skybox_dir;
        } else if (discription == "name") {
            ss >> skybox_name;
        } else if (discription == "SH") {
            env_map_type = 1;
        }else if (discription == "}") {
            break;
        } else {
            std::cerr << "Bad line: " << line << std::endl;
        }
    }
}