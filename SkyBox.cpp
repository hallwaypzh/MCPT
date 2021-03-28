#include "SkyBox.h"

#include "stb_image.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <glm/gtx/component_wise.hpp>

SkyBox::SkyBox(std::string dir, std::string name) {
    std::string path = dir + name;
    this->name = name;
    int w, h, c;
    unsigned char * tex_data = stbi_load(path.c_str(), &w, &h, &c, 0);
    
    size = (int) std::min(w / 4, h / 3);
    
    //int xs[6] = {1, 0, 1, 2, 3, 1};
    //int ys[6] = {0, 1, 1, 1, 1, 2};
    int xs[6] = {0, 2, 1, 1, 3, 1};
    int ys[6] = {1, 1, 2, 0, 1, 1};
    //std::string names[] = {"y+.jpg", "x-.jpg", "z+.jpg", "x+.jpg", "z-.jpg", "y-.jpg"};
    for (int i = 0; i < 6; i++) {
        env_map[i].resize(size * size, glm::vec3(0.f));
    }
    for (int i = 0; i < 6; i++) {
        int x1 = xs[i] * size;
        int y1 = ys[i] * size;
        //std::cout << names[i] << " ";
        //std::cout << y1 << ":" << y1 + size << " " << x1 << ":" << x1 + size << std::endl;
        for (int row = y1; row < y1 + size; row++) {
            for (int col = x1; col < x1 + size; col++) {
                glm::vec3 color;
                for (int k = 0; k < c; k++) {
                    //std::cout << k << std::endl;
                    color[k] = tex_data[((row * w) + col) * c + k]  / 255.f;
                }
                env_map[i][((row - y1)* size + (col -  x1))] =  color;                
            }
        }
        //DEBUG, ouput each plane of the box as jpg
        /*
        unsigned char* data = new unsigned char[size * size * 3];
        for (int row = 0; row < size; row++) {
            for (int col = 0; col < size; col++) {
                data[(row * size + col) * 3 + 0] = (unsigned char) (env_map[i][row * size + col].r * 255.f);
                data[(row * size + col) * 3 + 1] = (unsigned char) (env_map[i][row * size + col].g * 255.f);
                data[(row * size + col) * 3 + 2] = (unsigned char) (env_map[i][row * size + col].b * 255.f);                
            }
        }        
        stbi_write_jpg(names[i].c_str(), size, size, 3, data, 100);
        std::cout << "save img " << names[i] << std::endl;
        delete [] data;
        std::cout << i << std::endl;
        */
    }
    stbi_image_free(tex_data);    
    

    glm::vec3 eye(0.f, 0.f, 0.f);
    glm::vec3 centers[] = {
        glm::vec3(1.f, 0.f, 0.f), glm::vec3(-1.f, 0.f, 0.f), glm::vec3(0.f, -1.f, 0.f),
        glm::vec3(0.f, 1.f, 0.f), glm::vec3(0.f, 0.f, -1.f), glm::vec3(0.f, 0.f, 1.f) 
    };
    glm::vec3 ups[] = {
        glm::vec3(0.f, 1.f, 0.f), glm::vec3(0.f, 1.f, 0.f), glm::vec3(0.f, 0.f, 1.f),
        glm::vec3(0.f, 0.f, -1.f), glm::vec3(0.f, 1.f, 0.f), glm::vec3(0.f, 1.f, 0.f)
    };
    for (int i = 0; i < 6; i++) {
        cameras[i] = new Camera(eye, centers[i], ups[i], M_PI / 2, 1.f, size);
    }
    std::cout << "Created A Skybox:)\n";
}


glm::vec3 SkyBox::sampleEnvLight(glm::vec3 d) {
    d = glm::normalize(d);
    float val = glm::compMax(glm::abs(d));
    int id = -1;
    float u, v;
    glm::vec3 d1 = d / val;
    for (int i = 0; i < 3; i++) {
        if (fabs(d[i]) == val) {
            id = 2 * i + (d[i] > 0);
            //std::cout << id << std::endl;
            break;
        } 
    }
    glm::vec2 uv = cameras[id]->getScreenCoord(d);
    //nearest neighbor
    int i = (int) (uv.x * size);
    int j = (int) (uv.y * size);
    return env_map[id][(j * size + i)];
}