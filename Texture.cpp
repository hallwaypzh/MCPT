#include "Texture.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "stbi_image_write.h"

#include <cmath>

Texture::Texture(std::string dir, std::string name) {
    std::string path = dir + name;
    this->name = name;
    unsigned char *tex_data = stbi_load(path.c_str(), &w, &h, &c, 3);
    if (tex_data == NULL) {
        std::cerr << "fail to load texture: " << path << std::endl;
    }
    c = 3;
    data.resize(w * h);
    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
            glm::vec3 pixel(0.f);
            for (int k = 0; k < c; k++) {
                pixel[k] = (float)(tex_data[(j * w + i) * c + k]) / 255.f;
            }
            data[j * w + i] = (pixel);
        }
    }
    stbi_image_free(tex_data);
}

glm::vec3 Texture::getColor(int i, int j) const {
    // i, j in uv space
    // data is stored in the order of screen space
    j = h - 1 - j;
    i = w - 1 - i;
    return data[(j * w) + i];
}

glm::vec3 Texture::sampleColor(float u, float v) const {
    if (w == 0 || h == 0) {
        return glm::vec3(1.f);
    }
    //bi-linear interpolation
    //wrap the texture map around if necessary
    u = u - floorf(u);
    v = v - floorf(v);
    u *= w;
    v *= h;
    int i = (int) u;
    int j = (int) v;
    if ((u < 0.5) || (u > w - 0.5) || (v < 0.5) || (v > h - 0.5)) {
       return getColor(i, j);
    }
    glm::vec3 result(0.f);
    glm::vec3 c1, c2, c3, c4;
    int x1 = (int) (u - .5f);
    int y1 = (int) (v - .5f);
    int x2 = (int) (u + .5f);
    int y2 = (int) (v + .5f);
    float t = (u - (x1 + .5f));
    float s = (v - (y1 + .5f));
    c1 = getColor(x1, y1);
    c2 = getColor(x2, y1);
    c3 = getColor(x1, y2);
    c4 = getColor(x2, y2);
    return (1 - t) * (1 - s) * c1 + (1 - t) * s * c3 + t * (1 - s) * c2 + s * t * c3;
}