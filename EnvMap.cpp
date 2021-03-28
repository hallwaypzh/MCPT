#include "EnvMap.h"

#include "stb_image.h"

EnvMap::EnvMap(std::string dir, std::string name) {
    std::string path = dir + name;
    stbi_hdr_to_ldr_gamma(2.7f);
    unsigned char *data = stbi_load(path.c_str(), &W, &H, &C, 0);
    if (data == NULL) {
        std::cerr << "Fail to load environment map:(\n";
    }
    texels.resize(W * H, glm::vec3(0.f));
    Ls.resize(W * H, 0.f);
    glm::vec3 rgb2l(0.2126, 0.7152, 0.0722);
    float I = 0.f;
    float invH = M_PI / H;
    float invW = 2 * M_PI / W;
    for (int j = 0; j < H; j++) {
        float sinTheta = sinf(M_PI * (j + 0.5) * 1.f / H);
        for (int i = 0; i < W; i++) {
            texels[j * W + i].r = (float)(data[(j * W + i) * C + 0] / 255.f);
            texels[j * W + i].g = (float)(data[(j * W + i) * C + 1] / 255.f);
            texels[j * W + i].b = (float)(data[(j * W + i) * C + 2] / 255.f);
        }
    }
    //convert environment map into spherical harmonic lighting coefficients 
    for (int i = 0; i < 3; i++) {
        sh[i].resize(9, 0.f);
    }
       
    for (int j = 0; j < H; j++) {
        float theta = (j + 0.5f) * invH;
        float sinTheta = sinf(theta);
        float cosTheta = cosf(theta);
        for (int i = 0; i < W; i++) {
            float phi = (i + 0.5f) * invW;
            glm::vec3 d(sinTheta * cosf(phi), sinTheta * sinf(phi), cosTheta);
            std::vector<float> ys = Y(d);
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 9; l++) {
                    sh[k][l] += texels[j * W + i][k] * ys[l] * invH * invW;
                }
            }
        }
    }
}


glm::vec3 EnvMap::sampleLi(glm::vec3 d) {
    //assume d is already normalized
    glm::vec3 color(0.f);
    std::vector<float> ys = Y(d);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 9; j++) {
            color[i] += sh[i][j] * ys[j];
        }
    }
    return color;
}

glm::vec3 EnvMap::sampleLi1(glm::vec3 d) {
    d = glm::normalize(d);
    float sinTheta = sqrt(1.f - d.y * d.y);
    float cosPhi = d.z / sinTheta;
    float theta = acos(d.y);
    float phi = atan(d.x / d.z) + M_PI;
    float v = H * theta / M_PI;
    float u = W * phi / (2 * M_PI);
    assert(0 <= u <= W);
    assert(0 <= v <= H);
    int j = (int) v;
    int i = (int) u;
    return texels[(j * W + i)];
}