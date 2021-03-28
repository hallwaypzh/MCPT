#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <chrono>
#include <thread>
#include <mutex>

#include "TriangleMeshes.h"
#include "Triangle.h"
#include "Camera.h"
#include "KDTree.h"
#include "AABB.h"
#include "Scene.h"
#include "PathTracer.h"
#include "Sampler.h"
#include "SceneParser.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stbi_image_write.h"

void output_depthmap(std::string &filename, std::vector<float> depth_buffer, int H, int W);
void output_jpg(std::string name, std::vector<glm::vec3> img, int W, int H);

int main(int argc, char* argv[]) {
    if (argc < 2) {
      std::cerr << "Need to specify scene name!";
      return 0;
    } 
    std::string scene_name(argv[1]);
    SceneParser sceneParser(scene_name);
    PathTracer tracer(1024);
    
    int total = 0;
    int rpp = sceneParser.rayperpixel;
    int ss = sceneParser.supersampling;
    int W = sceneParser.width;
    int H = sceneParser.height;
    std::vector<glm::vec3> img(W * H, glm::vec3(0.f, 0.f, 0.f));
    Sampler sampler;
    std::mutex total_mutex;
    float aspect_ratio = W * 1.f / H;
    Camera camera(sceneParser.position, sceneParser.center, sceneParser.up, sceneParser.fovy, aspect_ratio, H);
    glm::mat4x4 view_mat = camera.getViewMat();
    TriangleMeshes triangleMeshes(sceneParser.mesh_dir, sceneParser.mesh_name);
    triangleMeshes.applyTransform(view_mat);
    Scene scene(triangleMeshes);

    glm::mat4x4 inv_view_mat = glm::inverse(view_mat);
    scene.setCam2World(inv_view_mat);

    if (sceneParser.has_skybox) {
      scene.createSkyBox(sceneParser.skybox_dir, sceneParser.skybox_name, sceneParser.env_map_type);
    }
    #pragma omp parallel for
    for (int j = 0; j < H; j++) {
      std::cout << "Starting to process " << j << "th row!:)\n";
        #pragma omp parallel for
        for (int i = 0; i < W; i++) {
            //jitters
            glm::vec2 u = sampler.uniformSample2D();
            std::vector<Ray> rays = camera.castRay(i, j, ss, u);
            for (auto ray : rays) {
              for (int k = 0; k < rpp; k++) {
                img[(j * W + i)] += tracer.Li(ray, scene, 0);
              }
             
            }
            img[(j * W + i)] /= (ss);
            img[(j * W + i)] /= rpp;
        }
        std::lock_guard<std::mutex> guard(total_mutex);
        total += 1;
        std::cout << "Finish to process " << j << "th row!\n";
        std::cout << total * 1.f / H * 100.f << std::endl;
    }
    output_jpg(sceneParser.output_name, img, W, H);    
    return 0;

}

void output_jpg(std::string name, std::vector<glm::vec3> img, int W, int H) {
  unsigned char *data = NULL;
  data = new unsigned char [W * H * 3];
  for (int j = 0; j < H; j++) {
    for (int i = 0; i < W; i++) {
      data[3 * (j * W + i) + 0] = (unsigned char) (std::min(img[(j * W + i)].r, 1.f) * 255.f);
      data[3 * (j * W + i) + 1] = (unsigned char) (std::min(img[(j * W + i)].g, 1.f) * 255.f);
      data[3 * (j * W + i) + 2] = (unsigned char) (std::min(img[(j * W + i)].b, 1.f) * 255.f);
      /*
      if (data[3 * (j * W + i) + 0] == 255) {
        std::cout << j << " " << i << std::endl;
        std::cout << (std::min(img[(j * W + i)].r, 1.f) * 255.f) << std::endl;
      }
      */
    }
  }
  stbi_write_jpg(name.c_str(), W, H, 3, data, 100);
}
