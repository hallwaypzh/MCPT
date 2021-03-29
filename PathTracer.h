#pragma once

#ifndef __PATHTRACER_H
#define __PATHTRACER_H

#include <glm/glm.hpp>

#include "Scene.h"
#include "Hit.h"
#include "Ray.h"
#include "Sampler.h"

class PathTracer {
    int max_depth;
    Sampler sampler;
    float threshold;
public:
    PathTracer(int depth) : max_depth(depth) { threshold = 0.5; }

    glm::vec3 directLight(const Ray &r, Scene &scene, Hit &hit, int sr) {
        glm::vec3 color(0.f);
        
        float thresh = scene.getLightNums() * 1.f / (scene.getLightNums() + 1.f);
        for (int k = 0; k < sr; k++) {
            float rv = sampler.uniformSample1D();
            if (rv < thresh) {
                LightSample light = scene.sampleLight(sampler);
                Ray shadow = hit.castRay(light.lp);
                if (glm::dot(hit.gn, shadow.d) < 0)
                    continue;
                if (!scene.intersect(shadow)) {
                    //glm::vec3 wi;
                    glm::vec3 f = hit.F(shadow.d, -r.d);
                    color += glm::dot(shadow.d, hit.norm) * f * light.getRadiance(hit.hitpoint) / thresh;
                }
            } else {
                glm::vec3 wi;
                float pdf;
                glm::vec3 f;
                Hit hit1;
                hit.sampleBRDF(sampler, -r.d, wi, f, pdf);
                Ray shadow(hit.hitpoint, wi);
                float t;
                if (!scene.intersect(shadow, t, hit1)) {
                    color += glm::dot(shadow.d, hit.norm) * f * scene.sampleSkybox(shadow) / (pdf  * (1 - thresh));
                } else {
                    color += glm::dot(wi, hit.norm) * f * hit1.material->le / (pdf * (1 - thresh));
                }
            }
        }
        /*
        if (scene.getLightNums() > 0) {
            for (int k = 0; k < sr; k++) {
                LightSample light = scene.sampleLight(sampler);
                Ray shadow = hit.castRay(light.lp);
                if (glm::dot(hit.gn, shadow.d) < 0)
                    continue;
                if (!scene.intersect(shadow)) {
                    //glm::vec3 wi;
                    glm::vec3 f = hit.F(shadow.d, -r.d);
                    color += glm::dot(shadow.d, hit.norm) * f * light.getRadiance(hit.hitpoint);
                }
            }
        } else {
            for (int k = 0; k < sr; k++) {
                glm::vec3 wi;
                float pdf;
                glm::vec3 f;
                hit.sampleBRDF(sampler, -r.d, wi, f, pdf);
                Ray shadow(hit.hitpoint, wi);
                if (!scene.intersect(shadow)) {
                    color += glm::dot(shadow.d, hit.norm) * f * scene.sampleSkybox(shadow) / pdf;
                }
            }
        }
        */
        return color / (1.f * sr);
    }

    glm::vec3 Li(const Ray &r, Scene &scene, int depth) {
        if (depth > max_depth)
            return glm::vec3(0.f);
        if (depth > 3) {
            float rr = sampler.uniformSample1D();
            if (rr > threshold) {
                return glm::vec3(0.f);
            }
        }

        float t = FLT_MAX;
        Hit hit;
        if (!scene.intersect(r, t, hit)) {
            return scene.sampleSkybox(r);
        }
        
        if (hit.material->is_light) {
            return hit.material->le;
        }
        
        //direct lighting
        glm::vec3 direct = directLight(r, scene, hit, 1);
        //indirect lighting
        glm::vec3 wi;
        float pdf;
        glm::vec3 f;
        hit.sampleBRDF(sampler, -r.d, wi, f, pdf);
        Ray ray2(hit.hitpoint, wi);
        glm::vec3 indirect = Li(ray2, scene, depth + 1);
        //return glm::clamp(direct + glm::dot(wi, hit.norm) * f * indirect / pdf, glm::vec3(0.f, 0.f, 0.f), glm::vec3(1.f, 1.f, 1.f));
        glm::vec3 lo =  hit.material->le + direct + glm::dot(wi, hit.norm) * f * indirect / pdf;
        if (depth > 3)
            lo = lo / threshold;
        return glm::clamp(lo, glm::vec3(0.f), glm::vec3(4.f));
    }
};

#endif