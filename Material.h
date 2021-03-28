#pragma once

#ifndef __MATERIAL_H
#define __MATERIAL_H

#include <string>
#include <glm/glm.hpp>

#include <glm/gtx/component_wise.hpp>

#include "Texture.h"
#include "Sampler.h"

struct Material {
    std::string name;

    glm::vec3 ka;
    glm::vec3 ks;
    glm::vec3 kd;
    glm::vec3 le;
    int illum;
    float shininess;
    float ior;
    bool is_light;
    Texture diffuse_tex;
    float kappa;
    float eplison;
    float invPi;
    int material_type; //0 for lambertian, 1 for phong lobe, 2 for glossy

    Material(std::string &name, glm::vec3 &ka, glm::vec3 &kd, glm::vec3 &ks,
             float shininess, float ior, int illum) {
        this->name = name;
        this->ka = ka;
        this->ks = ks;
        this->kd = kd;
        this->shininess = shininess;
        this->illum = illum;
        this->ior = ior;
        this->is_light = false;
        this->le = glm::vec3(0.f);
        this->diffuse_tex = Texture();
        // assume the material is not ideal black
        // e.g kd = 0 and ks = 0
        float rho1 = glm::length(kd);
        float rho2 = glm::length(ks);
        rho1 *= rho1;
        rho2 *= rho2;
        kappa = rho1 / (rho1 + rho2);
        eplison = 1e-6;
        invPi = 1.f / M_PI;
        if (glm::compMax(ks) < eplison)
            material_type = 0;
        else if (glm::compMax(kd) < eplison)
            material_type = 1;
        else
            material_type = 3;
    }

    void setLe(glm::vec3 le) {
        this->le = le;
        is_light = true;
    }

    glm::vec3 getBRDF(const glm::vec3 &wi, const glm::vec3 &wo, const glm::vec3 &norm) {
        if (dot(wo, norm) < 1e-6)
            return glm::vec3(0.f, 0.f, 0.f);
        glm::vec3 r = glm::normalize(2 * glm::dot(wi, norm) * norm  - wi);
        return 1.f / (float) M_PI * kd; //+ (shininess + 2) / (float) (2 * M_PI) * std::pow(std::max(0.f, glm::dot(r, wo)), shininess) * ks;
    }


    // sampling reference
    // http://www.pbr-book.org/3ed-2018/Monte_Carlo_Integration/2D_Sampling_with_Multidimensional_Transformations.html
    void uniformSampleHemisphere(const glm::vec2 &u, const glm::vec3 &wo, glm::vec3 &wi, float &pdf) {
        // assume wi, wo in local coordinate
        // u[i] is a r.v. with the distribution of U(0, 1)
        float r = std::max(1.f - u[0] * u[0], 0.f);
        wi.x = cosf(2 * M_PI * u[1]) * sqrt(r);
        wi.y = sinf(2 * M_PI * u[1]) * sqrt(r);
        wi.z = u[0];
        pdf = 0.5f / M_PI;
    }

    void cosineWeighedSampleHemisphere(const glm::vec2 &u, glm::vec3 &w, float &pdf) {
        // for diffuse term
        float r = std::max(0.f, 1 - u[0]);
        w.x = cosf(2 * M_PI * u[1]) * sqrt(r);
        w.y = sinf(2 * M_PI * u[1]) * sqrt(r);
        w.z = sqrt(u[0]);
        pdf = w.z / M_PI;
    }

    void cosineNWeighedSampleHemisphere(const glm::vec2 &u, glm::vec3 &w, float &pdf) {
        // for specular term
        float cos_theta = pow(u[0], 1.f / (1 + shininess));
        float sin_theta = sqrt(std::max(0.f, 1 - cos_theta * cos_theta));
        w.x = cosf(2 * M_PI * u[1]) * sin_theta;
        w.y = sinf(2 * M_PI * u[1]) * sin_theta;
        w.z = cos_theta;
        pdf = (shininess + 1) * pow(cos_theta, shininess) / (2 * M_PI);
    }


    float PDF(const glm::vec3 &wi, const glm::vec3 &wo) {
        // Multi importance sampling
        float pdf;
        float pdf1 = 0.f, pdf2 = 0.f;
        glm::vec3 wh = glm::normalize(wi + wo);
        
        pdf1 = (shininess + 1) * pow(wh.z, shininess) / (2 * M_PI);
        pdf2 = wi.z * invPi;
        if (material_type == 0) {
            assert(pdf2 > eplison);
            return pdf2;
        } else if (material_type == 1) {
            assert(pdf1 > eplison);
            return pdf1;
        } else {
            assert(pdf1 + pdf2 > eplison);
            return 0.5  * pdf1 + 0.5 * pdf2;
        }
    }

    void sampleF(Sampler &sampler, const glm::vec3 &wo, glm::vec3 &wi, glm::vec3 &f, float &pdf) {
        glm::vec2 u = sampler.uniformSample2D();
        // multi importance sampling
        if (material_type == 1) {
            // kd ---> 0
            glm::vec3 H;
            cosineNWeighedSampleHemisphere(u, H, pdf);
            wi = glm::normalize(2.f * H - wo);
            //double C = (shininess + 2) * (shininess + 4) / (8 * M_PI * (pow(2, - shininess / 2.f) + shininess));
            //f = (float) (pdf * 2 *  M_PI / (shininess + 1.f) *  C) * ks;
            //f = (shininess + 8) * invPi / 8 * ks * std::max((float) pow(H.z, shininess), 0.f);
            f = F(wi, wo);
        } else if (material_type == 0) {
            // ks ---> 0
            cosineWeighedSampleHemisphere(u, wi, pdf);
            //f = invPi * kd;
            f = F(wi, wo);
        } else {
            float t = sampler.uniformSample1D();
            float pdf1 = 0;
            float pdf2 = 0;
            glm::vec3 wi1, wi2;
            glm::vec3 R(-wo.x, -wo.y, wo.z);

            cosineWeighedSampleHemisphere(u, wi1, pdf1);

            cosineNWeighedSampleHemisphere(u, wi2, pdf2);
            wi2 = constructFrame(R) * wi2;
            if (t < 0.5) {
                wi = wi1;
            } else {
                wi = wi2;
            }
            f = F(wi, wo);
            pdf = (wi.z  > 0) ? (pdf1 + pdf2) * 0.5 : 1.f;
            /*
            if (t < 0.5) {
                cosineWeighedSampleHemisphere(u, wi, pdf);
                //f = invPi * kd;
                //pdf *= kappa;
            } else {
                glm::vec3 H, R;
                //cosineNWeighedSampleHemisphere(u, H, pdf);
                //H = glm::normalize(H);
                //wi = glm::normalize(2.f * H - wo);
                R = glm::vec3(-wo.x, -wo.y, wo.z);
                //in the coordinate frame of R
                cosineNWeighedSampleHemisphere(u, wi, pdf);
                wi = constructFrame(R) * wi;

                assert(wi.z > 0);
                //double C = (shininess + 2) * (shininess + 4) / (8 * M_PI * (pow(2, - shininess / 2.f) + shininess));
                //f = (float) (pdf * 2 *  M_PI / (shininess + 1.f) *  C) * ks;
                //pdf *= (1 - kappa);
                //f = (shininess + 8) * invPi / 8 * ks * std::max((float) pow(H.z, shininess), 0.f);
            }
            f = F(wi, wo);
            pdf = PDF(wi, wo);
            */
        }
    }
    void sampleF(Sampler &sampler, const glm::vec3 &wo, glm::vec3 &wi, glm::vec3 &f, float &pdf, glm::vec2 uv) {
        glm::vec2 u = sampler.uniformSample2D();
        // multi importance sampling
        if (material_type == 1) {
            glm::vec3 H;
            cosineNWeighedSampleHemisphere(u, H, pdf);
            wi = glm::normalize(2.f * H - wo);
        } else if (material_type == 0) {
            cosineWeighedSampleHemisphere(u, wi, pdf);
        } else {
            float t = sampler.uniformSample1D();
            float pdf1 = 0;
            float pdf2 = 0;
            glm::vec3 wi1, wi2;
            glm::vec3 R(-wo.x, -wo.y, wo.z);
            cosineWeighedSampleHemisphere(u, wi1, pdf1);
            cosineNWeighedSampleHemisphere(u, wi2, pdf2);
            wi2 = constructFrame(R) * wi2;
            if (t < 0.5) {
                wi = wi1;
            } else {
                wi = wi2;
            }
            
            pdf = (wi.z  > 0) ? (pdf1 + pdf2) * 0.5 : 1.f;
        }
        f = F(wi, wo, uv);
    }
    void sampleDiffuse(Sampler &sampler, glm::vec3 &wi, glm::vec3 &f, float pdf) {
        glm::vec2 u = sampler.uniformSample2D();
        cosineWeighedSampleHemisphere(u, wi, pdf);
        f = kd * invPi;
    }

    glm::vec3 F(const glm::vec3 &wi, const glm::vec3 &wo) {
        // evaluate fr(wi, wo)
        glm::vec3 h = glm::normalize(wi + wo);
        glm::vec3 R = glm::vec3(-wo.x, -wo.y, wo.z);
        if (wi.z < 0 || wi.z * wo.z < 0)
            return glm::vec3(0.f);
        //return this->kd * invPi + (shininess + 8) * invPi / 8 * ks * std::max((float) pow(h.z, shininess), 0.f);
        //return this->kd + ks * std::max((float) pow(h.z, shininess), 0.f);
        return  this->kd * invPi + (shininess + 2) * invPi / 2 * ks * std::max((float) pow(glm::dot(R, wi), shininess), 0.f);
        //return  this->kd +  ks * std::max((float) pow(glm::dot(R, wi), shininess), 0.f);
    }

    glm::vec3 F(const glm::vec3 &wi, const glm::vec3 &wo, const glm::vec2 &uv) {
        glm::vec3 h = glm::normalize(wi + wo);
        glm::vec3 R = glm::vec3(-wo.x, -wo.y, wo.z);
        if (wi.z < 0 || wi.z * wo.z < 0)
            return glm::vec3(0.f);
        return  this->kd * this->diffuse_tex.sampleColor(uv.x, uv.y) * invPi + (shininess + 2) * invPi / 2 * ks * std::max((float) pow(glm::dot(R, wi), shininess), 0.f);
    }

    glm::mat3x3 constructFrame(glm::vec3 &z) {
        glm::vec3 N = glm::normalize(z);
        glm::vec3 T, B;
        if (std::fabs(N.x) > std::fabs(N.y))
            T = glm::normalize(glm::vec3(N.z, 0.f, -N.x));
        else
            T = glm::normalize(glm::vec3(0.f, N.z, -N.y));
        B = glm::cross(N, T);
        return glm::mat3x3(T, B, N);
    }
};

#endif