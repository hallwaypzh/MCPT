/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.
    This file is part of pbrt.
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:
    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

//  my kdtree (both construction as well as intersection with ray) is modified from PBR's[Physically Based Rendering] implementation
//  http://www.pbr-book.org/3ed-2018/Primitives_and_Intersection_Acceleration/Kd-Tree_Accelerator.html


#pragma once
#ifndef __KDTREE_H
#define __KDTREE_H

#include <vector>
#include <algorithm>
#include <memory>

#include <glm/glm.hpp>

#include "Triangle.h"
#include "Ray.h"
#include "AABB.h"
#include "Hit.h"


struct KdAccelNode;
struct BoundEdge;

class KdTreeAccel {
  public:
    // KdTreeAccel Public Methods
    KdTreeAccel(std::vector<Triangle> &p,
                int isectCost = 80, int traversalCost = 1,
                float emptyBonus = 0.5, int maxPrims = 1, int maxDepth = -1);
    AABB WorldBound() const { return bounds; }
    ~KdTreeAccel();
    //bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
    bool IntersectP(const Ray &ray) const;
    bool Intersect(const Ray &ray, float &t) const;
    bool Intersect(const Ray &ray, float &t, Hit &hit) const;
  private:
    // KdTreeAccel Private Methods
    void buildTree(int nodeNum, const AABB &bounds,
                   const std::vector<AABB> &primBounds, int *primNums,
                   int nprims, int depth,
                   const std::unique_ptr<BoundEdge[]> edges[3], int *prims0,
                   int *prims1, int badRefines = 0);

    // KdTreeAccel Private Data
    const int isectCost, traversalCost, maxPrims;
    const float emptyBonus;
    std::vector<Triangle> primitives;
    std::vector<int> primitiveIndices;
    KdAccelNode *nodes;
    int nAllocedNodes, nextFreeNode;
    AABB bounds;
};

struct KdToDo {
    const KdAccelNode *node;
    float tMin, tMax;
};
#endif