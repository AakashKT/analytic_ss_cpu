
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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_SCENE_H
#define PBRT_CORE_SCENE_H

// core/scene.h*
#include "pbrt.h"
#include "geometry.h"
#include "primitive.h"
#include "light.h"

#include "lights/diffuse.h"
#include "shapes/triangle.h"
#include "core/geometry.h"
#include "accelerators/bvh.h"

namespace pbrt {

// Scene Declarations
class Scene {
  public:
    // Scene Public Methods
    Scene(std::shared_ptr<Primitive> aggregate,
          const std::vector<std::shared_ptr<Light>> &lights)
        : lights(lights), aggregate(aggregate) {
        // Scene Constructor Implementation
        worldBound = aggregate->WorldBound();
        
        std::vector<Vector3f> previousCgs;
        for (const auto &light : lights) {
            light->Preprocess(*this);
            if (light->flags & (int)LightFlags::Infinite)
                infiniteLights.push_back(light);
            
            if(light->flags & (int)LightFlags::Area) {
                const std::shared_ptr<DiffuseAreaLight> &areaLight = std::static_pointer_cast<DiffuseAreaLight>(light);
                const std::shared_ptr<Triangle> &shape = std::static_pointer_cast<Triangle>(areaLight->shape);
                const std::shared_ptr<TriangleMesh> &lightMesh = std::static_pointer_cast<TriangleMesh>(shape->mesh);

                bool found = false;
                for(auto cg : previousCgs)
                  if(cg == lightMesh->cg)
                    found = true;

                if(!found) {
                    previousCgs.push_back(lightMesh->cg);
                    this->areaLightShapes.push_back(shape);
                    this->areaLightMeshes.push_back(lightMesh);
                    this->areaLightRadiance.push_back(areaLight->Lemit);
                }
            }
        }

        const std::shared_ptr<BVHAccel> &agg = std::static_pointer_cast<BVHAccel>(this->aggregate);
        int numBlockers = agg->primitives.size();
        previousCgs.clear();
        for (int j=0; j < numBlockers; ++j) {
          const std::shared_ptr<GeometricPrimitive> &blockerPri = std::static_pointer_cast<GeometricPrimitive>(agg->primitives[j]);
          const std::shared_ptr<Triangle> &blockerTri = std::static_pointer_cast<Triangle>(blockerPri->shape);
          const std::shared_ptr<TriangleMesh> &blockerMesh = std::static_pointer_cast<TriangleMesh>(blockerTri->mesh);

          if(blockerMesh->isOccluder) {

            bool found = false;
            for(auto cg : previousCgs)
              if(cg == blockerMesh->cg)
                found = true;
            
            if(!found) {
              this->blockerMeshes.push_back(blockerMesh);
              previousCgs.push_back(blockerMesh->cg);
            }
          }
        }
    }
    const Bounds3f &WorldBound() const { return worldBound; }
    bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
    bool IntersectP(const Ray &ray) const;
    bool IntersectTr(Ray ray, Sampler &sampler, SurfaceInteraction *isect,
                     Spectrum *transmittance) const;

    // Scene Public Data
    std::vector<std::shared_ptr<Light>> lights;
    std::vector<std::shared_ptr<TriangleMesh>> areaLightMeshes;
    std::vector<std::shared_ptr<Triangle>> areaLightShapes;
    std::vector<std::shared_ptr<TriangleMesh>> blockerMeshes;
    std::vector<Spectrum> areaLightRadiance;
    // Store infinite light sources separately for cases where we only want
    // to loop over them.
    std::vector<std::shared_ptr<Light>> infiniteLights;

    std::shared_ptr<Primitive> aggregate;
    
  private:
    // Scene Private Data
    Bounds3f worldBound;
};

}  // namespace pbrt

#endif  // PBRT_CORE_SCENE_H
