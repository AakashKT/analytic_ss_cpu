#include "integrators/ratio.h"
#include "bssrdf.h"
#include "camera.h"
#include "film.h"
#include "interaction.h"
#include "paramset.h"
#include "scene.h"
#include "stats.h"

#include "lights/diffuse.h"
#include "shapes/triangle.h"

namespace pbrt {

Ratio::Ratio(const Bounds2i &pixelBounds, std::shared_ptr<const Camera> camera, std::shared_ptr<Sampler> sampler, bool visibility)
    : SamplerIntegrator(camera, sampler, pixelBounds), visibility(visibility) {}

void Ratio::Preprocess(const Scene &scene, Sampler &sampler) {
	// Compute number of samples to use for each light
	for (const auto &light : scene.lights)
		nLightSamples.push_back(sampler.RoundCount(light->nSamples));

	// Request samples for sampling all lights
	for (int i = 0; i < 1; ++i) {
		for (size_t j = 0; j < scene.lights.size(); ++j) {
			sampler.Request2DArray(nLightSamples[j]);
			sampler.Request2DArray(nLightSamples[j]);
		}
	}
}

Spectrum Ratio::LiWrite(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, std::ofstream &file, int depth) const {
	
	return this->Li(ray, scene, sampler, arena, depth);
}

Spectrum Ratio::Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth) const {

	ProfilePhase p(Prof::SamplerIntegratorLi);
    Spectrum L(0.f);
    // Find closest ray intersection or return background radiance
    SurfaceInteraction isect;
    if (!scene.Intersect(ray, &isect)) {
        for (const auto &light : scene.lights) L += light->Le(ray);
        return L;
    }

    // Compute scattering functions for surface interaction
    isect.ComputeScatteringFunctions(ray, arena);
    if (!isect.bsdf)
        return Li(isect.SpawnRay(ray.d), scene, sampler, arena, depth);
    Vector3f wo = isect.wo;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);
    if (scene.lights.size() > 0) {
        // Compute direct lighting for _DirectLightingIntegrator_ integrator
        if (visibility == true)
            L += UniformSampleAllLights(isect, scene, arena, sampler,
                                        nLightSamples);
        else
            L += UniformSampleAllLightsNoV(isect, scene, arena, sampler,
                                        nLightSamples);
    }
    if (depth + 1 < 1) {
        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, isect, scene, sampler, arena, depth);
        L += SpecularTransmit(ray, isect, scene, sampler, arena, depth);
    }
    return L;
}

Ratio *CreateRatio(const ParamSet &params,
                                     std::shared_ptr<Sampler> sampler,
                                     std::shared_ptr<const Camera> camera) {

	Bounds2i pixelBounds = camera->film->GetSampleBounds();

	bool visibility = params.FindOneBool("visibility", true);

    return new Ratio(pixelBounds, camera, sampler, visibility);
}

}