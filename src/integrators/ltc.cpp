#include "integrators/ltc.h"
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

LTC::LTC(const Bounds2i &pixelBounds, std::shared_ptr<const Camera> camera, std::shared_ptr<Sampler> sampler)
    : SamplerIntegrator(camera, sampler, pixelBounds) {}

void LTC::Preprocess(const Scene &scene, Sampler &sampler) {
    
}

Float LTC::integrateEdge(Vector3f v1, Vector3f v2) const {
	Float t1 = acos( Clamp( Dot(v1, v2), -1, 1 ) );
	Vector3f cp = Cross(v1, v2);
	Float t2 = cp.z;

	if(t1 > 0.0001)
		return t2*t1/sin(t1);
	else
		return t2;
}

Float LTC::integrateEdges(std::vector<std::pair<Vector3f, Vector3f>> &edges, Vector3f cg) const {
	Float irradiance = 0.0;

	for(auto edge : edges) {
		auto v1 = edge.first;
		auto v2 = edge.second;
		auto cp = Cross(v1-cg, v2-cg);

		if(Dot(cp, -cg) >= 0)
			irradiance += integrateEdge(v2, v1);
		else
			irradiance += integrateEdge(v1, v2);
	}

	return irradiance;
}

Vector3f LTC::applyLocalTransform(SurfaceInteraction &si, std::vector<std::pair<Vector3f, Vector3f>> &edges, Vector3f cg) const {
	Vector3f misoc1 = si.bsdf->WorldToLocal(si.wo);
	misoc1.z = 0.f;
	misoc1 = misoc1 / misoc1.Length();
	Vector3f misoc3(0.f, 0.f, 1.f);
	Vector3f misoc2 = Cross(misoc3, misoc1);
	misoc2 = misoc2 / misoc2.Length(); 

	for(int i=0; i<edges.size(); i++) {
		edges[i].first = si.bsdf->WorldToLocal(edges[i].first);
		Vector3f v1 = edges[i].first;
		edges[i].first.x = Dot(misoc1, v1);
		edges[i].first.y = Dot(misoc2, v1);
		edges[i].first.z = Dot(misoc3, v1);
		edges[i].first = edges[i].first / edges[i].first.Length();

		edges[i].second = si.bsdf->WorldToLocal(edges[i].second);
		Vector3f v2 = edges[i].second;
		edges[i].second.x = Dot(misoc1, v2);
		edges[i].second.y = Dot(misoc2, v2);
		edges[i].second.z = Dot(misoc3, v2);
		edges[i].second = edges[i].second / edges[i].second.Length();
	}

	Vector3f cg_ = si.bsdf->WorldToLocal(cg);
	cg.x = Dot(misoc1, cg_);
	cg.y = Dot(misoc2, cg_);
	cg.z = Dot(misoc3, cg_);
	cg = cg / cg.Length();
	return cg;
}

Vector3f LTC::applyLTC(SurfaceInteraction &si, std::vector<std::pair<Vector3f, Vector3f>> &edges, Vector3f cg, Float alpha_, Float *amp) const {
	Vector3f wo_local = si.bsdf->WorldToLocal(si.wo);
	Float theta = SphericalTheta(wo_local);
	
	Float theta_idx_f = theta * 63.0 * 2.0 / Pi;
	Float alpha_idx_f = alpha_ * 63.0;

	int theta_idx1 = floor(theta_idx_f);
	int theta_idx2 = ceil(theta_idx_f);
	int alpha_idx1 = floor(alpha_idx_f);
	int alpha_idx2 = ceil(alpha_idx_f);

	Float tw1 = theta_idx_f - theta_idx1;
	Float tw2 = 1.0 - tw1;
	Float aw1 = alpha_idx_f - alpha_idx1;
	Float aw2 = 1.0 - aw1;
		
	Vector3f c1;
	Vector3f c2;
	Vector3f c3;

	c1.x = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*64][0] + tw1 * tabMinv[alpha_idx1+theta_idx2*64][0]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*64][0] + tw1 * tabMinv[alpha_idx2+theta_idx2*64][0]);
	c1.y = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*64][3] + tw1 * tabMinv[alpha_idx1+theta_idx2*64][3]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*64][3] + tw1 * tabMinv[alpha_idx2+theta_idx2*64][3]);
	c1.z = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*64][6] + tw1 * tabMinv[alpha_idx1+theta_idx2*64][6]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*64][6] + tw1 * tabMinv[alpha_idx2+theta_idx2*64][6]);

	c2.x = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*64][1] + tw1 * tabMinv[alpha_idx1+theta_idx2*64][1]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*64][1] + tw1 * tabMinv[alpha_idx2+theta_idx2*64][1]);
	c2.y = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*64][4] + tw1 * tabMinv[alpha_idx1+theta_idx2*64][4]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*64][4] + tw1 * tabMinv[alpha_idx2+theta_idx2*64][4]);
	c2.z = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*64][7] + tw1 * tabMinv[alpha_idx1+theta_idx2*64][7]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*64][7] + tw1 * tabMinv[alpha_idx2+theta_idx2*64][7]);

	c3.x = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*64][2] + tw1 * tabMinv[alpha_idx1+theta_idx2*64][2]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*64][2] + tw1 * tabMinv[alpha_idx2+theta_idx2*64][2]);
	c3.y = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*64][5] + tw1 * tabMinv[alpha_idx1+theta_idx2*64][5]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*64][5] + tw1 * tabMinv[alpha_idx2+theta_idx2*64][5]);
	c3.z = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*64][8] + tw1 * tabMinv[alpha_idx1+theta_idx2*64][8]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*64][8] + tw1 * tabMinv[alpha_idx2+theta_idx2*64][8]);

	*amp = aw2 * (tw2 * tabAmplitude[alpha_idx1+theta_idx1*64] + tw1 * tabAmplitude[alpha_idx1+theta_idx2*64]) + aw1 * (tw2 * tabAmplitude[alpha_idx2+theta_idx1*64] + tw1 * tabAmplitude[alpha_idx2+theta_idx2*64]);

	for(int i=0; i<edges.size(); i++) {
		Vector3f v1 = edges[i].first;
		edges[i].first.x = Dot(c1, v1);
		edges[i].first.y = Dot(c2, v1);
		edges[i].first.z = Dot(c3, v1);
		if(edges[i].first.Length() != 0)
			edges[i].first = edges[i].first / edges[i].first.Length();

		Vector3f v2 = edges[i].second;
		edges[i].second.x = Dot(c1, v2);
		edges[i].second.y = Dot(c2, v2);
		edges[i].second.z = Dot(c3, v2);
		if(edges[i].second.Length() != 0)
			edges[i].second = edges[i].second / edges[i].second.Length();
	}

	Vector3f cg_(cg);
	cg.x = Dot(c1, cg_);
	cg.y = Dot(c2, cg_);
	cg.z = Dot(c3, cg_);
	if(cg.Length() != 0)
		cg = cg / cg.Length();
	return cg;
}

Vector3f LTC::projectToUnitSphere(std::vector<std::pair<Vector3f, Vector3f>> &edges, Vector3f cg) const {
	for(int i=0; i<edges.size(); i++) {
		edges[i].first = edges[i].first / edges[i].first.Length();
		edges[i].second = edges[i].second / edges[i].second.Length();
	}

	return cg / cg.Length();
}

Spectrum LTC::LiWrite(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, std::ofstream &file, int depth) const {
	
	return this->Li(ray, scene, sampler, arena, depth);
}

Spectrum LTC::Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth) const {

	ProfilePhase p(Prof::SamplerIntegratorLi);
	Spectrum L(0.0);

	SurfaceInteraction isect;
	bool foundIntersection = scene.Intersect(ray, &isect);
	
	if(foundIntersection) {
		Spectrum lightRadiance = isect.Le(-ray.d);
		isect.ComputeScatteringFunctions(ray, arena, true);

		if(lightRadiance.y() != 0.0f || !isect.bsdf) {
			L += lightRadiance;
		}
		else {
			float alpha = 0.f;
			auto diffuse = isect.getParams(isect, &alpha);

			int numLights = scene.areaLightMeshes.size();
			for (size_t j = 0; j < numLights; ++j) {
				float irradiance = 0.0;
				std::vector<std::pair<Vector3f, Vector3f>> edges;
				const std::shared_ptr<TriangleMesh> &lightMesh = scene.areaLightMeshes[j];

				for(int i=0; i<lightMesh->vertexIndices.size(); i+=3) {
					edges.clear();

					Normal3f fn = lightMesh->faceNormals[i/3];

					Vector3f v1 = Vector3f( lightMesh->p[ lightMesh->vertexIndices[i] ] );
					Vector3f v2 = Vector3f( lightMesh->p[ lightMesh->vertexIndices[i+1] ] );
					Vector3f v3 = Vector3f( lightMesh->p[ lightMesh->vertexIndices[i+2] ] );

					v1 = v1 - Vector3f(isect.p);
					v2 = v2 - Vector3f(isect.p);
					v3 = v3 - Vector3f(isect.p);

					Vector3f ecg = (v1+v2+v3) / 3.0;

					if(Dot(-ecg, fn) < 0)
						continue;

					edges.push_back(std::pair<Vector3f, Vector3f>(v1, v2));
					edges.push_back(std::pair<Vector3f, Vector3f>(v2, v3));
					edges.push_back(std::pair<Vector3f, Vector3f>(v3, v1));

					Float amplitude = 0.f;
					ecg = this->projectToUnitSphere(edges, ecg);
					ecg = this->applyLocalTransform(isect, edges, ecg);
					// ecg = this->projectToUnitSphere(edges, ecg);
					// amplitude = 1.0;
					ecg = this->applyLTC(isect, edges, ecg, alpha, &amplitude);
					
					irradiance += amplitude * integrateEdges(edges, ecg);
				}

				irradiance = irradiance * 0.5 * InvPi;
				if(irradiance >= 0.f)
					L += scene.areaLightRadiance[j] *  diffuse * irradiance;
			}
		}
	}

	return L;
}

LTC *CreateLTC(const ParamSet &params,
                                     std::shared_ptr<Sampler> sampler,
                                     std::shared_ptr<const Camera> camera) {

	Bounds2i pixelBounds = camera->film->GetSampleBounds();

    return new LTC(pixelBounds, camera, sampler);
}

}