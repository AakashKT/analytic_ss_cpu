#include "integrators/ltc+silhouette.h"
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

LTCSilhouette::LTCSilhouette(const Bounds2i &pixelBounds, std::shared_ptr<const Camera> camera, std::shared_ptr<Sampler> sampler)
    : SamplerIntegrator(camera, sampler, pixelBounds) {}

void LTCSilhouette::Preprocess(const Scene &scene, Sampler &sampler) {
    
}

Float LTCSilhouette::integrateEdge(Vector3f v1, Vector3f v2) const {
	Float t1 = acos( Clamp( Dot(v1, v2), -1, 1 ) );
	Vector3f cp = Cross(v1, v2);
	Float t2 = cp.z;

	if(t1 > 0.0001)
		return t2*t1/sin(t1);
	else
		return t2;
}

Float LTCSilhouette::integrateEdges(std::vector<std::pair<Vector3f, Vector3f>> &edges, Vector3f cg) const {
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

	return std::abs(irradiance);
}

Vector3f LTCSilhouette::clipToHorizon(SurfaceInteraction &si, std::vector<std::pair<Vector3f, Vector3f>> &edges, Vector3f cg) const {
	if(edges.size() == 0)
		return cg;
	
	std::vector<std::pair<Vector3f, Vector3f>> ogEdges(edges);
	std::vector<Vector3f> tempEdge;
	edges.clear();

	Vector3f h1(0.f, 1.f, 0.f);
	Vector3f h2(1.f, 0.f, 0.f);
	Vector3f n1(0.f, 0.f, 1.f);

	for(auto edge : ogEdges) {
		Float d1 = Dot(edge.first, n1);
		Float d2 = Dot(edge.second, n1);
		if(d1 <= 0.f && d2 <= 0.f)
			continue;
		else if(d1 > 0.f && d2 > 0.f)
			edges.push_back(edge);
		else {
			std::vector<Vector3f> newEdge;
			Vector3f n2 = Cross(edge.first, edge.second);

			std::vector<Vector3f> currentList;
			currentList.push_back(edge.first);
			currentList.push_back(edge.second);
	
			for(auto vertex : currentList) {
				Float dp = Dot(vertex, n1);

				if(dp > 0.f)
					newEdge.push_back(vertex);
				else if(dp == 0) {
					newEdge.push_back(vertex);
					tempEdge.push_back(vertex);
				}
				else {
					Vector3f l = Cross(n1, n2);
					l.z = 0.f;
					Vector3f i1 = l / l.Length();
					Vector3f i2 = -i1;

					if(Dot(i1, (currentList[0]+currentList[1])/2.0) >= 0.f) {
						newEdge.push_back(i1);
						tempEdge.push_back(i1);
					}
					else {
						newEdge.push_back(i2);
						tempEdge.push_back(i2);
					}
				}
			}

			if(tempEdge.size() >= 2) {
				std::pair<Vector3f, Vector3f> temp;
				temp.first = tempEdge[0];
				temp.second = tempEdge[1];
				edges.push_back(temp);
				
				tempEdge.clear();
			}
			
			std::pair<Vector3f, Vector3f> temp;
			temp.first = newEdge[0];
			temp.second = newEdge[1];
			edges.push_back(temp);
		}
	}

	cg = Vector3f(0.f, 0.f, 0.f);
	if(edges.size() == 0)
		return cg;

	for(auto edge : edges)
		cg += edge.first;
	cg /= (Float)edges.size();

	return cg;
}

Vector3f LTCSilhouette::applyLocalTransform(SurfaceInteraction &si, std::vector<std::pair<Vector3f, Vector3f>> &edges, Vector3f cg) const {
	Vector3f misoc1 = si.bsdf->WorldToLocal(si.wo);
	misoc1.z = 0.f;
	misoc1 = misoc1 / misoc1.Length();
	Vector3f misoc3(0.f, 0.f, 1.f);
	Vector3f misoc2 = Cross(misoc3, misoc1);
	misoc2 = misoc2 / misoc2.Length(); 

	std::vector<bool> processed;

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

		processed.push_back(false);
	}

	Vector3f cg_ = si.bsdf->WorldToLocal(cg);
	cg.x = Dot(misoc1, cg_);
	cg.y = Dot(misoc2, cg_);
	cg.z = Dot(misoc3, cg_);
	cg = cg / cg.Length();

	if(edges.size() == 0)
		return cg;

	// Order clockwise w.r.t. cg
	for(int i=0; i<edges.size(); i++) {
		Vector3f v1 = edges[i].first;
		Vector3f v2 = edges[i].second;
		Vector3f cp = Cross(v1-cg, v2-cg);
		// cp = cp / cp.Length();

		if(Dot(-cp, cg) >= 0) {
			edges[i].first = v2;
			edges[i].second = v1;
		}
	}

	// Make continuous 
	std::vector<std::pair<Vector3f, Vector3f>> ogEdges(edges);
	edges.clear();

	auto currentEdge = ogEdges[0];
	processed[0] = true;
	edges.push_back(currentEdge);

	while(true) {
		if(edges.size() == ogEdges.size())
			break;

		Float minimum = -1.f;
		int minimumIdx = -1;
		for(int i=0; i<ogEdges.size(); i++) {
			if(processed[i])
				continue;

			auto edge = ogEdges[i];

			Vector3f v1 = currentEdge.second - cg;
			v1 = v1 / v1.Length();
			Vector3f v2 = edge.first - cg;
			v2 = v2 / v2.Length();
			Float angle = acos( Clamp( Dot(v1, v2), -1, 1 ) );

			if(minimum == -1.f) {
				minimum = angle;
				minimumIdx = i;
			}
			else if(angle <= minimum) {
				minimum = angle;
				minimumIdx = i;
			}
		}

		auto nextEdge = ogEdges[minimumIdx];
		nextEdge.first = currentEdge.second;
		edges.push_back(nextEdge);
		processed[minimumIdx] = true;
		currentEdge = nextEdge;
	}

	edges[edges.size()-1].second = edges[0].first;

	return cg;
}

Vector3f LTCSilhouette::applyLTC(SurfaceInteraction &si, std::vector<std::pair<Vector3f, Vector3f>> &edges, Vector3f cg, Float alpha_, Float *amp) const {
	if(edges.size() == 0)
		return cg;
		
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
		edges[i].first = edges[i].first / edges[i].first.Length();

		Vector3f v2 = edges[i].second;
		edges[i].second.x = Dot(c1, v2);
		edges[i].second.y = Dot(c2, v2);
		edges[i].second.z = Dot(c3, v2);
		edges[i].second = edges[i].second / edges[i].second.Length();
	}

	Vector3f cg_(cg);
	cg.x = Dot(c1, cg_);
	cg.y = Dot(c2, cg_);
	cg.z = Dot(c3, cg_);
	cg = cg / cg.Length();
	return cg;
}

Vector3f LTCSilhouette::projectToUnitSphere(std::vector<std::pair<Vector3f, Vector3f>> &edges, Vector3f cg) const {
	for(int i=0; i<edges.size(); i++) {
		edges[i].first = edges[i].first / edges[i].first.Length();
		edges[i].second = edges[i].second / edges[i].second.Length();
	}

	return cg / cg.Length();
}

Spectrum LTCSilhouette::LiWrite(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, std::ofstream &file, int depth) const {

	return this->Li(ray, scene, sampler, arena, depth);
}

Spectrum LTCSilhouette::Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth) const {

	ProfilePhase p(Prof::SamplerIntegratorLi);
	Spectrum L(0.0);

	SurfaceInteraction isect;
	bool foundIntersection = scene.Intersect(ray, &isect);
	
	if(!foundIntersection) {

	}
	else if(foundIntersection) {
		Spectrum lightRadiance = isect.Le(-ray.d);
		isect.ComputeScatteringFunctions(ray, arena, true);

		auto shadingPoint = Vector3f(isect.p) + 1e-3 * Vector3f(isect.n);

		if(lightRadiance.y() != 0.0f || !isect.bsdf) {
			L += lightRadiance;
		}
		else {
			std::vector< std::vector<std::pair<Vector3f, Vector3f>> > lightsEdges;
			std::vector<Vector3f> lightsCg;
			int numLights = scene.areaLightMeshes.size();
			for (size_t j = 0; j < numLights; ++j) {
				const std::shared_ptr<TriangleMesh> &lightMesh = scene.areaLightMeshes[j];
				Vector3f avgNormal(0.f, 0.f, 0.f);
				auto silhouetteEdges = lightMesh->computeSilhouetteEdges(shadingPoint, avgNormal);

				lightsEdges.push_back(silhouetteEdges);
				lightsCg.push_back(scene.areaLightMeshes[j]->cg - shadingPoint);
			}

			float alpha = 0.f;
			auto diffuse = isect.getParams(isect, &alpha);

			Spectrum irradiance = 0.0;
			int idx = 0;
			Float amplitude = 0.f;
			for(auto silhouette : lightsEdges) {			
				auto cg = this->projectToUnitSphere(silhouette, lightsCg[idx]);

				cg = this->applyLocalTransform(isect, silhouette, cg);
				cg = this->clipToHorizon(isect, silhouette, cg);
				cg = this->applyLTC(isect, silhouette, cg, alpha, &amplitude);
				irradiance += scene.areaLightRadiance[idx] * amplitude * integrateEdges(silhouette, cg);

				idx++;
			}
			irradiance = irradiance * 0.5 * InvPi;
			if(irradiance.y() >= 0.f)
				L += diffuse * irradiance;
		}
	}

	return L;
}

LTCSilhouette *CreateLTCSilhouette(const ParamSet &params,
                                     std::shared_ptr<Sampler> sampler,
                                     std::shared_ptr<const Camera> camera) {

	Bounds2i pixelBounds = camera->film->GetSampleBounds();

    return new LTCSilhouette(pixelBounds, camera, sampler);
}

}