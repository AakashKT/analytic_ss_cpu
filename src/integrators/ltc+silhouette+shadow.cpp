#include "integrators/ltc+silhouette+shadow.h"
#include "bssrdf.h"
#include "camera.h"
#include "film.h"
#include "interaction.h"
#include "paramset.h"
#include "scene.h"
#include "stats.h"

#include "lights/diffuse.h"
#include "shapes/triangle.h"
#include "core/geometry.h"
#include "accelerators/bvh.h"

namespace pbrt {

STAT_COUNTER("Integrator/Zero Blockers", zeroBlockers);
STAT_COUNTER("Integrator/One Blockers", oneBlockers);
STAT_COUNTER("Integrator/Two Blockers", twoBlockers);
STAT_COUNTER("Integrator/Three Blockers", threeBlockers);

float sgn(float x) {
	if(x < 0)
		return -1.f;
	else
		return 1.f;
}

LTCSilhouetteShadow::LTCSilhouetteShadow(const Bounds2i &pixelBounds, std::shared_ptr<const Camera> camera, std::shared_ptr<Sampler> sampler)
    : SamplerIntegrator(camera, sampler, pixelBounds) {}

void LTCSilhouetteShadow::Preprocess(const Scene &scene, Sampler &sampler) {
    
}

Float LTCSilhouetteShadow::integrateEdge(Vector3f v1, Vector3f v2) const {
	Float t1 = acos(Clamp(Dot(v1, v2), -1, 1));
	Vector3f cp = Cross(v1, v2);
	Float t2 = cp.z;

	if(t1 > 0.0001)
		return t2*t1/sin(t1);
	else
		return t2;

	// float x = Dot(v1, v2);
    // float y = std::abs(x);

    // float a = 0.8543985 + (0.4965155 + 0.0145206*y)*y;
    // float b = 3.4175940 + (4.1616724 + y)*y;
    // float v = a / b;

    // float theta_sintheta = (x > 0.0) ? v : 0.5*(1.0 / sqrt(std::max(1.0 - x*x, 1e-7)) ) - v;

	// auto rval = Cross(v1, v2)*theta_sintheta;

	// return rval.z;
}

Float LTCSilhouetteShadow::integrateEdges(std::vector<std::pair<Vector3f, Vector3f>> &edges) const {
	Float irradiance = 0.0;

	for(auto edge : edges) {
		auto v1 = edge.first;
		auto v2 = edge.second;

		irradiance += integrateEdge(v2, v1);
	}
	
	return std::abs(irradiance);
	// return irradiance;
}

void LTCSilhouetteShadow::clipToHorizon(SurfaceInteraction &si, std::vector<std::pair<Vector3f, Vector3f>> &edges) const {

	if(edges.size() == 0)
		return;

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
			
			bool firstBelow = false;
			int idx = 0;
			for(auto vertex : currentList) {
				Float dp = Dot(vertex, n1);

				if(dp > 0.f)
					newEdge.push_back(vertex);
				else if(dp == 0.f) {
					newEdge.push_back(vertex);
					tempEdge.push_back(vertex);

					if (idx == 0)
						firstBelow = true;
					else
						firstBelow = false;
				}
				else {
					if (idx == 0)
						firstBelow = true;
					else
						firstBelow = false;

					Vector3f l = Cross(n1, n2);
					l.z = 0.f;
					Vector3f i1 = l / l.Length();
					Vector3f i2 = -i1;

					Float th = acos(Dot(currentList[0], currentList[1]));

					Float th1 = acos(Dot(i1, currentList[0]));
					Float th2 = acos(Dot(i1, currentList[1]));

					Float d1 = std::abs((th1+th2) - th);

					th1 = acos(Dot(i2, currentList[0]));
					th2 = acos(Dot(i2, currentList[1]));

					Float d2 = std::abs((th1+th2) - th);

					if(d1 <= 1e-10) {
						newEdge.push_back(i1);
						tempEdge.push_back(i1);
					}
					else if(d2 <= 1e-10) {
						newEdge.push_back(i2);
						tempEdge.push_back(i2);
					}
					else {
						if(d1 <= d2) {
							newEdge.push_back(i1);
							tempEdge.push_back(i1);
						}
						else {
							newEdge.push_back(i2);
							tempEdge.push_back(i2);
						}
					}
				}
				idx++;
			}

			if (firstBelow) {
				if (tempEdge.size() >= 2) {
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
			else {
				std::pair<Vector3f, Vector3f> temp;
				temp.first = newEdge[0];
				temp.second = newEdge[1];
				edges.push_back(temp);

				if (tempEdge.size() >= 2) {
					std::pair<Vector3f, Vector3f> temp;
					temp.first = tempEdge[1];
					temp.second = tempEdge[0];

					edges.push_back(temp);
					tempEdge.clear();
				}
			}
		}
	}

	for(int i=0; i<edges.size(); i++) {
		if(edges[i].first.z == 0.f)
			edges[i].first.z = 1e-3;
		
		if(edges[i].second.z == 0.f)
			edges[i].second.z = 1e-3;
	}

	if (edges.size() <= 2)
		edges.clear();
}

Vector3f LTCSilhouetteShadow::applyLocalTransform(SurfaceInteraction &si, std::vector<std::pair<Vector3f, Vector3f>> &edges, Vector3f &avgNormal, Vector3f cg) const {
	
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

		edges[i].second = si.bsdf->WorldToLocal(edges[i].second);
		Vector3f v2 = edges[i].second;
		edges[i].second.x = Dot(misoc1, v2);
		edges[i].second.y = Dot(misoc2, v2);
		edges[i].second.z = Dot(misoc3, v2);
	}

	Vector3f cg_ = si.bsdf->WorldToLocal(cg);
	cg.x = Dot(misoc1, cg_);
	cg.y = Dot(misoc2, cg_);
	cg.z = Dot(misoc3, cg_);

	Vector3f avgNormal_ = si.bsdf->WorldToLocal(avgNormal);
	avgNormal.x = Dot(misoc1, avgNormal_);
	avgNormal.y = Dot(misoc2, avgNormal_);
	avgNormal.z = Dot(misoc3, avgNormal_);
	
	return cg;
}

void LTCSilhouetteShadow::applyLTC(SurfaceInteraction &si, std::vector<std::pair<Vector3f, Vector3f>> &edges, Float alpha_, Float *amp) const {
	if(edges.size() == 0)
		return;

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

	int ub = 64;
		
	Vector3f c1;
	Vector3f c2;
	Vector3f c3;

	c1.x = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*ub][0] + tw1 * tabMinv[alpha_idx1+theta_idx2*ub][0]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*ub][0] + tw1 * tabMinv[alpha_idx2+theta_idx2*ub][0]);
	c1.y = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*ub][3] + tw1 * tabMinv[alpha_idx1+theta_idx2*ub][3]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*ub][3] + tw1 * tabMinv[alpha_idx2+theta_idx2*ub][3]);
	c1.z = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*ub][6] + tw1 * tabMinv[alpha_idx1+theta_idx2*ub][6]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*ub][6] + tw1 * tabMinv[alpha_idx2+theta_idx2*ub][6]);

	c2.x = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*ub][1] + tw1 * tabMinv[alpha_idx1+theta_idx2*ub][1]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*ub][1] + tw1 * tabMinv[alpha_idx2+theta_idx2*ub][1]);
	c2.y = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*ub][4] + tw1 * tabMinv[alpha_idx1+theta_idx2*ub][4]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*ub][4] + tw1 * tabMinv[alpha_idx2+theta_idx2*ub][4]);
	c2.z = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*ub][7] + tw1 * tabMinv[alpha_idx1+theta_idx2*ub][7]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*ub][7] + tw1 * tabMinv[alpha_idx2+theta_idx2*ub][7]);

	c3.x = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*ub][2] + tw1 * tabMinv[alpha_idx1+theta_idx2*ub][2]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*ub][2] + tw1 * tabMinv[alpha_idx2+theta_idx2*ub][2]);
	c3.y = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*ub][5] + tw1 * tabMinv[alpha_idx1+theta_idx2*ub][5]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*ub][5] + tw1 * tabMinv[alpha_idx2+theta_idx2*ub][5]);
	c3.z = aw2 * (tw2 * tabMinv[alpha_idx1+theta_idx1*ub][8] + tw1 * tabMinv[alpha_idx1+theta_idx2*ub][8]) + aw1 * (tw2 * tabMinv[alpha_idx2+theta_idx1*ub][8] + tw1 * tabMinv[alpha_idx2+theta_idx2*ub][8]);

	*amp = aw2 * (tw2 * tabAmplitude[alpha_idx1+theta_idx1*ub] + tw1 * tabAmplitude[alpha_idx1+theta_idx2*ub]) + aw1 * (tw2 * tabAmplitude[alpha_idx2+theta_idx1*ub] + tw1 * tabAmplitude[alpha_idx2+theta_idx2*ub]);
	
	for(int i=0; i<edges.size(); i++) {
		Vector3f v1 = edges[i].first;
		edges[i].first.x = Dot(c1, v1);
		edges[i].first.y = Dot(c2, v1);
		edges[i].first.z = Dot(c3, v1);

		Vector3f v2 = edges[i].second;
		edges[i].second.x = Dot(c1, v2);
		edges[i].second.y = Dot(c2, v2);
		edges[i].second.z = Dot(c3, v2);
	}
}

void LTCSilhouetteShadow::projectToUnitSphere(std::vector<std::pair<Vector3f, Vector3f>> &edges) const {
	for(int i=0; i<edges.size(); i++) {
		edges[i].first = edges[i].first / edges[i].first.Length();
		edges[i].second = edges[i].second / edges[i].second.Length();

		if(SphericalTheta(edges[i].first) < Pi/2.0+1e-3 && SphericalTheta(edges[i].first) > Pi/2.0-1e-3) {
			edges[i].first = SphericalDirection(sin(Pi/2-1e-3), cos(Pi/2-1e-3), SphericalPhi(edges[i].first));
		}

		if(SphericalTheta(edges[i].second) < Pi/2.0+1e-3 && SphericalTheta(edges[i].second) > Pi/2.0-1e-3) {
			edges[i].second = SphericalDirection(sin(Pi/2-1e-3), cos(Pi/2-1e-3), SphericalPhi(edges[i].second));
		}
	}
}

void LTCSilhouetteShadow::unprojectXY(SurfaceInteraction &si, std::vector<std::pair<Vector3f, Vector3f>> &edges, Vector3f lookAt, bool deltaCorrection, Float delta) const {

	Vector3f N(0.f, 0.f, 1.f);
	Vector3f L(lookAt);
	L = L / L.Length();
	if(Dot(L, N) >= 0.95)
		N = Vector3f(0.f, 1.f, 0.f);

	Vector3f S = Cross(N, L);
	S = S / S.Length();
	Vector3f Ud = Cross(L, S);
	Ud = Ud / Ud.Length();

	Vector3f c1;
	c1.x = S.x;
	c1.y = Ud.x;
	c1.z = L.x;

	Vector3f c2;
	c2.x = S.y;
	c2.y = Ud.y;
	c2.z = L.y;

	Vector3f c3;
	c3.x = S.z;
	c3.y = Ud.z;
	c3.z = L.z;

	for(int i=0; i<edges.size(); i++) {
		Vector3f v1 = edges[i].first;
		v1.x *= v1.z;
		v1.y *= v1.z;
		if(deltaCorrection)
			v1.z = v1.z - 1.f - std::abs(delta);
		else
			v1.z = v1.z - 0.f;
		edges[i].first.x = Dot(c1, v1);
		edges[i].first.y = Dot(c2, v1);
		edges[i].first.z = Dot(c3, v1);

		Vector3f v2 = edges[i].second;
		v2.x *= v2.z;
		v2.y *= v2.z;
		if(deltaCorrection)
			v2.z = v2.z - 1.f - std::abs(delta);
		else
			v2.z = v2.z - 0.f;
		edges[i].second.x = Dot(c1, v2);
		edges[i].second.y = Dot(c2, v2);
		edges[i].second.z = Dot(c3, v2);
	}
}

Float LTCSilhouetteShadow::projectXY(SurfaceInteraction &si, std::vector<std::pair<Vector3f, Vector3f>> &edges, Vector3f lookAt, bool deltaCorrection) const {

	Vector3f N(0.f, 0.f, 1.f);
	Vector3f L(lookAt);
	L = L / L.Length();
	if(Dot(L, N) >= 0.95)
		N = Vector3f(0.f, 1.f, 0.f);

	Vector3f S = Cross(N, L);
	S = S / S.Length();
	Vector3f Ud = Cross(L, S);
	Ud = Ud / Ud.Length();

	Vector3f c1 = S;
	Vector3f c2 = Ud;
	Vector3f c3 = L;

	Float delta = 0.f;
	for(int i=0; i<edges.size(); i++) {
		Vector3f v1 = edges[i].first;
		edges[i].first.x = Dot(c1, v1);
		edges[i].first.y = Dot(c2, v1);
		edges[i].first.z = Dot(c3, v1);

		Vector3f v2 = edges[i].second;
		edges[i].second.x = Dot(c1, v2);
		edges[i].second.y = Dot(c2, v2);
		edges[i].second.z = Dot(c3, v2);

		delta = std::min(delta, edges[i].first.z);
		delta = std::min(delta, edges[i].second.z);
	}

	if(!deltaCorrection)
		delta = 0.f;

	float randomScale = rand()/Float(RAND_MAX/2.0) - 1;
	for(int i=0; i<edges.size(); i++) {
		if(deltaCorrection)
			edges[i].first.z += std::abs(delta) + 1.f;
		else
			edges[i].first.z += 0.f;	
		edges[i].first.x /= edges[i].first.z;
		edges[i].first.y /= edges[i].first.z;
		edges[i].first.x *= ( randomScale * 1e-3 + 1.0);
		edges[i].first.y *= ( randomScale * 1e-3 + 1.0);

		if(deltaCorrection)
			edges[i].second.z += std::abs(delta) + 1.f;
		else
			edges[i].second.z += 0.f;
		edges[i].second.x /= edges[i].second.z;
		edges[i].second.y /= edges[i].second.z;
		edges[i].second.x *= ( randomScale * 1e-3 + 1.0);
		edges[i].second.y *= ( randomScale * 1e-3 + 1.0);
	}

	return delta;
}

void LTCSilhouetteShadow::sortVertices(SurfaceInteraction &si, std::vector<std::pair<Vector3f, Vector3f>> &edges, Vector3f cg) const {

	if(edges.size() == 0)
		return;

	std::vector<std::pair<Vector3f, Vector3f>> ogEdges(edges);
	edges.clear();

	std::vector<bool> processed;
	for(int i=0; i<ogEdges.size(); i++)
		processed.push_back(false);

	Float delta = this->projectXY(si, ogEdges, cg, false);

	auto currentEdge = ogEdges[0];
	processed[0] = true;
	edges.push_back(currentEdge);

	while(true) {
		if(edges.size() == ogEdges.size())
			break;

		Float minimum = -1.f;
		int minimumIdx = -1;
		bool isReverse = false;
		for(int i=0; i<ogEdges.size(); i++) {
			if(processed[i])
				continue;

			auto edge = ogEdges[i];

			Vector2f v(currentEdge.second.x, currentEdge.second.y);

			Vector2f v1(edge.first.x, edge.first.y);
			Vector2f v2(edge.second.x, edge.second.y);

			Vector2f diff = v1-v;
			// Float dist1 = diff.Length();
			Float dist1 = sqrt( diff.x*diff.x + diff.y*diff.y );
			diff = v2-v;
			// Float dist2 = diff.Length();
			Float dist2 = sqrt( diff.x*diff.x + diff.y*diff.y );

			if(minimum == -1.f) {
				if(dist1 < dist2) {
					isReverse = false;
					minimum = dist1;
				}
				else {
					isReverse = true;
					minimum = dist2;
				}
				minimumIdx = i;
			}
			else {
				if(dist1 <= dist2) {
					if(dist1 <= minimum) {
						isReverse = false;
						minimum = dist1;
						minimumIdx = i;
					}
				}
				else if(dist2 <= dist1) {
					if(dist2 <= minimum) {
						isReverse = true;
						minimum = dist2;
						minimumIdx = i;
					}
				}
			}
		}

		auto nextEdge = ogEdges[minimumIdx];
		if(isReverse) {
			nextEdge.first = ogEdges[minimumIdx].second;
			nextEdge.second = ogEdges[minimumIdx].first;
		}

		nextEdge.first = currentEdge.second;
		edges.push_back(nextEdge);
		processed[minimumIdx] = true;

		currentEdge = nextEdge;
	}

	edges[edges.size()-1].second = edges[0].first;

	this->unprojectXY(si, edges, cg, false, delta);
}

bool LTCSilhouetteShadow::isVisible(SurfaceInteraction &si, std::vector<std::pair<Vector3f, Vector3f>> &lightEdges, Float lightRadius, Vector3f lightsAvgNormal, 
										std::vector<std::pair<Vector3f, Vector3f>> &blockerEdges, Float blockerRadius) const {
	Vector3f lcg(0.f, 0.f, 0.f);
	for(auto edge : lightEdges) {
		lcg += edge.first;
	}
	lcg /= (Float) lightEdges.size();
	lcg = -lcg;

	Vector3f bcg(0.f, 0.f, 0.f);
	for(auto edge : blockerEdges) {
		bcg += edge.first;
	}
	bcg /= (Float) blockerEdges.size();

	auto corners = computeFarPlaneCorners(si, lightEdges, lightRadius, -lcg);
	Vector3f tr = corners[0];
	Vector3f tl = corners[1];
	Vector3f bl = corners[2];
	Vector3f br = corners[3];

	// Remove blockers that lie beyond far plane 
	Vector3f farPlaneVec = bcg + lcg;
	farPlaneVec = farPlaneVec / farPlaneVec.Length();

	Vector3f topPlaneNormal = Normalize(Cross(tr, tl-tr));
	Vector3f bottomPlaneNormal = Normalize(Cross(br, br-bl));
	Vector3f rightPlaneNormal = Normalize(Cross(br, tr-br));
	Vector3f leftPlaneNormal = Normalize(Cross(bl, bl-tl));

	bool farPlane = Dot(lightsAvgNormal, farPlaneVec) >= 0 ? true : false;

	Float temp = 0.f;
	if(farPlane) {
		temp = Dot(bcg, topPlaneNormal);
		bool topPlane = ( temp > 0 || (temp <= 0 && std::abs(temp) <= blockerRadius) ) ? true : false;
		return true;
		if(topPlane) {
			temp = Dot(bcg, bottomPlaneNormal);
			bool bottomPlane = ( temp > 0 || (temp <= 0 && std::abs(temp) <= blockerRadius) ) ? true : false;
			// return true;
			if(bottomPlane) {
				temp = Dot(bcg, rightPlaneNormal);
				bool rightPlane = ( temp > 0 || (temp <= 0 && std::abs(temp) <= blockerRadius) ) ? true : false;
				if(rightPlane) {
					temp = Dot(bcg, leftPlaneNormal);
					bool leftPlane = ( temp > 0 || (temp <= 0 && std::abs(temp) <= blockerRadius) ) ? true : false;
					// return true;
					if(leftPlane) {
						return true;
					}
				}
			}
		}
	}

	return false;

	// bool topPlane = ( Dot(bcg, topPlaneNormal) >= 0 || Dot(bcg, topPlaneNormal) < blockerRadius ) ? true : false;
	// bool bottomPlane = ( Dot(bcg, bottomPlaneNormal) >= 0 || Dot(bcg, bottomPlaneNormal) < blockerRadius ) ? true : false;
	// bool rightPlane = ( Dot(bcg, rightPlaneNormal) >= 0 || Dot(bcg, rightPlaneNormal) < blockerRadius ) ? true : false;
	// bool leftPlane = ( Dot(bcg, leftPlaneNormal) >= 0 || Dot(bcg, leftPlaneNormal) < blockerRadius ) ? true : false;

	// return farPlane;
	// return farPlane && topPlane && bottomPlane && rightPlane && leftPlane;
}

std::vector<Vector3f> LTCSilhouetteShadow::computeFarPlaneCorners(SurfaceInteraction &si, std::vector<std::pair<Vector3f, Vector3f>> edges, Float radius, Vector3f cg) const {
	// Returns Top right, Top Left, Bottom left, Bottom Right
	std::vector<Vector3f> corners(4);

	Vector3f N(0.f, 0.f, 1.f);
	Vector3f L(cg);
	L = L / L.Length();
	if(Dot(L, N) >= 0.95)
		N = Vector3f(0.f, 1.f, 0.f);

	Vector3f S = Cross(N, L);
	S = S / S.Length();
	Vector3f Ud = Cross(L, S);
	Ud = Ud / Ud.Length();

	corners[0] = cg + radius * Ud + radius * S;
	corners[1] = cg + radius * Ud - radius * S;
	corners[2] = cg - radius * Ud - radius * S;
	corners[3] = cg - radius * Ud + radius * S;

	return corners;

	// Vector3f cg(0.f, 0.f, 0.f);
	// for(auto edge : edges) {
	// 	cg += edge.first;
	// }
	// cg /= (Float) edges.size();

	// // Camera look at toward CG ; Forward transform
	// Float delta = this->projectXY(si, edges, cg, true);

	// Float maxx, maxy, minx, miny, maxz, minz;
	// for(int i=0; i<edges.size(); i++) {
	// 	if(i == 0) {
	// 		maxx = edges[i].first.x;
	// 		maxy = edges[i].first.y;
	// 		maxz = edges[i].first.z;
	// 		minx = edges[i].first.x;
	// 		miny = edges[i].first.y; 
	// 		minz = edges[i].first.z; 
	// 	}
	// 	else {
	// 		maxx = std::max(maxx, edges[i].first.x);
	// 		maxy = std::max(maxy, edges[i].first.y);
	// 		maxz = std::max(maxz, edges[i].first.z);
	// 		minx = std::min(minx, edges[i].first.x);
	// 		miny = std::min(miny, edges[i].first.y);
	// 		minz = std::min(minz, edges[i].first.z);
	// 	}
	// }

	// std::vector<std::pair<Vector3f, Vector3f>> tempCorners;
	// tempCorners.push_back(std::pair<Vector3f, Vector3f>(Vector3f(maxx, maxy, minz), Vector3f(maxx, maxy, minz)));
	// tempCorners.push_back(std::pair<Vector3f, Vector3f>(Vector3f(minx, maxy, minz), Vector3f(minx, maxy, minz)));
	// tempCorners.push_back(std::pair<Vector3f, Vector3f>(Vector3f(minx, miny, minz), Vector3f(minx, miny, minz)));
	// tempCorners.push_back(std::pair<Vector3f, Vector3f>(Vector3f(maxx, miny, minz), Vector3f(maxx, miny, minz)));

	// // Inverse Transform
	// this->unprojectXY(si, tempCorners, cg, true, delta);

	// corners[0] = tempCorners[0].first;
	// corners[1] = tempCorners[1].first;
	// corners[2] = tempCorners[2].first;
	// corners[3] = tempCorners[3].first;

	// return corners;
}

Spectrum LTCSilhouetteShadow::LiWrite(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, std::ofstream &logFile, int depth) const {

	///* LOG */
	// logFile << "P " << std::to_string(sampler.currentPixel.x) << " " << std::to_string(sampler.currentPixel.y) << std::endl;
	// std::cout << "P " << std::to_string(sampler.currentPixel.x) << " " << std::to_string(sampler.currentPixel.y) << std::endl;
	///* ENDLOG */
				
	ProfilePhase p(Prof::SamplerIntegratorLi);
	Spectrum L(0.0);

	SurfaceInteraction isect;
	bool foundIntersection = scene.Intersect(ray, &isect);
	
	if(!foundIntersection) {
		// LOG		
		// logFile << "0" << std::endl;
		// logFile << "0" << std::endl;
		// logFile << "0" << std::endl;
		// logFile << "0" << std::endl;
		// logFile << "0" << std::endl;
		// logFile << "0" << std::endl;
		// ENDLOG
	}
	else if(foundIntersection) {
		Spectrum lightRadiance = isect.Le(-ray.d);
		isect.ComputeScatteringFunctions(ray, arena, true);

		Triangle *isectShape = (Triangle *) isect.shape;
    	const std::shared_ptr<TriangleMesh> &isectMesh = std::static_pointer_cast<TriangleMesh>(isectShape->mesh);

		if(lightRadiance.y() != 0.0f || isectMesh->isLight) {
			L += lightRadiance;

			// LOG		
			// logFile << "0" << std::endl;
			// logFile << "0" << std::endl;
			// logFile << "0" << std::endl;
			// logFile << "0" << std::endl;
			// logFile << "0" << std::endl;
			// logFile << "0" << std::endl;
			// ENDLOG
		}
		else {
			/* Initialize */
			Vector3f projectionVector(0.f, 0.f, 1.f);

			float alpha = 0.f;
			auto diffuse = isect.getParams(isect, &alpha);
			// alpha = alpha*alpha;

			Spectrum irradiance(0.0);
			int idx = 0;
			Float amplitude = 0.f;

			/* Find and store light source silhouettes
					Also sort the individual vertices in continus manner.
					Project to unit sphere and clip to horizon. */
			std::vector< std::vector<std::pair<Vector3f, Vector3f>> > lightsEdges, lightsEdgesGlobal;
			std::vector<Vector3f> lightsAvgNormal, lightsGlobalCg;
			std::vector<Spectrum> lightsEmit;
			std::vector<Float> lightsRadii;
			
			// LOG
			// std::vector< std::vector<std::pair<Vector3f, Vector3f>> > nonClipped, clipped;
			// LOG

			int numLights = scene.areaLightMeshes.size();
			for (int j=0; j < numLights; ++j) {
				const std::shared_ptr<TriangleMesh> &lightMesh = scene.areaLightMeshes[j];
				Vector3f avgNormal(0.f, 0.f, 0.f);
				auto silhouetteEdges = lightMesh->computeSilhouetteEdges(isect.p, avgNormal);
				auto gcg = lightMesh->cg - Vector3f(isect.p);

				auto cg = this->applyLocalTransform(isect, silhouetteEdges, avgNormal, gcg);
				auto globalEdges = silhouetteEdges;
				this->sortVertices(isect, silhouetteEdges, cg);
				this->projectToUnitSphere(silhouetteEdges);
				// std::vector<std::pair<Vector3f, Vector3f>> temp1(silhouetteEdges); //LOG
				this->clipToHorizon(isect, silhouetteEdges);
				// std::vector<std::pair<Vector3f, Vector3f>> temp2(silhouetteEdges); //LOG
				this->projectXY(isect, silhouetteEdges, projectionVector);

				if(silhouetteEdges.size() != 0) {
					lightsAvgNormal.push_back(avgNormal);
					lightsEdgesGlobal.push_back(globalEdges);
					lightsEdges.push_back(silhouetteEdges);
					lightsEmit.push_back(scene.areaLightRadiance[j]);
					lightsGlobalCg.push_back(gcg);

					lightsRadii.push_back(lightMesh->sphereRadius);

					// LOG
					// nonClipped.push_back(temp1);
					// clipped.push_back(temp2);
					// LOG
				}
			}

			/* Find and store blocker (occluder) silhouettes
					Also sort the individual vertices in continus manner.
					Project to unit sphere and clip to horizon. */
			std::vector< std::vector<std::pair<Vector3f, Vector3f>> > blockersEdges, blockersEdgesGlobal;
			std::vector<Vector3f> blockersAvgNormal, blockersGlobalCg;
			std::vector<Float> blockersRadii;

			// LOG
			// std::vector< std::vector<std::pair<Vector3f, Vector3f>> > occNonClipped, occClipped;
			// LOG

			int numBlockers = scene.blockerMeshes.size();
			for (int j=0; j < numBlockers; ++j) {
				const std::shared_ptr<TriangleMesh> &blockerMesh = scene.blockerMeshes[j];

				if(isectMesh->cg == blockerMesh->cg || !blockerMesh->isOccluder) {
					continue;
				}
				
				Vector3f avgNormal(0.f, 0.f, 0.f);
				auto silhouetteEdges = blockerMesh->computeSilhouetteEdges(isect.p, avgNormal);
				auto gcg = blockerMesh->cg - Vector3f(isect.p);

				auto cg = this->applyLocalTransform(isect, silhouetteEdges, avgNormal, gcg);
				auto globalEdges = silhouetteEdges;
				this->sortVertices(isect, silhouetteEdges, cg);
				this->projectToUnitSphere(silhouetteEdges);
				// std::vector<std::pair<Vector3f, Vector3f>> temp1(silhouetteEdges); // LOG
				this->clipToHorizon(isect, silhouetteEdges);
				// std::vector<std::pair<Vector3f, Vector3f>> temp2(silhouetteEdges); // LOG
				this->projectXY(isect, silhouetteEdges, projectionVector);

				if(silhouetteEdges.size() != 0) {
					blockersAvgNormal.push_back(avgNormal);
					blockersEdgesGlobal.push_back(globalEdges);
					blockersEdges.push_back(silhouetteEdges);
					blockersGlobalCg.push_back(gcg);

					blockersRadii.push_back(blockerMesh->sphereRadius);
					
					// LOG
					// occNonClipped.push_back(temp1);
					// occClipped.push_back(temp2);
					// LOG
				}
			}

			if(lightsEdges.size() == 0) {
				// LOG		
				// logFile << "0" << std::endl;
				// logFile << "0" << std::endl;
				// logFile << "0" << std::endl;
				// logFile << "0" << std::endl;
				// logFile << "0" << std::endl;
				// logFile << "0" << std::endl;
				// ENDLOG
			
				return L;
			}

			// LOG
			// std::vector< std::vector<std::pair<Vector3f, Vector3f>> > ghClip, ghClipLtc;
			// LOG

			int lidx = 0;
			for(auto light : lightsEdges) {
				Spectrum currentIrradiance(0.f);

				std::vector<std::pair<Vector3f, Vector3f>> currentLight(light);
				Vector3f currentLightCg(lightsGlobalCg[lidx]);
				
				if(blockersEdges.size() == 0) {
					/* Compute unoccluded radiance */
					std::vector<std::pair<Vector3f, Vector3f>> currentLight_(currentLight);

					this->unprojectXY(isect, currentLight_, projectionVector);
					this->projectToUnitSphere(currentLight_);
					// std::vector<std::pair<Vector3f, Vector3f>> temp1(currentLight_); // LOG
					this->applyLTC(isect, currentLight_, alpha, &amplitude);
					this->projectToUnitSphere(currentLight_);
					this->clipToHorizon(isect, currentLight_);
					this->projectToUnitSphere(currentLight_);
					// std::vector<std::pair<Vector3f, Vector3f>> temp2(currentLight_); // LOG

					currentIrradiance += lightsEmit[lidx] * amplitude * integrateEdges(currentLight_);

					// LOG
					// ghClip.push_back(temp1);
					// ghClipLtc.push_back(temp2);
					// LOG
				}
				else if(blockersEdges.size() != 0) {					
					/* Create current light source polygon */
					std::vector<PolyClip::Point2d> currentLightVertices;
					for(auto vertex : currentLight) {
						auto ver = PolyClip::Point2d(vertex.first.x, vertex.first.y);
						currentLightVertices.push_back(ver);
					}

					/* Extract verts for all blocker polygons */
					std::vector< std::vector<PolyClip::Point2d> > blockersPolygons;
					idx = 0;
					for(auto blocker : blockersEdges) {

						if(!this->isVisible(isect, lightsEdgesGlobal[lidx], lightsRadii[lidx], lightsAvgNormal[lidx], blockersEdgesGlobal[idx], blockersRadii[idx])) {
							idx++;
							continue;
						}
						
						std::vector<PolyClip::Point2d> blockerVertices;
						for(auto vertex : blocker) {
							auto ver = PolyClip::Point2d(vertex.first.x, vertex.first.y);							
							blockerVertices.push_back(ver);
						}
						
						blockersPolygons.push_back(blockerVertices);
						idx++;
					}

					if(blockersPolygons.size() == 0)
						zeroBlockers++;
					else if(blockersPolygons.size() == 1)
						oneBlockers++;
					else if(blockersPolygons.size() == 2)
						twoBlockers++;
					else if(blockersPolygons.size() == 3)
						threeBlockers++;

					if(blockersPolygons.size() != 0) {
						/* Union of occluder polygons */
						DLL *unionFirst = new DLL();
						DLL *unionLast = unionFirst;
						
						DLL *current = unionFirst;
						for(int i=0; i<blockersPolygons.size(); i++) {
							current->polygon = blockersPolygons[i];
							current->marker = true;
							
							DLL *newElem = new DLL();
							newElem->prev = current;
							
							current->next = newElem;
							current = newElem;
						}
						unionLast = current;

						// DLL *current1 = unionFirst;
						// while(true) {
						// 	if(current1->next == NULL)
						// 		break;
								
						// 	DLL *current2 = current1->next;
						// 	while(true) {
						// 		if(current1->polygon.size() != 0 && current2->polygon.size() != 0) {
						// 			PolyClip::Polygon c1(current1->polygon);
						// 			PolyClip::Polygon c2(current2->polygon);

						// 			/* GH Phase 1 */
						// 			PolyClip::PloygonOpration::DetectIntersection(c1, c2);

						// 			/* GH Phase 2, 3 */
						// 			int indicator = -1;
						// 			bool didIntersect = false;
						// 			std::vector<std::vector<PolyClip::Point2d>> possibleRes;
						// 			PolyClip::PloygonOpration::Mark(c1, c2, possibleRes, PolyClip::MarkUnion, &didIntersect, &indicator);

						// 			if(!didIntersect) {
						// 				if(indicator == 1) {
						// 					if(current1->prev != NULL)
						// 						current1->prev->next = current1->next;
												
						// 					if(current1->next != NULL)
						// 						current1->next->prev = current1->prev;

						// 					current1->polygon.clear();
						// 					break;
						// 				}
						// 				else if(indicator == 2) {
						// 					if(current2->prev != NULL)
						// 						current2->prev->next = current2->next;
												
						// 					if(current2->next != NULL)
						// 						current2->next->prev = current2->prev;

						// 					current2->polygon.clear();
						// 				}
						// 			}
						// 			else if(didIntersect) {
						// 				std::vector<std::vector<PolyClip::Point2d>> diffRes = PolyClip::PloygonOpration::ExtractUnionResults(c1);
										
						// 				// Remove last point ; because it is the same as the first point
						// 				// TODO: Does union having > 1 elements mean that the rest need to be subtracted?
						// 				// diffRes[0].pop_back();
						// 				current1->polygon = diffRes[0];

						// 				if(current2->prev != NULL)
						// 					current2->prev->next = current2->next;
											
						// 				if(current2->next != NULL)
						// 					current2->next->prev = current2->prev;

						// 				current2->polygon.clear();
						// 			}
						// 		}

						// 		if(current2->next == NULL)
						// 			break;
						// 		else
						// 			current2 = current2->next;
						// 	}

						// 	current1 = current1->next;
						// }

						/* Loop over all blocker polygons and take set difference with current light source polygon */
						int numElems = 1;
						DLL *firstElem = new DLL();
						DLL *lastElem = firstElem;
						firstElem->polygon = currentLightVertices;
						firstElem->marker = false;

						DLL *currentUnion = unionFirst;
						while(true) {
							if(currentUnion->polygon.size() != 0) {
								DLL *current = firstElem;
								while(true) {

									if(current->polygon.size() != 0 && !current->marker) {
										PolyClip::Polygon lightPolygon(current->polygon);
										PolyClip::Polygon blockerPolygon(currentUnion->polygon);

										/* GH Phase 1 */
										PolyClip::PloygonOpration::DetectIntersection(lightPolygon, blockerPolygon);

										/* GH Phase 2, 3 */
										int indicator = -1;
										bool didIntersect = false;
										std::vector<std::vector<PolyClip::Point2d>> possibleRes;
										PolyClip::PloygonOpration::Mark(lightPolygon, blockerPolygon, possibleRes, PolyClip::MarkDifferentiate, &didIntersect, &indicator);

										if(!didIntersect) {
											if(indicator == 1) {
												if(current->prev != NULL)
													current->prev->next = current->next;
												
												if(current->next != NULL)
													current->next->prev = current->prev;
												
												current->polygon.clear();

												numElems -= 1;
											}
											else if(indicator == 2) {
												DLL *newElem = new DLL();
												newElem->next = firstElem;
												newElem->polygon = currentUnion->polygon;
												newElem->marker = true;

												firstElem->prev = newElem;
												firstElem = newElem;

												numElems += 1;
											}
										}
										else if(didIntersect) {
											// std::vector<std::vector<PolyClip::Point2d>> diffRes;
											std::vector<std::vector<PolyClip::Point2d>> diffRes = PolyClip::PloygonOpration::ExtractDifferentiateResults(lightPolygon);

											if(current->prev != NULL)
												current->prev->next = current->next;
												
											if(current->next != NULL)
												current->next->prev = current->prev;
											
											current->polygon.clear();
											numElems -= 1;
											
											int pidx = 0;
											for(auto polyVerts : diffRes) {
												// polyVerts.pop_back();
												
												DLL *newElem = new DLL();
												newElem->next = firstElem;
												newElem->polygon = polyVerts;
												newElem->marker = false;

												firstElem->prev = newElem;
												firstElem = newElem;

												numElems += 1;
												pidx++;
											}
										}
									}

									if(current->next == NULL)
										break;
									else
										current = current->next;
								}
							}

							if(currentUnion->next == NULL)
								break;
							else
								currentUnion = currentUnion->next;
						}

						current = firstElem;

						// LOG
						// std::vector<std::pair<Vector3f, Vector3f>> temp1, temp2;
						// LOG

						while(true) {
							std::vector<std::pair<Vector3f, Vector3f>> intersectionEdges;

							auto polygon = current->polygon;

							if(polygon.size() != 0) {
								for (int j=0; j<polygon.size(); j++) {
									Vector3f v1(polygon[j].x_, polygon[j].y_, 1.0);
									Vector3f v2(polygon[(j+1)%polygon.size()].x_, polygon[(j+1)%polygon.size()].y_, 1.0);

									intersectionEdges.push_back(std::pair<Vector3f, Vector3f>(v1, v2));
								}

								this->unprojectXY(isect, intersectionEdges, projectionVector);
								this->projectToUnitSphere(intersectionEdges);
								// temp1 = intersectionEdges; // LOG
								this->applyLTC(isect, intersectionEdges, alpha, &amplitude);
								this->projectToUnitSphere(intersectionEdges);
								this->clipToHorizon(isect, intersectionEdges);
								this->projectToUnitSphere(intersectionEdges);
								// temp2 = intersectionEdges; // LOG

								if(current->marker)
									currentIrradiance = currentIrradiance - ( lightsEmit[lidx] * amplitude * integrateEdges(intersectionEdges) );
								else
									currentIrradiance = currentIrradiance + ( lightsEmit[lidx] * amplitude * integrateEdges(intersectionEdges) );
							}
							
							if(current->next == NULL)
								break;
							else
								current = current->next;
						}

						// LOG
						// ghClip.push_back(temp1);
						// ghClipLtc.push_back(temp2);
						// LOG
					}
					else {
						std::vector<std::pair<Vector3f, Vector3f>> currentLight_(currentLight);

						this->unprojectXY(isect, currentLight_, projectionVector);
						this->projectToUnitSphere(currentLight_);
						// std::vector<std::pair<Vector3f, Vector3f>> temp1(currentLight_); // LOG
						this->applyLTC(isect, currentLight_, alpha, &amplitude);
						this->projectToUnitSphere(currentLight_);
						this->clipToHorizon(isect, currentLight_);
						this->projectToUnitSphere(currentLight_);
						// std::vector<std::pair<Vector3f, Vector3f>> temp2(currentLight_); // LOG

						currentIrradiance += lightsEmit[lidx] * amplitude * integrateEdges(currentLight_);

						// LOG
						// ghClip.push_back(temp1);
						// ghClipLtc.push_back(temp2);
						// LOG
					}
				}

				/* Final radiance computation */
				irradiance = 0.5 * InvPi * currentIrradiance;
				if(irradiance.y() >= 0.f)
					L += diffuse * irradiance;

				lidx++;
			}

			// LOG
			// logFile << std::to_string(nonClipped.size()) << std::endl;
			// for(auto item : nonClipped) {
			// 	logFile << std::to_string(item.size()) << std::endl;

			// 	for(auto edge : item) {
			// 		logFile << SphericalTheta(edge.first) << " " << SphericalPhi(edge.first) << " ";
			// 	}
			// 	logFile << std::endl;
			// }

			// logFile << std::to_string(clipped.size()) << std::endl;
			// for(auto item : clipped) {
			// 	logFile << std::to_string(item.size()) << std::endl;

			// 	for(auto edge : item) {
			// 		logFile << SphericalTheta(edge.first) << " " << SphericalPhi(edge.first) << " ";
			// 	}
			// 	logFile << std::endl;
			// }

			// logFile << std::to_string(occNonClipped.size()) << std::endl;
			// for(auto item : occNonClipped) {
			// 	logFile << std::to_string(item.size()) << std::endl;

			// 	for(auto edge : item) {
			// 		logFile << SphericalTheta(edge.first) << " " << SphericalPhi(edge.first) << " ";
			// 	}
			// 	logFile << std::endl;
			// }

			// logFile << std::to_string(occClipped.size()) << std::endl;
			// for(auto item : occClipped) {
			// 	logFile << std::to_string(item.size()) << std::endl;

			// 	for(auto edge : item) {
			// 		logFile << SphericalTheta(edge.first) << " " << SphericalPhi(edge.first) << " ";
			// 	}
			// 	logFile << std::endl;
			// }

			// logFile << std::to_string(ghClip.size()) << std::endl;
			// for(auto item : ghClip) {
			// 	logFile << std::to_string(item.size()) << std::endl;

			// 	for(auto edge : item) {
			// 		logFile << SphericalTheta(edge.first) << " " << SphericalPhi(edge.first) << " ";
			// 	}
			// 	logFile << std::endl;
			// }

			// logFile << std::to_string(ghClipLtc.size()) << std::endl;
			// for(auto item : ghClipLtc) {
			// 	logFile << std::to_string(item.size()) << std::endl;

			// 	for(auto edge : item) {
			// 		logFile << SphericalTheta(edge.first) << " " << SphericalPhi(edge.first) << " ";
			// 	}
			// 	logFile << std::endl;
			// }
			// LOG
		}
	}

	return L;
}

Vector3f LTCSilhouetteShadow::gnomonicUnproject(std::vector<std::pair<Vector3f, Vector3f>> &edges, Vector3f cg_, Vector3f refp) 
		const {
		
	Float theta = cg_.y * Pi / 180.0;
	Float phi = cg_.x * Pi / 180.0;
	if(theta == 0.f)
		phi == 0.f;
	Vector3f cg(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), cos(theta));
	cg = cg / cg.Length();

	for(int i=0; i<edges.size(); i++) {
		// Edge point 1
		theta = edges[i].first.y * Pi / 180.0;
		phi = edges[i].first.x * Pi / 180.0;
		if(theta == 0.f)
			phi == 0.f;
		edges[i].first = Vector3f(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta));
		edges[i].first = edges[i].first / edges[i].first.Length();

		// Edge point 2
		theta = edges[i].second.y * Pi / 180.0;
		phi = edges[i].second.x * Pi / 180.0;
		if(theta == 0.f)
			phi == 0.f;
		edges[i].second = Vector3f(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta));
		edges[i].second = edges[i].second / edges[i].second.Length();
	}

	return cg;
	
	// Float thetas_ = atan2( sqrt( pow(refp.x, 2) + pow(refp.y, 2) ), refp.z );
	// Float phis_ = atan2(refp.y, refp.x);

	// Float phis = Pi/2.0 - thetas_;
	// Float lambds = phis_;

	// // Cg
	// Float x = cg_.x;
	// Float y = cg_.y;
	// Float rho = sqrt( pow(x, 2) + pow(y, 2) );
	// Float c = atan(rho);

	// Float num = y*sin(c)*cos(phis);
	// Float phip = 0.f;
	// if(num == 0)
	// 	phip = asin( cos(c)*sin(phis) );
	// else
	// 	phip = asin( cos(c)*sin(phis) + num/rho );
	
	// Float lambdp = lambds + atan2( x*sin(c), rho*cos(phis)*cos(c) - y*sin(phis)*sin(c) );

	// Float thetap_ = Pi/2.0 - phip;
	// Float phip_ = lambdp;

	// Vector3f cg(sin(thetap_)*cos(phip_), sin(thetap_)*sin(phip_), cos(thetap_));
	
	// for(int i=0; i<edges.size(); i++) {
	// 	// Edge point 1
	// 	x = edges[i].first.x;
	// 	y = edges[i].first.y;
	// 	rho = sqrt( pow(x, 2) + pow(y, 2) );
	// 	c = atan(rho);

	// 	num = y*sin(c)*cos(phis);
	// 	phip = 0.f;
	// 	if(num == 0)
	// 		phip = asin( cos(c)*sin(phis) );
	// 	else
	// 		phip = asin( cos(c)*sin(phis) + num/rho );
		
	// 	lambdp = lambds + atan2( x*sin(c), rho*cos(phis)*cos(c) - y*sin(phis)*sin(c) );

	// 	thetap_ = Pi/2.0 - phip;
	// 	phip_ = lambdp;

	// 	edges[i].first = Vector3f(sin(thetap_)*cos(phip_), sin(thetap_)*sin(phip_), cos(thetap_));

	// 	// Edge point 2
	// 	x = edges[i].second.x;
	// 	y = edges[i].second.y;
	// 	rho = sqrt( pow(x, 2) + pow(y, 2) );
	// 	c = atan(rho);

	// 	num = y*sin(c)*cos(phis);
	// 	phip = 0.f;
	// 	if(num == 0)
	// 		phip = asin( cos(c)*sin(phis) );
	// 	else
	// 		phip = asin( cos(c)*sin(phis) + num/rho );
		
	// 	lambdp = lambds + atan2( x*sin(c), rho*cos(phis)*cos(c) - y*sin(phis)*sin(c) );

	// 	thetap_ = Pi/2.0 - phip;
	// 	phip_ = lambdp;

	// 	edges[i].second = Vector3f(sin(thetap_)*cos(phip_), sin(thetap_)*sin(phip_), cos(thetap_));
	// }

	// return cg;
}

Vector3f LTCSilhouetteShadow::gnomonicProject(std::vector<std::pair<Vector3f, Vector3f>> &edges, Vector3f cg_, Vector3f refp) 
		const {
	
	Float cgx = SphericalPhi(cg_) * 180.0 / Pi;
	Float cgy = SphericalTheta(cg_) * 180.0 / Pi;
	if(cgy == 0.f)
		cgx = 0.f;
	Vector3f cg(cgx, cgy, 0.f);

	for(int i=0; i<edges.size(); i++) {
		// Edge point 1
		Float fx = SphericalPhi(edges[i].first) * 180.0 / Pi;
		Float fy = SphericalTheta(edges[i].first) * 180.0 / Pi;
		if(fy == 0.f)
			fx = 0.f;
		edges[i].first.x = fx;
		edges[i].first.y = fy;
		edges[i].first.z = 0.f;

		// Edge point 2
		fx = SphericalPhi(edges[i].second) * 180.0 / Pi;
		fy = SphericalTheta(edges[i].second) * 180.0 / Pi;
		if(fy == 0.f)
			fx = 0.f;
		edges[i].second.x = fx;
		edges[i].second.y = fy;
		edges[i].second.z = 0.f;
	}

	// Order clockwise w.r.t. cg
	// for(int i=0; i<edges.size(); i++) {
	// 	Vector3f v1 = edges[i].first;
	// 	Vector3f v2 = edges[i].second;
	// 	Vector3f cp = Cross(v1-cg, v2-cg);
	// 	// cp = cp / cp.Length();

	// 	if(Dot(cp, Vector3f(0.f, 0.f, 1.f)) < 0) {
	// 		edges[i].first = v2;
	// 		edges[i].second = v1;
	// 	}
	// }

	return cg;

	// Float thetap_ = atan2( sqrt( pow(refp.x, 2) + pow(refp.y, 2) ), refp.z );
	// Float phip_ = atan2(refp.y, refp.x);

	// Float phip = Pi/2.0 - thetap_;
	// Float lambdp = phip_;

	// // Cg
	// Vector3f cg(cg_);
	// Float thetas_ = atan2( sqrt( pow(cg.x, 2) + pow(cg.y, 2) ), cg.z );
	// Float phis_ = atan2(cg.y, cg.x);

	// Float phis = Pi/2.0 - thetas_;
	// Float lambds = phis_;

	// Float cosc = sin(phip)*sin(phis) + cos(phip)*cos(phis)*cos(lambds-lambdp);

	// cg.x = cos(phis)*sin(lambds-lambdp) / cosc;
	// cg.y = ( cos(phip)*sin(phis) - sin(phip)*cos(phis)*cos(lambds-lambdp) ) / cosc;
	// cg.z = 0.f;
	
	// for(int i=0; i<edges.size(); i++) {
	// 	// Edge point 1
	// 	Vector3f v = edges[i].first;
	// 	thetas_ = atan2( sqrt( pow(v.x, 2) + pow(v.y, 2) ), v.z );
	// 	phis_ = atan2(v.y, v.x);

	// 	phis = Pi/2.0 - thetas_;
	// 	lambds = phis_;

	// 	cosc = sin(phip)*sin(phis) + cos(phip)*cos(phis)*cos(lambds-lambdp);

	// 	edges[i].first.x = cos(phis)*sin(lambds-lambdp) / cosc;
	// 	edges[i].first.y = ( cos(phip)*sin(phis) - sin(phip)*cos(phis)*cos(lambds-lambdp) ) / cosc;
	// 	edges[i].first.z = 0.f;

	// 	// Edge point 2
	// 	v = edges[i].second;

	// 	thetas_ = atan2( sqrt( pow(v.x, 2) + pow(v.y, 2) ), v.z );
	// 	phis_ = atan2(v.y, v.x);

	// 	phis = Pi/2.0 - thetas_;
	// 	lambds = phis_;

	// 	cosc = sin(phip)*sin(phis) + cos(phip)*cos(phis)*cos(lambds-lambdp);

	// 	edges[i].second.x = cos(phis)*sin(lambds-lambdp) / cosc;
	// 	edges[i].second.y = ( cos(phip)*sin(phis) - sin(phip)*cos(phis)*cos(lambds-lambdp) ) / cosc;
	// 	edges[i].second.z = 0.f;
	// }

	// // // Order clockwise w.r.t. cg
	// for(int i=0; i<edges.size(); i++) {
	// 	Vector3f v1 = edges[i].first;
	// 	Vector3f v2 = edges[i].second;
	// 	Vector3f cp = Cross(v1-cg, v2-cg);
	// 	// cp = cp / cp.Length();

	// 	if(Dot(cp, Vector3f(0.f, 0.f, 1.f)) < 0) {
	// 		edges[i].first = v2;
	// 		edges[i].second = v1;
	// 	}
	// }

	// return cg;
}



LTCSilhouetteShadow *CreateLTCSilhouetteShadow(const ParamSet &params,
                                     std::shared_ptr<Sampler> sampler,
                                     std::shared_ptr<const Camera> camera) {

	Bounds2i pixelBounds = camera->film->GetSampleBounds();

    return new LTCSilhouetteShadow(pixelBounds, camera, sampler);
}

}