
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

// shapes/triangle.cpp*
#include "shapes/triangle.h"
#include "texture.h"
#include "textures/constant.h"
#include "paramset.h"
#include "sampling.h"
#include "efloat.h"
#include "ext/rply.h"
#include <array>

namespace pbrt {

STAT_PERCENT("Intersections/Ray-triangle intersection tests", nHits, nTests);

// Triangle Local Definitions
static void PlyErrorCallback(p_ply, const char *message) {
    Error("PLY writing error: %s", message);
}

// Triangle Method Definitions
STAT_RATIO("Scene/Triangles per triangle mesh", nTris, nMeshes);
TriangleMesh::TriangleMesh(
    const Transform &ObjectToWorld, int nTriangles, const int *vertexIndices,
    int nVertices, const Point3f *P, const Vector3f *S, const Normal3f *N,
    const Point2f *UV, const std::shared_ptr<Texture<Float>> &alphaMask,
    const std::shared_ptr<Texture<Float>> &shadowAlphaMask,
    const int *fIndices)
    : nTriangles(nTriangles),
      nVertices(nVertices),
      vertexIndices(vertexIndices, vertexIndices + 3 * nTriangles),
      alphaMask(alphaMask),
      shadowAlphaMask(shadowAlphaMask) {
    ++nMeshes;
    nTris += nTriangles;
    triMeshBytes += sizeof(*this) + this->vertexIndices.size() * sizeof(int) +
                    nVertices * (sizeof(*P) + (N ? sizeof(*N) : 0) +
                                 (S ? sizeof(*S) : 0) + (UV ? sizeof(*UV) : 0) +
                                 (fIndices ? sizeof(*fIndices) : 0));

    // Transform mesh vertices to world space
    p.reset(new Point3f[nVertices]);
    for (int i = 0; i < nVertices; ++i) p[i] = ObjectToWorld(P[i]);

    // Copy _UV_, _N_, and _S_ vertex data, if present
    if (UV) {
        uv.reset(new Point2f[nVertices]);
        memcpy(uv.get(), UV, nVertices * sizeof(Point2f));
    }
    if (N) {
        n.reset(new Normal3f[nVertices]);
        for (int i = 0; i < nVertices; ++i) n[i] = ObjectToWorld(N[i]);
    }
    if (S) {
        s.reset(new Vector3f[nVertices]);
        for (int i = 0; i < nVertices; ++i) s[i] = ObjectToWorld(S[i]);
    }

    if (fIndices)
        faceIndices = std::vector<int>(fIndices, fIndices + nTriangles);
    
    this->isLight = false;
    
    for(int i=0; i<this->vertexIndices.size(); i+=3) {
        int i1 = this->vertexIndices[i];
        int i2 = this->vertexIndices[i+1];
        int i3 = this->vertexIndices[i+2];

        Vector3f v1 = Vector3f( this->p[ this->vertexIndices[i] ] );
        Vector3f v2 = Vector3f( this->p[ this->vertexIndices[i+1] ] );
        Vector3f v3 = Vector3f( this->p[ this->vertexIndices[i+2] ] );
        Vector3f cp = Cross(v1-v2, v1-v3);
        cp = cp / cp.Length();

        Normal3f n1 = n[i1];
        Normal3f n2 = n[i2];
        Normal3f n3 = n[i3];

        Normal3f finalNormal = (n1+n2+n3)/3.0f;
        finalNormal /= finalNormal.Length();

        // if(Dot(finalNormal, cp) >= 0)
        //     finalNormal = Normal3f(cp);
        // else
        //     finalNormal = Normal3f(-cp);

        faceNormals.push_back(finalNormal);
    }

    // std::cout << std::to_string(faceNormals.size()) << std::endl;
    
    std::map<std::string, int> vds;
    for(int i=0; i<this->vertexIndices.size(); i+=3) {
        Vector3f v1 = Vector3f( this->p[ this->vertexIndices[i] ] );
        Vector3f v2 = Vector3f( this->p[ this->vertexIndices[i+1] ] );
        Vector3f v3 = Vector3f( this->p[ this->vertexIndices[i+2] ] );

        int v1idx, v2idx, v3idx;

        std::string v1k = std::to_string(v1.x) + "_" + std::to_string(v1.y) + "_" + std::to_string(v1.z);
        std::string v2k = std::to_string(v2.x) + "_" + std::to_string(v2.y) + "_" + std::to_string(v2.z);
        std::string v3k = std::to_string(v3.x) + "_" + std::to_string(v3.y) + "_" + std::to_string(v3.z);

        auto v1r = vds.find(v1k);
        auto v2r = vds.find(v2k);
        auto v3r = vds.find(v3k);

        if(v1r == vds.end()) {
            this->uniqueVertices.push_back(v1);

            vds[v1k] = this->uniqueVertices.size()-1;
            v1idx = vds[v1k];
        }
        else 
            v1idx = vds[v1k];
        
        if(v2r == vds.end()) {
            this->uniqueVertices.push_back(v2);

            vds[v2k] = this->uniqueVertices.size()-1;
            v2idx = vds[v2k];
        }
        else 
            v2idx = vds[v2k];

        if(v3r == vds.end()) {
            this->uniqueVertices.push_back(v3);

            vds[v3k] = this->uniqueVertices.size()-1;
            v3idx = vds[v3k];
        }
        else 
            v3idx = vds[v3k];
        
        std::vector<int> temp;
        temp.push_back(v1idx);
        temp.push_back(v2idx);
        temp.push_back(v3idx);

        uniqueFaceIndices.push_back(temp);
    }

    // std::cout << std::to_string(uniqueFaceIndices.size()) << std::endl;

    this->cg = Vector3f(0.f, 0.f, 0.f);
    for (int i = 0; i < this->uniqueVertices.size(); ++i) {
        this->cg += Vector3f(this->uniqueVertices[i]);
    }
    this->cg /= Float(this->uniqueVertices.size());

    this->sphereRadius = 0.f;
    for(int i = 0; i < this->uniqueVertices.size(); ++i) {
        Vector3f temp = this->uniqueVertices[i] - this->cg;
        this->sphereRadius = std::max(this->sphereRadius, temp.Length());
    }
    // this->sphereRadius += 1.0;

    // std::cout << this->cg << std::endl;
    // std::cout << "-------" << std::endl;

    std::map<std::string, int> ds;
    int face_index = 0;
    for(int i=0; i<this->uniqueFaceIndices.size(); i++) {
        int i1 = this->uniqueFaceIndices[i][0];
        int i2 = this->uniqueFaceIndices[i][1];
        int i3 = this->uniqueFaceIndices[i][2];

        int eidx1, eidx2, eidx3;

        std::string k1 = std::to_string(i1) + "_" + std::to_string(i2);
        std::string k1_ = std::to_string(i2) + "_" + std::to_string(i1);
        auto k1r = ds.find(k1);
        auto k1_r = ds.find(k1_);
        if(k1r != ds.end() || k1_r != ds.end()) {
            int dsidx = -1;
            if(k1r != ds.end())
                dsidx = ds[k1];
            else if(k1_r != ds.end())
                dsidx = ds[k1_];

            this->edgeIndices[dsidx][3] = face_index;
            eidx1 = dsidx;
        }
        else if(k1r == ds.end() && k1_r == ds.end()) {
            std::vector<int> temp;
            temp.push_back(i1);
            temp.push_back(i2);
            temp.push_back(face_index);
            temp.push_back(-1);

            this->edgeIndices.push_back(temp);

            ds[k1] = this->edgeIndices.size()-1;
            eidx1 = ds[k1];
        }

        std::string k2 = std::to_string(i2) + "_" + std::to_string(i3);
        std::string k2_ = std::to_string(i3) + "_" + std::to_string(i2);
        auto k2r = ds.find(k2);
        auto k2_r = ds.find(k2_);
        if(k2r != ds.end() || k2_r != ds.end()) {
            int dsidx = -1;
            if(k2r != ds.end())
                dsidx = ds[k2];
            else if(k2_r != ds.end())
                dsidx = ds[k2_];

            this->edgeIndices[dsidx][3] = face_index;
            eidx2 = dsidx;
        }
        else if(k2r == ds.end() && k2_r == ds.end()) {
            std::vector<int> temp;
            temp.push_back(i2);
            temp.push_back(i3);
            temp.push_back(face_index);
            temp.push_back(-1);

            this->edgeIndices.push_back(temp);

            ds[k2] = this->edgeIndices.size()-1;
            eidx2 = ds[k2];
        }

        std::string k3 = std::to_string(i3) + "_" + std::to_string(i1);
        std::string k3_ = std::to_string(i1) + "_" + std::to_string(i3);
        auto k3r = ds.find(k3);
        auto k3_r = ds.find(k3_);
        if(k3r != ds.end() || k3_r != ds.end()) {
            int dsidx = -1;
            if(k3r != ds.end())
                dsidx = ds[k3];
            else if(k3_r != ds.end())
                dsidx = ds[k3_];

            this->edgeIndices[dsidx][3] = face_index;
            eidx3 = dsidx;
        }
        else if(k3r == ds.end() && k3_r == ds.end()) {
            std::vector<int> temp;
            temp.push_back(i3);
            temp.push_back(i1);
            temp.push_back(face_index);
            temp.push_back(-1);

            this->edgeIndices.push_back(temp);

            ds[k3] = this->edgeIndices.size()-1;
            eidx3 = ds[k3];
        }

        face_index++;
    }

    // for(auto temp : this->edgeIndices) {
    //     for(auto temp1 : temp) {
    //         std::cout << std::to_string(temp1) << " ";
    //     }
    //     std::cout << std::endl;
    // }
}

std::vector< std::pair<Vector3f, Vector3f> > TriangleMesh::getTriangle(Point3f sp, Vector3f &avgNormal, Vector3f &gcg, int idx) {
    Vector3f shadingPoint = Vector3f(sp);
    std::vector< std::pair<Vector3f, Vector3f> > edges;
    avgNormal = Normalize( Vector3f(this->faceNormals[idx]) );

    Vector3f v1 = Vector3f(this->uniqueVertices[ this->uniqueFaceIndices[idx][0] ]);
    Vector3f v2 = Vector3f(this->uniqueVertices[ this->uniqueFaceIndices[idx][1] ]);
    Vector3f v3 = Vector3f(this->uniqueVertices[ this->uniqueFaceIndices[idx][2] ]);

    gcg = (v1+v2+v3) / 3.0;
    gcg = gcg - shadingPoint;

    Float cnd = Dot(-Normalize(gcg), avgNormal); 
    if(cnd >= 0) {
        edges.push_back(std::pair<Vector3f, Vector3f>(v1-shadingPoint, v2-shadingPoint));
        edges.push_back(std::pair<Vector3f, Vector3f>(v2-shadingPoint, v3-shadingPoint));
        edges.push_back(std::pair<Vector3f, Vector3f>(v3-shadingPoint, v1-shadingPoint));

        // Vector3f rand1(rand()/Float(RAND_MAX), rand()/Float(RAND_MAX), rand()/Float(RAND_MAX));
        // rand1 *= 1e-2;
        // Vector3f rand2(rand()/Float(RAND_MAX), rand()/Float(RAND_MAX), rand()/Float(RAND_MAX));
        // rand2 *= 1e-2;
        // Vector3f rand3(rand()/Float(RAND_MAX), rand()/Float(RAND_MAX), rand()/Float(RAND_MAX));
        // rand3 *= 1e-2;

        // edges.push_back(std::pair<Vector3f, Vector3f>(v1-shadingPoint+rand1, v2-shadingPoint+rand1));
        // edges.push_back(std::pair<Vector3f, Vector3f>(v2-shadingPoint+rand2, v3-shadingPoint+rand2));
        // edges.push_back(std::pair<Vector3f, Vector3f>(v3-shadingPoint+rand3, v1-shadingPoint+rand3));
    }

    return edges;
}

std::vector< std::pair<Vector3f, Vector3f> > TriangleMesh::computeSilhouetteEdges(Point3f sp, Vector3f &avgNormal) {
    std::vector< std::pair<Vector3f, Vector3f> > edges;
    avgNormal = Vector3f(0.f, 0.f, 0.f);

    int eidx = 0;
    Vector3f shadingPoint = Vector3f(sp);
    for(auto edge : this->edgeIndices) {
        std::pair<Vector3f, Vector3f> sil;

        Vector3f v1 = Vector3f(this->uniqueVertices[ edge[0] ]);
        Vector3f v2 = Vector3f(this->uniqueVertices[ edge[1] ]);
        Vector3f cg = Vector3f((v1+v2)/2.0);

        if(edge[2] == -1) {
            Normal3f fn = this->faceNormals[ edge[3] ];
            Vector3f view = Vector3f(shadingPoint - cg);
            view = view / view.Length();

            Float cnd = Dot(view, fn); 
            if(cnd >= 0) {
                sil.first = v1 - shadingPoint;
                sil.second = v2 - shadingPoint;

                avgNormal += Vector3f(fn);
                edges.push_back(sil);
            }
        }
        else if(edge[3] == -1) {
            Normal3f fn = this->faceNormals[ edge[2] ];
            Vector3f view = Vector3f(shadingPoint - cg);
            view = view / view.Length();

            Float cnd = Dot(view, fn); 
            if(cnd >= 0) {
                sil.first = v1 - shadingPoint;
                sil.second = v2 - shadingPoint;

                avgNormal += Vector3f(fn);
                edges.push_back(sil);
            }
        }
        else {
            Normal3f fn1 = this->faceNormals[ edge[2] ];
            Normal3f fn2 = this->faceNormals[ edge[3] ];

            Vector3f view = Vector3f(shadingPoint - cg);
            view = view / view.Length();

            Float cnd1 = Dot(view, fn1); 
            Float cnd2 = Dot(view, fn2); 

            if((cnd1 >= 0 && cnd2 < 0) || (cnd1 < 0 && cnd2 >= 0)) {
                sil.first = v1 - shadingPoint;
                sil.second = v2 - shadingPoint;

                if((cnd1 >= 0 && cnd2 < 0))
                    avgNormal += Vector3f(fn1);
                else if((cnd1 < 0 && cnd2 >= 0))
                    avgNormal += Vector3f(fn2);
                
                edges.push_back(sil);
            }
        }

        eidx++;
    }

    if(edges.size() != 0) {
        avgNormal /= avgNormal.Length();
    }

    return edges;
}

std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(
    const Transform *ObjectToWorld, const Transform *WorldToObject,
    bool reverseOrientation, int nTriangles, const int *vertexIndices,
    int nVertices, const Point3f *p, const Vector3f *s, const Normal3f *n,
    const Point2f *uv, const std::shared_ptr<Texture<Float>> &alphaMask,
    const std::shared_ptr<Texture<Float>> &shadowAlphaMask,
    const int *faceIndices) {
    std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>(
        *ObjectToWorld, nTriangles, vertexIndices, nVertices, p, s, n, uv,
        alphaMask, shadowAlphaMask, faceIndices);
    std::vector<std::shared_ptr<Shape>> tris;
    tris.reserve(nTriangles);
    for (int i = 0; i < nTriangles; ++i)
        tris.push_back(std::make_shared<Triangle>(ObjectToWorld, WorldToObject,
                                                  reverseOrientation, mesh, i));
    return tris;
}

std::vector<std::shared_ptr<Shape>> CreateTriangleMeshWithParams( 
    const Transform *ObjectToWorld, const Transform *WorldToObject,
    bool reverseOrientation, int nTriangles, const int *vertexIndices,
    int nVertices, const Point3f *p, const Vector3f *s, const Normal3f *n,
    const Point2f *uv, const std::shared_ptr<Texture<Float>> &alphaMask,
    const std::shared_ptr<Texture<Float>> &shadowAlphaMask,
    const int *faceIndices, const ParamSet &params) {
    std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>(
        *ObjectToWorld, nTriangles, vertexIndices, nVertices, p, s, n, uv,
        alphaMask, shadowAlphaMask, faceIndices);

    mesh->isOccluder = params.FindOneBool("occluder", false);
    mesh->isLight = params.FindOneBool("lightsource", false);

    std::vector<std::shared_ptr<Shape>> tris;
    tris.reserve(nTriangles);
    for (int i = 0; i < nTriangles; ++i)
        tris.push_back(std::make_shared<Triangle>(ObjectToWorld, WorldToObject,
                                                  reverseOrientation, mesh, i));
    return tris;
}

bool WritePlyFile(const std::string &filename, int nTriangles,
                  const int *vertexIndices, int nVertices, const Point3f *P,
                  const Vector3f *S, const Normal3f *N, const Point2f *UV,
                  const int *faceIndices) {
    p_ply plyFile =
        ply_create(filename.c_str(), PLY_DEFAULT, PlyErrorCallback, 0, nullptr);
    if (plyFile == nullptr)
        return false;

    ply_add_element(plyFile, "vertex", nVertices);
    ply_add_scalar_property(plyFile, "x", PLY_FLOAT);
    ply_add_scalar_property(plyFile, "y", PLY_FLOAT);
    ply_add_scalar_property(plyFile, "z", PLY_FLOAT);
    if (N) {
        ply_add_scalar_property(plyFile, "nx", PLY_FLOAT);
        ply_add_scalar_property(plyFile, "ny", PLY_FLOAT);
        ply_add_scalar_property(plyFile, "nz", PLY_FLOAT);
    }
    if (UV) {
        ply_add_scalar_property(plyFile, "u", PLY_FLOAT);
        ply_add_scalar_property(plyFile, "v", PLY_FLOAT);
    }
    if (S)
        Warning("%s: PLY mesh will be missing tangent vectors \"S\".",
                filename.c_str());

    ply_add_element(plyFile, "face", nTriangles);
    ply_add_list_property(plyFile, "vertex_indices", PLY_UINT8, PLY_INT);
    if (faceIndices)
        ply_add_scalar_property(plyFile, "face_indices", PLY_INT);
    ply_write_header(plyFile);

    for (int i = 0; i < nVertices; ++i) {
        ply_write(plyFile, P[i].x);
        ply_write(plyFile, P[i].y);
        ply_write(plyFile, P[i].z);
        if (N) {
            ply_write(plyFile, N[i].x);
            ply_write(plyFile, N[i].y);
            ply_write(plyFile, N[i].z);
        }
        if (UV) {
            ply_write(plyFile, UV[i].x);
            ply_write(plyFile, UV[i].y);
        }
    }

    for (int i = 0; i < nTriangles; ++i) {
        ply_write(plyFile, 3);
        ply_write(plyFile, vertexIndices[3 * i]);
        ply_write(plyFile, vertexIndices[3 * i + 1]);
        ply_write(plyFile, vertexIndices[3 * i + 2]);
        if (faceIndices)
            ply_write(plyFile, faceIndices[i]);
    }
    ply_close(plyFile);
    return true;
}

Bounds3f Triangle::ObjectBound() const {
    // Get triangle vertices in _p0_, _p1_, and _p2_
    const Point3f &p0 = mesh->p[v[0]];
    const Point3f &p1 = mesh->p[v[1]];
    const Point3f &p2 = mesh->p[v[2]];
    return Union(Bounds3f((*WorldToObject)(p0), (*WorldToObject)(p1)),
                 (*WorldToObject)(p2));
}

Bounds3f Triangle::WorldBound() const {
    // Get triangle vertices in _p0_, _p1_, and _p2_
    const Point3f &p0 = mesh->p[v[0]];
    const Point3f &p1 = mesh->p[v[1]];
    const Point3f &p2 = mesh->p[v[2]];
    return Union(Bounds3f(p0, p1), p2);
}

bool Triangle::Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                         bool testAlphaTexture) const {
    ProfilePhase p(Prof::TriIntersect);
    ++nTests;
    // Get triangle vertices in _p0_, _p1_, and _p2_
    const Point3f &p0 = mesh->p[v[0]];
    const Point3f &p1 = mesh->p[v[1]];
    const Point3f &p2 = mesh->p[v[2]];

    // Perform ray--triangle intersection test

    // Transform triangle vertices to ray coordinate space

    // Translate vertices based on ray origin
    Point3f p0t = p0 - Vector3f(ray.o);
    Point3f p1t = p1 - Vector3f(ray.o);
    Point3f p2t = p2 - Vector3f(ray.o);

    // Permute components of triangle vertices and ray direction
    int kz = MaxDimension(Abs(ray.d));
    int kx = kz + 1;
    if (kx == 3) kx = 0;
    int ky = kx + 1;
    if (ky == 3) ky = 0;
    Vector3f d = Permute(ray.d, kx, ky, kz);
    p0t = Permute(p0t, kx, ky, kz);
    p1t = Permute(p1t, kx, ky, kz);
    p2t = Permute(p2t, kx, ky, kz);

    // Apply shear transformation to translated vertex positions
    Float Sx = -d.x / d.z;
    Float Sy = -d.y / d.z;
    Float Sz = 1.f / d.z;
    p0t.x += Sx * p0t.z;
    p0t.y += Sy * p0t.z;
    p1t.x += Sx * p1t.z;
    p1t.y += Sy * p1t.z;
    p2t.x += Sx * p2t.z;
    p2t.y += Sy * p2t.z;

    // Compute edge function coefficients _e0_, _e1_, and _e2_
    Float e0 = p1t.x * p2t.y - p1t.y * p2t.x;
    Float e1 = p2t.x * p0t.y - p2t.y * p0t.x;
    Float e2 = p0t.x * p1t.y - p0t.y * p1t.x;

    // Fall back to double precision test at triangle edges
    if (sizeof(Float) == sizeof(float) &&
        (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
        double p2txp1ty = (double)p2t.x * (double)p1t.y;
        double p2typ1tx = (double)p2t.y * (double)p1t.x;
        e0 = (float)(p2typ1tx - p2txp1ty);
        double p0txp2ty = (double)p0t.x * (double)p2t.y;
        double p0typ2tx = (double)p0t.y * (double)p2t.x;
        e1 = (float)(p0typ2tx - p0txp2ty);
        double p1txp0ty = (double)p1t.x * (double)p0t.y;
        double p1typ0tx = (double)p1t.y * (double)p0t.x;
        e2 = (float)(p1typ0tx - p1txp0ty);
    }

    // Perform triangle edge and determinant tests
    if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
        return false;
    Float det = e0 + e1 + e2;
    if (det == 0) return false;

    // Compute scaled hit distance to triangle and test against ray $t$ range
    p0t.z *= Sz;
    p1t.z *= Sz;
    p2t.z *= Sz;
    Float tScaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
    if (det < 0 && (tScaled >= 0 || tScaled < ray.tMax * det))
        return false;
    else if (det > 0 && (tScaled <= 0 || tScaled > ray.tMax * det))
        return false;

    // Compute barycentric coordinates and $t$ value for triangle intersection
    Float invDet = 1 / det;
    Float b0 = e0 * invDet;
    Float b1 = e1 * invDet;
    Float b2 = e2 * invDet;
    Float t = tScaled * invDet;

    // Ensure that computed triangle $t$ is conservatively greater than zero

    // Compute $\delta_z$ term for triangle $t$ error bounds
    Float maxZt = MaxComponent(Abs(Vector3f(p0t.z, p1t.z, p2t.z)));
    Float deltaZ = gamma(3) * maxZt;

    // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
    Float maxXt = MaxComponent(Abs(Vector3f(p0t.x, p1t.x, p2t.x)));
    Float maxYt = MaxComponent(Abs(Vector3f(p0t.y, p1t.y, p2t.y)));
    Float deltaX = gamma(5) * (maxXt + maxZt);
    Float deltaY = gamma(5) * (maxYt + maxZt);

    // Compute $\delta_e$ term for triangle $t$ error bounds
    Float deltaE =
        2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

    // Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
    Float maxE = MaxComponent(Abs(Vector3f(e0, e1, e2)));
    Float deltaT = 3 *
                   (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
                   std::abs(invDet);
    if (t <= deltaT) return false;

    // Compute triangle partial derivatives
    Vector3f dpdu, dpdv;
    Point2f uv[3];
    GetUVs(uv);

    // Compute deltas for triangle partial derivatives
    Vector2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
    Vector3f dp02 = p0 - p2, dp12 = p1 - p2;
    Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
    bool degenerateUV = std::abs(determinant) < 1e-8;
    if (!degenerateUV) {
        Float invdet = 1 / determinant;
        dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
        dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
    }
    if (degenerateUV || Cross(dpdu, dpdv).LengthSquared() == 0) {
        // Handle zero determinant for triangle partial derivative matrix
        Vector3f ng = Cross(p2 - p0, p1 - p0);
        if (ng.LengthSquared() == 0)
            // The triangle is actually degenerate; the intersection is
            // bogus.
            return false;

        CoordinateSystem(Normalize(ng), &dpdu, &dpdv);
    }

    // Compute error bounds for triangle intersection
    Float xAbsSum =
        (std::abs(b0 * p0.x) + std::abs(b1 * p1.x) + std::abs(b2 * p2.x));
    Float yAbsSum =
        (std::abs(b0 * p0.y) + std::abs(b1 * p1.y) + std::abs(b2 * p2.y));
    Float zAbsSum =
        (std::abs(b0 * p0.z) + std::abs(b1 * p1.z) + std::abs(b2 * p2.z));
    Vector3f pError = gamma(7) * Vector3f(xAbsSum, yAbsSum, zAbsSum);

    // Interpolate $(u,v)$ parametric coordinates and hit point
    Point3f pHit = b0 * p0 + b1 * p1 + b2 * p2;
    Point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];

    // Test intersection against alpha texture, if present
    if (testAlphaTexture && mesh->alphaMask) {
        SurfaceInteraction isectLocal(pHit, Vector3f(0, 0, 0), uvHit, -ray.d,
                                      dpdu, dpdv, Normal3f(0, 0, 0),
                                      Normal3f(0, 0, 0), ray.time, this);
        if (mesh->alphaMask->Evaluate(isectLocal) == 0) return false;
    }

    // Fill in _SurfaceInteraction_ from triangle hit
    *isect = SurfaceInteraction(pHit, pError, uvHit, -ray.d, dpdu, dpdv,
                                Normal3f(0, 0, 0), Normal3f(0, 0, 0), ray.time,
                                this, faceIndex);

    // Override surface normal in _isect_ for triangle
    isect->n = isect->shading.n = Normal3f(Normalize(Cross(dp02, dp12)));
    if (reverseOrientation ^ transformSwapsHandedness)
        isect->n = isect->shading.n = -isect->n;

    if (mesh->n || mesh->s) {
        // Initialize _Triangle_ shading geometry

        // Compute shading normal _ns_ for triangle
        Normal3f ns;
        if (mesh->n) {
            ns = (b0 * mesh->n[v[0]] + b1 * mesh->n[v[1]] + b2 * mesh->n[v[2]]);
            if (ns.LengthSquared() > 0)
                ns = Normalize(ns);
            else
                ns = isect->n;
        } else
            ns = isect->n;

        // Compute shading tangent _ss_ for triangle
        Vector3f ss;
        if (mesh->s) {
            ss = (b0 * mesh->s[v[0]] + b1 * mesh->s[v[1]] + b2 * mesh->s[v[2]]);
            if (ss.LengthSquared() > 0)
                ss = Normalize(ss);
            else
                ss = Normalize(isect->dpdu);
        } else
            ss = Normalize(isect->dpdu);

        // Compute shading bitangent _ts_ for triangle and adjust _ss_
        Vector3f ts = Cross(ss, ns);
        if (ts.LengthSquared() > 0.f) {
            ts = Normalize(ts);
            ss = Cross(ts, ns);
        } else
            CoordinateSystem((Vector3f)ns, &ss, &ts);

        // Compute $\dndu$ and $\dndv$ for triangle shading geometry
        Normal3f dndu, dndv;
        if (mesh->n) {
            // Compute deltas for triangle partial derivatives of normal
            Vector2f duv02 = uv[0] - uv[2];
            Vector2f duv12 = uv[1] - uv[2];
            Normal3f dn1 = mesh->n[v[0]] - mesh->n[v[2]];
            Normal3f dn2 = mesh->n[v[1]] - mesh->n[v[2]];
            Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
            bool degenerateUV = std::abs(determinant) < 1e-8;
            if (degenerateUV) {
                // We can still compute dndu and dndv, with respect to the
                // same arbitrary coordinate system we use to compute dpdu
                // and dpdv when this happens. It's important to do this
                // (rather than giving up) so that ray differentials for
                // rays reflected from triangles with degenerate
                // parameterizations are still reasonable.
                Vector3f dn = Cross(Vector3f(mesh->n[v[2]] - mesh->n[v[0]]),
                                    Vector3f(mesh->n[v[1]] - mesh->n[v[0]]));
                if (dn.LengthSquared() == 0)
                    dndu = dndv = Normal3f(0, 0, 0);
                else {
                    Vector3f dnu, dnv;
                    CoordinateSystem(dn, &dnu, &dnv);
                    dndu = Normal3f(dnu);
                    dndv = Normal3f(dnv);
                }
            } else {
                Float invDet = 1 / determinant;
                dndu = (duv12[1] * dn1 - duv02[1] * dn2) * invDet;
                dndv = (-duv12[0] * dn1 + duv02[0] * dn2) * invDet;
            }
        } else
            dndu = dndv = Normal3f(0, 0, 0);
        if (reverseOrientation) ts = -ts;
        isect->SetShadingGeometry(ss, ts, dndu, dndv, true);
    }

    *tHit = t;
    ++nHits;
    return true;
}

bool Triangle::IntersectP(const Ray &ray, bool testAlphaTexture) const {
    ProfilePhase p(Prof::TriIntersectP);
    ++nTests;
    // Get triangle vertices in _p0_, _p1_, and _p2_
    const Point3f &p0 = mesh->p[v[0]];
    const Point3f &p1 = mesh->p[v[1]];
    const Point3f &p2 = mesh->p[v[2]];

    // Perform ray--triangle intersection test

    // Transform triangle vertices to ray coordinate space

    // Translate vertices based on ray origin
    Point3f p0t = p0 - Vector3f(ray.o);
    Point3f p1t = p1 - Vector3f(ray.o);
    Point3f p2t = p2 - Vector3f(ray.o);

    // Permute components of triangle vertices and ray direction
    int kz = MaxDimension(Abs(ray.d));
    int kx = kz + 1;
    if (kx == 3) kx = 0;
    int ky = kx + 1;
    if (ky == 3) ky = 0;
    Vector3f d = Permute(ray.d, kx, ky, kz);
    p0t = Permute(p0t, kx, ky, kz);
    p1t = Permute(p1t, kx, ky, kz);
    p2t = Permute(p2t, kx, ky, kz);

    // Apply shear transformation to translated vertex positions
    Float Sx = -d.x / d.z;
    Float Sy = -d.y / d.z;
    Float Sz = 1.f / d.z;
    p0t.x += Sx * p0t.z;
    p0t.y += Sy * p0t.z;
    p1t.x += Sx * p1t.z;
    p1t.y += Sy * p1t.z;
    p2t.x += Sx * p2t.z;
    p2t.y += Sy * p2t.z;

    // Compute edge function coefficients _e0_, _e1_, and _e2_
    Float e0 = p1t.x * p2t.y - p1t.y * p2t.x;
    Float e1 = p2t.x * p0t.y - p2t.y * p0t.x;
    Float e2 = p0t.x * p1t.y - p0t.y * p1t.x;

    // Fall back to double precision test at triangle edges
    if (sizeof(Float) == sizeof(float) &&
        (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
        double p2txp1ty = (double)p2t.x * (double)p1t.y;
        double p2typ1tx = (double)p2t.y * (double)p1t.x;
        e0 = (float)(p2typ1tx - p2txp1ty);
        double p0txp2ty = (double)p0t.x * (double)p2t.y;
        double p0typ2tx = (double)p0t.y * (double)p2t.x;
        e1 = (float)(p0typ2tx - p0txp2ty);
        double p1txp0ty = (double)p1t.x * (double)p0t.y;
        double p1typ0tx = (double)p1t.y * (double)p0t.x;
        e2 = (float)(p1typ0tx - p1txp0ty);
    }

    // Perform triangle edge and determinant tests
    if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
        return false;
    Float det = e0 + e1 + e2;
    if (det == 0) return false;

    // Compute scaled hit distance to triangle and test against ray $t$ range
    p0t.z *= Sz;
    p1t.z *= Sz;
    p2t.z *= Sz;
    Float tScaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
    if (det < 0 && (tScaled >= 0 || tScaled < ray.tMax * det))
        return false;
    else if (det > 0 && (tScaled <= 0 || tScaled > ray.tMax * det))
        return false;

    // Compute barycentric coordinates and $t$ value for triangle intersection
    Float invDet = 1 / det;
    Float b0 = e0 * invDet;
    Float b1 = e1 * invDet;
    Float b2 = e2 * invDet;
    Float t = tScaled * invDet;

    // Ensure that computed triangle $t$ is conservatively greater than zero

    // Compute $\delta_z$ term for triangle $t$ error bounds
    Float maxZt = MaxComponent(Abs(Vector3f(p0t.z, p1t.z, p2t.z)));
    Float deltaZ = gamma(3) * maxZt;

    // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
    Float maxXt = MaxComponent(Abs(Vector3f(p0t.x, p1t.x, p2t.x)));
    Float maxYt = MaxComponent(Abs(Vector3f(p0t.y, p1t.y, p2t.y)));
    Float deltaX = gamma(5) * (maxXt + maxZt);
    Float deltaY = gamma(5) * (maxYt + maxZt);

    // Compute $\delta_e$ term for triangle $t$ error bounds
    Float deltaE =
        2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

    // Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
    Float maxE = MaxComponent(Abs(Vector3f(e0, e1, e2)));
    Float deltaT = 3 *
                   (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
                   std::abs(invDet);
    if (t <= deltaT) return false;

    // Test shadow ray intersection against alpha texture, if present
    if (testAlphaTexture && (mesh->alphaMask || mesh->shadowAlphaMask)) {
        // Compute triangle partial derivatives
        Vector3f dpdu, dpdv;
        Point2f uv[3];
        GetUVs(uv);

        // Compute deltas for triangle partial derivatives
        Vector2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
        Vector3f dp02 = p0 - p2, dp12 = p1 - p2;
        Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
        bool degenerateUV = std::abs(determinant) < 1e-8;
        if (!degenerateUV) {
            Float invdet = 1 / determinant;
            dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
            dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
        }
        if (degenerateUV || Cross(dpdu, dpdv).LengthSquared() == 0) {
            // Handle zero determinant for triangle partial derivative matrix
            Vector3f ng = Cross(p2 - p0, p1 - p0);
            if (ng.LengthSquared() == 0)
                // The triangle is actually degenerate; the intersection is
                // bogus.
                return false;

            CoordinateSystem(Normalize(Cross(p2 - p0, p1 - p0)), &dpdu, &dpdv);
        }

        // Interpolate $(u,v)$ parametric coordinates and hit point
        Point3f pHit = b0 * p0 + b1 * p1 + b2 * p2;
        Point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];
        SurfaceInteraction isectLocal(pHit, Vector3f(0, 0, 0), uvHit, -ray.d,
                                      dpdu, dpdv, Normal3f(0, 0, 0),
                                      Normal3f(0, 0, 0), ray.time, this);
        if (mesh->alphaMask && mesh->alphaMask->Evaluate(isectLocal) == 0)
            return false;
        if (mesh->shadowAlphaMask &&
            mesh->shadowAlphaMask->Evaluate(isectLocal) == 0)
            return false;
    }
    ++nHits;
    return true;
}

Float Triangle::Area() const {
    // Get triangle vertices in _p0_, _p1_, and _p2_
    const Point3f &p0 = mesh->p[v[0]];
    const Point3f &p1 = mesh->p[v[1]];
    const Point3f &p2 = mesh->p[v[2]];
    return 0.5 * Cross(p1 - p0, p2 - p0).Length();
}

Interaction Triangle::Sample(const Point2f &u, Float *pdf) const {
    Point2f b = UniformSampleTriangle(u);
    // Get triangle vertices in _p0_, _p1_, and _p2_
    const Point3f &p0 = mesh->p[v[0]];
    const Point3f &p1 = mesh->p[v[1]];
    const Point3f &p2 = mesh->p[v[2]];
    Interaction it;
    it.p = b[0] * p0 + b[1] * p1 + (1 - b[0] - b[1]) * p2;
    // Compute surface normal for sampled point on triangle
    it.n = Normalize(Normal3f(Cross(p1 - p0, p2 - p0)));
    // Ensure correct orientation of the geometric normal; follow the same
    // approach as was used in Triangle::Intersect().
    if (mesh->n) {
        Normal3f ns(b[0] * mesh->n[v[0]] + b[1] * mesh->n[v[1]] +
                    (1 - b[0] - b[1]) * mesh->n[v[2]]);
        it.n = Faceforward(it.n, ns);
    } else if (reverseOrientation ^ transformSwapsHandedness)
        it.n *= -1;

    // Compute error bounds for sampled point on triangle
    Point3f pAbsSum =
        Abs(b[0] * p0) + Abs(b[1] * p1) + Abs((1 - b[0] - b[1]) * p2);
    it.pError = gamma(6) * Vector3f(pAbsSum.x, pAbsSum.y, pAbsSum.z);
    *pdf = 1 / Area();
    return it;
}

Float Triangle::SolidAngle(const Point3f &p, int nSamples) const {
    // Project the vertices into the unit sphere around p.
    std::array<Vector3f, 3> pSphere = {
        Normalize(mesh->p[v[0]] - p), Normalize(mesh->p[v[1]] - p),
        Normalize(mesh->p[v[2]] - p)
    };

    // http://math.stackexchange.com/questions/9819/area-of-a-spherical-triangle
    // Girard's theorem: surface area of a spherical triangle on a unit
    // sphere is the 'excess angle' alpha+beta+gamma-pi, where
    // alpha/beta/gamma are the interior angles at the vertices.
    //
    // Given three vertices on the sphere, a, b, c, then we can compute,
    // for example, the angle c->a->b by
    //
    // cos theta =  Dot(Cross(c, a), Cross(b, a)) /
    //              (Length(Cross(c, a)) * Length(Cross(b, a))).
    //
    Vector3f cross01 = (Cross(pSphere[0], pSphere[1]));
    Vector3f cross12 = (Cross(pSphere[1], pSphere[2]));
    Vector3f cross20 = (Cross(pSphere[2], pSphere[0]));

    // Some of these vectors may be degenerate. In this case, we don't want
    // to normalize them so that we don't hit an assert. This is fine,
    // since the corresponding dot products below will be zero.
    if (cross01.LengthSquared() > 0) cross01 = Normalize(cross01);
    if (cross12.LengthSquared() > 0) cross12 = Normalize(cross12);
    if (cross20.LengthSquared() > 0) cross20 = Normalize(cross20);

    // We only need to do three cross products to evaluate the angles at
    // all three vertices, though, since we can take advantage of the fact
    // that Cross(a, b) = -Cross(b, a).
    return std::abs(
        std::acos(Clamp(Dot(cross01, -cross12), -1, 1)) +
        std::acos(Clamp(Dot(cross12, -cross20), -1, 1)) +
        std::acos(Clamp(Dot(cross20, -cross01), -1, 1)) - Pi);
}

std::vector<std::shared_ptr<Shape>> CreateTriangleMeshShape(
    const Transform *o2w, const Transform *w2o, bool reverseOrientation,
    const ParamSet &params,
    std::map<std::string, std::shared_ptr<Texture<Float>>> *floatTextures) {
    
    int nvi, npi, nuvi, nsi, nni;
    const int *vi = params.FindInt("indices", &nvi);
    const Point3f *P = params.FindPoint3f("P", &npi);
    const Point2f *uvs = params.FindPoint2f("uv", &nuvi);
    if (!uvs) uvs = params.FindPoint2f("st", &nuvi);
    std::vector<Point2f> tempUVs;
    if (!uvs) {
        const Float *fuv = params.FindFloat("uv", &nuvi);
        if (!fuv) fuv = params.FindFloat("st", &nuvi);
        if (fuv) {
            nuvi /= 2;
            tempUVs.reserve(nuvi);
            for (int i = 0; i < nuvi; ++i)
                tempUVs.push_back(Point2f(fuv[2 * i], fuv[2 * i + 1]));
            uvs = &tempUVs[0];
        }
    }
    if (uvs) {
        if (nuvi < npi) {
            Error(
                "Not enough of \"uv\"s for triangle mesh.  Expected %d, "
                "found %d.  Discarding.",
                npi, nuvi);
            uvs = nullptr;
        } else if (nuvi > npi)
            Warning(
                "More \"uv\"s provided than will be used for triangle "
                "mesh.  (%d expcted, %d found)",
                npi, nuvi);
    }
    if (!vi) {
        Error(
            "Vertex indices \"indices\" not provided with triangle mesh shape");
        return std::vector<std::shared_ptr<Shape>>();
    }
    if (!P) {
        Error("Vertex positions \"P\" not provided with triangle mesh shape");
        return std::vector<std::shared_ptr<Shape>>();
    }
    const Vector3f *S = params.FindVector3f("S", &nsi);
    if (S && nsi != npi) {
        Error("Number of \"S\"s for triangle mesh must match \"P\"s");
        S = nullptr;
    }
    const Normal3f *N = params.FindNormal3f("N", &nni);
    if (N && nni != npi) {
        Error("Number of \"N\"s for triangle mesh must match \"P\"s");
        N = nullptr;
    }
    for (int i = 0; i < nvi; ++i)
        if (vi[i] >= npi) {
            Error(
                "trianglemesh has out of-bounds vertex index %d (%d \"P\" "
                "values were given",
                vi[i], npi);
            return std::vector<std::shared_ptr<Shape>>();
        }

    int nfi;
    const int *faceIndices = params.FindInt("faceIndices", &nfi);
    if (faceIndices && nfi != nvi / 3) {
        Error("Number of face indices, %d, doesn't match number of faces, %d",
              nfi, nvi / 3);
        faceIndices = nullptr;
    }

    std::shared_ptr<Texture<Float>> alphaTex;
    std::string alphaTexName = params.FindTexture("alpha");
    if (alphaTexName != "") {
        if (floatTextures->find(alphaTexName) != floatTextures->end())
            alphaTex = (*floatTextures)[alphaTexName];
        else
            Error("Couldn't find float texture \"%s\" for \"alpha\" parameter",
                  alphaTexName.c_str());
    } else if (params.FindOneFloat("alpha", 1.f) == 0.f)
        alphaTex.reset(new ConstantTexture<Float>(0.f));

    std::shared_ptr<Texture<Float>> shadowAlphaTex;
    std::string shadowAlphaTexName = params.FindTexture("shadowalpha");
    if (shadowAlphaTexName != "") {
        if (floatTextures->find(shadowAlphaTexName) != floatTextures->end())
            shadowAlphaTex = (*floatTextures)[shadowAlphaTexName];
        else
            Error(
                "Couldn't find float texture \"%s\" for \"shadowalpha\" "
                "parameter",
                shadowAlphaTexName.c_str());
    } else if (params.FindOneFloat("shadowalpha", 1.f) == 0.f)
        shadowAlphaTex.reset(new ConstantTexture<Float>(0.f));

    return CreateTriangleMesh(o2w, w2o, reverseOrientation, nvi / 3, vi, npi, P,
                              S, N, uvs, alphaTex, shadowAlphaTex, faceIndices);
}

}  // namespace pbrt
