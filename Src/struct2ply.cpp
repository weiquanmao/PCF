#include "struct2ply.h"
#include "plyhead.h"
#include "GeometryObject.h"
#include "ultitool.h"

#include <wrap/io_trimesh/export.h>
#include <wrap/io_trimesh/import.h>
#include <vcg/space/index/kdtree/kdtree.h>
#include <vcg/complex/algorithms/point_sampling.h>

#include <windows.h> 
#include <fstream>
#include <iostream>
#include <io.h>

using namespace std;

const unsigned int  nColorChannel = 3;
const unsigned char Color_Gray[nColorChannel] = { 100, 100, 100 };
const unsigned char Color_Solid[nColorChannel] = { 255,   0,   0 };
const unsigned char Color_Solid_Cube[nColorChannel] = { 255,   0,   0 };
const unsigned char Color_Solid_Cylinder[nColorChannel] = { 255,   0,   0 };
const unsigned char Color_Noise[nColorChannel] = { 0, 255,   0 };
const unsigned char Color_Plane[10][nColorChannel] =
{
    { 215,  55, 101 },
    { 255, 245,  55 },
    { 97,  54, 130 },
    { 194,  85,  38 },
    { 0, 192, 210 },
    { 0, 178, 111 },
    { 250, 202,  87 },
    { 255,  95,  61 },
    { 128,  88, 189 },
    { 149,  85,  66 }
};

const double _PI = 3.141592653;
const double _CirN = 120;

double CalArea(const ObjPatch *patch)
{
    double s = 0;
    if (patch != 0) {
        if (patch->type() == Patch_Circle) {
            ObjCircle *circle = (ObjCircle *)patch;
            s = _PI*circle->m_radius*circle->m_radius;
        }
        if (patch->type() == Patch_Rectangle) {
            ObjRect *rect = (ObjRect *)patch;
            s = rect->width() * rect->height();
        }
    }
    return s;
}
double CalArea(const ObjSolid *solid)
{
    double s = 0;
    if (solid != 0) {
        if (solid->type() == Solid_Cube) {
            ObjCube *cube = (ObjCube *)solid;
            double x = cube->DimX();
            double y = cube->DimY();
            double z = cube->DimZ();
            s = 2 * (x*y + x*z + y*z);
        }
        if (solid->type() == Solid_Cylinder) {
            ObjCylinder *cyl = (ObjCylinder *)solid;
            double r = cyl->m_radius;
            double l = cyl->m_length;
            s = 2 * (_PI*r*r) + l*(2 * _PI*r);
        }
    }
    return s;
}
double CalArea(const ObjSet &os)
{
    double S = 0.0;
    for (int i = 0; i < os.m_PlaneList.size(); ++i)
        S += CalArea(os.m_PlaneList[i]);
    for (int i = 0; i < os.m_SolidList.size(); ++i)
        S += CalArea(os.m_SolidList[i]);
    return S;
}

vcg::Color4b IdxColor(int colorcode){
    vcg::Color4b c(
        Color_Plane[colorcode % 10][0],
        Color_Plane[colorcode % 10][1],
        Color_Plane[colorcode % 10][2],
        255);
    return c;
}

int SampleTri(
    vector<vcg::Point3f> &vert,
    vector<vcg::Point3f> &norm,
    vector<vcg::Color4b> &color,
    const vcg::Point3f &p1,
    const vcg::Point3f &p2,
    const vcg::Point3f &p3,
    const vcg::Point3f &n,
    const vcg::Color4b &c,
    const int k)
{
    int count_n = 0;
    for (int i = 0; i < k; ++i, ++count_n) {
        vcg::Point3f _p;
        vcg::Color4b _c;
        tri_random(_p, _c, p1, p2, p3, c, c, c);
        vert.push_back(_p);
        norm.push_back(n);
        color.push_back(_c);
    }
    return count_n;
}

int SamplePatch(
    const ObjPatch *patch,
    vector<vcg::Point3f> &vert, 
    vector<vcg::Point3f> &norm,
    vector<vcg::Color4b> &color,
    const int k)
{
    int count_n = 0;
    if (patch != 0) {
        vcg::Point3f n = patch->m_N;
        vcg::Color4b c = IdxColor(patch->m_index);
        if (patch->type() == Patch_Circle) {
            ObjCircle *circle = (ObjCircle *)patch;
            double _a = 2*_PI / _CirN;
            int _n = k*1.0 / _CirN + 0.5;
            double r = circle->m_radius;
            vcg::Point3f o, u, v;
            o = circle->m_O;
            vcg::GetUV(n, u, v);
            vcg::Point3f r1 = v*r;
            for (int i = 1; i < _CirN; i++) {
                double rang = i*_a;
                vcg::Point3f r2 = (u*sin(rang) + v*cos(rang))*r;
                count_n += SampleTri(
                    vert, norm, color,
                    o, o+r1, o+r2,
                    n, c, _n);
                r1 = r2;
            }
        }
        if (patch->type() == Patch_Rectangle) {
            ObjRect *rect = (ObjRect *)patch;
            int _n = k / 2.0 + 0.5;
            vcg::Point3f p0 = rect->m_O;
            vcg::Point3f p1 = rect->m_O + rect->m_AX;
            vcg::Point3f p2 = rect->m_O + rect->m_AY;
            vcg::Point3f p3 = rect->m_O + rect->m_AX + rect->m_AY;
            count_n += SampleTri(
                vert, norm, color,
                p0, p1, p3,
                n, c, _n);
            count_n += SampleTri(
                vert, norm, color,
                p0, p2, p3,
                n, c, _n);
        }
    }
    return count_n;
}

int SampleSolid(
    const ObjSolid *solid,
    vector<vcg::Point3f> &vert,
    vector<vcg::Point3f> &norm,
    vector<vcg::Color4b> &color,
    const int k)
{
    int count_n = 0;
    if (solid != 0) {
        if (solid->type() == Solid_Cube) {
            ObjCube *cube = (ObjCube *)solid;
            double x = cube->DimX();
            double y = cube->DimY();
            double z = cube->DimZ();
            double s1 = x*y;
            double s2 = x*z;
            double s3 = y*z;
            double s = (s1 + s2 + s3);
            double k1 = s1 / s;
            double k2 = s2 / s;
            double k3 = s3 / s;
            vcg::Point3f pA0 = cube->m_O;
            vcg::Point3f pA1 = cube->m_O + cube->m_AX;
            vcg::Point3f pA2 = cube->m_O + cube->m_AY;
            vcg::Point3f pA3 = cube->m_O + cube->m_AX + cube->m_AY;
            vcg::Point3f pB0 = pA0 + cube->m_AZ;
            vcg::Point3f pB1 = pA1 + cube->m_AZ;
            vcg::Point3f pB2 = pA2 + cube->m_AZ;
            vcg::Point3f pB3 = pA3 + cube->m_AZ;
            vcg::Point3f n1 = cube->m_AX; n1.Normalize();
            vcg::Point3f n2 = cube->m_AY; n2.Normalize();
            vcg::Point3f n3 = cube->m_AZ; n3.Normalize();
            vcg::Point3f n11 = -n1;
            vcg::Point3f n12 = -n2;
            vcg::Point3f n13 = -n3;
            vcg::Color4b c = vcg::Color4b(255, 0, 0, 255);
            count_n += SampleTri(
                vert, norm, color,
                pA0, pA1, pA3,
                n13, c, k*k1 / 4.0);
            count_n += SampleTri(
                vert, norm, color,
                pA0, pA2, pA3,
                n13, c, k*k1 / 4.0);

            count_n += SampleTri(
                vert, norm, color,
                pB0, pB1, pB3,
                n3, c, k*k1 / 4.0);
            count_n += SampleTri(
                vert, norm, color,
                pB0, pB2, pB3,
                n3, c, k*k1 / 4.0);

            count_n += SampleTri(
                vert, norm, color,
                pA0, pA1, pB1,
                n12, c, k*k2 / 4.0);
            count_n += SampleTri(
                vert, norm, color,
                pA0, pB0, pB1,
                n12, c, k*k2 / 4.0);

            count_n += SampleTri(
                vert, norm, color,
                pA0, pA2, pB2,
                n11, c, k*k3 / 4.0);
            count_n += SampleTri(
                vert, norm, color,
                pA0, pB0, pB2,
                n11, c, k*k3 / 4.0);

            count_n += SampleTri(
                vert, norm, color,
                pA3, pA1, pB1,
                n1, c, k*k3 / 4.0);
            count_n += SampleTri(
                vert, norm, color,
                pA3, pB3, pB1,
                n1, c, k*k3 / 4.0);

            count_n += SampleTri(
                vert, norm, color,
                pA3, pA2, pB2,
                n2, c, k*k2 / 2.0);
            count_n += SampleTri(
                vert, norm, color,
                pA3, pB3, pB2,
                n2, c, k*k2 / 4.0);
        }
        if (solid->type() == Solid_Cylinder) {
           ObjCylinder *cyl = (ObjCylinder *)solid;
           vcg::Color4b c = vcg::Color4b(0, 0, 255, 255);

           double r = cyl->m_radius;
           double l = cyl->m_length;

           vcg::Point3f n1 = cyl->m_N;
           vcg::Point3f n2 = -cyl->m_N;
           vcg::Point3f o1 = cyl->m_O + cyl->m_N*(l*0.5);
           vcg::Point3f o2 = cyl->m_O - cyl->m_N*(l*0.5);

           double s = CalArea(cyl);
           double s_c = _PI*r*r;
           double s_l = s - 2 * s_c;
           int k1 = k*(s_c / s) / _CirN+0.5;
           int k2 = k*(s_l / s) / (_CirN*2) + 0.5;
           

           vcg::Point3f u, v;
           vcg::GetUV(n1, u, v);

           vcg::Point3f r1 = v*r;
           double _a = 2*_PI / _CirN;

           for (int i = 1; i < _CirN; i++) {
               double rang = i*_a;
               vcg::Point3f r2 = (u*sin(rang) + v*cos(rang))*r;
               count_n += SampleTri(
                   vert, norm, color,
                   o1, o1+r1, o1+r2,
                   n1, c, k1);
               count_n += SampleTri(
                   vert, norm, color,
                   o2, o2 + r1, o2 + r2,
                   n2, c, k1);
               vcg::Point3f _n = tri_norm(o1 + r1, o2 + r1, o2 + r2);
               count_n += SampleTri(
                   vert, norm, color,
                   o1 + r1, o1 + r2, o2 + r2,
                   _n, c, k2);
               count_n += SampleTri(
                   vert, norm, color,
                   o1 + r1, o2 + r1, o2 + r2,
                   _n, c, k2);
               r1 = r2;
           }
        }
    }
    return count_n;
}

void ReSample(
    const std::string &PlyFile,
    const unsigned int number_sample)
{
    MeshModel *SrtMM = 0;
    int mask = 0;
    SrtMM = new MeshModel(PlyFile.c_str(), "");
    int ret = vcg::tri::io::Importer<CMeshO>::Open(SrtMM->cm, PlyFile.c_str(), mask);

    MeshModel *DstMM = new MeshModel("", "");
    float radius = SurfaceSampling<CMeshO, MeshSampler<CMeshO>>::ComputePoissonDiskRadius(SrtMM->cm, number_sample);
    MeshSampler<CMeshO> mps(DstMM->cm);
    SurfaceSampling<CMeshO, MeshSampler<CMeshO>>::PoissonDiskParam pp;
    SurfaceSampling<CMeshO, MeshSampler<CMeshO>>::PoissonDiskPruningByNumber(mps, SrtMM->cm, number_sample, radius, pp, 0.005);
    mask = DstMM->mask();
    ret = vcg::tri::io::Exporter<CMeshO>::Save(DstMM->cm, PlyFile.c_str(), mask);
    delete SrtMM;
    delete DstMM;
}
bool struct2ply(
    const std::string &StructFile,
    const unsigned int number_sample,
    const std::string &PlyFile)
{
    if (StructFile.empty())
        return false;

    ObjSet *objset = LoadObjSet(StructFile.c_str());
    if (objset == 0)
        return false;


    string outply = PlyFile;
    if (outply.empty())
        outply = StructFile.substr(0, StructFile.rfind('.')) + ".ply";
    
    double TotalS = CalArea(*objset);

	cout
		<< endl
		<< "============================================\n"
		<< "[ -----          Struct2Ply          ----- ]\n"
		<< "[>] " << StructFile << " --> " << outply
		<< endl;
	
	double r = number_sample*1.0 / TotalS;
	vector<vcg::Point3f> vert;
	vector<vcg::Point3f> norm;
	vector<vcg::Color4b> color;
    for (int i = 0; i < objset->m_PlaneList.size(); ++i) {
        int k = CalArea(objset->m_PlaneList[i])*r + 0.5;
        SamplePatch(objset->m_PlaneList[i], vert, norm, color, k);
    }
    for (int i = 0; i < objset->m_SolidList.size(); ++i) {
        int k = CalArea(objset->m_SolidList[i])*r + 0.5;
        SampleSolid(objset->m_SolidList[i], vert, norm, color, k);
    }
       
	cout << ">> # Number of Points :" << vert.size() << endl;
	if (savePly(outply.c_str(), vert, norm, color) != vert.size()) {
		cout << ">> [Error] : Failed to Save Ply File." << endl;
		return false;
    }
    else {
        //cout << ">> Tp ReSample..." << endl;
        //ReSample(outply, number_sample);
        cout << ">> Done." << endl;
    }

	return true;
}

int savePly(
	const string &plyPath,
	const vector<vcg::Point3f> &V,
	const vector<vcg::Point3f> &N,
	const vector<vcg::Color4b> &C)
{
	ofstream fout(plyPath, ios::out);
	if (fout.is_open()) {
		fout.precision(6);

		bool hasNorm = !N.empty();
		bool hasColor = !C.empty();

		if ((hasNorm  && V.size() != N.size()) ||
			(hasColor && V.size() != C.size()))
			return -2;

		int Count = V.size();
		char plyHead[1024];
		if (hasNorm && hasColor) sprintf(plyHead, PlyHeadPNC, Count);
		else if (hasNorm)       sprintf(plyHead, PlyHeadPN, Count);
		else if (hasColor)       sprintf(plyHead, PlyHeadPC, Count);
		else                     sprintf(plyHead, PlyHeadP, Count);
		fout << plyHead;

		for (int i = 0; i < Count; ++i) {
			fout << V[i].X() << " " << V[i].Y() << " " << V[i].Z();
			if (hasNorm)
				fout << " " << N[i].X() << " " << N[i].Y() << " " << N[i].Z();
			if (hasColor)
				fout << " " << (unsigned int)C[i].V(0) << " " << (unsigned int)C[i].V(1) << " " << (unsigned int)C[i].V(2);
			fout << endl;
		}
		fout.close();
		return Count;
	}
	return -1;
}

int PlyPts(const std::string &plyFile)
{
    MeshModel *SrtMM = 0;
    int mask = 0;
    SrtMM = new MeshModel(plyFile.c_str(), "");
    int ret = vcg::tri::io::Importer<CMeshO>::Open(SrtMM->cm, plyFile.c_str(), mask);

    int k = SrtMM->cm.vn;

    delete SrtMM;
    return k;
}