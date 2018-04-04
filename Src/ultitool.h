#ifndef _ULTITOOL_H_FILE_
#define _ULTITOOL_H_FILE_

#include "MeshDoc.h"
#include <vcg/complex/algorithms/point_sampling.h>
#include <random>

double tri_area(const vcg::Point3f &p1, const vcg::Point3f &p2, const vcg::Point3f &p3);
vcg::Point3f tri_norm(const vcg::Point3f &p1, const vcg::Point3f &p2, const vcg::Point3f &p3);
void tri_random(
    vcg::Point3f &p, vcg::Color4b &c,
    const vcg::Point3f &p1, const vcg::Point3f &p2, const vcg::Point3f &p3,
    const vcg::Color4b &c1, const vcg::Color4b &c2, const vcg::Color4b &c3
);
void RandomTranslate(
	vcg::Point3f &p, const double t, 
	vcg::math::MarsenneTwisterRNG &rnd,
	std::default_random_engine &grnd);
void RandomRotate(
	vcg::Point3f &n, const double r, 
	vcg::math::MarsenneTwisterRNG &rnd,
	std::default_random_engine &grnd);

#endif // !_ULTITOOL_H_FILE_

