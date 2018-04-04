#include "ultitool.h"
#include "plyHead.h"

#include <vcg/complex/algorithms/point_sampling.h>

double tri_area(const vcg::Point3f &p1, const vcg::Point3f &p2, const vcg::Point3f &p3)
{
    vcg::Point3f d1 = p1 - p3;
    vcg::Point3f d2 = p2 - p3;
    vcg::Point3f s = d1 ^ d2;
    return s.Norm() / 2.0;
}

vcg::Point3f tri_norm(const vcg::Point3f &p1, const vcg::Point3f &p2, const vcg::Point3f &p3)
{
    vcg::Point3f d1 = p2 - p1;
    vcg::Point3f d2 = p3 - p2;
    vcg::Point3f s = d1 ^ d2;
    return s.Normalize();
}
void tri_random(
    vcg::Point3f &p, vcg::Color4b &c,
    const vcg::Point3f &p1, const vcg::Point3f &p2, const vcg::Point3f &p3,
    const vcg::Color4b &c1, const vcg::Color4b &c2, const vcg::Color4b &c3
)
{
    vcg::Point3f u = vcg::tri::SurfaceSampling<CMeshO>::RandomBarycentric();
    p = p1*u[0] + p2*u[1] + p3*u[2];
    c.lerp(c1, c2, c3, u);
}

void RandomTranslate(
	vcg::Point3f &p, const double t, 
	vcg::math::MarsenneTwisterRNG &rnd,
	std::default_random_engine &grnd)
{
	std::normal_distribution<double> distribution(0.0,1.0);
	vcg::Point3f direction = vcg::math::GeneratePointInUnitBallUniform<float>(rnd);
	direction.normalized();
	double l = t*distribution(grnd);

	p += direction*l;
}
void RandomRotate(
	vcg::Point3f &n, const double r,
	vcg::math::MarsenneTwisterRNG &rnd,
	std::default_random_engine &grnd)
{
	std::normal_distribution<double> distribution(1, 1);
	vcg::Point3f direction = vcg::math::GeneratePointInUnitBallUniform<float>(rnd);
	direction.normalized();
	double ang = r*distribution(grnd)*3.141592653/180.0;
	
	vcg::Point3f n1 = direction * (n*direction);
	vcg::Point3f n2 = n - n1;
	vcg::Point3f n3 = direction ^ n2;
	vcg::Point3f n4 = (n2*cos(ang)) + (n3*sin(ang));

	n = n1 + n4;
}