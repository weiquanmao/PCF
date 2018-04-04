#ifndef _OBJ2PLY_H_FILE_
#define _OBJ2PLY_H_FILE_

#include "MeshDoc.h"

bool obj2ply(
	const std::string &inputObjFolder,
	const std::string &outputPlyFolder,
	const unsigned int number_sample,
	const unsigned int noise_loc_dis = 0U, 
	const unsigned int noise_norm_ang = 0U);
int savePly(
	const std::string &plyPath,
	const std::vector<vcg::Point3f> &V, 
	const std::vector<vcg::Point3f> &N, 
	const std::vector<vcg::Color4b> &C);
bool addPlyNoise(
	const std::string &plyPath, 
	const unsigned int noise_loc_dis,
	const unsigned int noise_norm_ang);

#endif // !_OBJ2PLY_H_FILE_

