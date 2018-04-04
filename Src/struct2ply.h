#ifndef _OBJ2PLY_H_FILE_
#define _OBJ2PLY_H_FILE_

#include "MeshDoc.h"

int PlyPts(const std::string &plyFile);
bool struct2ply(
	const std::string &StructFile,
	const unsigned int number_sample,
    const std::string &PlyFile = "");
int savePly(
	const std::string &plyPath,
	const std::vector<vcg::Point3f> &V, 
	const std::vector<vcg::Point3f> &N, 
	const std::vector<vcg::Color4b> &C);

#endif // !_OBJ2PLY_H_FILE_

