#include "PointCloudFit.h"
#include "PointCloudFitUtil.h"

#include <eigenlib/Eigen/Core>
#include <eigenlib/Eigen/Eigenvalues>


std::vector<ObjCube*> PCFit::DetectCubeFromPlanes(std::vector<ObjPlane*> &planes)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
    std::vector<ObjCube*> cubes;

    const double TRDis = Threshold_PRDis;
    const double TAng = Threshold_PRAng;
    const double TIoU = Threshold_PRIoU;

	// A. Find Correlate Planes
    std::vector< std::vector<ObjPlane*> > CubeFaces;
    int nGroups = CubeFaceInferring(CubeFaces, planes, TRDis, TAng, TIoU, true);
    
	// B. Estimate Cube
    for (int i = 0; i < nGroups; ++i) {
        ObjCube *Cube = CubeMeasure(CubeFaces.at(i), TAng);
        cubes.push_back(Cube);
    }
	// C. Attach Planes To Cube
    int nAdded = AttachToCube(planes, CubeFaces, cubes, TAng, true);
   
	// D. Set Label And Remove Surface Planes
	CMeshO::PerVertexAttributeHandle<PtType> type_hi = 
        vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
    std::vector<int> IdCubePlanes;
    for (int i = 0; i < CubeFaces.size(); ++i)
        for (int j = 0; j < CubeFaces[i].size(); ++j) {
            IdCubePlanes.push_back(CubeFaces[i][j]->m_PlaneIndex);
            delete CubeFaces[i][j];
        }

	for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); vi++) {
		if (!(*vi).IsD()) {
			int code = type_hi[vi];
            if (code & Pt_OnPlane == 0)
                continue;
            for (int k = 0; k < IdCubePlanes.size(); ++k) {
                if (code == IdCubePlanes.at(k)) {
                    type_hi[vi] = Pt_OnCube;
                    break;
                }
            }
		}
	}

	return cubes;
}