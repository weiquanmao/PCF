#include "PointCloudFit.h"
#include "PointCloudFitUtil.h"

#include <eigenlib/Eigen/Core>
#include <eigenlib/Eigen/Eigenvalues>


std::vector<ObjCube*> PCFit::DetectCubeFromPlanes(std::vector<ObjPatch*> &planes)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
    std::vector<ObjCube*> cubes;

    const double TRDis = Threshold_PRDis;
    const double TAng = Threshold_PRAng;
    const double TIoU = Threshold_PRIoU;
    // A. Split Rect and Circle
    std::vector<ObjRect*> rects;
    std::vector<ObjCircle*> circles;
    for (int i = 0; i < planes.size(); ++i) {
        ObjPatch *patch = planes.at(i);
        if (patch->type() == Patch_Rectangle)
            rects.push_back((ObjRect*)patch);
        if (patch->type() == Patch_Circle)
            circles.push_back((ObjCircle*)patch);
    }

	// B. Find Correlate Planes
    std::vector< std::vector<ObjRect*> > CubeFaces;
    int nGroups = CubeFaceInferring(CubeFaces, rects, TRDis, TAng, TIoU, true);
    
	// C. Estimate Cube
    for (int i = 0; i < nGroups; ++i) {
        ObjCube *Cube = new ObjCube(Pt_OnCube + i);
        CubeMeasure(Cube, CubeFaces.at(i), TAng);
        cubes.push_back(Cube);
    }
	// D. Attach Planes To Cube
    // int nAdded = AttachToCube(rects, CubeFaces, cubes, TAng, true);
   
	// E. Set Label
	CMeshO::PerVertexAttributeHandle<PtType> type_hi = 
        vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
    std::vector<int> IdCubePlanes;
    for (int i = 0; i < CubeFaces.size(); ++i)
        for (int j = 0; j < CubeFaces[i].size(); ++j) {
            IdCubePlanes.push_back(CubeFaces[i][j]->m_index);
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

    // F. Return planes
    std::vector<ObjPatch*> _planes;
    for (int i = 0; i < rects.size(); ++i)
        _planes.push_back(rects.at(i));
    for (int i = 0; i < circles.size(); ++i)
        _planes.push_back(circles.at(i));
    planes.swap(_planes);


	return cubes;
}