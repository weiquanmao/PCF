#include "PointCloudFit.h"
#include "PointCloudFitUtil.h"

std::vector<ObjCube*> PCFit::DetectCubeFromPlanes(std::vector<ObjPatch*> &patches)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
    std::vector<ObjCube*> cubes;

    const double TRDis = Threshold_PRDis;
    const double TAng = Threshold_PRAng;
    const double TIoU = Threshold_PRIoU;
    const double TDis = m_refa*Threshold_DisToSurface;

    // A. Split Rect and Circle
    std::vector<ObjRect*> rects;
    std::vector<ObjCircle*> circles;
    for (int i = 0; i < patches.size(); ++i) {
        ObjPatch *onePatch = patches.at(i);
        if (onePatch->type() == Patch_Rectangle)
            rects.push_back((ObjRect*)onePatch);
        if (onePatch->type() == Patch_Circle)
            circles.push_back((ObjCircle*)onePatch);
    }
    flog("    >> Split [ %d ] patch(es) with [ %d ] rectangle(s) and [ %d ] circles\n",
        patches.size(), rects.size(), circles.size());

	// B. Cube Face Inferring
    flog("    >> Infer cubes for [ %d ] planes ... \n", rects.size());
    std::vector< std::vector<ObjRect*> > CubeFaces;
    int nGroups = CubeFaceInferring(CubeFaces, rects, TRDis, TAng, TIoU, true);
    
	// C. Estimate Cube 
    _ResetObjCode(Pt_OnCube);
    for (int i = 0; i < nGroups; ++i) {
        flog("    >> Estimate [ No.%d ] cube ... \n", i+1);
        ObjCube *oneCube = new ObjCube(_GetObjCode(Pt_OnCube));
        CubeMeasure(CubeFaces.at(i), oneCube, TAng);
        cubes.push_back(oneCube);
    }
	// D. Attach Planes To Cube
    int nAdded = AttachToCube(mesh, rects, CubeFaces, cubes, TAng, TDis, true);
   
	// E. Set Label And Delete Patches Memory
    CMeshO::PerVertexAttributeHandle<PtType> type_hi =
        vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);

    std::vector<int> CubeFaceId;
    for (int i = 0; i < CubeFaces.size(); ++i)
        for (int j = 0; j < CubeFaces[i].size(); ++j) {
            CubeFaceId.push_back(CubeFaces[i][j]->m_index);
            delete CubeFaces[i][j];
        }
	for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); vi++) {      
        if (!(*vi).IsD() && ((type_hi[vi] & Pt_OnPlane) == Pt_OnPlane)) {
            for (int k = 0; k < CubeFaceId.size(); ++k) {
                if (type_hi[vi] == CubeFaceId.at(k)) {
                    type_hi[vi] = Pt_OnCube;
                    break;
                }
            }
        }
    }

    // F. Return Remained Patches
    std::vector<ObjPatch*> _patches;
    for (int i = 0; i < rects.size(); ++i)
        _patches.push_back(rects.at(i));
    for (int i = 0; i < circles.size(); ++i)
        _patches.push_back(circles.at(i));
    patches.swap(_patches);


	return cubes;
}
