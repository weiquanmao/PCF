#include "PointCloudFit.h"
#include "PointCloudFitUtil.h"


ObjCylinder* PCFit::DetectCylinderSymAxis()
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
    CMeshO::PerVertexAttributeHandle<PtType> type_hi =
        vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, _MyPtAttri);

	// -- Get Point List And Normal List
    std::vector<int> indexList;
    std::vector<vcg::Point3f> pointList;
    std::vector<vcg::Point3f> normList;
    vcg::Point3f center = GetPointList(indexList, pointList, normList, true);
    assert(normList.size() == pointList.size());

    ObjCylinder *cylinder = new ObjCylinder();
    vcg::Point3f PO, N;

	// -- Detect Symmetric Axis
    if (!DetectSymAxis(pointList, normList, m_refa, PO, N)) {
        delete cylinder;
        return 0;
    }

    cylinder->m_pO = PO;
    cylinder->m_N = N;

	// -- Estimate Radiu & Identify Surface Points
    std::vector<int> OnClyList;
    AttachToCylinder(OnClyList, *cylinder, pointList, normList, m_refa, 0.0);
    if (OnClyList.size() < pointList.size()*Threshold_NPtsCyl) {
        delete cylinder;
        return 0;
    }
    
	for (int iter = 0; iter < 0; iter++) {
        PO = cylinder->m_pO;
        N = cylinder->m_N;

		std::vector<vcg::Point3f> pointListChecked;
		std::vector<vcg::Point3f> directionListChecked;
		for (int i = 0; i<OnClyList.size(); i++) {
			int index = OnClyList.at(i);
			pointListChecked.push_back(pointList.at(index));
			directionListChecked.push_back(normList.at(index));
		}
		if(!DetectSymAxis(pointListChecked, directionListChecked, m_refa, PO, N)) {
            delete cylinder;
            return 0;
        }

        AttachToCylinder(OnClyList, *cylinder, pointList, normList, m_refa, cylinder->m_radius);
        if (OnClyList.size() < pointList.size()*Threshold_NPtsCyl) {
            delete cylinder;
            return 0;
        }

		if (90 - abs(90 - vcg::Angle(N, cylinder->m_N)*_R2D) < 2) {
            delete cylinder;
            return 0;
        }

		//const int NOnCly = OnClyList.size();
		//double *data = new double[3 * NOnCly];
		//double *X = data;
		//double *Y = data + NOnCly;
		//double *Z = data + 2 * NOnCly;

		//for (int i = 0; i<NOnCly; i++)
		//{
		//	int index = OnClyList.at(i);
		//	X[i] = pointList.at(index).X();
		//	Y[i] = pointList.at(index).Y();
		//	Z[i] = pointList.at(index).Z();
		//}
		//int row = 3, col = NOnCly;
		//double *PC, V[3];//V[3] USELESS
		//PC = new double[9];
		//int ret = PCA(data, row, col, PC, V);
		//delete[] data;
		//if (ret == -1)
		//	return false;
		//vcg::Point3f refN = vcg::Point3f(PC[0], PC[3], PC[6]);
		//if ((90 - abs(90 - vcg::Angle(refN, N)*_R2D)) < 5)
		//	N = refN;

		//delete[] PC;
	}

    
	for (int i = 0; i<OnClyList.size(); i++)
	{
		int index = indexList.at(OnClyList.at(i));
		type_hi[index] = Pt_OnCylinder;
	}

    cylinder->m_pO += center;
	return cylinder;
}
