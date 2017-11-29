#include "PointCloudFit.h"
#include "PointCloudFitUtil.h"
#include "gco/GCoptimization.h"
#include <wrap/io_trimesh/io_mask.h>

ObjCylinder* PCFit::DetectCylinderSymAxis()
{
    ObjCylinder *cylinder = 0;
	CMeshO &mesh = m_meshDoc.mesh->cm;

    // Check Norm
    CMeshO::PerVertexAttributeHandle<PtType> type_hi =
        vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
    if (!m_meshDoc.mesh->hasDataMask(vcg::tri::io::Mask::IOM_VERTNORMAL)) {
        flog("\n\n[=DetectCylinderSymAxis=]: [ ): ] Normals are needed to detected cylinder. \n");
        return 0;
    }

	// -- Get Point List And Normal List
    std::vector<int> indexList;
    std::vector<vcg::Point3f> pointList;
    std::vector<vcg::Point3f> normList;
    vcg::Point3f center = GetPointList(indexList, pointList, normList, true);
    assert(normList.size() == pointList.size());
    if (pointList.size() < mesh.vn*Threshold_NPtsCylinder)
        return 0;

	// -- Detect Symmetric Axis
    vcg::Point3f PO, N;
    if (!DetectSymAxis(pointList, normList, PO, N, m_refa))
        return 0;

    _ResetObjCode(Pt_OnCylinder);
    cylinder = new ObjCylinder(_GetObjCode(Pt_OnCylinder));
    cylinder->m_O = PO;
    cylinder->m_N = N;

	// -- Estimate Radius & Identify Surface Points
    std::vector<int> OnClyList;
    AttachToCylinder(OnClyList, *cylinder, pointList, normList, m_refa, 0.0);
    if (OnClyList.size() < pointList.size()*Threshold_NPtsCylinder) {
        delete cylinder;
        return 0;
    }
    
	for (int iter = 0; iter < 0; iter++) {
        PO = cylinder->m_O;
        N = cylinder->m_N;

		std::vector<vcg::Point3f> pointListChecked;
		std::vector<vcg::Point3f> directionListChecked;
		for (int i = 0; i<OnClyList.size(); i++) {
			int index = OnClyList.at(i);
			pointListChecked.push_back(pointList.at(index));
			directionListChecked.push_back(normList.at(index));
		}
		if(!DetectSymAxis(pointListChecked, directionListChecked, PO, N, m_refa)) {
            delete cylinder;
            return 0;
        }

        AttachToCylinder(OnClyList, *cylinder, pointList, normList, m_refa, cylinder->m_radius);
        if (OnClyList.size() < pointList.size()*Threshold_NPtsCylinder) {
            delete cylinder;
            return 0;
        }

		if ( CheckAng00(VCGAngle(N, cylinder->m_N)) < 2) {
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
		//if (CheckAng00(VCGAngle(refN, N)) < 5)
		//	N = refN;

		//delete[] PC;
	}
    
    // -- Set Labels
	for (int i = 0; i<OnClyList.size(); i++)
	{
		int index = indexList.at(OnClyList.at(i));
		type_hi[index] = Pt_OnCylinder;
	}

    // -- Move Back
    cylinder->m_O += center;
    
	return cylinder;
}

std::vector<ObjCylinder*> PCFit::DetectCylinderGCO(const int expCylinderNum, const int iteration)
{
    
    CMeshO &mesh = m_meshDoc.mesh->cm;
    std::vector<ObjCylinder*> clyinders;

    // -- Check Norm
    if (!m_meshDoc.mesh->hasDataMask(vcg::tri::io::Mask::IOM_VERTNORMAL)) {
        flog("\n\n[=DetectCylinderGCO=]: [ ): ] Normals are needed to detected cylinder. \n");
        return clyinders;
    }

    // -- Get Normalized Point List (Moved So That the Center is [0,0])
    std::vector<int> indexList;
    std::vector<vcg::Point3f> pointList;
    std::vector<vcg::Point3f> normList;
    vcg::Point3f center = GetPointList(indexList, pointList, normList, true);
    assert(normList.size() == pointList.size());
    if (pointList.size() < mesh.vn*Threshold_NPtsCylinder)
        return clyinders;

    // -- Detect Init Cylinders
    std::vector<ObjCylinder*> cylCandidates;
    const double TDis = m_refa * Threshold_DisToSurface;
    const double TAng = Threshold_AngToSurface;
    const double inlierRatio = Threshold_NPtsCylinder;
    const int maxN = expCylinderNum;
    vcg::Box3f Box = mesh.bbox; Box.Translate(-center);   
    double maxRatio = DetectCylinderRansac(pointList, normList, cylCandidates, TDis, TAng, maxN, inlierRatio, &Box);
    if (cylCandidates.empty())
        return cylCandidates;


    // -- Fit by GCO
    // E(f)       = Sigma_p{Dp(lp)} + lambda*Sigma_(pg:N){Vpg(lq,pq)} + labelEnergy*|L| .
    // Dp(lp)     = ||p-lp||_2 / a + [¡Ï(p,lp)/ang]^2 , if lp != 0;
    //            = noiseEnergy                      , otherwise (i.e. lp = 0).
    // Vpq(lp,lq) = w_{p,q}*d_{lp,lq}.
    // d_{lp,lq}  = 1 , if lp != lq;
    //            = 0 , otherwise (i.e. lp = lq).
    // w_{p,q}    = lambda * exp{ - (||p-q||_2) ^ 2 / 2*delta^2} .
    // numNeighbors := KNN numbers of neighbor system N, i.e. |N| .
	const int TInlier = inlierRatio*pointList.size();

    const int maxLoop = 3;
    const int maxGCOIteration = iteration > 0 ? iteration : 10;
    const int numNeighbors = 5;
    const int NoiseEnergy = Threshold_DisToSurface*2;
    const int LabelEnergy = TInlier*NoiseEnergy;
    const int lambda = NoiseEnergy;
    const int delta = NoiseEnergy;

    
    
    int *gcoResult = new int[pointList.size()];
    try {
        for (int _iter = 0; _iter < maxLoop; ++_iter) {
			flog("    >> [ No.%d ] loop for MCF-GCO ...\n", _iter+1);
            int numSite = pointList.size();
            int numLabel = cylCandidates.size() + 1;
            GCoptimizationGeneralGraph *gco = new GCoptimizationGeneralGraph(numSite, numLabel);
            MPFGCOCost gcoCost =
                MPFGCOGeneratCost(cylCandidates, pointList, normList, m_refa, NoiseEnergy, LabelEnergy);
            MPFGCONeighbors gcoNei =
                MPFGCOParseNeighbors(mesh, indexList, lambda, delta, numNeighbors, m_refa);
            // -- Set [Data Energy]
            gco->setDataCost(gcoCost.dataCost);
            // -- Set [Smooth Energy]
            gco->setSmoothCost(gcoCost.smoothCost);
            // -- Set [Label Energy]
            gco->setLabelCost(gcoCost.labelCost);
            // -- Set Neighbors System
            gco->setAllNeighbors(gcoNei.neighborsCounts, gcoNei.neighborsIndexes, gcoNei.neighborsWeights);

            gco->setLabelOrder(true);
            gco->setVerbosity(1);
            gco->expansion(maxGCOIteration);

            // -- Get Result <& Re-Estimate>
            for (int i = 0; i < numSite; i++)
                gcoResult[i] = gco->whatLabel(i);
            GCOReEstimat(cylCandidates, pointList, gcoResult, TInlier);

            // -- Check Inlier Ratio
            std::vector<int> inliers;
            for (int i = 0; i < cylCandidates.size(); ++i) {
                int _inliers = CylinderInliers(*cylCandidates.at(i), pointList, normList, TDis, TAng);
                if (_inliers < TInlier) {
                    flog("    >> Quit cylinder [ Id.%d ] with [ %d > %d ] points ...\n", cylCandidates[i]->m_index, _inliers, TInlier);
                    delete cylCandidates[i];
                    cylCandidates.erase(cylCandidates.begin() + i);
                    i--;
                }
                else
                    inliers.push_back(_inliers);
            }
            assert(inliers.size() == CylinderInliers.size());

            // -- Merge The Closer
            for (int i = 0; i < cylCandidates.size()-1; ++i) {
                for (int j = i + 1; j < cylCandidates.size(); ++j) {
                    if (CloseCylinders(*cylCandidates.at(i), *cylCandidates.at(j))) {
                        if (inliers.at(i) < inliers.at(j)) {
                            std::swap(cylCandidates.at(i), cylCandidates.at(j));
                            std::swap(inliers.at(i), inliers.at(j));
                        }
                        flog("    >> Remove cylinder [ Id.%d |%d ] as it is similar to[ Id.%d | %d ]  ...\n",
                            cylCandidates[i]->m_index, inliers[i],
                            cylCandidates[j]->m_index, inliers[j]);
                        cylCandidates.erase(cylCandidates.begin() + j);
                        inliers.erase(inliers.begin() + j);
                        j--;                       
                    }
                }
            }

            // -- Cleaning Up GCO Memory
            gcoCost.memRelease();
            gcoNei.memRelease();

            // -- Continue OR Break
            if (cylCandidates.empty())
                break;
        }
    }
    catch (GCException e) {
        e.Report();
    }
	// -- Cleaning Up GCO Memory
	indexList.clear();
	pointList.clear();
	normList.clear();
	delete[] gcoResult;


	// -- Move Back
	for (int i = 0; i<cylCandidates.size(); ++i)
		cylCandidates.at(i)->m_O += center;

	// -- Attach
	clyinders.swap(cylCandidates);
	resetType(Pt_OnPlane);
	AttachToCylinder(mesh, clyinders, TDis, TAng, inlierRatio);
   
    return clyinders;
}
