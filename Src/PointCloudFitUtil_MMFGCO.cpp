#include "PointCloudFitUtil.h"
#include <vcg/space/index/kdtree/kdtree.h>


////////////////////////////////////////////////
// ------- Multi-Model Fitting with GCO -------
////////////////////////////////////////////////
// MPFGCO : Multi-Plane Fitting with GCO
// MCFGCO : Multi-Cylinder Fitting with GCO
// MMFGCO : Multi-Model Fitting with GCO
// GCO    : Graph Cut Optimization (http://vision.csd.uwo.ca/code/)
// GTE    : Geometric Tools Engine (https://www.geometrictools.com/)

MPFGCONeighbors MPFGCOParseNeighbors(
    CMeshO &mesh,
    const std::vector<int> &ptIndex,
    const int lambda, const double delta,
    const int numNeighbors,
    const double unit_a)
{
    const int NVert = mesh.vert.size();
    const int NPts = ptIndex.size();

    int *NeiCount = new int[NPts];
    int *_NeiIndex = new int[NPts*numNeighbors];
    int **NeiIndex = new int*[NPts];
    int *_NeiWeight = new int[NPts*numNeighbors];
    int **NeiWeight = new int*[NPts];
    for (int i = 0; i < NPts; ++i) {
        NeiCount[i] = 0;
        NeiIndex[i] = _NeiIndex + i*numNeighbors;
        NeiWeight[i] = _NeiWeight + i*numNeighbors;
    }

    //--------------

    int *indexMap = new int[NVert];
    for (int i = 0; i < NVert; ++i)
        indexMap[i] = -1;
    for (int i = 0; i < NPts; ++i) {
        assert(ptIndex[i] < NVert);
        indexMap[ptIndex[i]] = i;
    }
    vcg::VertexConstDataWrapper<CMeshO> ww(mesh);
    vcg::KdTree<float> KDTree(ww);
    vcg::KdTree<float>::PriorityQueue queue;
    const int knnNum = numNeighbors * 2 > 5 ? numNeighbors * 2 : 5;
    const double _r = -1.0 / 2 * (unit_a*unit_a)*(delta*delta);
    CMeshO::VertexIterator vi = mesh.vert.begin();
    for (int i = 0; i < NPts; ++i) {
        vcg::Point3f p = (vi + ptIndex[i])->cP();
        KDTree.doQueryK(p, knnNum, queue);
        int neiNumQuery = queue.getNofElements();
        int neiNum = 0;
        for (int k = 0; k < neiNumQuery; k++) {
            int neightId = queue.getIndex(k);
            if ((vi + neightId)->IsD() ||
                indexMap[neightId] == -1)
                continue;
            // Vpq(lp,lq) = w_{p,q}*d_{lp,lq}
            // d_{lp,lq} = 1 if lp != lq, otherwise 0; 
            // (see implement of smooth cost in function MPFGCOGeneratCost())
            // w_{p,q} = lambda * exp{ - (||p-q||_2/a) ^ 2 / 2*delta^2}
            double w_pq = vcg::SquaredDistance(p, (vi + neightId)->cP());
            w_pq = lambda *exp(w_pq*_r);

            assert(indexMap[neightId] != -1);
            NeiIndex[i][neiNum] = indexMap[neightId];
            NeiWeight[i][neiNum] = int(w_pq + 0.5);
            neiNum++;
            if (neiNum >= numNeighbors)
                break;
        }
        NeiCount[i] = neiNum;
    }
    delete[] indexMap;

    //--------------

    MPFGCONeighbors gcoNei;
    gcoNei.numNeighbor = numNeighbors;
    gcoNei.numSite = NPts;
    gcoNei.neighborsCounts = NeiCount;
    gcoNei.neighborsIndexes = NeiIndex;
    gcoNei.neighborsWeights = NeiWeight;

#if defined(_ReportOut_)
    reportMat<int>(_NeiWeight, NPts, numNeighbors, "../~NeighborWeight~.txt");
#endif

    return gcoNei;
}

MPFGCOCost MPFGCOGeneratCost(
    const std::vector<vcg::Plane3f> &planes,
    const std::vector<vcg::Point3f> &points,
    const std::vector<vcg::Point3f> &norms,
    const double unit_a,
    const int cost_noise,
    const int cost_label)
{
    const bool bHasNorm = norms.empty() ? false : true;
    if (bHasNorm)
        assert(norms.size() == points.size());

    const double angCost = 15.0;
    const double angr = 1.0 / (angCost*angCost);

    const int NPts = points.size();
    const int NPlane = planes.size();
    const int NLabel = NPlane + 1;
    // Data Energy
    // Dp(lp) = ||p-lp||_2 / a + [¡Ï(p,lp)/ang]^2,
    // |p-lp||_2 / a : distance form p to lp (the plane) in unit a.
    int *DataCost = new int[NPts *NLabel];
    for (int i = 0; i < NPts; ++i)
        DataCost[i*NLabel] = cost_noise;

    for (int i = 1; i < NLabel; ++i) {
        vcg::Plane3f ple = planes[i - 1];
        // Dis
        for (int j = 0; j < NPts; ++j) {
            double d = abs(vcg::SignedDistancePlanePoint(ple, points.at(j)));
            DataCost[j*NLabel + i] = int(d / unit_a + 0.5);
        }
        // Ang
        if (bHasNorm) {
            for (int j = 0; j < NPts; ++j) {
                double ang = CheckAng00(VCGAngle(ple.Direction(), norms[j]));
                DataCost[j*NLabel + i] += int(ang*ang*angr);
            }
        }
    }
    // Smooth Energy
    // Vpq(lp,lq) = w_{p,q}*d_{lp,lq}
    // d_{lp,lq} = 1 if lp != lq, otherwise 0; 
    // w_{p,q} = lambda * exp{ - (||p-q||_2/a) ^ 2 / 2*delta^2}
    // (see implement of neighbor system in function MPFGCOParseNeighbors())
    int *SmoothCost = new int[NLabel*NLabel];
    for (int i = 0; i < NLabel; ++i) {
        for (int j = 0; j < NLabel; ++j)
            SmoothCost[i*NLabel + j] = 1;
        SmoothCost[i*NLabel + i] = 0;
    }

    MPFGCOCost gcoCost;
    gcoCost.numLabel = NLabel;
    gcoCost.numSite = NPts;
    gcoCost.dataCost = DataCost;
    gcoCost.smoothCost = SmoothCost;
    gcoCost.labelCost = cost_label;

#if defined(_ReportOut_)
    reportMat<int>(DataCost, NPts, NLabel, "../~CostData~.txt");
    reportMat<int>(SmoothCost, NLabel, NLabel, "../~CostSmooth~.txt");
#endif

    return gcoCost;
}
MPFGCOCost MPFGCOGeneratCost(
    const std::vector<ObjCylinder*> &cylinders,
    const std::vector<vcg::Point3f> &points,
    const std::vector<vcg::Point3f> &norms,
    const double unit_a,
    const int cost_noise,
    const int cost_label)
{
    const bool bHasNorm = norms.empty() ? false : true;
    if (bHasNorm)
        assert(norms.size() == points.size());

    const double angCost = 15.0;
    const double angr = 1.0 / angCost;

    const int NPts = points.size();
    const int NCylinders = cylinders.size();
    const int NLabel = NCylinders + 1;
    // Data Energy
    // Dp(lp) = ||p-lp||_2 / a + [¡Ï(p,lp)/ang],
    // |p-lp||_2 / a : distance form p to lp (the plane) in unit a.
    int *DataCost = new int[NPts *NLabel];
    for (int i = 0; i < NPts; ++i)
        DataCost[i*NLabel] = cost_noise;

    for (int i = 1; i < NLabel; ++i) {
        ObjCylinder *cyl = cylinders[i - 1];
        // Dis
        for (int j = 0; j < NPts; ++j) {
            double d = abs(SignedDistanceCylinderPoint(*cyl, points.at(j)));
            DataCost[j*NLabel + i] = int(d / unit_a + 0.5);
        }
        // Ang
        if (bHasNorm) {
            for (int j = 0; j < NPts; ++j) {
                double ang = AngCylinderPoint(*cyl, points.at(j), norms.at(j));
                DataCost[j*NLabel + i] += int(ang*angr);
            }
        }
    }
    // Smooth Energy
    // Vpq(lp,lq) = w_{p,q}*d_{lp,lq}
    // d_{lp,lq} = 1 if lp != lq, otherwise 0; 
    // w_{p,q} = lambda * exp{ - (||p-q||_2/a) ^ 2 / 2*delta^2}
    // (see implement of neighbor system in function MPFGCOParseNeighbors())
    int *SmoothCost = new int[NLabel*NLabel];
    for (int i = 0; i < NLabel; ++i) {
        for (int j = 0; j < NLabel; ++j)
            SmoothCost[i*NLabel + j] = 1;
        SmoothCost[i*NLabel + i] = 0;
    }

    MPFGCOCost gcoCost;
    gcoCost.numLabel = NLabel;
    gcoCost.numSite = NPts;
    gcoCost.dataCost = DataCost;
    gcoCost.smoothCost = SmoothCost;
    gcoCost.labelCost = cost_label;

#if defined(_ReportOut_)
    reportMat<int>(DataCost, NPts, NLabel, "../~CostData~.txt");
    reportMat<int>(SmoothCost, NLabel, NLabel, "../~CostSmooth~.txt");
#endif

    return gcoCost;
}

std::vector<double> GCOReEstimat(
    std::vector<vcg::Plane3f> &planes,
    const std::vector<vcg::Point3f> &pointList,
    const int *labels, const unsigned int TInlier)
{
    std::vector<double> errors;

    const int NPlane = planes.size();
    const int NPoint = pointList.size();

    std::vector<vcg::Plane3f> optPlanes;
    std::vector<std::vector<int>> planeVerList;
    planeVerList.resize(NPlane + 1);
    std::vector<int> removeIdx;
    for (int i = 0; i < NPoint; ++i) {
        int label = labels[i];
        planeVerList[label].push_back(i);
    }
    
    for (int k = 1; k < NPlane + 1; k++) {
        if (planeVerList[k].size() <= TInlier) {
            flog("    >> Quit the [ No.%d ] plane with [ %d ] points ...\n", k, planeVerList[k].size());
            continue;
        }
        flog("    >> Fine Fit the [ No.%d ] plane ...\n", k);
        vcg::Plane3f plane;
        double err = FinePlane(pointList, planeVerList[k], plane);
        errors.push_back(err);
        optPlanes.push_back(plane);
    }

    planes.swap(optPlanes);

    return errors;
}

std::vector<double> GCOReEstimat(
    std::vector<ObjCylinder*> &cylinders,
    const std::vector<vcg::Point3f> &pointList,
    const int *labels, const unsigned int TInlier)
{
    std::vector<double> errors;

    const int NCylinder = cylinders.size();
    const int NPoint = pointList.size();

    std::vector<ObjCylinder*> optCylinders;
    std::vector<std::vector<int>> cylinderVerList;
    cylinderVerList.resize(NCylinder + 1);
    for (int i = 0; i < NPoint; ++i) {
        int label = labels[i];
        cylinderVerList[label].push_back(i);
    }
    _ResetObjCode(Pt_OnCylinder);
    for (int k = 1; k < NCylinder + 1; k++) {
        if (cylinderVerList[k].size() <= TInlier) {
            flog("    >> Quit the [ No.%d ] cylinder with [ %d < %d ] points ...\n", k, cylinderVerList[k].size(), TInlier);
            continue;
        }
        flog("    >> Fine Fit the [ No.%d ] cylinder with [ %d > %d ] points ...\n", k, cylinderVerList[k].size(), TInlier);
        double err = 0;
        ObjCylinder* cylinder = FineCylinder(pointList, cylinderVerList[k], err);
        errors.push_back(err);
        optCylinders.push_back(cylinder);
    }

    cylinders.swap(optCylinders);
    for (int i = 0; i < optCylinders.size(); ++i)
        delete optCylinders.at(i);

    return errors;
}