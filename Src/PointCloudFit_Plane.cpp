#include "PointCloudFit.h"
#include "PointCloudFitUtil.h"
#include "gco/GCoptimization.h"

std::vector<ObjPlane*> PCFit::DetectPlanesHT(const int expPN)
{
    CMeshO &mesh = m_meshDoc.mesh->cm;

    // -- Get Normalized Point List (Moved So That the Center is [0,0])
    std::vector<int> indexList;
    std::vector<vcg::Point3f> pointList;
    std::vector<vcg::Point3f> normList;
    vcg::Point3f center = GetPointList(indexList, pointList, normList, true);

	// -- Calculate intercept
	double candicate1 = mesh.bbox.Diag() * sqrt(3.0) / 2.0;
	double candicate2 = mesh.bbox.DimX() + mesh.bbox.DimY() + mesh.bbox.DimZ();
	double intercept = (candicate1 < candicate2) ? candicate1 : candicate2;

    // -- Calculate Thresholds
	const double _planeDisThreshold = m_refa*Threshold_DisToPlane;
	const double _planeAngThreshold = Threshold_AngToPlane;
	const int _THard = fmax(300, pointList.size()*0.01);

    // -- Detect Planes
    std::vector<ObjPlane*> planes;
    std::vector<vcg::Point4f> infPlanes; // Useless

    infPlanes = DetectHTPlanes(
        pointList, normList,
        intercept, m_refa, Precision_HT,
        _planeDisThreshold, _planeAngThreshold,
        Threshold_NPtsPlane, _THard, expPN,
        0,
        &planes, &mesh, &indexList);

    // -- Move Back
    for (int i = 0; i<planes.size(); ++i) {
        planes.at(i)->m_pO += center;
    }

	indexList.clear();
	pointList.clear();
    normList.clear();

	return planes;
}

std::vector<ObjPlane*> PCFit::DetectPlanesGCO(const int expPN, const int iteration)
{ 
    std::vector<ObjPlane*> planes;
    CMeshO &mesh = m_meshDoc.mesh->cm;

    // -- Get Normalized Point List (Moved So That the Center is [0,0])
    std::vector<int> indexList;
    std::vector<vcg::Point3f> pointList;
    std::vector<vcg::Point3f> normList;
    vcg::Point3f center = GetPointList(indexList, pointList, normList, true);

    // -- Calculate intercept
    double candicate1 = mesh.bbox.Diag() * sqrt(3.0) / 2.0;
    double candicate2 = mesh.bbox.DimX() + mesh.bbox.DimY() + mesh.bbox.DimZ();
    double intercept = (candicate1 < candicate2) ? candicate1 : candicate2;

    // -- Calculate Thresholds
    const double _planeDisThreshold = m_refa*Threshold_DisToPlane;
    const double _planeAngThreshold = Threshold_AngToPlane;
    const int _THard = fmax(300, pointList.size()*0.01);

    // -- Detect Init Planes
    std::vector<vcg::Point4f> infPlanes;
    std::vector<double> errors;
    infPlanes = DetectHTPlanes(
        pointList, normList,
        intercept, m_refa, Precision_HT,
        _planeDisThreshold, _planeAngThreshold,
        //0.0, 0, expPN,
        Threshold_NPtsPlane, _THard, expPN,
        &errors);

    // -- Fit by GCO
    // E(f)       = Sigma_p{Dp(lp)} + lambda*Sigma_(pg:N){Vpg(lq,pq)} + labelEnergy*|L| .
    // Dp(lp)     = ||p-lp||_2 / a + [¡Ï(p,lp)/ang]^2 , if lp != 0;
    //            = noiseEnergy                      , otherwise (i.e. lp = 0).
    // Vpq(lp,lq) = w_{p,q}*d_{lp,lq}.
    // d_{lp,lq}  = 1 , if lp != lq;
    //            = 0 , otherwise (i.e. lp = lq).
    // w_{p,q}    = lambda * exp{ - (||p-q||_2) ^ 2 / 2*delta^2} .
    // numNeighbors := KNN numbers of neighbor system N, i.e. |N| .
    int NoiseEnergy = 4;
    int LabelEnergy = 200;
    int lambda = 20;
    int delta = 10;

    int numNeighbors = 7;

    int *gcoResult = new int[pointList.size()];
	try{
        int numSite = pointList.size();
        int numLabel = infPlanes.size() + 1;
		GCoptimizationGeneralGraph *gco = new GCoptimizationGeneralGraph(numSite, numLabel);
        MPFGCOCost gcoCost = 
            MPFGCOGeneratCost(infPlanes, pointList, normList, m_refa, NoiseEnergy, LabelEnergy);
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
        gco->expansion(iteration);
        
        // -- Get Result <& Re-Estimate>
        for (int i = 0; i < numSite; i++)
            gcoResult[i] = gco->whatLabel(i);
        // errors = GCOReEstimat(infPlanes, pointList, gcoResult);

        // -- Cleaning Up
        gcoCost.memRelease();
        gcoNei.memRelease();
	}
	catch (GCException e){
		e.Report();
	}
    planes = ExtractPlanes(mesh, indexList, pointList, infPlanes.size(), gcoResult);


    // -- Move Back
    for (int i = 0; i<planes.size(); ++i) {
        planes.at(i)->m_pO += center;
    }
    
    
    indexList.clear();
    pointList.clear();
    normList.clear();
    delete[] gcoResult;

    return planes;
}