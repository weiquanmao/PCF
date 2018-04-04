#include "PointCloudFit.h"
#include "PCA/PCA.h"

#include <wrap/io_trimesh/io_mask.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/space/index/kdtree/kdtree.h>

int LabelDeleteAsNoise(CMeshO &mesh)
{
    int nNoise = 0;
    CMeshO::PerVertexAttributeHandle<PtType> type_hi =
        vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
    for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi) {
        if ((*vi).IsD()) {
            nNoise++;
            type_hi[vi] = Pt_Noise;
        }
    }
    return nNoise;
}
int RegionGrow(
    CMeshO &mesh, std::vector<std::vector<int>> &clusterIdx,
    const int stepn, const double dis)
{
    assert(dis >= 0);
    assert(stepn >= 3);

    // 0. Bulid KD-Tree
    vcg::VertexConstDataWrapper<CMeshO> ww(mesh);
    vcg::KdTree<float> KDTree(ww);
    vcg::KdTree<float>::PriorityQueue queue;

    // 2. Iteration
    std::vector<std::vector<int>> _clusters;
    int clusterID = 0;    // No. of Clusters
    int maxCluster = 0;	 // Index of the Largest Cluster

    CMeshO::PerVertexAttributeHandle<PtType> type_hi =
        vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
    int idx = 0;
    for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi, ++idx) {
        if (!(vi->IsD()) && type_hi[vi] == Pt_Undefined) {
            // Seed Point
            clusterID++;
            std::vector<int> _cluster;
            std::vector<CMeshO::VertexIterator> _clusterDump;
            _clusterDump.push_back(vi);           
            _cluster.push_back(idx);
            type_hi[vi] = clusterID;
            // Grow by KNN
            while (!_clusterDump.empty()) {
                CMeshO::VertexIterator curSeed = *(_clusterDump.end() - 1);
                _clusterDump.pop_back();
#if 1
                KDTree.doQueryK(curSeed->cP(), stepn, queue);
                int neighbours = queue.getNofElements();
                for (int k = 0; k < neighbours; k++) {
                    int nborIdx = queue.getIndex(k);
                    CMeshO::VertexIterator nbor = mesh.vert.begin() + nborIdx;
                    if (!(nbor->IsD()) && type_hi[nbor] == Pt_Undefined &&
                        Distance(curSeed->cP(), nbor->cP()) < dis) {
                        _clusterDump.push_back(nbor);
                        _cluster.push_back(nborIdx);
                        type_hi[nbor] = clusterID;
                    }
                }
#else // Query All at once
                std::vector<unsigned int> neiPtIdx;
                std::vector<float> neiPtSaureDis;
                KDTree.doQueryDist(curSeed->cP(), dis*dis, neiPtIdx, neiPtSaureDis);
                int neighbours = neiPtIdx.size();
                for (int k = 0; k < neighbours; k++) {
                    CMeshO::VertexIterator nbor = mesh.vert.begin() + neiPtIdx[k];
                    if (!(nbor->IsD()) && type_hi[nbor] == Pt_Undefined) {
                        type_hi[nbor] = clusterID;
                        _clusterDump.push_back(nbor);
                        _cluster.push_back(neiPtIdx[k]);
                    }
                }
#endif
            }
            _clusters.push_back(_cluster);
            if (_cluster.size() > _clusters.at(maxCluster).size())
                maxCluster = clusterID-1;
        }
    }
    clusterIdx.swap(_clusters);

    return maxCluster;
}
double MeshSizeAlongN(CMeshO &mesh, const vcg::Point3f n)
{
    double minDis = DBL_MAX;
    double maxDis = -DBL_MAX;
    vcg::Point3f normalN = n;
    normalN.Normalize();
    for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi)
    {
        if ((*vi).IsD())
            continue;
        double proLoc = (*vi).cP()*normalN;
        if (proLoc > maxDis) maxDis = proLoc;
        if (proLoc < minDis) minDis = proLoc;
    }
    return abs(maxDis - minDis);
}

//-----------------------------------

vcg::Point3f PCFit::GetPointList(
    std::vector<int> &indexList,
    std::vector<vcg::Point3f> &pointList,
    std::vector<vcg::Point3f> &normList,
    const bool moveToCenter)
{
    CMeshO &mesh = m_meshDoc.mesh->cm;
    bool normalSupportted = m_meshDoc.mesh->hasDataMask(vcg::tri::io::Mask::IOM_VERTNORMAL);

    vcg::tri::UpdateBounding<CMeshO>::Box(mesh);
    vcg::Point3f center(0.0, 0.0, 0.0);
    if (moveToCenter) {
        center = mesh.bbox.Center();
        flog(
            "    [--GetBorder--]: #nPts-%d \n"
            "      | #MinPt  : < %7.3f, %7.3f, %7.3f > \n"
            "      | #MaxPt  : < %7.3f, %7.3f, %7.3f > \n"
            "      | #Center : < %7.3f, %7.3f, %7.3f > \n"
            "    [--GetBorder--]: Done in %.4f seconds.\n",
            mesh.VN(),
            mesh.bbox.min.X(), mesh.bbox.min.Y(), mesh.bbox.min.Z(),
            mesh.bbox.max.X(), mesh.bbox.max.Y(), mesh.bbox.max.Z(),
            center.X(), center.Y(), center.X());
    }
    indexList.clear();
    pointList.clear();
    normList.clear();
    // -- Get Point List And Normal List (If Exist)   
    indexList.reserve(mesh.vn);
    pointList.reserve(mesh.vn);
    if (normalSupportted)
        normList.reserve(mesh.vn);
    CMeshO::PerVertexAttributeHandle<PtType> type_hi =
        vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
    int index = 0;
    for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++index, ++vi)
    {
        if (type_hi[vi] == Pt_Undefined && !(vi->IsD())) {
            indexList.push_back(index);
            pointList.push_back((*vi).cP() - center);
            if (normalSupportted)
                normList.push_back((*vi).cN());
        }
    }
    if (moveToCenter)
        flog("    >> Got #%d Normalized Pts.\n", pointList.size());
    else
        flog("    >> Got #%d Pts.\n", pointList.size());
    return center;
}


int PCFit::DeNoiseKNN()
{
	CMeshO &mesh = m_meshDoc.mesh->cm;	

	QTime time;
	time.start();
	flog(
		"    [--DeNoiseKNN--]: #nPts-%d \n"
		"      | #NIteration : %d \n"
		"      | #NNeighbord : %d \n"
		"      | #DisRatio   : %.3f \n"
		"      | #Gap        : ",
		mesh.vn, DeNoise_MaxIteration, DeNoise_KNNNeighbors, DeNoise_DisRatioOfOutlier);
    //----[[

	// 1. Build Kd-tree
	vcg::VertexConstDataWrapper<CMeshO> ww(mesh);
	vcg::KdTree<float> KDTree(ww);
	vcg::KdTree<float>::PriorityQueue queue;

	// 2. Iterations
	const int nIter = DeNoise_MaxIteration <= 0 ? 100 : DeNoise_MaxIteration;
	for (int i = 0; i<nIter; i++) {
		bool hasOutlier = false;

		// 1)Update Threshold
		vcg::tri::UpdateBounding<CMeshO>::Box(mesh);
        double Gap = mesh.bbox.Dim().V(mesh.bbox.MinDim())*DeNoise_DisRatioOfOutlier;
		//----
		flog("[%d]-%.3f ", i + 1, Gap);
		//----

		// 2)KNN Check for Each Point
		for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi) {
			if (!(*vi).IsD()) {
				KDTree.doQueryK((*vi).cP(), DeNoise_KNNNeighbors, queue);
				int neighbours = queue.getNofElements();
				float avgDist = 0;
				for (int i = 0; i < neighbours; i++) {
					int neightId = queue.getIndex(i);
					avgDist += Distance((*vi).cP(), mesh.vert[neightId].cP());
				}
				if (avgDist / neighbours > Gap) {
					vcg::tri::Allocator<CMeshO>::DeleteVertex(mesh, *vi);
					hasOutlier = true;
				}
			}
		}

		if (!hasOutlier)
			break;
	}

	// 3. Delete Noise Pts
    int nNoise = LabelDeleteAsNoise(mesh);

	// 4. Update Box Size
	vcg::tri::UpdateBounding<CMeshO>::Box(mesh);

	//----]]
	flog("\n    [--DeNoiseKNN--]: Done, %d outlier(s) were found in %.4f seconds.\n",
		nNoise, time.elapsed()/1000.0);

	return nNoise;
}
int PCFit::DeNoiseRegGrw()
{
	CMeshO &mesh = m_meshDoc.mesh->cm;

	QTime time;
	time.start();
	flog(
		"    [--DeNoiseRegGrw--]: #nPts-%d \n"
		"      | #NIteration : %d \n"
        "      | #KNNStep    : %d \n"
		"      | #DisRatio   : %.3f \n"
		"      | #Gap        : ",
		mesh.vn,
		DeNoise_MaxIteration, DeNoise_GrowNeighbors, DeNoise_DisRatioOfOutlier);
    //----[[

    // 1. Iteration Grow
    const int nStep = DeNoise_GrowNeighbors;
	const int nIter = DeNoise_MaxIteration <= 0 ? 100 : DeNoise_MaxIteration;
	CMeshO::PerVertexAttributeHandle<PtType> type_hi = 
		vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
	for (int _iter = 0; _iter<nIter; _iter++) {
		// 1)Update Threshold
		vcg::tri::UpdateBounding<CMeshO>::Box(mesh);
        double Gap = mesh.bbox.Dim().V(mesh.bbox.MinDim())*DeNoise_DisRatioOfOutlier;

        // 2)Region Growing
        std::vector<std::vector<int>> clusterIdx;
		int maxCluster = RegionGrow(mesh, clusterIdx, nStep, Gap);
        CMeshO::VertexIterator vi = mesh.vert.begin();
        for (int i = 0; i < clusterIdx.size(); ++i) {
            for (int j = 0; j < clusterIdx[i].size(); ++j)
                type_hi[clusterIdx[i][j]] = Pt_Undefined;
            if (i != maxCluster) {
                for (int j = 0; j < clusterIdx[i].size(); ++j)
                    vcg::tri::Allocator<CMeshO>::DeleteVertex(mesh, *(vi + clusterIdx[i][j]));
            }            
        }

		//----
		flog("[%d]-%.3f-<%d in %d> ", _iter+1, Gap, clusterIdx[maxCluster].size(), clusterIdx.size());
		//----

		if (clusterIdx.size() == 1)
			break;
	}
    // 2. Delete Noise Pts
    int nNoise = LabelDeleteAsNoise(mesh);

    // 3. Update Box Size
	vcg::tri::UpdateBounding<CMeshO>::Box(mesh);

	//----]]
	flog("\n    [--DeNoiseRegGrw--]: Done, found %d outlier(s) in %.4f seconds.\n",
		nNoise, time.elapsed() / 1000.0);

	return nNoise;
}

bool PCFit::PCADimensions(	
	std::vector<vcg::Point3f> &PDirections,
    vcg::Point3f &PSize)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;

	QTime time;
	time.start();
    //----[[
    // 1. Call PCA
	const int N = mesh.vn;
	int row = 3, col = N;
	double *data = new double[3 * N];
	int Count = 0;
	for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi)
	{
		if ((*vi).IsD())
			continue;
		data[Count] = (*vi).cP().X();
		data[N + Count] = (*vi).cP().Y();
		data[2 * N + Count] = (*vi).cP().Z();
		++Count;
	}
	assert(Count == N);
	double PC[9], V[3];//V[3] useless
	int ret = PCA(data, row, col, PC, V);
	delete[] data;
	if (ret == -1) {
        //----
		flog("    [--PCADimensions--] PCA failed. [--PCADimensions--]\n");
        //----
		doFailure();
	}
    
    // 2. Normalize And Assign
	vcg::Point3f NX = vcg::Point3f(PC[0], PC[3], PC[6]);
	vcg::Point3f NY = vcg::Point3f(PC[1], PC[4], PC[7]);
	vcg::Point3f NZ = vcg::Point3f(PC[2], PC[5], PC[8]);
	NX.Normalize();
	NY.Normalize();
	NZ.Normalize();

	PSize.X() = MeshSizeAlongN(mesh,NX);
	PSize.Y() = MeshSizeAlongN(mesh,NY);
	PSize.Z() = MeshSizeAlongN(mesh,NZ);

	PDirections.clear();
	PDirections.push_back(NX);
	PDirections.push_back(NY);
	PDirections.push_back(NZ);

    //----]]
	flog(
		"    [--PCADimensions--] #bPts - %d/%d \n"
		"       | -- Principal direction vectors, -- \n"
        "       | -- col by decreasing order.     -- \n"
		"       |  %7.4f %7.4f %7.4f - < %7.4f > \n"
		"       |  %7.4f %7.4f %7.4f - < %7.4f > \n"
		"       |  %7.4f %7.4f %7.4f - < %7.4f > \n"
		"    [--PCADimensions--]: Done in %.4f seconds.\n",
		Count, N,
		NX.X(), NY.X(), NZ.X(), PSize.X(),
		NX.Y(), NY.Y(), NZ.Y(), PSize.Y(),
		NX.Z(), NY.Z(), NZ.Z(), PSize.Z(),
		time.elapsed() / 1000.0);

	return true;
}
double PCFit::Roughness()
{
	CMeshO &mesh = m_meshDoc.mesh->cm;

	QTime time;
	time.start();

    //----[[

    // 1. Build KD-Tree
	vcg::VertexConstDataWrapper<CMeshO> ww(mesh);
	vcg::KdTree<float> KDTree(ww);
	vcg::KdTree<float>::PriorityQueue queue;

    // 2. Fit at Each Point
    // Key Parameter
    // rafa = a*miu + b*std	
    const int knn = 30;
	const int a = 1.0;
	const int b = 3.0;
	const int N = mesh.vn;
	
	int Count = 0;
	double R = 0.0;
	for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi)
	{
		if ((*vi).IsD())
			continue;

        // 1) Query KNN
		KDTree.doQueryK((*vi).cP(), knn, queue);
		int neighbours = queue.getNofElements();
		
        // 2) Fit KNN
		std::vector<vcg::Point3f> Pts;
		vcg::Plane3f ple;
		for (int i = 0; i < neighbours; i++) {
			int neightId = queue.getIndex(i);
			Pts.push_back(mesh.vert[neightId].cP());
		}
		vcg::FitPlaneToPointSet(Pts, ple);

        // 3) Cal u and v
		double miu = 0.0;
		double sigma = 0.0;
		for (int i = 0; i < neighbours; ++i) {
			float d = vcg::SignedDistancePlanePoint(ple, Pts[i]);
			miu += abs(d);
			sigma += d*d;
		}
		R += a * (miu / neighbours) + b * sqrt(sigma / neighbours);
		++Count;
	}
	assert(Count == N);
	double roughness = R / Count;

    //----]]

	flog(
		"    [--Roughness--] #bPts - %d/%d \n"
		"       | RoughnessA - < %7.4f > \n"
		"    [--Roughness--]: Done in %.4f seconds.\n",
		Count, N, roughness, time.elapsed() / 1000.0);

	return roughness;
}
