#include "PointCloudFit.h"
#include "PCA.h"

#include <QTime>
#include <wrap/io_trimesh/io_mask.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/space/index/kdtree/kdtree.h>
#include <vcg/space/fitting3.h>

int PCFit::DeNoiseKNN()
{
	CMeshO &mesh = m_meshDoc.mesh->cm;	

	QTime time;
	time.start();
	printf(
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
		double DX = mesh.bbox.DimX();
		double DY = mesh.bbox.DimY();
		double DZ = mesh.bbox.DimZ();
		double Gap = DX<DY ? (DX<DZ ? DX : DZ) : (DY<DZ ? DY : DZ);
		Gap *= DeNoise_DisRatioOfOutlier;

		//----
		printf("[%d]-%.3f ", i + 1, Gap);
		//----

		// 2)Point Check
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
	int nNoise = 0;
	CMeshO::PerVertexAttributeHandle<SatePtType> type_hi = 
		vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<SatePtType>(mesh, _MySatePtAttri);
	for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi) {
		if ((*vi).IsD()) {
			nNoise++;
			type_hi[vi] = Pt_Noise;
		}
	}

	// 4. Update Box Size
	vcg::tri::UpdateBounding<CMeshO>::Box(mesh);

	//----]]
	printf(
		"\n    [--DeNoiseKNN--]: Done, %d outlier(s) were found in %.4f seconds.\n",
		nNoise, time.elapsed()/1000.0);

	return nNoise;
}

int RegionGrow(
	CMeshO &mesh, 
    const int stepn, const double dis, 
    std::vector<int> &nums, int *maxID, long *maxN)
{
	assert(dis >= 0);
	assert(stepn >= 3);

	// 0. Bulid KD-Tree
	vcg::VertexConstDataWrapper<CMeshO> ww(mesh);
	vcg::KdTree<float> KDTree(ww);
	vcg::KdTree<float>::PriorityQueue queue;
	
	// 2. Iteration
	nums.clear();
	int nCluster = 0;    // No. of Clusters
	long maxCount = 0;   // Size of The Largest Cluster
	int maxCluster = 0;	 // Index of the Largest Cluster
	CMeshO::PerVertexAttributeHandle<SatePtType> type_hi =
		vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<SatePtType>(mesh, _MySatePtAttri);
	for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi) {
		if ( !(vi->IsD()) && type_hi[vi] == Pt_Undefined) {
			// Seed Point
			nCluster++;
			long count = 1;
			type_hi[vi] = nCluster;
			std::vector<CMeshO::VertexIterator> dump;
			dump.push_back(vi);
			// Grow by KNN
			while (!dump.empty()) {
				CMeshO::VertexIterator curSeed = *(dump.end() - 1);
				dump.pop_back();
#if 1
				KDTree.doQueryK(curSeed->cP(), stepn, queue);
				int neighbours = queue.getNofElements();

				for (int k = 0; k < neighbours; k++) {
					CMeshO::VertexIterator nbor = mesh.vert.begin() + queue.getIndex(k);
					if (!(nbor->IsD()) &&
						type_hi[nbor] == Pt_Undefined &&
						Distance(curSeed->cP(), nbor->cP()) < dis) {
						type_hi[nbor] = nCluster;
						dump.push_back(nbor);
						count++;
					}
				}
#else // Query All at once
                std::vector<unsigned int> neiPtIdx;
                std::vector<float> neiPtSaureDis;
                KDTree.doQueryDist(curSeed->cP(), dis*dis, neiPtIdx, neiPtSaureDis);
                int neighbours = neiPtIdx.size();
                for (int k = 0; k < neighbours; k++) {
                    CMeshO::VertexIterator nbor = mesh.vert.begin() + neiPtIdx[k];
                    if ( !(nbor->IsD()) && type_hi[nbor] == Pt_Undefined ) {
                        type_hi[nbor] = nCluster;
                        dump.push_back(nbor);
                        count++;
                    }
                }
#endif
			}
			nums.push_back(count);
			if (count > maxCount) {
				maxCount = count;
				maxCluster = nCluster;
			}
		}
	}

	if (maxID != 0)
		*maxID = maxCluster;
	if (maxN != 0)
		*maxN = maxCount;

	return nCluster;
}
int PCFit::DeNoiseRegGrw()
{
	CMeshO &mesh = m_meshDoc.mesh->cm;

	QTime time;
	time.start();
	printf(
		"    [--DeNoiseRegGrw--]: #nPts-%d \n"
		"      | #NIteration : %d \n"
        "      | #KNNStep    : %d \n"
		"      | #DisRatio   : %.3f \n"
		"      | #Gap        : ",
		mesh.vn,
		DeNoise_MaxIteration, DeNoise_GrowNeighbors, DeNoise_DisRatioOfOutlier);
    //----[[
    const int nStep = DeNoise_GrowNeighbors;
	const int nIter = DeNoise_MaxIteration <= 0 ? 100 : DeNoise_MaxIteration;
	CMeshO::PerVertexAttributeHandle<SatePtType> type_hi = 
		vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<SatePtType>(mesh, _MySatePtAttri);
	for (int _iter = 0; _iter<nIter; _iter++) {

		// 1)Update Threshold
		vcg::tri::UpdateBounding<CMeshO>::Box(mesh);
		double DX = mesh.bbox.DimX();
		double DY = mesh.bbox.DimY();
		double DZ = mesh.bbox.DimZ();
		double Gap = DX<DY ? (DX<DZ ? DX : DZ) : (DY<DZ ? DY : DZ);
		Gap *= DeNoise_DisRatioOfOutlier;
		// 2)Region Growing
		std::vector<int> Nums;
		long maxN = 0;
		int maxCluster = 0;
		RegionGrow(mesh, nStep, Gap, Nums, &maxCluster, &maxN);
		for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi)
		{
			if ((type_hi[vi] & 0x0FF) > 0) {
				if (type_hi[vi] != maxCluster)
					vcg::tri::Allocator<CMeshO>::DeleteVertex(mesh, *vi);
				type_hi[vi] = Pt_Undefined;
			}
		}

		//----
		printf("[%d]-%.3f-<%d in %d> ", _iter+1, Gap, maxN, Nums.size());
		//----

		if (Nums.size() == 1)
			break;
	}

	int nNoise = 0;
	for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi) {
		if ((*vi).IsD()) {
			nNoise++;
			type_hi[vi] = Pt_Noise;
		}
	}
	vcg::tri::UpdateBounding<CMeshO>::Box(mesh);

	//----]]
	printf(
		"\n    [--DeNoiseRegGrw--]: Done, found %d outlier(s) in %.4f seconds.\n",
		nNoise, time.elapsed() / 1000.0);

	return nNoise;
}

double PCFit::GetMeshSizeAlongN(const vcg::Point3f n)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
	double minDis = DBL_MAX;
	double maxDis = -DBL_MAX;
	CMeshO::VertexIterator vi;
	int cnt = 0;
	for (vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi)
	{
		if ((*vi).IsD())
			continue;
		double proLoc = (*vi).cP().dot(n);
		if (proLoc > maxDis) maxDis = proLoc;
		if (proLoc < minDis) minDis = proLoc;
		cnt++;
	}
	return abs(maxDis - minDis);
}
bool PCFit::PCADimensionAna(
	vcg::Point3f &PSize,
	std::vector<vcg::Point3f> &PDirections,
	bool leftNoisePts)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;

	QTime time;
	time.start();

	const int N = leftNoisePts ? mesh.vn : mesh.vert.size();
	int row = 3, col = N;
	double *data = new double[3 * N];
	int Count = 0;
	for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi)
	{
		if (leftNoisePts && (*vi).IsD())
			continue;
		data[Count] = (*vi).cP().X();
		data[N + Count] = (*vi).cP().Y();
		data[2 * N + Count] = (*vi).cP().Z();
		++Count;
	}
	assert(Count == N);
	double *PC, V[3];//V[3] useless
	PC = new double[9];
	int ret = PCA(data, row, col, PC, V);
	delete[] data;
	if (ret == -1) {
		printf("    [--PCADirection--] PCA failed. [--PCADirection--]\n");
		doFailure();
	}

	vcg::Point3f NX = vcg::Point3f(PC[0], PC[3], PC[6]);
	vcg::Point3f NY = vcg::Point3f(PC[1], PC[4], PC[7]);
	vcg::Point3f NZ = vcg::Point3f(PC[2], PC[5], PC[8]);
	NX.Normalize();
	NY.Normalize();
	NZ.Normalize();

	delete[] PC;

	PSize.X() = GetMeshSizeAlongN(NX);
	PSize.Y() = GetMeshSizeAlongN(NY);
	PSize.Z() = GetMeshSizeAlongN(NZ);

	PDirections.clear();
	PDirections.push_back(NX);
	PDirections.push_back(NY);
	PDirections.push_back(NZ);

	printf(
		"    [--PCADirection--] #bPts - %d/%d \n"
		"       ( Principal direction vectors, col by decreasing order )\n"
		"       |  %7.4f %7.4f %7.4f - <%7.4f> \n"
		"       |  %7.4f %7.4f %7.4f - <%7.4f> \n"
		"       |  %7.4f %7.4f %7.4f - <%7.4f> \n"
		"    [--PCADirection--]: Done in %.4f seconds.\n",
		Count, N,
		NX.X(), NY.X(), NZ.X(), PSize.X(),
		NX.Y(), NY.Y(), NZ.Y(), PSize.Y(),
		NX.Z(), NY.Z(), NZ.Z(), PSize.Z(),
		time.elapsed() / 1000.0);

	return true;
}

double PCFit::RoughnessAna(bool leftNoisePts)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;

	QTime time;
	time.start();

	vcg::VertexConstDataWrapper<CMeshO> ww(mesh);
	vcg::KdTree<float> KDTree(ww);
	vcg::KdTree<float>::PriorityQueue queue;

    // Key Parameter
	const int knn = 20;
	const int a = 1.0;
	const int b = 1.0;
	const int N = leftNoisePts ? mesh.vn : mesh.vert.size();
	
	int Count = 0;
	double R = 0;
	for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi)
	{
		if (leftNoisePts && (*vi).IsD())
			continue;
		KDTree.doQueryK((*vi).cP(), knn, queue);
		int neighbours = queue.getNofElements();
		
		std::vector<vcg::Point3f> Pts;
		vcg::Plane3f ple;
		for (int i = 0; i < neighbours; i++) {
			int neightId = queue.getIndex(i);
			Pts.push_back(mesh.vert[neightId].cP());
		}
		vcg::FitPlaneToPointSet(Pts, ple);
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
	printf(
		"    [--RoughnessAna--] #bPts - %d/%d \n"
		"       | RoughnessA - <%7.4f> \n"
		"    [--RoughnessAna--]: Done in %.4f seconds.\n",
		Count, N, roughness, time.elapsed() / 1000.0);

	return roughness;
}


vcg::Point3f PCFit::GetPointList(
    std::vector<int> &indexList,
    std::vector<vcg::Point3f> &pointList,
    std::vector<vcg::Point3f> &normList,
    const bool bNormalize)
{
    CMeshO &mesh = m_meshDoc.mesh->cm;
    bool normalSupportted = m_meshDoc.mesh->hasDataMask(vcg::tri::io::Mask::IOM_VERTNORMAL);

    vcg::tri::UpdateBounding<CMeshO>::Box(mesh);
    vcg::Point3f center(0.0, 0.0, 0.0);
    if (bNormalize) {
        center = mesh.bbox.Center();
        printf(
            "\n"
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
    CMeshO::PerVertexAttributeHandle<SatePtType> type_hi =
        vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<SatePtType>(mesh, _MySatePtAttri);
    CMeshO::VertexIterator vi;
    int index;
    for (index = 0, vi = mesh.vert.begin(); vi != mesh.vert.end(); ++index, ++vi)
    {
        if (type_hi[vi] == Pt_Undefined && !(vi->IsD())) {
            indexList.push_back(index);
            pointList.push_back((*vi).cP() - center);
            if (normalSupportted)
                normList.push_back((*vi).cN());
        }
    }
    if (bNormalize)
        printf("    >> Got #%d Normalized Pts.\n", pointList.size());
    else
        printf("    >> Got #%d Pts.\n", pointList.size());
    return center;
}
