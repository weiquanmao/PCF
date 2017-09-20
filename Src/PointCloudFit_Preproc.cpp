#include "PointCloudFit.h"
#include "PCA.h"

#include <QTime>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/space/index/kdtree/kdtree.h>

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
		mesh.vn, DeNoise_MaxIteration, DeNoise_MaxNeighbors, DeNoise_DisRatioOfOutlier);
	//[[----

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
				KDTree.doQueryK((*vi).cP(), DeNoise_MaxNeighbors, queue);
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
	const double dis, const int stepn,
	std::vector<int> &nums,
	int *maxID, long *maxN)
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
			type_hi[vi] = Pt_Reserve + nCluster;
			std::vector<CMeshO::VertexIterator> dump;
			dump.push_back(vi);
			// Grow by KNN
			while (!dump.empty()) {
				CMeshO::VertexIterator curSeed = *(dump.end() - 1);
				dump.pop_back();
				KDTree.doQueryK(curSeed->cP(), stepn, queue);
				int neighbours = queue.getNofElements();

				for (int k = 0; k < neighbours; k++) {
					CMeshO::VertexIterator nbor = mesh.vert.begin() + queue.getIndex(k);
					if (!(nbor->IsD()) &&
						type_hi[nbor] == Pt_Undefined &&
						Distance(curSeed->cP(), nbor->cP()) < dis) {
						type_hi[nbor] = Pt_Reserve + nCluster;
						dump.push_back(nbor);
						count++;
					}
				}
			}
			nums.push_back(count);
			if (count > maxCount) {
				maxCount = count;
				maxCluster = Pt_Reserve + nCluster;
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
		"      | #NNeighbord : %d \n"
		"      | #DisRatio   : %.3f \n"
		"      | #Gap        : ",
		mesh.vn,
		DeNoise_MaxIteration, DeNoise_GrowNeighbors, DeNoise_DisRatioOfOutlier);
	//[[----

	const int _K = DeNoise_GrowNeighbors;
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
		RegionGrow(mesh, Gap, _K, Nums, &maxCluster, &maxN);
		for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi)
		{
			if (type_hi[vi] > Pt_Reserve) {
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