#include "PointCloudFit.h"
#include "PCA.h"
#ifdef _USE_OPENMP_
#include <omp.h>
#endif

#include <QTime>
#include <wrap/io_trimesh/io_mask.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/space/fitting3.h>

const vcg::Point4f _NON_PLANE = vcg::Point4f(0.0, 0.0, 0.0, 0.0);

// [1] Detect Planes
enum FixedAxis {
	FixedAxis_X = 0,
	FixedAxis_Y = 1,
	FixedAxis_Z = 2
};
int HoughPlaneOne(
	vcg::Point4f &Plane,
	const FixedAxis fix,
	const std::vector<vcg::Point3f> &PointList,
	const std::vector<vcg::Point3f> &NormList,
	const double _intercept,
	const double _a,
	const double _s)
{
	QTime time;
	time.start();


	Plane.SetZero();

	const bool bHasNorm = (NormList.size() > 0) ? true : false;
	if (bHasNorm)
		assert(PointList.size() == NormList.size());

	int index1, index2, index3;
	switch (fix)
	{
		case FixedAxis_Z: index1 = 2; index2 = 0; index3 = 1; break; //ZXY
		case FixedAxis_Y: index1 = 1; index2 = 2; index3 = 0; break; //YZX
		case FixedAxis_X:
		default:		  index1 = 0; index2 = 1; index3 = 2; break; //XYZ
	}

	const double scale1 = _s;
	const double scale2 = _a;
	int A = 2 * ceil(1 / scale1) + 1;
	int B = 2 * ceil(1 / scale1) + 1;
	int C = 2 * ceil(_intercept / scale2 + 1) + 1;
	unsigned int **facPlane = new unsigned int*[C];
	unsigned int *houghBuf = new unsigned int[A*B*C];
	memset(houghBuf, 0, sizeof(unsigned int)*A*B*C);
	for (int i = 0; i<C; i++)
		facPlane[i] = houghBuf + i*A*B;

	int count = 0;
	int a, b, c;
	double fa, fb, fc;
	double x, y, z;
	int num = PointList.size();
	auto iterP = PointList.begin();
	auto iterN = NormList.begin();
	for (iterP = PointList.begin(); iterP != PointList.end(); iterP++)
	{
		x = (*iterP).V(index1);
		y = (*iterP).V(index2);
		z = (*iterP).V(index3);

		double a_begin = -1.0;
		double a_end = 1.0;
		double b_begin = -1.0;
		double b_end = 1.0;
		if (bHasNorm)
		{
			double nx = (*iterN).V(index1);
			double ny = (*iterN).V(index2);
			double nz = (*iterN).V(index3);
			ny = ny / nx;
			nz = nz / nx;
			if (abs(ny) - 1.0>0.3 || abs(nz) - 1.0>0.3)
				continue;
			a_begin = ny - 0.3;
			a_end = ny + 0.3;
			b_begin = nz - 0.3;
			b_end = nz + 0.3;

			iterN++;
		}
		if (a_begin<-1) a_begin = -1;
		if (a_end>1) a_end = 1;
		if (b_begin<-1) b_begin = -1;
		if (b_end>1) b_end = 1;

		for (fa = a_begin; fa<a_end; fa += scale1)
		{
			for (fb = b_begin; fb<b_end; fb += scale1)
			{
				fc = -x - fa*y - fb*z;
				c = fc / scale2 + C / 2;
				if (c >= 0 && c<C)
				{
					a = (fa / scale1 + A / 2);
					b = (fb / scale1 + B / 2);
					(facPlane[c])[a*A + b] += 1;
					count++;
				}	
			}
		}
	}

	int Temp[4] = { 0, 0, 0, 0 };
	int maxVal = (facPlane[0])[0];

	unsigned *ptemp;
	for (c = 0; c<C; c++) {
		ptemp = facPlane[c];
		for (a = 0; a<A; a++) {
			for (b = 0; b<B; b++) {
				if (ptemp[b]>maxVal)
				{
					Temp[0] = a;
					Temp[1] = b;
					Temp[2] = c;
					Temp[3] = ptemp[b];
					maxVal = Temp[3];
				}
			}
			ptemp += A;
		}
	}
	delete[] houghBuf;
	delete[] facPlane;

	vcg::Point4f ret;
	double pa = (Temp[0] - A / 2)*scale1;
	double pb = (Temp[1] - B / 2)*scale1;
	double pc = (Temp[2] - C / 2)*scale2;
	int vote = Temp[3];

	switch (fix)
	{
		case FixedAxis_Z: Plane = vcg::Point4f(pa, pb, 1.0, pc); break; //ZXY
		case FixedAxis_Y: Plane = vcg::Point4f(pb, 1.0, pa, pc);; break; //YZX
		case FixedAxis_X:
		default:		  Plane = vcg::Point4f(1.0, pa, pb, pc);; break; //XYZ
	}

	printf(
		"      | [#Time-%7.4f]-[%d-Seted] : %7.4f %7.4f %7.4f %77.4f  #nPts-< %d >\n",
		time.elapsed() / 1000.0, fix+1, 
		Plane.X(), Plane.Y(), Plane.Z(), Plane.W(),
		vote);

	return vote;
}
int HoughPlane(
	vcg::Point4f &Plane,
	const std::vector<vcg::Point3f> &PointList,
	const std::vector<vcg::Point3f> &NormList,
	const double _intercept,
	const double _a,
	const double _s)
{//检测 + 拟合平面
	QTime time;
	time.start();
	printf(
		"    [--Plane_HT--]: #nPts-%d \n",
		PointList.size());

	vcg::Point4f P[3];
	int N[3];
#pragma omp parallel for
	for (int i = 0; i < 3; i++) {
		N[i] = HoughPlaneOne(P[i], FixedAxis(i), PointList, NormList, _intercept, _a, _s);
	}

	int retIdx = 0;
	if (N[0]>N[1])
		retIdx = (N[0]>N[2]) ? 0 : 2;
	else
		retIdx = (N[1]>N[2]) ? 1 : 2;
	
	printf(
		"      | [+][Checked] : [%d] #nPts-< %d > \n"
		"    [--Plane_HT--]: Done in %.4f seconds \n",
		retIdx + 1, N[retIdx], time.elapsed() / 1000.0);

	Plane = P[retIdx];
	return N[retIdx];
}

// [2] Surface Points Verification
inline vcg::Point4f P3TOP4(const vcg::Point3f P3) {
	vcg::Point4f P4;
	P4.X() = P3.X();
	P4.Y() = P3.Y();
	P4.Z() = P3.Z();
	P4.W() = 1.0;
	return P4;
}
std::vector<int> AttachToPlane(
	const std::vector<vcg::Point3f> &PointList,
	const std::vector<vcg::Point3f> &NormList,
	const vcg::Point4f &plane,
	const double _TDis, const double _TAng)
{//判断是否属于平面点
	QTime time;
	time.start();

	const bool bHasNorm = (NormList.size() > 0) ? true : false;
	if (bHasNorm)
		assert(PointList.size() == NormList.size());

	// Check Distance
	const double TDis = _TDis * sqrt(plane.V(0)*plane.V(0) + plane.V(1)*plane.V(1) + plane.V(2)*plane.V(2));
	std::vector<int> OnPlaneList;
	for (int i = 0; i<PointList.size(); i++) {
		vcg::Point4f pt = P3TOP4(PointList[i]);
		double d = abs(pt*plane);
		if (d<TDis)
			OnPlaneList.push_back(i);
	}

	// Check Norm
	if (bHasNorm && !OnPlaneList.empty())
	{
		vcg::Point3f NP = vcg::Point3f(plane.V(0), plane.V(1), plane.V(2));
		NP.Normalize();
		std::vector<int> ReReMoved;
		for (int i = 0; i<OnPlaneList.size(); ++i) {
			vcg::Point3f npt = NormList.at(OnPlaneList.at(i));
			npt.Normalize();
			double ang = 90 - abs(90 - vcg::AngleN(NP, npt)*_R2D);
			if (ang>_TAng)
				ReReMoved.push_back(i);
		}

		if (!ReReMoved.empty()) {
			for (int i = ReReMoved.size() - 1; i >= 0; --i)
				OnPlaneList.erase(OnPlaneList.begin() + ReReMoved.at(i));
		}
		ReReMoved.clear();
	}

	printf(
		"      [--Attach--]: Attach points to planes...\n"
		"        | #Threshold_Dis : %.4f\n"
		"        | #Threshold_Ang : %.4f-[%d]\n"
		"        | #nPts-OnPlane  : %d \n"
		"      [--Attach--]: Done in %.4f seconds. \n",
		TDis, _TAng, bHasNorm,
		OnPlaneList.size(), time.elapsed() / 1000.0);

	return OnPlaneList;
}

// [3] Coplanar Separation
void PicMaxRegion(
	const std::vector<vcg::Point3f> &PointList,
	std::vector<int> &index,
	const double _2TDis)
{
	QTime time;
	time.start();

	if (PointList.empty())
		return;

	std::vector<int> *tempLoc = new std::vector<int>;
	std::vector<int> *LocFinal = new std::vector<int>;

	std::vector<vcg::Point3f> seedStack;
	vcg::Point3f seed;
	std::vector<int> recorder;

	while (!index.empty())
	{
		//每一类
		seedStack.clear();
		seedStack.push_back(PointList.at(index.at(0)));
		tempLoc->clear();
		tempLoc->push_back(index.at(0));

		index.erase(index.begin());
		//区域生长
		while (!seedStack.empty())
		{
			seed = seedStack.at(0);
			seedStack.erase(seedStack.begin());
			recorder.clear();
			for (int i = 0; i<index.size(); i++)
			{
				int p = index.at(i);
				if (vcg::Distance(seed, PointList.at(p))<_2TDis)
				{
					seedStack.push_back(PointList.at(p));
					tempLoc->push_back(p);
					recorder.push_back(i);
				}
			}
			for (int i = recorder.size() - 1; i >= 0; i--)
				index.erase(index.begin() + recorder.at(i));
		}
		if (tempLoc->size() > LocFinal->size())
		{
			std::vector<int> *temp;
			temp = LocFinal;
			LocFinal = tempLoc;
			tempLoc = temp;
		}
	}
	for (auto iter : *LocFinal)
		index.push_back(iter);
	LocFinal->clear();
	tempLoc->clear();
	delete LocFinal;
	delete tempLoc;
	recorder.clear();
	std::sort(index.begin(), index.end());

	printf(
		"      [--MaxRegion--]: #Pts-%d\n"
		"        | #Threshold_Dis : %.4f\n"
		"      [--MaxRegion--]: Done in %.4f seconds. \n",
		PointList.size(), _2TDis, index.size(), time.elapsed() / 1000.0);
}

// [4] LS Fit
float EvalPlane(vcg::Plane3f &pl, std::vector<vcg::Point3f> posVec)
{
	float off = 0;
	for (size_t i = 0; i < posVec.size(); ++i) {
		float d = vcg::SignedDistancePlanePoint(pl, posVec[i]);
		off += d*d;
	}

	off = sqrt(off / float(posVec.size()));
	return off;
}
double FineFit(
	const std::vector<vcg::Point3f> &pointList,
	const std::vector<int> &planeVerList,
	vcg::Point4f &plane)
{
	QTime time;
	time.start();

	std::vector<vcg::Point3f> ExactVec;
	std::vector<float> WeightVec;
	vcg::Plane3f ple;
	for (int j = 0; j<planeVerList.size(); ++j)
	{
		vcg::Point3f p = pointList.at(planeVerList.at(j));
		ExactVec.push_back(p);
	}

	vcg::FitPlaneToPointSet(ExactVec, ple);
	plane.X() = ple.Direction().X();
	plane.Y() = ple.Direction().Y();
	plane.Z() = ple.Direction().Z();
	plane.W() = -ple.Offset();
	double fitError = EvalPlane(ple, ExactVec);

	printf(
		"      [--Fit_LS--]: #Pts-%d\n"
		"        | #Plane    : < %.4f, %.4f, %.4f, %.4f> \n"
		"        | #FitError : %.4f\n"
		"      [--Fit_LS--]: Done in %.4f seconds. \n",
		planeVerList.size(),
		plane.X(), plane.Y(), plane.Z(), plane.W(),
		fitError, time.elapsed() / 1000.0);

	return fitError;
}

// [5] Get Minimum-Bounding-Rectangle
double PatchDimensions(
	const std::vector<vcg::Point3f> & pts,
	const vcg::Point3f &nx,
	const vcg::Point3f &ny,
	vcg::Point3f &Pt_O, vcg::Point3f &Pt_Dx, vcg::Point3f &Pt_Dy,
	double &ex, double &ey
)
{
	const int NPP = pts.size();

	std::vector<std::pair<double, int>> lX;
	std::vector<std::pair<double, int>> lY;
	for (int i = 0; i<NPP; i++) {
		vcg::Point3f pp = pts.at(i);
		double px = pp*nx;
		double py = pp*ny;
		lX.push_back(std::pair<double, int>(px, i));
		lY.push_back(std::pair<double, int>(py, i));
	}
	std::sort(lX.begin(), lX.end());
	std::sort(lY.begin(), lY.end());

	vcg::Point3f MinX = pts.at(lX.at(0).second);
	vcg::Point3f MinY = pts.at(lY.at(0).second);
	double lx = lX.at(NPP - 1).first - lX.at(0).first;
	double ly = lY.at(NPP - 1).first - lY.at(0).first;
	Pt_Dx = nx*lx;
	Pt_Dy = ny*ly;
	double MinY2X = MinY*nx;
	double MinX2Y = MinX*ny;
	vcg::Point3f O1 = MinY - nx*(MinY2X - lX.at(0).first);
	vcg::Point3f O2 = MinX - ny*(MinX2Y - lY.at(0).first);
	Pt_O = (O1 + O2) / 2.0;
	ex = EIConfidence(lX);
	ey = EIConfidence(lY);

	return lx*ly;
}
void GetMBR(
	CMeshO &mesh,
	const vcg::Point4f &Plane,
	Sailboard *ABord,
	const std::vector<int> &PlaneVerList,	
	std::vector<int> &IndexList,
	std::vector<vcg::Point3f> &PointList,
	std::vector<vcg::Point3f> &NormList
)
{//提出数据+分析平面
	QTime time;
	time.start();

	if (ABord == 0 ||
		Plane == _NON_PLANE ||
		PlaneVerList.empty())
		return;

	assert(IndexList.size() == PointList.size());
	const bool bHasNorm = (NormList.size() > 0) ? true : false;
	if (bHasNorm)
		assert(IndexList.size() == NormList.size());


	vcg::Point3f NP = vcg::Point3f(Plane.V(0), Plane.V(1), Plane.V(2));
	if (Plane.V(3)>0) NP = -NP;
	double offset = abs(Plane.V(3)) / NP.Norm();
	NP.Normalize();
	ABord->m_N = NP;


	vcg::Plane3f VCGPLane = vcg::Plane3f();
	VCGPLane.SetDirection(NP);
	VCGPLane.SetOffset(offset);
	CMeshO::PerVertexAttributeHandle<SatePtType> type_hi = 
		vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<SatePtType>(mesh, _MySatePtAttri);
	std::vector<vcg::Point3f> OnPPt;

	const int PlaneCode = ABord->m_PlaneIndex;
	for (int i = PlaneVerList.size() - 1; i >= 0; --i)
	{
		int index = PlaneVerList.at(i);
		type_hi[IndexList.at(index)] = PlaneCode;
		IndexList.erase(IndexList.begin() + index);

		vcg::Point3f pp = VCGPLane.Projection(PointList.at(index));

		OnPPt.push_back(pp);
		PointList.erase(PointList.begin() + index);
	}
	if (bHasNorm)
	{
		for (int i = PlaneVerList.size() - 1; i >= 0; --i) {
			int index = PlaneVerList.at(i);
			NormList.erase(NormList.begin() + index);
		}
	}

	if (OnPPt.empty())
		return;

	const int N = OnPPt.size();
	double *data = new double[3 * N];
	for (int i = 0; i<N; ++i)
	{
		data[i] = OnPPt.at(i).X();
		data[N + i] = OnPPt.at(i).Y();
		data[2 * N + i] = OnPPt.at(i).Z();
	}
	int row = 3, col = N;
	double *PC, V[3];//V[3]无用
	PC = new double[9];
	int ret = PCA(data, row, col, PC, V);
	delete[] data;
	if (ret == -1)
		return;

	vcg::Point3f NX_0 = vcg::Point3f(PC[0], PC[3], PC[6]);
	vcg::Point3f NY_0 = vcg::Point3f(PC[1], PC[4], PC[7]);
	NY_0 = NP^NX_0;
	NX_0 = NY_0^NP;
	NX_0.Normalize();
	NY_0.Normalize();
	delete[] PC;

	vcg::Point3f O, DX, DY;
	double confidenceX, confidenceY;
	double SA = PatchDimensions(OnPPt, NX_0, NY_0, O, DX, DY, confidenceX, confidenceY);
	for (int i = -29; i < 30; i += 2) { // 并没想想中的慢,挺快的
		double k = tan(i *_D2R);
		vcg::Point3f nx = NX_0 + NY_0*k;
		nx.Normalize();
		vcg::Point3f ny = NP^nx;
		vcg::Point3f o, dx, dy;
		double ex, ey;
		double sa = PatchDimensions(OnPPt, nx, ny, o, dx, dy, ex, ey);
		if (sa < SA) {
			SA = sa;
			O = o; DX = dx; DY = dy;
			confidenceX = ex; confidenceY = ey;
		}
	}

	ABord->m_pO = O;
	ABord->m_dX = DX;
	ABord->m_dY = DY;
	ABord->m_sizeConfidence = vcg::Point3f(confidenceX, confidenceY, 0);
	//
	printf(
		"      [--PatchDim--]: #Pts-%d\n"
		"        | #Loc-PCA_X   : < %7.4f, %7.4f, %7.4f> \n"
		"        | #Loc-PCA_Y   : < %7.4f, %7.4f, %7.4f> \n"
		"        | #Loc-Adjed_O : < %7.4f, %7.4f, %7.4f> \n"
		"        | #Loc-Adjed_X : < %7.4f, %7.4f, %7.4f> - [%7.4f]\n"
		"        | #Loc-Adjed_Y : < %7.4f, %7.4f, %7.4f> - [%7.4f]\n"
		"      [--PatchDim--]: Done in %.4f seconds. \n",
		OnPPt.size(),
		NX_0.X(), NX_0.Y(), NX_0.Z(), NY_0.X(), NY_0.Y(), NY_0.Z(),
		O.X(), O.Y(), O.Z(),
		DX.X(), DX.Y(), DX.Z(), confidenceX,
		DY.X(), DY.Y(), DY.Z(), confidenceY,
		time.elapsed() / 1000.0);

	OnPPt.clear();
}



std::vector<Sailboard*> PCFit::DetectPlanes(const int expPN)
{
	// -- Move
	std::vector<Sailboard*> sailbords;
	CMeshO &mesh = m_meshDoc.mesh->cm;
	bool normalSupportted = m_meshDoc.mesh->hasDataMask(vcg::tri::io::Mask::IOM_VERTNORMAL);

	vcg::tri::UpdateBounding<CMeshO>::Box(mesh);
	vcg::Point3f center = mesh.bbox.Center();

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


	// -- Get Point List And Normal List (If Exist)
	std::vector<int> indexList;
	std::vector<vcg::Point3f> pointList;
	std::vector<vcg::Point3f> directionList;
	indexList.reserve(mesh.vn);
	pointList.reserve(mesh.vn);
	if (normalSupportted)
		directionList.reserve(mesh.vn);
	CMeshO::PerVertexAttributeHandle<SatePtType> type_hi =
		vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<SatePtType>(mesh, _MySatePtAttri);
	CMeshO::VertexIterator vi;
	int index;
	for (index = 0, vi = mesh.vert.begin(); vi != mesh.vert.end(); ++index, ++vi)
	{
		if (type_hi[vi] == Pt_Undefined) {
			indexList.push_back(index);
			pointList.push_back((*vi).cP() - center);
			if (normalSupportted)
				directionList.push_back((*vi).cN());
		}
	}
	printf("    >> Got and moved #%d Pts.\n", pointList.size());

	// -- Calculate intercept
	double candicate1 = mesh.bbox.Diag() * sqrt(3.0) / 2.0;
	double candicate2 = mesh.bbox.DimX() + mesh.bbox.DimY() + mesh.bbox.DimZ();
	double intercept = candicate1 < candicate2 ? candicate1 : candicate2;

	const double _planeDisThreshold = m_refa*Threshold_DisToPlane;
	const double _planeAngThreshold = m_refa*Threshold_AngToPlane;
	const int _planeNThreshold = pointList.size()*Threshold_NPtsPlane;
	int planeNum = 0;
	for (planeNum = 0; planeNum<expPN; planeNum++)
	{
		// HT Detection
		vcg::Point4f plane;
#if 0
		double _planeNT = _planeNThreshold;
#else
		double _planeNT = pointList.size()*Threshold_NPtsPlane;
#endif
		if (normalSupportted)
			_planeNT *= 0.5;
		printf("    >> Detecting the [ No.%d ] plane with #Minimum - [ %d ] ...\n", planeNum + 1, _planeNT);
		int NP = HoughPlane(plane, pointList, directionList, intercept, m_refa, Precision_HT); // Center At (0,0,0)
		if (plane == _NON_PLANE || NP <= _planeNT)
			break;

		// Surface Points Verification
		std::vector<int> planeVerList = AttachToPlane(pointList, directionList, plane, _planeDisThreshold, _planeAngThreshold);

		// Coplanar Separation
		PicMaxRegion(pointList, planeVerList, _planeDisThreshold*2.0);
		if (planeVerList.size() <= _planeNT)
			break;

		// Least Squre Fit
		Sailboard *oneBord = new Sailboard(Pt_OnPlane + planeNum);
		double err = FineFit(pointList, planeVerList, plane);
		// planeVerList = AttachToPlane(pointList,directionList,plane,_planeDThreshold); // AGAIN!
		oneBord->m_varN = err;

		// Get Minimum-Bounding-Rectangle
		GetMBR(mesh, plane, oneBord, planeVerList, indexList, pointList, directionList);

		sailbords.push_back(oneBord);
		planeVerList.clear();
	}
	indexList.clear();
	pointList.clear();
	directionList.clear();

	for (int i = 0; i<sailbords.size(); ++i) {
		sailbords.at(i)->m_pO += center;
	}

	return sailbords;
}