#include "PointCloudFit.h"
#include "PCA.h"
#include <fstream>
#include <vcg/complex/algorithms/update/bounding.h>

bool DetectSysAxis(
	const std::vector<vcg::Point3f> &PointList,
	const std::vector<vcg::Point3f> &NormList,
	const double _a,
	vcg::Point3f &PO, vcg::Point3f &NN)
{
	const int PNum = PointList.size();
	assert(NormList.size() == PNum);


	vcg::Box3f box;
	box.Set(PointList.at(0));
	for (int i = 1; i<PNum; i++)
		box.Add(PointList.at(i));
	const vcg::Point3f lo = box.min;
	const vcg::Point3f hi = box.max;

	const int A = ceil(box.DimX() / _a);
	const int B = ceil(box.DimY() / _a);
	const int C = ceil(box.DimZ() / _a);

	int **facPlane = new int*[C];
	int *houghBuf = new int[A*B*C];
	memset(houghBuf, 0, sizeof(int)*A*B*C);
	for (int i = 0; i<C; i++)
		facPlane[i] = houghBuf + i*A*B;

	//霍夫变换
	int X, Y, Z;
	double x, y, z;
	double x0, y0, z0;
	double nx, ny, nz;
	for (int i = 0; i<PNum; i++)
	{
		vcg::Point3f pp = PointList.at(i);
		vcg::Point3f nn = NormList.at(i);
		x0 = pp.X();
		y0 = pp.Y();
		z0 = pp.Z();
		nx = nn.X();
		ny = nn.Y();
		nz = nn.Z();

		for (X = 0; X<A; X++)
		{
			x = lo.X() + X*_a;
			y = (x - x0) / nx*ny + y0;
			z = (x - x0) / nx*nz + z0;
			Y = int((y - lo.Y()) / _a + 0.5);
			Z = int((z - lo.Z()) / _a + 0.5);
			if (Y >= 0 && Y<B && Z >= 0 && Z<C)
				(facPlane[Z])[Y*A + X] += 1;
		}
	}

	//寻找最大值
	int max = (facPlane[0])[0];
	for (int c = 0; c<C; c++)
	{
		int *ptemp = facPlane[c];
		for (int b = 0; b<B; b++)
		{
			for (int a = 0; a<A; a++)
				if (ptemp[a]>max)
					max = ptemp[a];
			ptemp += A;
		}
	}
	int cut_value = 5;

	std::ofstream out("a.txt", std::ios::out);
	std::vector<double> vX, vY, vZ;
	std::vector<int> weight;
	int count = 0;

	double ux, uy, uz;
	ux = uy = uz = 0;

	for (int c = 0; c<C; c++)
	{
		int *ptemp = facPlane[c];
		for (int b = 0; b<B; b++)
		{
			for (int a = 0; a<A; a++)
			{
				if (ptemp[a]>cut_value)
				{
					x = lo.X() + a*_a;	y = lo.Y() + b*_a;	z = lo.Z() + c*_a;
					vX.push_back(x);	vY.push_back(y);	vZ.push_back(z);
					ux += x*ptemp[a];	uy += y*ptemp[a];	uz += z*ptemp[a];
					weight.push_back(ptemp[a]);
					count += ptemp[a];

					out << x << "\t" << y << "\t" << z << "\t" << ptemp[a] << std::endl;/////////
				}
			}
			ptemp += A;
		}
	}
	delete[] houghBuf;

	PO.X() = ux / (count*1.0);
	PO.Y() = uy / (count*1.0);
	PO.Z() = uz / (count*1.0);

	int row = 3;
	int col = count;
	double *data = new double[3 * col];
	double *pX = data;
	double *pY = data + col;
	double *pZ = data + 2 * col;
	int index = 0;
	for (int i = 0; i<vX.size(); i++)
	{
		for (int j = 0; j<weight[i]; j++)
		{
			pX[index] = vX[i];
			pY[index] = vY[i];
			pZ[index] = vZ[i];

			index++;
		}
	}

	vX.clear();
	vY.clear();
	vZ.clear();

	double *PC, V[3];//V[3]无用
	PC = new double[9];
	int ret;
	ret = PCA(data, row, col, PC, V);
	delete[] data;
	if (ret == -1)
		return false;

	NN.X() = PC[0];
	NN.Y() = PC[3];
	NN.Z() = PC[6];
	delete[] PC;

	return true;
}
typedef struct {
	int index;
	double r;
	double l;
} CylinderPointInf;
bool CPIByL(CylinderPointInf a, CylinderPointInf b)
{
	return a.l<b.l;
}
double Clus3(const std::vector<double> &data)
{
	int n_begin = data.size()*0.2;
	int n_end = data.size()*0.8;
	double u1, u2, u3;
	double v1, v2, v3;
	double n1, n2, n3;
	double V = *(data.end() - 1) - *data.begin();
	V = V*V*data.size();
	double tempV;
	double U;
	if (n_end - n_begin<4)
	{
		U = 0;
		for (int k = n_begin; k <= n_end; k++)
			U += data[k];
		U = U / (1.0*(n_end - n_begin + 1));
		return U;
	}
	//cout<<n_begin<<"\t"<<n_end<<endl;
	for (int i = n_begin + 1; i<n_end - 1; i++)
	{
		for (int j = i; j<n_end - 1; j++)
		{
			//cout<<i<<"\t"<<j<<endl;
			u1 = u2 = u3 = 0;
			n1 = n2 = n3 = 0;
			for (int k = n_begin; k<i; k++)
			{
				u1 += data[k];
				n1++;
			}
			for (int k = i; k <= j; k++)
			{
				u2 += data[k];
				n2++;
			}
			for (int k = j + 1; k <= n_end; k++)
			{
				u3 += data[k];
				n3++;
			}
			u1 = u1 / (n1*1.0);
			u2 = u2 / (n2*1.0);
			u3 = u3 / (n3*1.0);

			v1 = v2 = v3 = 0;
			for (int k = n_begin; k<i; k++)
			{
				v1 = (data[k] - u1)*(data[k] - u1);
			}
			for (int k = i; k <= j; k++)
			{
				v2 = (data[k] - u2)*(data[k] - u2);
			}
			for (int k = j + 1; k <= n_end; k++)
			{
				v3 = (data[k] - u3)*(data[k] - u3);
			}

			tempV = v1 + v2 + v3;
			if (tempV<V)
			{
				U = u2;
				V = tempV;
			}
		}
	}

	return U;
}
double EstRadious(const std::vector<double> &RList, double &weight)
{
	const double cutRatio = 0.2;
	double sum = 0.0;
	double sum2 = 0.0;
	int N = 0;
	for (int i = RList.size()*cutRatio; i<RList.size()*(1.0 - 2 * cutRatio); ++i) {
		double r = RList.at(i);
		sum += r;
		sum2 += r*r;
		N++;
	}
	double meanR = sum / (N*1.0);
	double var = sqrt(sum2 / (N*1.0) - meanR*meanR);

	weight = 1.0 - var / meanR;
	return meanR;
}
std::vector<int> AttachToCylinder(
	const std::vector<vcg::Point3f> &PointList,
	const std::vector<vcg::Point3f> &NormList,
	const vcg::Point3f &PO,
	const vcg::Point3f &NN,
	const double _a, const double RT,
	double &R, double &L,
	double &Rw, double &Lw,
	vcg::Point3f &NewPO)
{
	std::ofstream out("R.txt", std::ios::out);

	bool bRT = RT > 0 ? true : false;
	double distanceThreshold = bRT ? RT * 0.5 : _a * 20;

	// 1. 初步筛选表面点（记录索引、半径和投影位置）
	std::vector<CylinderPointInf> CPI;
	double MeanR = 0.0;
	for (int i = 0; i<PointList.size(); i++)
	{
		vcg::Point3f pp = PointList.at(i);
		vcg::Point3f np = NormList.at(i);

		vcg::Point3f nn = NN^np;
		nn.Normalize();
		double distance = abs(nn*(pp - PO));
		double intercept = NN*(pp - PO);
		double radius = sqrt((pp - PO).SquaredNorm() - intercept*intercept);
		if (distance<distanceThreshold &&
			(!bRT || (radius<RT*1.5 && radius>RT*0.5)))
		{
			CylinderPointInf temp;
			temp.index = i;
			temp.l = intercept;
			temp.r = radius;
			CPI.push_back(temp);
			MeanR += radius;
		}
	}
	MeanR = MeanR / CPI.size();

	// 2. 沿轴向排序
	std::sort(CPI.begin(), CPI.end(), CPIByL);

	// 3. 取最大段（剔除两端异常点）
	double maxGap = _a*0.5;
	int n_begin;
	int n_end;
	double maxLength = 0;
	for (int i = 0; i<CPI.size();)
	{
		int j = i;
		for (; j<CPI.size() - 1; j++)
		{
			if (CPI[j + 1].l - CPI[j].l>maxGap)
				break;
		}
		if (CPI[j].l - CPI[i].l>maxLength)
		{
			n_begin = i;
			n_end = j;
			maxLength = CPI[j].l - CPI[i].l;
		}
		i = j + 1;
	}

	// 4. 分层验证（二次验证）
	double lstep = _a*0.5;
	int index_begin, index_end;
	index_begin = index_end = 0;
	double l_begin = CPI[0].l;
	double l_end = l_begin + lstep;

	std::vector<int> FOnCylList;
	std::vector<double> FmeanR_temp;
	std::vector<int> FnCount;
	double ML = 0;
	double FminL, FmaxL;

	std::vector<int> OnCylList;
	std::vector<double> meanR_temp;
	std::vector<int> nCount;
	double minL, maxL;
	minL = CPI[CPI.size() - 1].l;
	maxL = CPI[0].l;
	for (int i = n_begin; i <= n_end; i++)
	{
		if ((CPI[i].l) > l_end || i == CPI.size() - 1)
		{
			std::vector<double> R_temp;
			for (int j = index_begin; j <= index_end; j++)
				R_temp.push_back(CPI[j].r);

			std::sort(R_temp.begin(), R_temp.end());
			double uR = Clus3(R_temp);
			out << l_begin << "\t" << uR << std::endl;
			if (uR<1.5*MeanR && uR>0.8*MeanR) {

				meanR_temp.push_back(uR);
				nCount.push_back(R_temp.size());
				R_temp.clear();

				for (int j = index_begin; j <= index_end; j++)
				{
					if (CPI[j].r < 1.2 * uR &&
						CPI[j].r > 0.8 * uR) {
						OnCylList.push_back(CPI[j].index);
						if (CPI[j].l < minL) minL = CPI[j].l;
						if (CPI[j].l > maxL) maxL = CPI[j].l;
					}
				}
			}
			else if (maxL - minL > ML) {
				FOnCylList = OnCylList;
				FmeanR_temp = meanR_temp;
				FnCount = nCount;
				ML = maxL - minL;
				FmaxL = maxL;
				FminL = minL;
				minL = maxL;
				maxL = minL - 1.0;
				OnCylList.clear();
				meanR_temp.clear();
				nCount.clear();
			}
			//---------------------
			index_begin = index_end = i;
			l_begin = CPI[i].l;
			l_end = l_begin + lstep;
		}
		else
			index_end = i;
	}

	Lw = EIConfidence(FnCount);
	R = EstRadious(FmeanR_temp, Rw);
	L = FmaxL - FminL;
	NewPO = PO + NN*(FmaxL + FminL) / 2.0;


	meanR_temp.clear();
	CPI.clear();
	nCount.clear();


	return FOnCylList;
}


MainBodyCylinder* PCFit::DetectCylinder()
{
	CMeshO &mesh = m_meshDoc.mesh->cm;

	vcg::tri::UpdateBounding<CMeshO>::Box(mesh);
	vcg::Point3f center = mesh.bbox.Center();

	// -- Get Point List And Normal List
	std::vector<int> indexList;
	std::vector<vcg::Point3f> pointList;
	std::vector<vcg::Point3f> directionList;
	indexList.reserve(mesh.vn);
	pointList.reserve(mesh.vn);
	directionList.reserve(mesh.vn);
	CMeshO::PerVertexAttributeHandle<SatePtType> type_hi = vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<SatePtType>(mesh, _MySatePtAttri);
	CMeshO::VertexIterator vi;
	int index;
	for (index = 0, vi = mesh.vert.begin(); vi != mesh.vert.end(); ++index, ++vi)
	{
		if (type_hi[vi] == Pt_Undefined) {
			indexList.push_back(index);
			pointList.push_back((*vi).cP() - center);
			directionList.push_back((*vi).cN());
		}
	}


	//

	// 1.估计对称轴
	vcg::Point3f PO, N;
	if (!DetectSysAxis(pointList, directionList, m_refa, PO, N)) {
		indexList.clear();
		pointList.clear();
		directionList.clear();
		return 0;
	}


	// 2.主要是估个半径+初步确定表面点
	double R, L;
	double RWeight, LWeight;
	vcg::Point3f UpdatePO;
	std::vector<int> OnClyList = AttachToCylinder(pointList, directionList, PO, N.Normalize(), m_refa, 0, R, L, RWeight, LWeight, UpdatePO);
	if (OnClyList.size() < pointList.size()*Threshold_NPtsCyl) {
		indexList.clear();
		pointList.clear();
		directionList.clear();
		return 0;
	}
	PO = UpdatePO;

	for (int iter = 0; iter < 0; iter++) {
		vcg::Point3f PO_old = PO;
		vcg::Point3f N_old = N;
		std::vector<vcg::Point3f> pointListChecked;
		std::vector<vcg::Point3f> directionListChecked;
		for (int i = 0; i<OnClyList.size(); i++) {
			int index = OnClyList.at(i);
			pointListChecked.push_back(pointList.at(index));
			directionListChecked.push_back(directionList.at(index));
		}
		bool bAxis = DetectSysAxis(pointListChecked, directionListChecked, m_refa, PO, N);
		pointListChecked.clear();
		directionListChecked.clear();
		if (!bAxis)
			break;

		// 4.再次确定点
		double RT = R;
		OnClyList = AttachToCylinder(pointList, directionList, PO, N, m_refa, RT, R, L, RWeight, LWeight, UpdatePO);

		PO = UpdatePO;
		if (90 - abs(90 - vcg::Angle(N, N_old)*_R2D) < 2)
			break;
		// 5.用表面点估计轴向
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
		//double *PC, V[3];//V[3]无用
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

	if (OnClyList.size() < pointList.size()*Threshold_NPtsCyl) {
		indexList.clear();
		pointList.clear();
		directionList.clear();
		return 0;
	}

	for (int i = 0; i<OnClyList.size(); i++)
	{
		int index = indexList.at(OnClyList.at(i));
		type_hi[index] = Pt_OnMBCylinder;
	}

	indexList.clear();
	pointList.clear();
	directionList.clear();

	MainBodyCylinder *MainCyl = new MainBodyCylinder();
	MainCyl->m_length = L;
	MainCyl->m_radius = R;
	MainCyl->m_length_weight = LWeight;
	MainCyl->m_radius_weight = RWeight;
	MainCyl->m_pO = PO + center;
	MainCyl->m_N = N;

	return MainCyl;
}
