#include "PointCloudFitUtil.h"
#include "gte/Mathematics/GteApprCylinder3.h"
#include "tool/PCA.h"

#include <random>
#include <list>


bool DetectSymAxis(
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    vcg::Point3f &PO, vcg::Point3f &NN,
    const double _a)
{
    const int PNum = PointList.size();
    assert(NormList.size() == PNum);

    QTime time;
    time.start();


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

    // -- Voting 
    for (int i = 0; i<PNum; i++)
    {
        vcg::Point3f pp = PointList.at(i);
        vcg::Point3f nn = NormList.at(i);
        double px = pp.X();
        double py = pp.Y();
        double pz = pp.Z();
        double nx = nn.X();
        double ny = nn.Y();
        double nz = nn.Z();
#ifdef _USE_OPENMP_
#pragma omp parallel for
#endif // !_USE_OPENMP_
        for (int X = 0; X<A; X++)
        {
            double x = lo.X() + X*_a - px;
            double y = x / nx*ny + py;
            double z = x / nx*nz + pz;
            int Y = int((y - lo.Y()) / _a + 0.5);
            int Z = int((z - lo.Z()) / _a + 0.5);
            if (Y >= 0 && Y<B && Z >= 0 && Z<C)
                (facPlane[Z])[Y*A + X] += 1;
        }
    }

    // -- Find Max Vote
    int maxVote = (facPlane[0])[0];
    for (int c = 0; c<C; c++) {
        int *ptemp = facPlane[c];
        for (int b = 0; b<B; b++) {
            for (int a = 0; a<A; a++)
                if (ptemp[a]>maxVote)
                    maxVote = ptemp[a];
            ptemp += A;
        }
    }
#if 1
    const int TVote = maxVote * 0.5;
#else
    const int TVote = 5;
#endif


    
    // -- Pick Candicate Points (& Calculate Centroid)
    std::vector<double> vX, vY, vZ;
    std::vector<int> weight;
    double ux = 0.0;
    double uy = 0.0;
    double uz = 0.0;
    int count = 0;
#if defined(_ReportOut_)
    std::ofstream _out("../~Axis~.txt", std::ios::out);
#endif
    for (int c = 0; c<C; c++) {
        int *ptemp = facPlane[c];
        for (int b = 0; b<B; b++) {
            for (int a = 0; a<A; a++) {
                if (ptemp[a]>TVote)
                {
                    double x = lo.X() + a*_a;  vX.push_back(x);  ux += x*ptemp[a];
                    double y = lo.Y() + b*_a;  vY.push_back(y);  uy += y*ptemp[a];
                    double z = lo.Z() + c*_a;  vZ.push_back(z);  uz += z*ptemp[a];                  		
                    weight.push_back(ptemp[a]);
                    count += ptemp[a];
#if defined(_ReportOut_)
                    _out << x << "\t" << y << "\t" << z << "\t" << ptemp[a] << std::endl;
#endif                    
                }
            }
            ptemp += A;
        }
    }
#if defined(_ReportOut_)
    _out.close();
#endif
    delete[] houghBuf;


    // -- [Assign Location]
    PO.X() = ux / (count*1.0);
    PO.Y() = uy / (count*1.0);
    PO.Z() = uz / (count*1.0);

    // -- Estimate Norm
    int row = 3;
    int col = count;
    double *data = new double[3 * col];
    double *pX = data;
    double *pY = data + col;
    double *pZ = data + 2 * col;
    int index = 0;
    for (int i = 0; i<vX.size(); i++)
    {
        for (int j = 0; j<weight[i]; j++, index++)
        {
            pX[index] = vX[i];
            pY[index] = vY[i];
            pZ[index] = vZ[i];
        }
    }

    vX.clear();
    vY.clear();
    vZ.clear();

    double PC[9], V[3];//V[3] IS USELESS
    int ret = PCA(data, row, col, PC, V);
    delete[] data;
    if (ret == -1) {
        flog("    [--DetectSymAxis--] PCA failed. [--DetectSymAxis--]\n");
        doFailure();
        return false;
    }

    // -- [Assign Norm]
    NN.X() = PC[0];
    NN.Y() = PC[3];
    NN.Z() = PC[6];

    flog(
        "      [--DetectSymAxis--]: #nPts-%d...\n"
        "        | #RefA : %.4f\n"
        "        | #TVote: %d vs [ %d ]\n"
        "        | #Axis : < %7.4f, %7.4f, %7.4f > \n"
        "        | #Loc  : < %7.4f, %7.4f, %7.4f > \n"
        "      [--DetectSymAxis--]: Done in %.4f seconds. \n",
        PNum, _a, TVote, maxVote,
        PO.X(), PO.Y(), PO.Z(),
        NN.X(), NN.Y(), NN.Z(),
        time.elapsed() / 1000.0);

    return true;
}


bool CmpCPIByL(const CylinderPointInfo &a, const CylinderPointInfo &b)
{
    return a.l<b.l;
}
template<class T>
void StdVar(const std::vector<T> &data, const int begin, const int end, double &miu, double &var)
{
    assert(begin >= 0 && begin < data.size());
    assert(end >= 0 && end < data.size());
    assert(begin <= end);


    double u = 0.0;
    double v = 0.0;
    int n = (end - begin + 1);
    for (int k = begin; k <= end; k++) {
        T val = data[k];
        u += val;
        v += val*val;
    }
    u = u / n;
    v = v / n - u*u;
    miu = u;
    var = v;
}
template<class T>
double RobustMean(const std::vector<T> &data)
{
    std::vector<T> orderData = data;
    std::sort(orderData.begin(), orderData.end());

    int n_begin = orderData.size()*0.2;
    int n_end = orderData.size()*0.8;
    double u1, u2, u3;
    double v1, v2, v3;

    double MaxVar = *(orderData.end() - 1) - *orderData.begin();
    MaxVar = MaxVar*MaxVar*orderData.size();

    double U = 0;
    if (n_end - n_begin<4)
    {
        for (int k = n_begin; k <= n_end; k++)
            U += orderData[k];
        U = U / (1.0*(n_end - n_begin + 1));
    }
    else {
        for (int i = n_begin + 1; i < n_end - 1; i++)
        {
            for (int j = i; j < n_end - 1; j++)
            {
                StdVar(orderData, n_begin, i - 1, u1, v1);
                StdVar(orderData, i, j, u2, v2);
                StdVar(orderData, j + 1, n_end, u3, v3);
                double SumV = v1 + v2 + v3;
                if (SumV < MaxVar)
                {
                    U = u2;
                    MaxVar = SumV;
                }
            }
        }
    }
    return U;
}
double EstRadius(const std::vector<double> &RList)
{
    const double cutRatio = 0.2;
    const int n1 = RList.size()*cutRatio;
    const int n2 = RList.size()*(1.0 - 2 * cutRatio);
    double sum = 0.0;
    int N = n2-n1+1;    
    for (int i = n1; i<n2; ++i) {
        double r = RList.at(i);
        sum += r;
    }
    double meanR = sum / N;

    return meanR;
}
void AttachToCylinder(
    std::vector<int> &PtOnCylinder,
    ObjCylinder &cylinder,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const double _a, const double TR)
{
    assert(PointList.size() == NormList.size());
    QTime time;
    time.start();

    const bool bTR = TR > 0 ? true : false;
    const double TDis = bTR ? TR * 0.5 : _a * 20;
    const double TRMax = TR * 1.5;
    const double TRMin = TR * 0.5;
    const double TLayerContinueMax = 1.5;
    const double TLayerContinueMin = 0.8;
    const double TLayerPtRMax = 1.2;
    const double TLayerPtRMin = 0.8;

    const vcg::Point3f NN = cylinder.m_N;
    const vcg::Point3f PO = cylinder.m_O;


    flog(
        "      [--AttachToCylinder--]: #nPts-%d...\n"
        "        | #Axis       : < %7.4f, %7.4f, %7.4f > \n"
        "        | #Loc        : < %7.4f, %7.4f, %7.4f > \n",
        PointList.size(),
        PO.X(), PO.Y(), PO.Z(),
        NN.X(), NN.Y(), NN.Z());

   
    // 1. Initial
    std::vector<CylinderPointInfo> CPI;
    CPI.reserve(PointList.size());
    double MeanR = 0.0;
    for (int i = 0; i<PointList.size(); i++)
    {
        vcg::Point3f pp = PointList.at(i);
        vcg::Point3f np = NormList.at(i);

        vcg::Point3f nn = NN^np;
        vcg::Point3f dp = pp - PO;
        nn.Normalize();
        double distance = abs(nn*dp);
        double intercept = NN*dp;
        double radius = sqrt(dp.SquaredNorm() - intercept*intercept);
        if (distance<TDis && (!bTR || (radius<TRMax && radius>TRMin)))
        {
            CPI.push_back(CylinderPointInfo(i,radius,intercept));
            MeanR += radius;
        }
    }
    MeanR = MeanR / CPI.size();

    // 2. Sort
    std::sort(CPI.begin(), CPI.end(), CmpCPIByL);

    // 3. The Max Continuous Part
#if 1
    const double maxGap = _a;
    int MaxPartBegin, MaxPartEnd;
    double MaxPartLength = 0;
    for (int i = 0; i<CPI.size();)
    {
        int j = i;
        for (; j<CPI.size() - 1; j++)
        {
            if (CPI[j + 1].l - CPI[j].l>maxGap)
                break;
        }
        if (CPI[j].l - CPI[i].l>MaxPartLength)
        {
            MaxPartBegin = i;
            MaxPartEnd = j;
            MaxPartLength = CPI[j].l - CPI[i].l;
        }
        i = j + 1;
    }
#else
    int MaxPartBegin = 0;
    int MaxPartEnd = CPI.size()-1;
    double MaxPartLength = *(CPI.end()-1) - *(CPI.begin());
#endif

    flog(
        "        | ---------INIT---------\n"
        "        | #nPts       : %d\n"
        "        | #Thresholds : %7.4f [ %7.4f, %7.4f ]\n"
        "        | #Reslut_R   : %7.4f\n"
        "        | #Reslut_L   : [ %d : %7.4f : %d] -> %7.4f\n"
        "        | ---------INIT---------\n",
        CPI.size(),
        TDis, TRMin, TRMax,
        MeanR,
        MaxPartBegin, maxGap, MaxPartEnd, MaxPartLength);

    // 4. Identify at Each Layer
    const double TLayerContinueRMax = MeanR * TLayerContinueMax;
    const double TLayerContinueRMin = MeanR * TLayerContinueMin;

    std::vector<int> PtOnCylList;
    std::vector<double> PtMeanR;
    std::vector<int> PtCount;
    double lBegin, lEnd;

    std::vector<int> _PtOnCylList;
    std::vector<double> _PtMeanR;
    std::vector<int> _PtCount;
    int minIndex = MaxPartBegin;
    int maxIndex = MaxPartBegin;
    double minL = CPI[MaxPartBegin].l;
    double maxL = CPI[MaxPartBegin].l;
    double maxLength = maxL-minL;

    const double layerGap = _a*0.5;
    double layerLBegin = CPI[MaxPartBegin].l;
    double  layerLEnd = layerLBegin + layerGap;

#if defined(_ReportOut_)
    std::ofstream _out("../~R~.txt", std::ios::out);
#endif
    for (int i = MaxPartBegin; i <= MaxPartEnd; i++)
    {
        if ((CPI[i].l) > layerLEnd || i == MaxPartEnd)
        {
            if ((CPI[i].l) > layerLEnd && i == MaxPartEnd)
                maxIndex = i;

            std::vector<double> LayerR;
            for (int k = minIndex; k <= maxIndex; k++)
                LayerR.push_back(CPI[k].r);
            double uR = RobustMean(LayerR);

            bool bContinue = uR<TLayerContinueRMax && uR>TLayerContinueRMin;
            const double _TRLMax = TLayerPtRMax * uR;
            const double _TRLMin = TLayerPtRMin * uR;
            if (bContinue) {
                _PtMeanR.push_back(uR);
                _PtCount.push_back(LayerR.size());
                for (int k = minIndex; k <= maxIndex; k++) {
                    if (CPI[k].r < _TRLMax && CPI[k].r > _TRLMin) {
                        _PtOnCylList.push_back(CPI[k].index);
                        maxL = CPI[k].l;
                    }
                }
            }
            if (!bContinue || maxIndex == MaxPartEnd) {
                if (maxL - minL > maxLength) {
                    PtOnCylList.swap(_PtOnCylList);
                    PtMeanR.swap(_PtMeanR);
                    PtCount.swap(_PtCount);
                    lBegin = minL;
                    lEnd = maxL;
                    
                    _PtOnCylList.clear();
                    _PtMeanR.clear();
                    _PtCount.clear();

                    maxLength = maxL - minL;
                    if (i != MaxPartEnd)
                        minL = maxL = CPI[i + 1].l;

                    flog("        | #CircleGap  : %d [ %7.4f, %7.4f ] -> %7.4f \n",
                        PtOnCylList.size(), lBegin, lEnd, maxLength);
                }
            }
            //---------------------
            minIndex = maxIndex = i;
            layerLBegin = CPI[i].l;
            layerLEnd = layerLBegin + layerGap;

#if defined(_ReportOut_)
            _out << (layerLBegin + layerLEnd) / 2.0 << "\t" << uR << std::endl;
#endif
        }
        else
            maxIndex = i;
    }
#if defined(_ReportOut_)
    _out.close();
#endif

    cylinder.m_radius = EstRadius(PtMeanR);
    cylinder.m_length = lEnd - lBegin;
    cylinder.m_O = PO + NN*(lEnd + lBegin) / 2.0;

    PtOnCylinder.swap(PtOnCylList);

    flog(
        "        | #FReslut_R  : %7.4f\n"
        "        | #FReslut_L  : [ %7.4f, %7.4f ] -> %7.4f\n"
        "      [--AttachToCylinder--]: Done in %.4f seconds. \n",
        cylinder.m_radius, lBegin, lEnd, cylinder.m_length,
        time.elapsed() / 1000.0);
}


//------------------------------------
double SignedDistanceCylinderPoint(
    const ObjCylinder &cyl, const vcg::Point3f &p)
{
    vcg::Point3f op = p - cyl.m_O;
    double l1 = (op * cyl.m_N) / (cyl.m_N.Norm());
    double l2 = op.Norm();
    double l = (l2*l2 - l1*l1);
    l = l < 0 ? 0 : sqrt(l);
    return l - cyl.m_radius;
}
double AngCylinderPoint(
    const ObjCylinder &cyl, const vcg::Point3f &p, const vcg::Point3f &n)
{
    vcg::Point3f op = p - cyl.m_O;
    double l = (op * cyl.m_N) / (cyl.m_N.SquaredNorm());
    vcg::Point3f np = op - cyl.m_N*l;
    return CheckAng00(VCGAngle(np, n));
}
bool CylinderInlier(
	const ObjCylinder &cyl,
	const vcg::Point3f &p, const double TDis,
	const bool containInner)
{
	bool isInlier = false;
	double d = SignedDistanceCylinderPoint(cyl, p);
	double l = (p - cyl.m_O)*cyl.m_N / cyl.m_N.Norm();

	if (containInner)
		isInlier = (d < TDis) && (abs(l) - cyl.m_length/2.0 < TDis);
	else
		isInlier = (abs(d) < TDis);
	return isInlier;
}
bool CylinderInlier(
	const ObjCylinder &cyl,
	const vcg::Point3f &p, const double TDis,
	const vcg::Point3f &n, const double TAng,
	const bool containInner)
{
	bool isInlier = false;
	double d = SignedDistanceCylinderPoint(cyl, p);
	double l = (p - cyl.m_O)*cyl.m_N / cyl.m_N.Norm();

	if (containInner) {
		if ((d < -TDis) && (abs(l) - cyl.m_length / 2.0 < -TDis))
			isInlier = true;
		else
			isInlier = 
			(d < TDis) && (abs(l) - cyl.m_length / 2.0 < TDis) &&
			(AngCylinderPoint(cyl, p, n) < TAng);
	}
	else
		isInlier = (abs(d) < TDis) && (AngCylinderPoint(cyl, p, n) < TAng);


	return isInlier;
}
int CylinderInliers(
    const ObjCylinder &cyl,
    const std::vector<vcg::Point3f> &pointList,
    const std::vector<vcg::Point3f> &normList,
    const double TDis, const double TAng,
	const bool containInner,
    std::vector<int> *inlierIdx)
{
    const bool bHasNorm = normList.empty() ? false : true;
    if (bHasNorm)
        assert(normList.size() == pointList.size());
	
    std::vector<int> inliers;
    inliers.reserve(pointList.size());
    for (int i = 0; i < pointList.size(); ++i) {
		bool isInlier = false;
		if (bHasNorm)
			isInlier = CylinderInlier(cyl, pointList[i], TDis, normList.at(i), TAng, containInner);
		else
			isInlier = CylinderInlier(cyl, pointList[i], TDis, containInner);
        if (isInlier)
            inliers.push_back(i);
    }
    int retN = inliers.size();
    if (inlierIdx != 0)
        inlierIdx->swap(inliers);
    return retN;
}
bool CloseCylinders(
    const ObjCylinder &cyl1,
    const ObjCylinder &cyl2,
    const double TRDis, const double TAng)
{
    double meanR = (cyl1.m_radius + cyl2.m_radius)/2.0;
    double dR = abs(cyl1.m_radius - cyl2.m_radius) / 2.0;
    if (dR/meanR > TRDis)
        return false;
    double ang = CheckAng00(VCGAngle(cyl1.m_N, cyl2.m_N));
    if (ang > TAng)
        return false;
    if (abs(SignedDistanceCylinderPoint(cyl1, cyl2.m_O)) / meanR > TRDis)
        return false;

    return true;
}
ObjCylinder *EstCylinderTwoPoint(
    const vcg::Point3f p1, const vcg::Point3f n1,
    const vcg::Point3f p2, const vcg::Point3f n2,
    const double TDisDeviation, const int TAngRequired)
{
    if (CheckAng00(VCGAngle(n1, n2)) < TAngRequired)
        return 0;

    vcg::Point3f N = n1 ^ n2; N.normalized();

    // x1 := n1, z1 := N
    vcg::Point3f y1 = N ^ n1; 
    double p2_in_y1 = (p2 - p1)*y1;
    double n2_in_y1 = n2*y1;
    double lambda2  = p2_in_y1 / n2_in_y1;
    vcg::Point3f a2 = p2 - (n2 * lambda2);
    double       r2 = (n2 * lambda2).Norm();
    // x2 := n2, z1 := N
    vcg::Point3f y2 = N ^ n2;
    double p1_in_y2 = (p1 - p2)*y2;
    double n1_in_y2 = n1*y2;
    double lambda1  = p1_in_y2 / n1_in_y2;
    vcg::Point3f a1 = p1 - (n1 * lambda1);
    double       r1 = (n1 * lambda1).Norm();

    if (abs(r1 - r2) / (r1 + r2) > TDisDeviation)
        return 0;

    
    vcg::Point3f O = (a1 + a2) / 2.0;
    double R = (r1 + r2) / 2.0;
    ObjCylinder *cyl = new ObjCylinder(0);
    cyl->m_O = O;
    cyl->m_N = N;
    cyl->m_radius = R;
    cyl->m_length = 1.0;

    return cyl;
}
ObjCylinder *FineCylinder(
    const std::vector<vcg::Point3f> &pointList,
    const std::vector<int> &cylVerList,
    double &err)
{
    QTime time;
    time.start();

    std::vector<gte::Vector3<double>> positions;
    positions.reserve(cylVerList.size());
    for (unsigned int i = 0; i < cylVerList.size(); ++i)
    {
        vcg::Point3f pt = pointList.at(cylVerList.at(i));
        gte::Vector3<double> data;
        data[0] = pt.X();
        data[1] = pt.Y();
        data[2] = pt.Z();
        positions.push_back(data);
    }
    unsigned int numThreads = std::thread::hardware_concurrency();
    gte::ApprCylinder3<double> fitter(numThreads, 1024, 512);

    unsigned int numVertices = static_cast<unsigned int>(positions.size());
    gte::Cylinder3<double> cylinder;
    double minError = fitter(numVertices, positions.data(), cylinder);

    ObjCylinder *cyl = new ObjCylinder(_GetObjCode(Pt_OnCylinder));
    cyl->m_O.X() = cylinder.axis.origin[0];
    cyl->m_O.Y() = cylinder.axis.origin[1];
    cyl->m_O.Z() = cylinder.axis.origin[2];
    cyl->m_N.X() = cylinder.axis.direction[0];
    cyl->m_N.Y() = cylinder.axis.direction[1];
    cyl->m_N.Z() = cylinder.axis.direction[2];
    cyl->m_radius = cylinder.radius;
    cyl->m_length = cylinder.height;
    err = minError;
    
    

    flog(
        "      [--Fit_Cylinder--]: #Pts-%d\n"
        "        | #ID        : %d\n"
        "        | #min error : %7.4f\n"
        "        | #Center    : < %7.4f, %7.4f, %7.4f > \n"
        "        | #direction : < %7.4f, %7.4f, %7.4f > \n"
        "        | #radius    : %7.4f\n"
        "        | #height    : %7.4f\n"
        "      [--Fit_Cylinder--]: Done in %.4f seconds. \n",
        cylVerList.size(), cyl->m_index, err,
        cyl->m_O.X(), cyl->m_O.Y(), cyl->m_O.Z(),
        cyl->m_N.X(), cyl->m_N.Y(), cyl->m_N.Z(),
        cyl->m_radius, cyl->m_length,
        time.elapsed() / 1000.0);

    return cyl;
}


typedef std::pair<ObjCylinder*, double> WCylinder;
int RansacTimes(const double Percentage, const double InlierRatio, const int MiniSupportNum) {
    // P := Percentage
    // I := InlierRatio
    // N := MiniSupportNum
    // P = 1 - (1-I^N)^T
    // >> T = log(1-P)/log(1-I^N)
    return log(1 - Percentage) / log(1 - pow(InlierRatio, MiniSupportNum)) + 0.5;
}
void AddLimitList(
    std::list<WCylinder> &cyls, const WCylinder &cyl, const int maxN)
{
    if (cyls.size() < maxN) {
        // Just Add One
        cyls.push_back(cyl);
    }
    else {
        // Add One And Pop The Last
        auto iter = cyls.begin();
        for (;iter != cyls.end();iter++) {
            if (iter->second < cyl.second)
                break;
        }
        if (iter != cyls.end()) {
            cyls.insert(iter, cyl);
            cyls.pop_back();
        }
    }
}
double DetectCylinderRansac(
    const std::vector<vcg::Point3f> &pointList,
    const std::vector<vcg::Point3f> &normList,
    std::vector<ObjCylinder*> &cylCandidates,
    const double TDis, const double TAng,
    const int maxN, const double inlierRatio,
    const vcg::Box3f *constriantBox)
{
    assert(pointList.size() == normList.size());
    assert(inlierRatio > 0.0 && inlierRatio < 1.0);
    QTime time;
    time.start();

    if (pointList.size() == 0)
        return -1.0;

    const int NPts = pointList.size();

    const double Percentage = 0.9;
    const double TDisDeviation = (1-0.6) / (1+0.6);
    const int TAngRequired = 15;
    const int SupportN = 2;
    const int TheoryIteration = RansacTimes(Percentage, inlierRatio*0.1, SupportN);
    const int MaxIteration = TheoryIteration < 10 ? 10 : (TheoryIteration > 1e4 ? 1e4 : TheoryIteration);
    const double MaxR = (constriantBox == 0) ? 0.0 : constriantBox->Dim().V(constriantBox->MinDim());

    std::list<WCylinder> cyls;
    double maxInlierRatio = 0.0;
    int iter = 0;
    std::default_random_engine randomEngine; // Keep the result same
    std::vector<int> idx;
    idx.reserve(NPts);
    for (int i = 0; i < NPts; ++i)
        idx.push_back(i);
    while (iter < MaxIteration) {
        iter++;       
        vcg::Point3f p[2], n[2];
        // Random Pick
        for (int i = 0; i < 2; ++i) {
            // Random distribution
            std::uniform_int_distribution<int> distribution(i, NPts - 1);
            int k = distribution(randomEngine);
            std::swap(idx[i], idx[k]);
            p[i] = pointList.at(idx[i]);
            n[i] = normList.at(idx[i]);
        }     
        // Estimate
        ObjCylinder *cyl = EstCylinderTwoPoint(p[0], n[0], p[1], n[1], TDisDeviation, TAngRequired);
        if (cyl == 0) // Fail to Est by Two Points
            continue;

        // Good Model Check
        if (constriantBox != 0 && (
            !constriantBox->IsIn(cyl->m_O) ||
            cyl->m_radius > MaxR / 2.0)) { // Not A Good Model For Loc & Size
            delete cyl;
            continue;
        }
        double inlierRatio_this = CylinderInliers(*cyl, pointList, normList, TDis, TAng) * 1.0 / NPts;
        if (inlierRatio_this < inlierRatio) { // Not A Good Model For Inlier Ratio
            delete cyl;
            continue;
        }

        // A Good Model
        AddLimitList(cyls, WCylinder(cyl, inlierRatio_this), maxN);

        // Update Param
        if (inlierRatio_this > maxInlierRatio) {
            maxInlierRatio = inlierRatio_this;
        }
    }

    std::vector<ObjCylinder*> _cylCandidates;
    _cylCandidates.reserve(cyls.size());
    for (auto iter = cyls.begin(); iter != cyls.end(); ++iter)
        _cylCandidates.push_back(iter->first);
    cylCandidates.swap(_cylCandidates);
    flog(
        "      [--DetectCylinderRansac--]: #nPts-%d...\n"
        "        | #ExpCand    : %d \n"
        "        | #ExpInlier  : %7.4f \n"
        "        | #CanNum     : %d \n"
        "        | #MaxInlier  : %7.4f \n"
        "        | #Iteration  : %d - [ %d | %d ] \n"
        "      [--DetectCylinderRansac--]: Done in %.4f seconds. \n",
        NPts, maxN, inlierRatio,
        cyls.size(), maxInlierRatio,
        iter, MaxIteration, TheoryIteration,
        time.elapsed() / 1000.0);    

    return maxInlierRatio;
}
double FineCylinderLength(
    CMeshO &mesh,
    ObjCylinder &cyl,
    std::vector<int> &inlierIdx)
{
    CMeshO::VertexIterator vi = mesh.vert.begin();
    std::vector<std::pair<double, int>> LList;
    vcg::Point3f N = cyl.m_N;
    N.normalized();
    for (int i = 0; i < inlierIdx.size(); ++i) {
        int idx = inlierIdx.at(i);
        vcg::Point3f &pt = (vi + idx)->cP();
        double l = (pt - cyl.m_O)*N;
        LList.push_back(std::pair<double, int>(l, idx));
    }
    std::sort(LList.begin(), LList.end());

    int idxBegin, idxEnd;
    PatchDimensionOne(LList, idxBegin, idxEnd);

    std::vector<int> _inlierIdx;
    _inlierIdx.reserve(idxEnd - idxBegin + 1);
    for (int i = idxBegin; i <= idxEnd; ++i)
        _inlierIdx.push_back(LList.at(i).second);

    double newLength =abs( LList.at(idxEnd).first - LList.at(idxBegin).first );
    cyl.m_length = newLength;
    cyl.m_O = cyl.m_O + N*(LList.at(idxBegin).first + LList.at(idxEnd).first) / 2.0;
    inlierIdx.swap(_inlierIdx);

    return newLength;
}
int AttachToCylinder(
	CMeshO &mesh,
	std::vector<ObjCylinder*> &cylinders,
	const double TDis, const double TAng,
	const double TInlierRatio)
{
	QTime time;
	time.start();
	CMeshO::PerVertexAttributeHandle<PtType> type_hi =
		vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
	int count = 0;
	for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); vi++) {
		if (!vi->IsD() && type_hi[vi] == Pt_Undefined)
			count++;
	}
	const int T = count*TInlierRatio;
	flog(
		"      [--AttachToCylinder--]: #NPts-%d...\n"
		"        | #InputCyl   : %d \n"
		"        | #TInlier    : %.4f \n",
		count,
		cylinders.size(), TInlierRatio);

	count = 0;
	for (int k = 0; k < cylinders.size(); ++k) {
		ObjCylinder *cyl = cylinders.at(k);
		std::vector<int> inlierIdx;
		int idx = 0;
		for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); vi++, idx++) {
			if (!vi->IsD() && type_hi[vi] == Pt_Undefined &&
				CylinderInlier(*cyl, vi->cP(), TDis, vi->cN(), TAng, false))
				inlierIdx.push_back(idx);
		}
		if (inlierIdx.size() < T) {
			flog("        | >> Quit [ Id.%d ] plane with inliers [ %d < %d ] \n", cyl->m_index, inlierIdx.size(), T);
			delete cyl;
			cylinders.erase(cylinders.begin() + k);
			k--;			
		}
		else {
            FineCylinderLength(mesh, *cyl, inlierIdx);
			for (int i = 0; i < inlierIdx.size(); ++i)
				type_hi[inlierIdx[i]] = Pt_OnCylinder;
			count += inlierIdx.size();
		}
	}
	flog(
		"        | #OutPutCyl  : %d \n"
		"        | #NPts       : %d \n"
		"      [--AttachToCylinder--]: Done in %.4f seconds. \n",
		cylinders.size(), count,
		time.elapsed() / 1000.0);

	return count;
}
