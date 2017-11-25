#include "PointCloudFitUtil.h"
#include "tool/PCA.h"

bool DetectSymAxis(
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const double _a,
    vcg::Point3f &PO, vcg::Point3f &NN
)
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

    
#if 0
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
    int TVote = maxVote * 0.5;
#else
    int TVote = 5;
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
    if (ret == -1)
        return false;

    // -- [Assign Norm]
    NN.X() = PC[0];
    NN.Y() = PC[3];
    NN.Z() = PC[6];

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
double EstRadius(const std::vector<double> &RList, double &weight)
{
    const double cutRatio = 0.2;
    const int n1 = RList.size()*cutRatio;
    const int n2 = RList.size()*(1.0 - 2 * cutRatio);
    double sum = 0.0;
    double sum2 = 0.0;
    int N = n2-n1+1;    
    for (int i = n1; i<n2; ++i) {
        double r = RList.at(i);
        sum += r;
        sum2 += r*r;
    }
    double meanR = sum / N;
    double var = sqrt(sum2 / N - meanR*meanR);

    weight = 1.0 - var / meanR;
    return meanR;
}
void AttachToCylinder(
    std::vector<int> &PtOnCylinder,
    ObjCylinder &cylinder,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const double _a, const double TR)
{
    const bool bTR = TR > 0 ? true : false;
    const double TDis = bTR ? TR * 0.5 : _a * 20;
    const double TRMax = TR * 1.5;
    const double TRMin = TR * 0.5;
    const double TLayerContinueMax = 1.5;
    const double TLayerContinueMin = 0.8;
    const double TLayerPtRMax = 1.2;
    const double TLayerPtRMin = 0.8;


    const vcg::Point3f NN = cylinder.m_N;
    const vcg::Point3f PO = cylinder.m_N;

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

    cylinder.m_EIConfL = EIConfidence(PtCount);
    cylinder.m_radius = EstRadius(PtMeanR, cylinder.m_EIConfR);
    cylinder.m_length = lEnd - lBegin;
    cylinder.m_O = PO + NN*(lEnd + lBegin) / 2.0;

    PtOnCylinder.swap(PtOnCylList);
}