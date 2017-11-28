#include "PointCloudFitUtil.h"
#include "tool/PCA.h"

// [1] Detect Planes
int HoughPlaneOne(
    vcg::Plane3f &Plane,
    const FixedAxis fix,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const double _intercept, const double _a, const double _s)
{
    QTime time;
    time.start();

    Plane = _NON_PLANE;

    const bool bHasNorm = (NormList.size() > 0) ? true : false;
    if (bHasNorm)
        assert(PointList.size() == NormList.size());

    int index1, index2, index3;
    switch (fix)
    {
    case FixedAxis_Z: index1 = 2; index2 = 0; index3 = 1; break; //ZXY
    case FixedAxis_Y: index1 = 1; index2 = 2; index3 = 0; break; //YZX
    case FixedAxis_X: index1 = 0; index2 = 1; index3 = 2; break; //XYZ
    default: return -1;
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

    // 1.0*x + fa*y + fb*z = fc;
    int a, b, c;
    double fa, fb, fc;
    double x, y, z;
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
            iterN++;
            ny = ny / nx;
            nz = nz / nx;
            if (abs(ny) - 1.0>0.3 || abs(nz) - 1.0>0.3)
                continue;
            a_begin = ny - 0.3;
            a_end = ny + 0.3;
            b_begin = nz - 0.3;
            b_end = nz + 0.3;

            if (a_begin<-1) a_begin = -1;
            if (a_end>1) a_end = 1;
            if (b_begin<-1) b_begin = -1;
            if (b_end>1) b_end = 1;
        }

        for (fa = a_begin; fa<a_end; fa += scale1) {
            for (fb = b_begin; fb<b_end; fb += scale1) {
                fc = -x - fa*y - fb*z;
                c = fc / scale2 + C / 2;
                if (c >= 0 && c<C)
                {
                    a = (fa / scale1 + A / 2);
                    b = (fb / scale1 + B / 2);
                    (facPlane[c])[a*A + b] += 1;
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

    
    double pa = (Temp[0] - A / 2)*scale1;
    double pb = (Temp[1] - B / 2)*scale1;
    double pc = (Temp[2] - C / 2)*scale2;
    int vote = Temp[3];

    vcg::Point3f N;
    switch (fix)
    {
    case FixedAxis_Z: N = vcg::Point3f( pa,  pb, 1.0); break; //ZXY
    case FixedAxis_Y: N = vcg::Point3f( pb, 1.0, pa ); break; //YZX
    case FixedAxis_X: N = vcg::Point3f(1.0,  pa, pb ); break; //XYZ
    default: return -1;
    }
    // N*p = Off
    Plane.Set(N, -pc);

    flog(
        "        | [#Time-%7.4f]-[%c-Seted]: %d-Pts \n"
        "        |  >> %7.4fX + %7.4fY + %7.4fZ + %7.4f = 0 |=> %7.4fX + %7.4fY + %7.4fZ = %7.4f \n",
        time.elapsed() / 1000.0, 
        fix == FixedAxis_X ? 'X' : (fix == FixedAxis_Y ? 'Y' : 'Z'),
        vote,
        N.X(), N.Y(), N.Z(), pc,
        Plane.Direction().X(), Plane.Direction().Y(), Plane.Direction().Z(), Plane.Offset()
        );

    return vote;
}
int HoughPlane(
    vcg::Plane3f &Plane,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const double _intercept, const double _a, const double _s)
{
    QTime time;
    time.start();
    flog(
        "      [--Plane_HT--]: #nPts-%d \n",
        PointList.size());

    vcg::Plane3f P[3];
    int N[3];
#ifdef _USE_OPENMP_
#pragma omp parallel for
#endif // !_USE_OPENMP_
    for (int i = 0; i < 3; i++) {
        N[i] = HoughPlaneOne(P[i], FixedAxis(i), PointList, NormList, _intercept, _a, _s);
    }

    int retIdx = 0;
    if (N[0]>N[1])
        retIdx = (N[0]>N[2]) ? 0 : 2;
    else
        retIdx = (N[1]>N[2]) ? 1 : 2;

    Plane = P[retIdx];

    flog(
        "        | [+][Checked] : [%d] #nPts-< %d > \n"
        "      [--Plane_HT--]: Done in %.4f seconds \n",
        retIdx + 1, N[retIdx], time.elapsed() / 1000.0);

    return N[retIdx];
}


// [2] Surface Points Verification
int AttachToPlane(
    std::vector<int> &PtIdx,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const vcg::Plane3f &Plane,
    const double _TDis, const double _TAng)
{
    QTime time;
    time.start();

    const bool bHasNorm = (NormList.size() > 0) ? true : false;
    if (bHasNorm)
        assert(PointList.size() == NormList.size());

    // Check Distance
    std::vector<int> OnPlaneIdxList;
    for (int i = 0; i<PointList.size(); i++) {
        double dis = vcg::SignedDistancePlanePoint(Plane, PointList[i]);
        if (abs(dis)<_TDis)
            OnPlaneIdxList.push_back(i);
    }

    // Check Norm
    if (bHasNorm && !OnPlaneIdxList.empty())
    {
        vcg::Point3f NP = Plane.Direction();
        std::vector<int> OnPlaneList_NormChecked;
        for (int i = 0; i<OnPlaneIdxList.size(); ++i) {
            vcg::Point3f npt = NormList.at(OnPlaneIdxList.at(i));
            npt.Normalize();
            double ang = 90 - abs(90 - R2D(vcg::AngleN(NP, npt)));
            if (ang<_TAng)
                OnPlaneList_NormChecked.push_back(OnPlaneIdxList.at(i));
        }

        OnPlaneIdxList.swap(OnPlaneList_NormChecked);
    }

    PtIdx.swap(OnPlaneIdxList);

    flog(
        "      [--AttachToPlane--]: Attach points to planes...\n"
        "        | #Threshold_Dis : %.4f\n"
        "        | #Threshold_Ang : %.4f-[%d]\n"
        "        | #nPts-OnPlane  : %d \n"
        "      [--AttachToPlane--]: Done in %.4f seconds. \n",
        _TDis, _TAng, bHasNorm,
        PtIdx.size(), time.elapsed() / 1000.0);
   
    return PtIdx.size();
}


// [3] Coplanar Separation
void PicMaxRegion(
    const std::vector<vcg::Point3f> &PointList,
    std::vector<int> &Index,
    const double _TDis)
{
    QTime time;
    time.start();

    if (PointList.empty())
        return;

    std::vector<int> tempLoc;
    std::vector<int> LocFinal;

    std::vector<vcg::Point3f> seedStack;
    vcg::Point3f seed;
    std::vector<int> recorder;
    const int NPtsAll = Index.size();
    while (!Index.empty())
    {
        // Each Cluster
        seedStack.clear();
        tempLoc.clear();
        seedStack.push_back(PointList.at(*(Index.end()-1)));       
        tempLoc.push_back(*(Index.end() - 1));

        Index.pop_back();
        // Growing
        while (!seedStack.empty())
        {
            seed = *(seedStack.end()-1);
            seedStack.pop_back();
            recorder.clear();
            for (int i = 0; i<Index.size(); i++)
            {
                int p = Index.at(i);
                if (vcg::Distance(seed, PointList.at(p))<_TDis)
                {
                    seedStack.push_back(PointList.at(p));
                    tempLoc.push_back(p);
                    recorder.push_back(i);
                }
            }
            for (int i = recorder.size() - 1; i >= 0; i--)
                Index.erase(Index.begin() + recorder.at(i));
        }
        if (tempLoc.size() > LocFinal.size())
            LocFinal.swap(tempLoc);
    }
    Index.swap(LocFinal);

    LocFinal.clear();
    tempLoc.clear();
    recorder.clear();
    std::sort(Index.begin(), Index.end());

    flog(
        "      [--MaxRegion--]: #Pts-%d\n"
        "        | #Threshold_Dis : %.4f\n"
        "        | #NPts_MaxReg   : %d\n"
        "      [--MaxRegion--]: Done in %.4f seconds. \n",
        NPtsAll, _TDis, Index.size(), time.elapsed() / 1000.0);
}


// [4] LS Fit
double FinePlane(
    const std::vector<vcg::Point3f> &pointList,
    const std::vector<int> &planeVerList,
    vcg::Plane3f &plane)
{
    QTime time;
    time.start();

    std::vector<vcg::Point3f> ExactVec;
    vcg::Plane3f ple;
    for (int j = 0; j<planeVerList.size(); ++j)
    {
        vcg::Point3f p = pointList.at(planeVerList.at(j));
        ExactVec.push_back(p);
    }

    vcg::FitPlaneToPointSet(ExactVec, ple);
    float fitError = 0;
    for (int i = 0; i < ExactVec.size(); ++i) {
        float d = vcg::SignedDistancePlanePoint(ple, ExactVec[i]);
        fitError += d*d;
    }
    fitError = sqrt(fitError / ExactVec.size());

    plane = ple;

    flog(
        "      [--Fit_LS--]: #Pts-%d\n"
        "        | #Plane    : %7.4fX + %7.4fY + %7.4fZ = %7.4f \n"
        "        | #FitError : %.4f\n"
        "      [--Fit_LS--]: Done in %.4f seconds. \n",
        planeVerList.size(),
        plane.Direction().X(), plane.Direction().Y(), plane.Direction().Z(), plane.Offset(),
        fitError, time.elapsed() / 1000.0);

    return fitError;
}

// [5] Get Minimum-Bounding-Rectangle
bool PatchDimensionOne(
    std::vector<std::pair<double, int>> &proList,
    int &idxBegin, int &idxEnd,
    double r)
{
    int n = proList.size();
    int p1 = 0;
    int p2 = n - 1;
    double l = abs(proList.at(p2).first - proList.at(p1).first);
    bool updated = true;
    while (updated) {
        updated = false;
        if (n <= 2)
            break;
        
        // r * (1.0/li) < n / l := r*l<li*n
        // remove this point

        // Move Left      
        while (p1 < p2) {
            double l1 = abs(proList.at(p1+1).first - proList.at(p1).first);           
            if (r*l < l1*n) {
                l -= l1;
                n--;
                p1++;
                updated = true;
            }
            else
                break;
        }
        // Move Right
        while (p1 < p2) {
            double l2 = abs(proList.at(p2).first - proList.at(p2-1).first);
            // r * (1.0/l1) < n / l := r*l<l1*n
            // remove this point
            if (r*l < l2*n) {
                l -= l2;
                n--;
                p2--;
                updated = true;
            }
            else
                break;
        }
    }

    if (n > proList.size() / 2) {
        idxBegin = p1;
        idxEnd = p2;
        return true;
    }
    else
        return false;
}
double PatchDimensions(
    const std::vector<vcg::Point3f> & pts,
    const vcg::Point3f &nx,
    const vcg::Point3f &ny,
    vcg::Point3f &Pt_O, vcg::Point3f &Pt_Dx, vcg::Point3f &Pt_Dy,
    double *ex, double *ey)
{
    const int NPP = pts.size();
    // Project
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

    // 
    int xBegin = 0;
    int xEnd = lX.size() - 1;
    int yBegin = 0;
    int yEnd = lY.size() - 1;
    PatchDimensionOne(lX, xBegin, xEnd);
    PatchDimensionOne(lY, yBegin, yEnd);

    // Get Patch
    vcg::Point3f MinX = pts.at(lX.at(xBegin).second);
    vcg::Point3f MinY = pts.at(lY.at(yBegin).second);
    double lx = lX.at(xEnd).first - lX.at(xBegin).first;
    double ly = lY.at(yEnd).first - lY.at(yBegin).first;
    Pt_Dx = nx*lx;
    Pt_Dy = ny*ly;
    double MinY2X = MinY*nx;
    double MinX2Y = MinX*ny;
    vcg::Point3f O1 = MinY - nx*(MinY2X - lX.at(xBegin).first);
    vcg::Point3f O2 = MinX - ny*(MinX2Y - lY.at(yBegin).first);
    Pt_O = (O1 + O2) / 2.0;
    if (ex != 0) {
        std::vector<double> list;
        for (int i = xBegin; i <= xEnd; ++i)
            list.push_back(lX.at(i).first);
        *ex = EIConfidence(list);
    }
    if (ey != 0) {
        std::vector<double> list;
        for (int i = yBegin; i <= yEnd; ++i)
            list.push_back(lY.at(i).first);
        *ey = EIConfidence(list);
    }
    return lx*ly;
}

ObjRect* ExtractMBR(
    CMeshO &mesh,
    const vcg::Plane3f &Plane,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<int> &IndexList,
    const std::vector<int> &PlaneVerList)
{
    ObjRect* oneRect = 0;

    QTime time;
    time.start();

    if (Plane == _NON_PLANE ||
        PlaneVerList.empty())
        return 0;

    assert(IndexList.size() == PointList.size());
    flog("      [--ExtractMBR--]: #Pts-%d\n", PlaneVerList.size());

    // -- Get Project Points & Set Labels  
    CMeshO::PerVertexAttributeHandle<PtType> type_hi =
        vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
    std::vector<vcg::Point3f> OnPPt;
    const int PlaneCode = _GetPlaneCode();
    for (int i = PlaneVerList.size() - 1; i >= 0; --i)
    {
        int index = PlaneVerList.at(i);

        // Set labels      
        type_hi[IndexList.at(index)] = PlaneCode;

        // Get project points
        vcg::Point3f pp = Plane.Projection(PointList.at(index));
        OnPPt.push_back(pp);
    }
  
    // -- Initial Direction By PCA
    const int N = OnPPt.size();
    double *data = new double[3 * N];
    for (int i = 0; i<N; ++i)
    {
        data[i] = OnPPt.at(i).X();
        data[N + i] = OnPPt.at(i).Y();
        data[2 * N + i] = OnPPt.at(i).Z();
    }
    int row = 3, col = N;
    double PC[9], V[3];//V[3] useless

    int ret = PCA(data, row, col, PC, V);
    delete[] data;
    assert(ret != -1);

    vcg::Point3f NP = Plane.Direction();
    vcg::Point3f NX_0 = vcg::Point3f(PC[0], PC[3], PC[6]);
    vcg::Point3f NY_0 = vcg::Point3f(PC[1], PC[4], PC[7]);

    NY_0 = NP^NX_0;
    NX_0 = NY_0^NP;
    NX_0.Normalize();
    NY_0.Normalize();
    flog(
        "        | #Loc-PCA_X   : < %7.4f, %7.4f, %7.4f > \n"
        "        | #Loc-PCA_Y   : < %7.4f, %7.4f, %7.4f > \n",
        NX_0.X(), NX_0.Y(), NX_0.Z(),
        NY_0.X(), NY_0.Y(), NY_0.Z()
    );

    // -- Refine Direction By Projection
    vcg::Point3f O, DX, DY;
    double confidenceX, confidenceY;
    double SA = PatchDimensions(OnPPt, NX_0, NY_0, O, DX, DY, &confidenceX, &confidenceY);
    for (int i = 1; i < 90; i += 1) { // Not too slow
        double k = tan(i *_D2R);
        vcg::Point3f nx = NX_0 + NY_0*k; nx.Normalize();       
        vcg::Point3f ny = NP ^nx;
        vcg::Point3f o, dx, dy;
        double ex, ey;
        double sa = PatchDimensions(OnPPt, nx, ny, o, dx, dy, &ex, &ey);
        if (sa < SA) {
            SA = sa;
            O = o; DX = dx; DY = dy;
            confidenceX = ex; confidenceY = ey;
            flog(
                "        | #RotAng-S    : < %d > -< %7.4f > \n"
                "        | #Loc-Adjed_O : < %7.4f, %7.4f, %7.4f > \n"
                "        | #Loc-Adjed_X : < %7.4f, %7.4f, %7.4f > - <%7.4f> - [%7.4f%]\n"
                "        | #Loc-Adjed_Y : < %7.4f, %7.4f, %7.4f > - <%7.4f> - [%7.4f%]\n",
                i, SA,
                O.X(), O.Y(), O.Z(),
                DX.X(), DX.Y(), DX.Z(), DX.Norm(), confidenceX * 100,
                DY.X(), DY.Y(), DY.Z(), DY.Norm(), confidenceY * 100);
        }
    }
    oneRect = new ObjRect(PlaneCode);
    oneRect->m_O = O;
    oneRect->m_N = NP;
    oneRect->m_AX = DX;
    oneRect->m_AY = DY;
    oneRect->m_EIConfX = confidenceX;
    oneRect->m_EIConfY = confidenceY;
    //
    flog(
        "        | #------------------------------------- \n"
        "        | #Loc-Adjed_O : < %7.4f, %7.4f, %7.4f > \n"
        "        | #Loc-Adjed_X : < %7.4f, %7.4f, %7.4f > - <%7.4f> - [%7.4f%]\n"
        "        | #Loc-Adjed_Y : < %7.4f, %7.4f, %7.4f > - <%7.4f> - [%7.4f%]\n"
        "        | #------------------------------------- \n"
        "      [--ExtractMBR--]: Done in %.4f seconds. \n",
        O.X(), O.Y(), O.Z(),
        DX.X(), DX.Y(), DX.Z(), DX.Norm(), confidenceX*100,
        DY.X(), DY.Y(), DY.Z(), DY.Norm(), confidenceY*100,
        time.elapsed() / 1000.0);
    OnPPt.clear();

    return oneRect;
}
ObjCircle* CircleCheck(
    const ObjRect *rect,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<int> &PlaneVerList)
{
    QTime time;
    time.start();

    if (rect == 0)
        return 0;

    const double T1 = 0.05;
    const double T2 = 0.9;

    const double w = rect->width();
    const double h = rect->height();
    const double r = (w + h) / 2;
    const double k1 = abs(w - h)*0.5 / r;
    if (k1 > T1)
        return 0;

    const vcg::Point3f &ro = rect->m_O;
    const vcg::Point3f &rx = rect->m_AX;
    const vcg::Point3f &ry = rect->m_AY;
    const double cx = w*0.5;
    const double cy = h*0.5;
    const double sx = 1.0 / w;
    const double sy = 1.0 / h;
    const double r2 = r*r*0.25;
    const double r4 = r*r*0.125;
    int count_r = 0;
    int count_c = 0;
    int count_ci = 0;
    for (int i = 0; i < PlaneVerList.size(); ++i) {
        const vcg::Point3f &p = PointList.at(PlaneVerList.at(i));
        double x = (p - ro)*rx*sx;
        double y = (p - ro)*ry*sy;
        if (x >= 0.0 && x <= w && y >= 0.0 && y <= h) {
            count_r++;
            double d = (x - cx)*(x - cx) + (y - cy)*(y - cy);
            if (d <= r2)
                count_c++;
            if (d <= r4)
                count_ci++;
        }
    }

    double k2 = count_c*1.0 / count_r;
    if (k2 < T2)
        return 0;

    double k3 = count_ci*1.0 / count_c;
    if (k3 > 0.6)
        return 0;

    ObjCircle *circle = new ObjCircle(rect->m_index);
    circle->m_N = rect->m_N;
    circle->m_O = ro + (rx + ry) / 2.0;
    circle->m_varN = rect->m_varN;
    circle->m_radius = r*0.5;


    flog(
        "      [--CircleCheck--]: < %7.4f, %7.4f >\n"
        "        | #|X-Y|/|X+Y| :  %7.4f \n"
        "        | #InnerRate   :  %7.4f \n"
        "        | #CircleRate  :  %7.4f \n"
        "        | #------------------------------------- \n"
        "      [--CircleCheck--]: Done in %.4f seconds. \n",
        T1,T2,k1,k2,k3,
        time.elapsed() / 1000.0);

    return circle;
}

int DetectHTPlanes(
    std::vector<vcg::Plane3f> &Planes,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const double _intercept, const double _a, const double _s,
    const double TDis, const double TAng,
    const double TNPtsRatio, const int TNPtsHard,
    const int ExpPlaneNum,
    std::vector<double> *errors,
    std::vector<ObjPatch*> *pPlanes,
    CMeshO *pMesh, std::vector<int> *pIndexList)
{
    assert(
        (pPlanes == 0 && pMesh == 0 && pIndexList == 0) ||
        (pPlanes != 0 && pMesh != 0 && pIndexList != 0)
    );
    const bool bHasNorm = (NormList.size() > 0) ? true : false;
    const bool bRetPlane = (pPlanes != 0) ? true : false;
    const bool bReportErr = (errors != 0) ? true : false;

    std::vector<int> indexList;
    if (bRetPlane)
        indexList = *pIndexList;
    std::vector<vcg::Point3f> pointList = PointList;
    std::vector<vcg::Point3f> normList = NormList;
      
    if (bHasNorm)
        assert(pointList.size() == normList.size());
    if (bReportErr)
        errors->clear();

    std::vector<vcg::Plane3f> planeVec;
    int planeNum = 0;
    _ResetPlaneCode();
    while (1)
    {
        if (pointList.empty())
            break;
        if (ExpPlaneNum >= 0 && planeNum >= ExpPlaneNum)
            break;


        // int _planeNT = pointList.size()*TNPtsRatio;
        int _planeNT = fmax(pointList.size()*TNPtsRatio, TNPtsHard);
        if (bHasNorm)
            _planeNT *= 0.5;

        flog("    >> Detecting the [ No.%d ] plane with #Minimum - [ %d ] ...\n", planeNum + 1, _planeNT);

        vcg::Plane3f plane;
        // HT Detection            
        int NP = HoughPlane(plane, pointList, normList, _intercept, _a, _s); // Center At (0,0,0)
        if (plane == _NON_PLANE || NP <= _planeNT)
            break;

        // Surface Points Verification
        std::vector<int> planeVerList;
        AttachToPlane(planeVerList, pointList, normList, plane, TDis, TAng);
        if (planeVerList.size() < _planeNT)
            break;

        // Coplanar Separation
        PicMaxRegion(pointList, planeVerList, TDis);
        if (planeVerList.size() < TNPtsHard) {
            flog(
                "      [ -- The number of the max region is too small,       -- ] \n"
                "      [ -- and this region result will be DISCARD.          -- ] \n"
                "      [ -- Detection of the [ No.%02d ] plane wil be retired. -- ] \n", 
                planeNum + 1);
        }
        else {
            // Least Squre Fit
            double err = FinePlane(pointList, planeVerList, plane);
            // planeVerList = AttachToPlane(pointList,directionList,plane,_planeDThreshold); // AGAIN!
            if (bReportErr)
                errors->push_back(err);

            if (bRetPlane) {

                // Extract the Minimum-Bounding-Rectangle
                ObjRect *oneRect = ExtractMBR(*pMesh, plane, pointList, indexList, planeVerList);
                if (oneRect != 0) {
                    oneRect->m_varN = err;
                    // Circle Check
                    ObjCircle *oneCircle = CircleCheck(oneRect, pointList, planeVerList);
                    if (oneCircle != 0) {
                        delete oneRect;
                        pPlanes->push_back(oneCircle);
                    }
                    else
                        pPlanes->push_back(oneRect);
                }

            }
            planeVec.push_back(plane);
            planeNum++;
        }

        // Remove Pts on Plane
        for (int i = planeVerList.size() - 1; i >= 0; --i) {
            int index = planeVerList.at(i);
            pointList.erase(pointList.begin() + index);
        }
        if (bHasNorm) {
            for (int i = planeVerList.size() - 1; i >= 0; --i) {
                int index = planeVerList.at(i);
                normList.erase(normList.begin() + index);
            }
        }
        if (bRetPlane) {
            for (int i = planeVerList.size() - 1; i >= 0; --i) {
                int index = planeVerList.at(i);
                indexList.erase(indexList.begin() + index);
            }
        }
        planeVerList.clear();

    }
    if (bRetPlane)
        indexList.clear();
    pointList.clear();
    normList.clear();

    Planes.swap(planeVec);
    return Planes.size();
}

int ExtractPatches(
    CMeshO &mesh,
    std::vector<ObjPatch*> &patches,
    const std::vector<int> &indexList,
    const std::vector<vcg::Point3f> &pointList,
    const int planeNUm, const int *labels)
{
    std::vector<ObjPatch*> _patches;

    const int NPlane = planeNUm;
    const int NPoint = pointList.size();

    std::vector<vcg::Point4f> optPlanes;
    std::vector<std::vector<int>> planeVerList;
    planeVerList.resize(NPlane + 1);
    for (int i = 0; i < NPoint; ++i) {
        int label = labels[i];
        planeVerList[label].push_back(i);
    }
    for (int k = 1; k < NPlane + 1; k++) {
        if (planeVerList[k].empty())
            continue;
        vcg::Plane3f plane;
        double err = FinePlane(pointList, planeVerList[k], plane);
        ObjRect *onePatch = ExtractMBR(mesh, plane, pointList, indexList, planeVerList[k]);
        if (onePatch != 0) {
            onePatch->m_varN = err;
            _patches.push_back(onePatch);
        }
    }
    patches.swap(_patches);
    return patches.size();
}
