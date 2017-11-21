#include "PointCloudFitUtil.h"

#include <random>

#include "PCA.h"
#include <wrap/io_trimesh/io_mask.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/space/index/kdtree/kdtree.h>
#include <vcg/space/fitting3.h>

// [1] Detect Planes
int HoughPlaneOne(
    vcg::Point4f &Plane,
    const FixedAxis fix,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const double _intercept, const double _a, const double _s)
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

    flog(
        "        | [#Time-%7.4f]-[%d-Seted] : %7.4f %7.4f %7.4f %7.4f  #nPts-< %d >\n",
        time.elapsed() / 1000.0,
        fix == FixedAxis_X ? 'X' : (fix == FixedAxis_Y ? 'Y' : 'Z'),
        Plane.X(), Plane.Y(), Plane.Z(), Plane.W(),
        vote);

    return vote;
}
int HoughPlane(
    vcg::Point4f &Plane,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const double _intercept, const double _a, const double _s)
{//¼ì²â + ÄâºÏÆ½Ãæ
    QTime time;
    time.start();
    flog(
        "      [--Plane_HT--]: #nPts-%d \n",
        PointList.size());

    vcg::Point4f P[3];
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

    flog(
        "        | [+][Checked] : [%d] #nPts-< %d > \n"
        "      [--Plane_HT--]: Done in %.4f seconds \n",
        retIdx + 1, N[retIdx], time.elapsed() / 1000.0);

    Plane = P[retIdx];
    return N[retIdx];
}


// [2] Surface Points Verification
int AttachToPlane(
    std::vector<int> &PtOnPlane,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const vcg::Point4f &plane,
    const double _TDis, const double _TAng)
{
    QTime time;
    time.start();

    const bool bHasNorm = (NormList.size() > 0) ? true : false;
    if (bHasNorm)
        assert(PointList.size() == NormList.size());

    // Check Distance
    const double TDis = _TDis * sqrt(plane.V(0)*plane.V(0) + plane.V(1)*plane.V(1) + plane.V(2)*plane.V(2));
    std::vector<int> OnPlaneList;
    for (int i = 0; i<PointList.size(); i++) {
        const vcg::Point3f &pt3 = PointList[i];
        vcg::Point4f pt4(pt3.X(), pt3.Y(), pt3.Z(), 1.0);
        if (abs(pt4*plane)<TDis)
            OnPlaneList.push_back(i);
    }

    // Check Norm
    if (bHasNorm && !OnPlaneList.empty())
    {
        vcg::Point3f NP = vcg::Point3f(plane.V(0), plane.V(1), plane.V(2));
        NP.Normalize();
        std::vector<int> ReMoved;
        for (int i = 0; i<OnPlaneList.size(); ++i) {
            vcg::Point3f npt = NormList.at(OnPlaneList.at(i));
            npt.Normalize();
            double ang = 90 - abs(90 - vcg::AngleN(NP, npt)*_R2D);
            if (ang>_TAng)
                ReMoved.push_back(i);
        }

        if (!ReMoved.empty()) {
            for (int i = ReMoved.size() - 1; i >= 0; --i)
                OnPlaneList.erase(OnPlaneList.begin() + ReMoved.at(i));
        }
        ReMoved.clear();
    }

    flog(
        "      [--Attach--]: Attach points to planes...\n"
        "        | #Threshold_Dis : %.4f\n"
        "        | #Threshold_Ang : %.4f-[%d]\n"
        "        | #nPts-OnPlane  : %d \n"
        "      [--Attach--]: Done in %.4f seconds. \n",
        TDis, _TAng, bHasNorm,
        OnPlaneList.size(), time.elapsed() / 1000.0);

    PtOnPlane.swap(OnPlaneList);
    return PtOnPlane.size();
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

    std::vector<int> tempLoc;
    std::vector<int> LocFinal;

    std::vector<vcg::Point3f> seedStack;
    vcg::Point3f seed;
    std::vector<int> recorder;
    const int NPtsAll = index.size();
    while (!index.empty())
    {
        // Each Cluster
        seedStack.clear();
        seedStack.push_back(PointList.at(index.at(0)));
        tempLoc.clear();
        tempLoc.push_back(index.at(0));

        index.erase(index.begin());
        // Growing
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
                    tempLoc.push_back(p);
                    recorder.push_back(i);
                }
            }
            for (int i = recorder.size() - 1; i >= 0; i--)
                index.erase(index.begin() + recorder.at(i));
        }
        if (tempLoc.size() > LocFinal.size())
            LocFinal.swap(tempLoc);
    }
    for (auto iter : LocFinal)
        index.push_back(iter);
    LocFinal.clear();
    tempLoc.clear();
    recorder.clear();
    std::sort(index.begin(), index.end());

    flog(
        "      [--MaxRegion--]: #Pts-%d\n"
        "        | #Threshold_Dis : %.4f\n"
        "        | #NPts_MaxReg   : %d\n"
        "      [--MaxRegion--]: Done in %.4f seconds. \n",
        NPtsAll, _2TDis, index.size(), time.elapsed() / 1000.0);
}


// [4] LS Fit
double FineFit(
    const std::vector<vcg::Point3f> &pointList,
    const std::vector<int> &planeVerList,
    vcg::Point4f &plane)
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
    plane.X() = ple.Direction().X();
    plane.Y() = ple.Direction().Y();
    plane.Z() = ple.Direction().Z();
    plane.W() = -ple.Offset();
    float fitError = 0;
    for (int i = 0; i < ExactVec.size(); ++i) {
        float d = vcg::SignedDistancePlanePoint(ple, ExactVec[i]);
        fitError += d*d;
    }
    fitError = sqrt(fitError / ExactVec.size());

    flog(
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
    double &ex, double &ey)
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

void ExtractMBR(
    CMeshO &mesh,
    ObjPlane &APlne,
    const vcg::Point4f &Plane,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<int> &IndexList,
    const std::vector<int> &PlaneVerList)
{
    QTime time;
    time.start();

    if (Plane == _NON_PLANE ||
        PlaneVerList.empty())
        return;

    assert(IndexList.size() == PointList.size());

    vcg::Point3f NP = vcg::Point3f(Plane.V(0), Plane.V(1), Plane.V(2));
    if (Plane.V(3)>0) NP = -NP;
    double offset = abs(Plane.V(3)) / NP.Norm();
    NP.Normalize();
    APlne.m_N = NP;

    vcg::Plane3f VCGPLane = vcg::Plane3f();
    VCGPLane.SetDirection(NP);
    VCGPLane.SetOffset(offset);
    CMeshO::PerVertexAttributeHandle<PtType> type_hi =
        vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, _MyPtAttri);
    std::vector<vcg::Point3f> OnPPt;

    const int PlaneCode = APlne.m_PlaneIndex;
    for (int i = PlaneVerList.size() - 1; i >= 0; --i)
    {
        int index = PlaneVerList.at(i);
        type_hi[IndexList.at(index)] = PlaneCode;

        vcg::Point3f pp = VCGPLane.Projection(PointList.at(index));
        OnPPt.push_back(pp);
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
    double *PC, V[3];//V[3] useless
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
    for (int i = -45; i < 45; i += 1) { // Not too slow
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

    APlne.m_pO = O;
    APlne.m_dX = DX;
    APlne.m_dY = DY;
    APlne.m_sizeConfidence = vcg::Point3f(confidenceX, confidenceY, 0);
    //
    flog(
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

std::vector<vcg::Point4f> DetectHTPlanes(
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const double _intercept, const double _a, const double _s,
    const double TDis, const double TAng,
    const double TNPtsRatio, const int TNPtsHard,
    const int ExpPlaneNum,
    std::vector<double> *errors,
    std::vector<ObjPlane*> *pPlanes,
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

    std::vector<vcg::Point4f> planeVec;
    int planeNum = 0;
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

        vcg::Point4f plane;
        // HT Detection            
        int NP = HoughPlane(plane, pointList, normList, _intercept, _a, _s); // Center At (0,0,0)
        if (plane == _NON_PLANE || NP <= _planeNT)
            break;

        // Surface Points Verification
        std::vector<int> planeVerList;
        AttachToPlane(planeVerList, pointList, normList, plane, TDis, TAng);

        // Coplanar Separation
        PicMaxRegion(pointList, planeVerList, TDis*2.0);
        if (planeVerList.size() <= _planeNT)
            break;

        // Least Squre Fit
        double err = FineFit(pointList, planeVerList, plane);
        // planeVerList = AttachToPlane(pointList,directionList,plane,_planeDThreshold); // AGAIN!
        if (bReportErr)
            errors->push_back(err);

        if (bRetPlane) {
            ObjPlane *onePle = new ObjPlane(Pt_OnPlane + planeNum);
            onePle->m_varN = err;

            // Extract the Minimum-Bounding-Rectangle
            ExtractMBR(*pMesh, *onePle, plane, pointList, indexList, planeVerList);

            pPlanes->push_back(onePle);
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

        planeVec.push_back(plane);
        planeNum++;
    }
    if (bRetPlane)
        indexList.clear();
    pointList.clear();
    normList.clear();

    return planeVec;
}

std::vector<ObjPlane*> ExtractPlanes(
    CMeshO &mesh,
    const std::vector<int> &indexList,
    const std::vector<vcg::Point3f> &pointList,
    const int planeNUm, const int *labels)
{
    std::vector<ObjPlane*> planes;

    const int NPlane = planeNUm;
    const int NPoint = pointList.size();

    std::vector<vcg::Point4f> optPlanes;
    std::vector<std::vector<int>> planeVerList;
    planeVerList.resize(NPlane + 1);
    for (int i = 0; i < NPoint; ++i) {
        int label = labels[i];
        planeVerList[label].push_back(i);
    }
    int PlaneID = Pt_OnPlane;
    for (int k = 1; k < NPlane + 1; k++) {
        if (planeVerList[k].empty())
            continue;
        ObjPlane *onePle = new ObjPlane(PlaneID++);

        vcg::Point4f plane;
        double err = FineFit(pointList, planeVerList[k], plane);
        onePle->m_varN = err;
        ExtractMBR(mesh, *onePle, plane, pointList, indexList, planeVerList[k]);

        planes.push_back(onePle);
    }
    return planes;
}

// Plane Fitting with GCO
// MPFGCO : Multi-Plane Fitting with GCO
MPFGCOCost MPFGCOGeneratCost(
    const std::vector<vcg::Point4f> &planes,
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
        vcg::Point4f ple = planes[i - 1];
        vcg::Point3f n(ple.X(), ple.Y(), ple.Z());
        double w = ple.W() / n.Norm();
        n.normalized();
        // Dis
        for (int j = 0; j < NPts; ++j) {
            double d = abs(points[j] * n + w);
            DataCost[j*NLabel + i] = int(d / unit_a + 0.5);
        }
        // Ang
        if (bHasNorm) {
            for (int j = 0; j < NPts; ++j) {
                double ang = 90 - abs(90 - vcg::AngleN(n, norms[j])*_R2D);
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
    gcoCost.numLabel   = NLabel;
    gcoCost.numSite    = NPts;
    gcoCost.dataCost   = DataCost;
    gcoCost.smoothCost = SmoothCost;
    gcoCost.labelCost  = cost_label;

#if defined(_ReportOut_)
    reportMat<int>(DataCost, NPts, NLabel, "../~CostData~.txt");
    reportMat<int>(SmoothCost, NLabel, NLabel, "../~CostSmooth~.txt");
#endif

    return gcoCost;
}

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
    const double _r = -1.0 / 2*(unit_a*unit_a)*(delta*delta);
    CMeshO::VertexIterator vi = mesh.vert.begin();
    for (int i = 0; i < NPts; ++i) {
        vcg::Point3f p = (vi + ptIndex[i])->cP();
        KDTree.doQueryK(p, knnNum, queue);
        int neiNumQuery = queue.getNofElements();
        int neiNum = 0;
        for (int k = 0; k < neiNumQuery; k++) {
            int neightId = queue.getIndex(k);
            if ((vi + neightId)->IsD())
                continue;
            // Vpq(lp,lq) = w_{p,q}*d_{lp,lq}
            // d_{lp,lq} = 1 if lp != lq, otherwise 0; 
            // (see implement of smooth cost in function MPFGCOGeneratCost())
            // w_{p,q} = lambda * exp{ - (||p-q||_2/a) ^ 2 / 2*delta^2}
            double w_pq = vcg::SquaredDistance(p, (vi+neightId)->cP());
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

std::vector<double> GCOReEstimat(
    std::vector<vcg::Point4f> &planes,
    const std::vector<vcg::Point3f> &pointList,
    const int *labels)
{
    std::vector<double> errors;

    const int NPlane = planes.size();
    const int NPoint = pointList.size();

    std::vector<vcg::Point4f> optPlanes;
    std::vector<std::vector<int>> planeVerList;
    planeVerList.resize(NPlane+1);
    for (int i = 0; i < NPoint; ++i) {
        int label = labels[i];
        planeVerList[label].push_back(i);
    }
    for (int k = 1; k < NPlane+1; k++) {
        if (planeVerList[k].empty())
            continue;

        vcg::Point4f plane;
        double err = FineFit(pointList, planeVerList[k], plane);
        errors.push_back(err);
        optPlanes.push_back(plane);
    }

    planes.swap(optPlanes);

    return errors;
}