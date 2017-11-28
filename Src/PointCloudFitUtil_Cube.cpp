#include "PointCloudFitUtil.h"
#include "tool/IOU.h"


bool IsParallel(const vcg::Point3f &L1, const vcg::Point3f &L2, const double AngTD)
{
    return CheckAng00(VCGAngle(L1, L2)) < abs(AngTD);
}
bool IsPerpendicular(const vcg::Point3f &L1, const vcg::Point3f &L2, const double AngTD)
{
    return CheckAng90(VCGAngle(L1, L2)) < abs(AngTD);
}
double PlaneIoU(const ObjRect *P1, const ObjRect *P2)
{
    double oX = (P2->m_O - P1->m_O) * P1->m_AX / P1->m_AX.SquaredNorm();
    double oY = (P2->m_O - P1->m_O) * P1->m_AY / P1->m_AY.SquaredNorm();

    double xX = P2->m_AX * P1->m_AX / P1->m_AX.SquaredNorm();
    double xY = P2->m_AX * P1->m_AY / P1->m_AY.SquaredNorm();

    double yX = P2->m_AY * P1->m_AX / P1->m_AX.SquaredNorm();
    double yY = P2->m_AY * P1->m_AY / P1->m_AY.SquaredNorm();

    IOU::Quad quad1(
        IOU::Point(0, 0),
        IOU::Point(0, 1),
        IOU::Point(1, 1),
        IOU::Point(1, 0));
    IOU::Quad quad2(
        IOU::Point(oX, oY),
        IOU::Point(oX + yX, oY + yY),
        IOU::Point(oX + yX + xX, oY + yY + xY),
        IOU::Point(oX + xX, oY + xY));
    if (quad2.isAntiClockWiseConvex())
        quad2.flipBack();

    double interArea = CalInterArea(quad1, quad2);

    double s1 = quad1.area();
    double s2 = quad2.area();
    assert(interArea >= 0.0);
    assert(s1 > 0.0 && s2 > 0.0);
    assert(interArea <= s1 && interArea <= s2);

    return interArea / (s1 + s2 - interArea);

}
void RemovePlanes(
    std::vector<ObjRect*> &planes,
    const std::vector<ObjRect*> &planesTobeRemoved,
    const bool memRelease)
{
    for (int i = 0; i < planesTobeRemoved.size(); ++i) {
        ObjRect* toBeRomeved = planesTobeRemoved.at(i);
        for (int j = 0; j < planes.size(); ++j) {
            if (planes.at(j) == toBeRomeved) {
                planes.erase(planes.begin() + j);
                if (memRelease)
                    delete toBeRomeved;
                break;
            }
        }
    }
}


// [1] Infer Cube Faces
bool IsOppoFaces(
    const ObjRect *P1, const ObjRect *P2,
    const double TAng, const double TIoU)
{
    // 1. Parallel Normals
    if (!IsParallel(P1->m_N, P2->m_N, TAng))
        return false;

    // 2. Parallel Edges
    if (!( 
        ( IsParallel(P1->m_AX, P2->m_AX, TAng) && IsParallel(P1->m_AY, P2->m_AY, TAng) ) ||
        ( IsParallel(P1->m_AX, P2->m_AY, TAng) && IsParallel(P1->m_AY, P2->m_AX, TAng))
        ) )
        return false;

    // 3. Face to Face
    double iou12 = PlaneIoU(P1, P2);
    double iou21 = PlaneIoU(P2, P1);
    if ( (iou12 + iou21)/2.0 < TIoU)
        return false;

    // All Passed
    return true;
}
bool IsAdjacencyFaces(
    const ObjRect *P1, const ObjRect *P2,
    const double TRDis, const double TAng,
    PlaneRelation *adjType)
{
    PlaneRelation _adjType = PlaneRelation_NoRetation;
    // 1. Perpendicular Normals
    if (!IsPerpendicular(P1->m_N, P2->m_N, TAng))
        return false;

    // 2. Parallel & Perpendicular Edges
    double angXX = VCGAngle(P1->m_AX, P2->m_AX);
    double angXY = VCGAngle(P1->m_AX, P2->m_AY);
    double angYY = VCGAngle(P1->m_AY, P2->m_AY);
    double angYX = VCGAngle(P1->m_AY, P2->m_AX);
    bool solu1 = 
        (CheckAng90(angXX) < TAng && CheckAng90(angXY) < TAng) && // X1 -| X2 and Y2
        (CheckAng00(angYY) < TAng || CheckAng00(angYX) < TAng);             // Y1 // X2 or Y2
    bool solu2 =
        (CheckAng90(angYX) < TAng && CheckAng90(angYY) < TAng) && // Y1 -| X2 and Y2
        (CheckAng00(angXY) < TAng || CheckAng00(angXX) < TAng);             // X1 // X2 or Y2
    if (!solu1 && !solu2)
        return false;

    // 3. Adjacency Faces
    vcg::Point3f O1 = P1->m_O + (P1->m_AX + P1->m_AY) / 2.0;
    vcg::Point3f O2 = P2->m_O + (P2->m_AX + P2->m_AY) / 2.0;
    vcg::Point3f O2O = O2 - O1;
    vcg::Point3f O2O_1 = O2O - (P1->m_N * (O2O*P1->m_N)) / P1->m_N.SquaredNorm();
    vcg::Point3f O2O_2 = O2O - (P2->m_N * (O2O*P2->m_N)) / P2->m_N.SquaredNorm();
    vcg::Point3f NN = P1->m_N^P2->m_N;
    if (!IsPerpendicular(O2O_1, NN, 30.0) ||
        !IsPerpendicular(O2O_2, NN, 30.0))
        return false;

    double Dis1 = abs(O2O*P1->m_N);
    double Dis2 = abs(O2O*P2->m_N);
    double l1, u1;
    double l2, u2;
    if (CheckAng90(VCGAngle(P1->m_AX, NN)) < CheckAng90(VCGAngle(P1->m_AY, NN))) {
        l1 = P1->width();
        u1 = P1->height();
    } else {
        u1 = P1->width();
        l1 = P1->height();
    }
    if (CheckAng90(VCGAngle(P2->m_AX, NN)) < CheckAng90(VCGAngle(P2->m_AY, NN))) {
        l2 = P2->width();
        u2 = P2->height();
    } else {
        u2 = P2->width();
        l2 = P2->height();
    }
    if (abs(u1 - u2) / (u1 + u2) > 0.5)
        return false;

    if (
        ((Dis1 - l2*0.5) / l2 > TRDis || (l2*0.5 - Dis1) / l2 > TRDis / 2.0) ||
        ((Dis2 - l1*0.5) / l1 > TRDis || (l1*0.5 - Dis2) / l1 > TRDis / 2.0)
        )
        return false;

    // All Passed
    if (adjType != 0) {
        vcg::Point3f O2O_proj = O2O - (P1->m_N * (O2O*P1->m_N)) / P1->m_N.SquaredNorm();
        if (CheckAng90(VCGAngle(O2O_proj, P1->m_AX)) < CheckAng90(VCGAngle(O2O_proj, P1->m_AY))) {
            // U or B
            if (O2O*P1->m_AY > 0)
                _adjType = PlaneRelation_Adjacency_B;
            else
                _adjType = PlaneRelation_Adjacency_U;
        }
        else {
            // R or L
            if (O2O*P1->m_AX > 0)
                _adjType = PlaneRelation_Adjacency_R;
            else
                _adjType = PlaneRelation_Adjacency_L;
        }
        *adjType = _adjType;
    }

    return true;
}
PlaneRelation EstPlaneRelation(
    const ObjRect *P1, const ObjRect *P2,
    const double TRDis, const double TAng, const double TIoU)
{
    PlaneRelation  relation = PlaneRelation_NoRetation;

    if (IsOppoFaces(P1, P2, TAng, TIoU))
        return PlaneRelation_AtOppo;
    else if (IsAdjacencyFaces(P1, P2, TRDis, TAng, &relation))
        return relation;

    return PlaneRelation_NoRetation;
}

bool BuildBox(
    ObjRect* cubeFace[6], ObjRect* Rect,
    const double TRDis, const double TAng, const double TIoU)
{
    //  [0]       [1]       [2]          [3]         [4]        [5]
    //   |         |         |            |           |          |
    //RefPlane  UpPlane  RightPlane  BottomPlane  LeftPlane  OppoPlane
    assert(
        cubeFace[0] != 0 &&
        Rect != 0
    );

    PlaneRelation P[6];
    for (int i = 0; i < 6; ++i) {
        if (cubeFace[i] == 0)
            P[i] = PlaneRelation_NoRetation;
        else {
            PlaneRelation Pi = EstPlaneRelation(cubeFace[i], Rect, TRDis, TAng, TIoU);
            if (Pi == PlaneRelation_NoRetation)
                return false;
            P[i] = Pi;
        }
    }

    int inloc = -1;
    if (P[0] == PlaneRelation_Adjacency_U &&
        P[1] == PlaneRelation_NoRetation &&
        P[2] != PlaneRelation_AtOppo &&
        (P[3] & 0x01) == 0 &&
        P[4] != PlaneRelation_AtOppo &&
        P[5] != PlaneRelation_AtOppo
        )
        inloc = 1;
    else if (P[0] == PlaneRelation_Adjacency_R &&
        P[1] != PlaneRelation_AtOppo &&
        P[2] == PlaneRelation_NoRetation &&
        P[3] != PlaneRelation_AtOppo &&
        (P[4] & 0x01) == 0 &&
        P[5] != PlaneRelation_AtOppo
        )
        inloc = 2;
    else if (P[0] == PlaneRelation_Adjacency_B &&
        (P[1] & 0x01) == 0 &&
        P[2] != PlaneRelation_AtOppo &&
        P[3] == PlaneRelation_NoRetation &&
        P[4] != PlaneRelation_AtOppo &&
        P[5] != PlaneRelation_AtOppo
        )
        inloc = 3;
    else if (P[0] == PlaneRelation_Adjacency_L &&
        P[1] != PlaneRelation_AtOppo &&
        (P[2] & 0x01) == 0 &&
        P[3] != PlaneRelation_AtOppo &&
        P[4] == PlaneRelation_NoRetation &&
        P[5] != PlaneRelation_AtOppo
        )
        inloc = 4;
    else if (P[0] == PlaneRelation_AtOppo &&
        P[1] != PlaneRelation_AtOppo &&
        P[2] != PlaneRelation_AtOppo &&
        P[3] != PlaneRelation_AtOppo &&
        P[4] != PlaneRelation_AtOppo &&
        P[5] == PlaneRelation_NoRetation
        )
        inloc = 5;
    else
        return false;
    // Check same side
    if (inloc == 5) {// Add an oppo
        for (int i = 1; i <= 4; ++i) {
            if (cubeFace[i] != 0 &&
                P[i] == EstPlaneRelation(cubeFace[i], cubeFace[0], TRDis, TAng, TIoU))
                return false;
        }
    }
    if (cubeFace[5]) {// Already an oppo
        if (EstPlaneRelation(Rect, cubeFace[0], TRDis, TAng, TIoU)
            == EstPlaneRelation(Rect, cubeFace[5], TRDis, TAng, TIoU))
            return false;
    }

    cubeFace[inloc] = Rect;
    return true;
}
std::vector<ObjRect*> CubeFaceInferringOne(
    const std::vector<ObjRect*> &Rects,
    const double TRDis, const double TAng, const double TIoU)
{
    ObjRect* tempList[6];
    std::vector<ObjRect*> finalOut;
    finalOut.reserve(6);
    for (int i = 0; i < Rects.size(); ++i) {
        // Clean
        for (int j = 0; j < 6; ++j)
            tempList[j] = 0;
        // Init
        tempList[0] = Rects.at(i);
        int pNum = 1;
        for (int j = 0; j < Rects.size(); j++) {
            if (j == i)
                continue;
            if ((i == 0 && j == 2) ||
                (i == 0 && j == 3))
                int stopMe = 0;
            if (BuildBox(tempList, Rects.at(j), TRDis, TAng, TIoU)) {
                pNum++;
                if (pNum == 6)
                    break;
            }
        }

        // Less than 2 patches
        if (pNum < 2)
            continue;
        // Have only 2 patches, but not in opposite relation.
        if (pNum == 2 &&
            !(
                (tempList[0] != 0 && tempList[5] != 0) ||
                (tempList[1] != 0 && tempList[3] != 0) ||
                (tempList[2] != 0 && tempList[4] != 0))
            )
            continue;

        if (pNum > finalOut.size()) {
            //------
            flog(
                "        |--------------------- \n"
                "        |");
            //------
            finalOut.clear();
            for (int j = 0; j < 6; ++j) {
                if (tempList[j] != 0) {
                    finalOut.push_back(tempList[j]);
                    //------
                    flog(" %d", tempList[j]->m_index);
                    //------
                } else {
                    //------
                    flog(" 000");
                    //------
                }
            }
            //------
            flog("\n");
            //------
            if (finalOut.size() == 6)
                break;
        }
    }
    
    return finalOut;
}
int CubeFaceInferring(
    std::vector< std::vector<ObjRect*> > &CubeFaces,
    std::vector<ObjRect*> &Rects,
    const double TRDis, const double TAng, const double TIoU,
    const bool remove)
{
    QTime time;
    time.start();
    flog(
        "      [--Cube_Infer--]: #Ple-%d\n"
        "        | #TRDis   : %.4f \n"
        "        | #TAng    : %.4f \n",
        Rects.size(),
        TRDis, TAng);

    std::vector<ObjRect*> _rects = Rects;
    std::vector< std::vector<ObjRect*> > _cubefaces;
    while (_rects.size() > 1) {
        std::vector<ObjRect*> faces = CubeFaceInferringOne(_rects, TRDis, TAng, TIoU);
        if (!faces.empty()) {
            _cubefaces.push_back(faces);
            RemovePlanes(_rects, faces, false);
        }
        else
            break;
    }

    if (remove) {
        for (int i=0; i<_cubefaces.size(); ++i)
            RemovePlanes(Rects, _cubefaces.at(i), false);
    }
    CubeFaces.swap(_cubefaces);

    flog(
        "        |===================== \n"
        "        | #Cube    : %d \n"
        "        | #LeftPle : %d \n"
        "      [--Cube_Infer--]: Done in %.4f seconds. \n",
        CubeFaces.size(), _rects.size(),
        time.elapsed() / 1000.0);

    return CubeFaces.size();
}


// [2] Estimate Cube
// [2.1] Estimate Orientation
typedef std::pair< vcg::Point3f, double > WPoint3f;
vcg::Point3f RobustOneD(
    const std::vector< WPoint3f > &FaceD,
    const std::vector< WPoint3f > &EdgeD)
{// ！！！同向单位向量输入！！！
    vcg::Point3f N;
    if (!FaceD.empty()) {
        if (FaceD.size() == 1)
            N = FaceD.at(0).first;
        else {
            double r1 = FaceD.at(0).second * FaceD.at(0).second;
            double r2 = FaceD.at(1).second * FaceD.at(1).second;
            double k1 = r2 / (r1 + r2);
            double k2 = r1 / (r1 + r2);
            N = (FaceD.at(0).first * k1) + (FaceD.at(1).first * k2);
        }
    }
    else {
        N.SetZero();
        double w = 0.0;
        for (int i = 0; i<EdgeD.size(); i++)
        {
            double r = EdgeD.at(i).second * EdgeD.at(i).second;
            N += EdgeD.at(i).first * r;
            w += r;
        }
        N = N / w;
    }
    N.Normalize();
    return N; 
        
    /*
    std::pair< vcg::Point3f, double > retDW;
    if (!FaceD.empty()) {
        if (FaceD.size() == 1)
        {
            retDW.first = FaceD.at(0).first;
            retDW.second = 1.0;
        }
        else
        {
            double r1 = FaceD.at(0).second * FaceD.at(0).second;
            double r2 = FaceD.at(1).second * FaceD.at(1).second;
            retDW.first = ((FaceD.at(0).first * r2) + (FaceD.at(1).first * r1)) / ((r1 + r2) / 2.0);
            retDW.second = 1.0 + FaceD.at(0).first*FaceD.at(1).first;
        }
    }
    else {
        vcg::Point3f d;
        double w;
        d.SetZero();
        w = 0.0;
        for (int i = 0; i<EdgeD.size(); i++)
        {
            d += EdgeD.at(i).first * EdgeD.at(i).second * EdgeD.at(i).second;
            if (EdgeD.at(i).second>w)
                w = EdgeD.at(i).second;
        }
        retDW.first = d;
        retDW.second = w;
    }
    retDW.first.Normalize();
    return retDW;
    */
              
}
void RobustDirections(
    const std::vector< WPoint3f > &FaceNX,
    const std::vector< WPoint3f > &FaceNY,
    const std::vector< WPoint3f > &FaceNZ,
    const std::vector< WPoint3f > &EdgeNX,
    const std::vector< WPoint3f > &EdgeNY,
    const std::vector< WPoint3f > &EdgeNZ,
    vcg::Point3f &NX, vcg::Point3f &NY, vcg::Point3f &NZ)
{
    vcg::Point3f nx = RobustOneD(FaceNX, EdgeNX);
    vcg::Point3f ny = RobustOneD(FaceNY, EdgeNY);
    vcg::Point3f nz = RobustOneD(FaceNZ, EdgeNZ);

    Eigen::Matrix3f A;
    A(0, 0) = nx.X(), A(0, 1) = nx.Y(), A(0, 2) = nx.Z();
    A(1, 0) = ny.X(), A(1, 1) = ny.Y(), A(1, 2) = ny.Z();
    A(2, 0) = nz.X(), A(2, 1) = nz.Y(), A(2, 2) = nz.Z();
    Eigen::JacobiSVD<Eigen::MatrixXf> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::Matrix3f V = svd.matrixV(), U = svd.matrixU();
    Eigen::Matrix3f S = U * Eigen::Matrix3f::Identity() * V.transpose(); // S = U^-1 * A * VT * -1  
    nx = vcg::Point3f(S(0, 0), S(0, 1), S(0, 2));
    ny = vcg::Point3f(S(1, 0), S(1, 1), S(1, 2));
    nz = vcg::Point3f(S(2, 0), S(2, 1), S(2, 2));

    NX = nx;
    NY = ny;
    NZ = nz;

    /*
    std::pair< vcg::Point3f,double > NXW = RobustOneD(FaceNX,EdgeNX);
    std::pair< vcg::Point3f,double > NYW = RobustOneD(FaceNY,EdgeNY);
    std::pair< vcg::Point3f,double > NZW = RobustOneD(FaceNZ,EdgeNZ);

    NX = NXW.first;
    NY = NYW.first;
    NZ = NZW.first;
    if (NXW.second > NYW.second && NXW.second > NZW.second) {
    if (NYW.second > NZW.second)
    {// 1. X>Y>Z
    NZ = NX^NY;
    NY = NZ^NX;
    }
    else
    {// 2. X>Z>Y
    NY = NZ^NX;
    NZ = NX^NY;
    }
    } else if (NYW.second > NZW.second) {
    if (NZW.second > NXW.second)
    {// 3. Y>Z>X
    NX = NY^NZ;
    NZ = NX^NY;
    }
    else
    {// 4. Y>X>Z
    NZ = NX^NY;
    NX = NY^NZ;
    }
    } else {
    if (NYW.second > NXW.second)
    {// 5. Z>Y>X
    NX = NY^NZ;
    NY = NZ^NX;
    }
    else
    {// 6. Z>X>Y
    NY = NZ^NX;
    NX = NY^NZ;
    }
    }
    NX.Normalize();
    NY.Normalize();
    NZ.Normalize();
    */
}
void RobustOrientation(
    const std::vector<ObjRect*> &CubePlanes,
    vcg::Point3f &NX, vcg::Point3f &NY, vcg::Point3f &NZ,
    const double TAng)
{
    QTime time;
    time.start();

    ObjRect *PP = CubePlanes.at(0);
    vcg::Point3f NXRef = PP->m_AX; NXRef.Normalize();
    vcg::Point3f NYRef = PP->m_AY; NYRef.Normalize();
    vcg::Point3f NZRef = PP->m_N; NZRef.Normalize();
    double ZVar = PP->m_varN;
    double XYWeight = std::min(PP->m_EIConfX, PP->m_EIConfY);
    std::vector< WPoint3f > FaceNX, FaceNY, FaceNZ; // 由面法向提供的指向信息
    std::vector< WPoint3f > EdgeNX, EdgeNY, EdgeNZ; // 由边方向提供的指向信息
    FaceNZ.push_back(WPoint3f(NZRef, ZVar));
    EdgeNX.push_back(WPoint3f(NXRef, XYWeight));
    EdgeNY.push_back(WPoint3f(NYRef, XYWeight));

    vcg::Point3f *pX, *pY, *pZ;
    std::vector< WPoint3f > *pXV, *pYV, *pZV;
    for (int i = 1; i<CubePlanes.size(); ++i)
    {
        PP = CubePlanes.at(i);
        vcg::Point3f nx = PP->m_AX;
        vcg::Point3f ny = PP->m_AY;
        vcg::Point3f nz = PP->m_N;
        double zvar = PP->m_varN;
        double xyweight = std::min(PP->m_EIConfX, PP->m_EIConfY);
        if (IsParallel(nz, NZRef, TAng))
        { // Z平行于Zref
            pZ = &NZRef; pZV = &FaceNZ;
            if (IsParallel(nx, NXRef, TAng)) { // X平行于XRef,Y平行于YRef
                pX = &NXRef; pXV = &EdgeNX;
                pY = &NYRef; pYV = &EdgeNY;
            } else { // Y平行于XRef,X平行于YRef
                pX = &NYRef; pXV = &EdgeNY;
                pY = &NXRef; pYV = &EdgeNX;
            }
        }
        else if(IsParallel(nz, NXRef, TAng))
        { // Z平行于XRef
            pZ = &NXRef; pZV = &FaceNX;
            if (IsParallel(nx, NZRef, TAng)) { // X平行于ZRef,Y平行于YRef
                pX = &NZRef; pXV = &EdgeNZ;
                pY = &NYRef; pYV = &EdgeNY;
            } else { // X平行于YRef,Y平行于ZRef
                pX = &NYRef; pXV = &EdgeNY;
                pY = &NZRef; pYV = &EdgeNZ;
            }
        }
        else {
            // Z平行于YRef
            pZ = &NYRef; pZV = &FaceNY;
            if (IsParallel(nx, NXRef, TAng)) { // X平行于XRef,Y平行于ZRef
                pX = &NXRef; pXV = &EdgeNX;
                pY = &NZRef; pYV = &EdgeNZ;
            } else { // X平行于ZRef,Y平行于XRef
                pX = &NZRef; pXV = &EdgeNZ;
                pY = &NXRef; pYV = &EdgeNX;
            }
        }

        if (*pX*nx<0) nx = -nx;
        if (*pY*ny<0) ny = -ny;
        if (*pZ*nz<0) nz = -nz;
        nx.Normalize();
        ny.Normalize();
        nz.Normalize();
        pZV->push_back(WPoint3f(nz, zvar));
        pXV->push_back(WPoint3f(nx, xyweight));
        pYV->push_back(WPoint3f(ny, xyweight));
    }
    assert((FaceNX.size() + EdgeNX.size()) == CubePlanes.size());
    assert((FaceNY.size() + EdgeNY.size()) == CubePlanes.size());
    assert((FaceNZ.size() + EdgeNZ.size()) == CubePlanes.size());
    RobustDirections(FaceNX, FaceNY, FaceNZ, EdgeNX, EdgeNY, EdgeNZ, NX, NY, NZ);

    flog(
        "      [--Cube_Orientation--]: #Ple-%d\n"
        "        | #TAng  : %.4f \n"
        "        | #NX    : < %.4f, %.4f, %.4f > \n"
        "        | #NY    : < %.4f, %.4f, %.4f > \n"
        "        | #NZ    : < %.4f, %.4f, %.4f > \n"
        "      [--Cube_Orientation--]: Done in %.4f seconds. \n",
        CubePlanes.size(), TAng,
        NX.X(), NX.Y(), NX.Z(),
        NY.X(), NY.Y(), NY.Z(),
        NZ.X(), NZ.Y(), NZ.Z(),
        time.elapsed() / 1000.0);
}


// [2.2] Estimate Dimenstion
typedef std::pair< double, double > WLoc;
typedef std::pair< WLoc, WLoc > WLocPair;
WLoc FaceInter(
    const vcg::Point3f &NRef,
    const vcg::Point3f &O,
    const vcg::Point3f &NX,
    const vcg::Point3f &NY)
{
    double sum = 0.0;
    double sum2 = 0.0;
    double loc = 0.0;
    loc = NRef*O;               sum += loc;	sum2 += loc*loc;
    loc = NRef*(O + NX);		sum += loc;	sum2 += loc*loc;
    loc = NRef*(O + NY);		sum += loc;	sum2 += loc*loc;
    loc = NRef*(O + NX + NY);	sum += loc;	sum2 += loc*loc;
    double locMiu = sum / 4.0;
    double locVar = sqrt(sum2 / 4.0 - locMiu*locMiu);
    return WLoc(locMiu, locVar);
}
WLoc LineProj(
    const vcg::Point3f &NRef,
    const vcg::Point3f &O,
    const vcg::Point3f &NMain,
    const vcg::Point3f &NSecned)
{
    double loc11 = NRef*O;
    double loc12 = NRef*(O + NMain);
    double loc21 = NRef*(O + NSecned);
    double loc22 = NRef*(O + NMain + NSecned);

    double loc1 = (loc11 + loc21) / 2.0;
    double loc2 = (loc12 + loc22) / 2.0;
    if (loc1>loc2)
        return WLoc(loc1, loc2);
    else
        return WLoc(loc2, loc1);
}
double RobustLengthOneD(
    const std::vector< WLoc > &FaceInter,
    std::vector< WLocPair > &EdgeInter,
    double &minProj, double &maxProj)
{
    double w = 0.0;

    if (FaceInter.size() >= 2) {
        if (FaceInter.at(0).first<FaceInter.at(1).first) {
            minProj = FaceInter.at(0).first;
            maxProj = FaceInter.at(1).first;
        }
        else {
            minProj = FaceInter.at(1).first;
            maxProj = FaceInter.at(0).first;
        }
        double ll = maxProj - minProj;
        w = (1.0 - FaceInter.at(0).second / ll) * (1.0 - FaceInter.at(1).second / ll);
    }
    else {
        int N = EdgeInter.size();
        std::vector<WLoc> upList, downList;
        for (int i = 0; i<N; i++) {
            upList.push_back(EdgeInter.at(i).first);
            downList.push_back(EdgeInter.at(i).second);
        }

        double up, upw, down, downw;
        //up = upw = down = downw = 0.0;
        //for (int i = 0; i < N; ++i) {
        //	up += upList.at(i).first;
        //	upw += upList.at(i).second;
        //	down += downList.at(i).first;
        //	downw += downList.at(i).second;
        //}
        //up = up / N;
        //upw = upw / N;
        //down = down / N;
        //downw = downw / N;

        std::sort(upList.begin(), upList.end());
        std::sort(downList.begin(), downList.end());
        if (N % 2 == 0) {
            int index1 = N / 2 - 1;
            int index2 = N / 2;
            up = (upList.at(index1).first + upList.at(index2).first) / 2.0;
            upw = (upList.at(index1).second + upList.at(index2).second) / 2.0;
            down = (downList.at(index1).first + downList.at(index2).first) / 2.0;
            downw = (downList.at(index1).second + downList.at(index2).second) / 2.0;
        }
        else {
            int index = N / 2;
            up = upList.at(index).first;
            upw = upList.at(index).second;
            down = downList.at(index).first;
            downw = downList.at(index).second;
        }

        if (FaceInter.size() != 0) {
            double fi = FaceInter.at(0).first;
            double var = FaceInter.at(0).second;
            if (abs(fi - up)<abs(fi - down)) {
                up = fi;
                upw = 1.0 - var / (up - down);
            }
            else {
                down = fi;
                downw = 1.0 - var / (up - down);
            }
        }
        minProj = down;
        maxProj = up;
        w = upw*downw;
    }
    return w;
}
void RobustDimention(
    const std::vector<ObjRect*> &CubePlanes,
    const vcg::Point3f &NX, const vcg::Point3f &NY, const vcg::Point3f &NZ,
    vcg::Point3f &OP, vcg::Point3f &SIZE, vcg::Point3f &WEIGHTS,
    const double TAng)
{
    QTime time;
    time.start();

    // <位置,权值>
    std::vector< WLoc > FaceInterX, FaceInterY, FaceInterZ;
    std::vector< WLocPair > EdgeInterX, EdgeInterY, EdgeInterZ;

    const vcg::Point3f *pNX, *pNY, *pNZ;
    std::vector< WLoc > *pFZ;
    std::vector< WLocPair > *pEX, *pEY;
    for (int i = 0; i<CubePlanes.size(); ++i)
    {
        ObjRect *PP = CubePlanes.at(i);
        vcg::Point3f op = PP->m_O;
        vcg::Point3f nx = PP->m_AX;
        vcg::Point3f ny = PP->m_AY;
        vcg::Point3f nz = PP->m_N;
        double xyweight = std::min(PP->m_EIConfX, PP->m_EIConfY);
        if (IsParallel(nz, NX, TAng))
        { // YZ平面平行
            pNZ = &NX; pFZ = &FaceInterX;
            if (IsParallel(nx, NY, TAng)) { // x//Y y//Z
                pNX = &NY; pEX = &EdgeInterY;
                pNY = &NZ; pEY = &EdgeInterZ;
            } else {// x//Z y//Y
                pNX = &NZ; pEX = &EdgeInterZ;
                pNY = &NY; pEY = &EdgeInterY;
            }
        }
        else if (IsParallel(nz, NY, TAng))
        { // XZ平面平行
            pNZ = &NY; pFZ = &FaceInterY;
            if (IsParallel(nx, NX, TAng)) { // x//X y//Z
                pNX = &NX; pEX = &EdgeInterX;
                pNY = &NZ; pEY = &EdgeInterZ;
            } else {// x//Z y//X
                pNX = &NZ; pEX = &EdgeInterZ;
                pNY = &NX; pEY = &EdgeInterX;
            }
        }
        else
        { // XY平面平行
            pNZ = &NZ; pFZ = &FaceInterZ;
            if (IsParallel(nx, NX, TAng)) { // x//X y//Y
                pNX = &NX; pEX = &EdgeInterX;
                pNY = &NY; pEY = &EdgeInterY;
            } else {// x//Y y//X
                pNX = &NY; pEX = &EdgeInterY;
                pNY = &NX; pEY = &EdgeInterX;
            }
        }
        WLoc IF = FaceInter(*pNZ, op, nx, ny);
        pFZ->push_back(IF);
        WLoc IX = LineProj(*pNX, op, nx, ny);
        pEX->push_back( WLocPair( WLoc(IX.first, xyweight), WLoc(IX.second, xyweight) ) );
        WLoc IY = LineProj(*pNY, op, ny, nx);
        pEY->push_back( WLocPair( WLoc(IY.first, xyweight), WLoc(IY.second, xyweight) ) );

    }
    assert((FaceInterX.size() + EdgeInterX.size()) == CubePlanes.size());
    assert((FaceInterY.size() + EdgeInterY.size()) == CubePlanes.size());
    assert((FaceInterZ.size() + EdgeInterZ.size()) == CubePlanes.size());

    double minx, miny, minz, maxx, maxy, maxz;
    double wx = RobustLengthOneD(FaceInterX, EdgeInterX, minx, maxx);
    double wy = RobustLengthOneD(FaceInterY, EdgeInterY, miny, maxy);
    double wz = RobustLengthOneD(FaceInterZ, EdgeInterZ, minz, maxz);
    SIZE = vcg::Point3f(maxx - minx, maxy - miny, maxz - minz);
    WEIGHTS = vcg::Point3f(wx, wy, wz);
    OP = NX*minx + NY*miny + NZ*minz;

    flog(
        "      [--Cube_Dimention--]: #Ple-%d\n"
        "        | #TAng  : %.4f \n"
        "        | #OP    : < %.4f, %.4f, %.4f > \n"
        "        | #NX    : %.4f - [%.4f%] - < %.4f, %.4f, %.4f > \n"
        "        | #NY    : %.4f - [%.4f%] - < %.4f, %.4f, %.4f > \n"
        "        | #NZ    : %.4f - [%.4f%] - < %.4f, %.4f, %.4f > \n"
        "      [--Cube_Dimention--]: Done in %.4f seconds. \n",
        CubePlanes.size(), TAng,
        OP.X(), OP.Y(), OP.Z(),
        SIZE.X(), WEIGHTS.X(), NX.X(), NX.Y(), NX.Z(),
        SIZE.Y(), WEIGHTS.Y(), NY.X(), NY.Y(), NY.Z(),
        SIZE.Z(), WEIGHTS.Z(), NZ.X(), NZ.Y(), NZ.Z(),
        time.elapsed() / 1000.0);
}

void CubeMeasure(    
    const std::vector<ObjRect*> &CubePlanes,
    ObjCube *cube, const double TAng)
{
    // 1.Direction
    vcg::Point3f NX, NY, NZ;
    RobustOrientation(CubePlanes, NX, NY, NZ, TAng);
    // 2.Dimention
    vcg::Point3f OP, size, sizeweights;
    RobustDimention(CubePlanes, NX, NY, NZ, OP, size, sizeweights, TAng);

    cube->m_AX = NX*size.X();
    cube->m_AY = NY*size.Y();
    cube->m_AZ = NZ*size.Z();
    cube->m_O = OP;
    cube->m_EIConfX = sizeweights.X();
    cube->m_EIConfY = sizeweights.Y();
    cube->m_EIConfZ = sizeweights.Z();
}

// [4] Merge And Cut
bool MergeToCube(
    CMeshO &mesh,
    const ObjRect *rect, const ObjCube *Cube,
    std::vector<ObjRect*> &planeSplit,
    const double TAng, const double TDis)
{
    if (rect == 0 || Cube == 0)
        return false;

    vcg::Point3f co = Cube->m_O;
    vcg::Point3f cx = Cube->m_AX;
    vcg::Point3f cy = Cube->m_AY;
    vcg::Point3f cz = Cube->m_AZ;

    vcg::Point3f po = rect->m_O;
    vcg::Point3f pn = rect->m_N;
    vcg::Point3f px = rect->m_AX;
    vcg::Point3f py = rect->m_AY;
    vcg::Point3f p_center = po + (px + py) / 2.0;

    vcg::Point3f _o = co;
    vcg::Point3f _n = cz;
    vcg::Point3f _u = cx;
    vcg::Point3f _v = cy;
    // 1. Check Close Parallel
    double _s = 0.0;
    if (IsParallel(pn, cx,TAng)) { // in Y_Z
        _s = Cube->DimX();
        _n = cx; _u = cy; _v = cz;       
    }
    else if (IsParallel(pn, cy, TAng)) { // in X_Z
        _s = Cube->DimY();
        _n = cy; _u = cz; _v = cx;
    }
    else if (IsParallel(pn, cz, TAng)) { // in X_Y
        _s = Cube->DimZ();
        _n = cz; _u = cx; _v = cy;
    }
    double _l = (p_center-co)*_n/_n.Norm();
    if (abs(_l / _s) < 0.2 &&
        abs(_l) < 2*TDis)// at near
        _o = co;
    else if (abs(_l / _s - 1.0) < 0.2 &&
        abs(_l - _s) < 2*TDis)// at far
        _o = co + _n;
    else
        return false;
    // 2. Check Ratio
    const int nCode = rect->m_index;
    CMeshO::PerVertexAttributeHandle<PtType> type_hi =
        vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
    std::vector<int> idxOnPlane, idxOnCube;
    std::vector<vcg::Point3f> ptsOnPlane, ptsOnCube;
    int idx = 0;
    double _d1 = _u.Norm();
    double _d2 = _v.Norm();
    double _d = 2 * TDis;
    double r1 = 1.0 / _d1;
    double r2 = 1.0 / _d2;   
    for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); vi++, idx++) {
        if (type_hi[vi] == nCode) { // Belong to plane
            vcg::Point3f &p = vi->cP();
            
            double x = ((p-_o)*_u)*r1;
            double y = ((p-_o)*_v)*r2;
            if (x >= -_d && x <= _d1+_d && y >= -_d && y <= _d2+_d) {
                idxOnCube.push_back(idx);
                ptsOnCube.push_back(p);
            }
            else {
                idxOnPlane.push_back(idx);
                ptsOnPlane.push_back(p);
            }               
        }
    }
    double ratio = ptsOnCube.size()*1.0 / (ptsOnCube.size()+ptsOnPlane.size());
    const double T1 = 0.8;
    if (ratio >= T1) {// Have Enogh Cover
        for (int i = ptsOnPlane.size() - 1; i >= 0; --i) {
            int index = idxOnPlane.at(i);
            type_hi[index] = Pt_Undefined;
        }
        return true;
    }
    else {// Need Cut
        while (1) {
            // Pick Max
            std::vector<int> idxList;
            for (int i = 0; i < ptsOnPlane.size(); ++i)
                idxList.push_back(i);
            PicMaxRegion(ptsOnPlane, idxList, TDis);
            if (idxList.size() < 300) {
                for (int i = idxList.size() - 1; i >= 0; --i) {
                    int index = idxList.at(i);
                    type_hi[idxOnPlane.at(index)] = Pt_Undefined;
                }
                break;
            }
            // Fit
            vcg::Plane3f onePlane;
            double err = FinePlane(ptsOnPlane, idxList, onePlane);

            // MBR
            int newCode = _GetPlaneCode();
            ObjRect *oneRect = ExtractMBR(mesh, onePlane, ptsOnPlane, idxOnPlane, idxList);
            if (oneRect != 0) {
                oneRect->m_varN = err;
                // A New Split
                planeSplit.push_back(oneRect);
                for (int i = idxList.size() - 1; i >= 0; --i) {
                    int index = idxList.at(i);
                    type_hi[idxOnPlane.at(index)] = newCode;
                    ptsOnPlane.erase(ptsOnPlane.begin() + index);
                    idxOnPlane.erase(idxOnPlane.begin() + index);
                }
            }
            else
                break;
        }
        return true;
    }

    return false;   
}
int AttachToCube(
    CMeshO &mesh,
    std::vector<ObjRect*> &planes,
    std::vector< std::vector<ObjRect*>> &CubeFaces,
    const std::vector<ObjCube*> &cubes,   
    const double TAng, const double TDis,
    const bool remove)
{
    std::vector<ObjRect*> planesAdded;
    std::vector<ObjRect*> planesSplit;
    for (int i = 0; i < planes.size(); i++) {
        ObjRect *ple = planes.at(i);
        std::vector<ObjRect*> split;
        flog("    >> Attach [ ID.%d ] ... \n", ple->m_index);
        for (int k = 0; k < cubes.size(); ++k) {          
            if (MergeToCube(mesh, ple, cubes.at(k), split, TAng, TDis)) {
                flog("    [ -- Attached to the [ No.%d ] cube with [ %d ] split(s) ... -- ]\n", k+1, split.size());
                CubeFaces.at(k).push_back(ple);
                planesAdded.push_back(ple);
                for (int r = 0; r < split.size(); ++r) {
                    planesSplit.push_back(split.at(r));
                }
                break;
            }
        }
    }

    if (remove)
        RemovePlanes(planes, planesAdded, false);
    for (int i = 0; i < planesSplit.size(); ++i) {
        planes.push_back(planesSplit.at(i));
    }
    return planesAdded.size();
}
