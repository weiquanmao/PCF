#ifndef _POINT_CLOUD_FIT_UTIL_H_FILE_
#include "PointCloudFit.h"
#include "utility/flog.h"
#ifdef _USE_OPENMP_
#include <omp.h>
#endif
#include <QTime>

#ifndef _PCF_ANG_
#define _PCF_ANG_
#define _D2R 0.017453292
#define _R2D 57.29577951
#define D2R(deg) (deg*_D2R)
#define R2D(rad) (rad*_R2D)

#define VCGAngle(N1, N2) R2D(abs(vcg::Angle(N1, N2)))
#define CheckAng00(ang) (90.0 - abs(90.0 - ang))
#define CheckAng90(ang) abs(90.0 - ang)
#endif


#define _ReportOut_

template<class T> static void reportMat(
    T* mat, const int row, const int col, 
    const char *file)
{
    std::ofstream fout(file);
    for (int i = 0; i <row; ++i) {
        for (int j = 0; j < col; ++j) {
            fout << mat[i*col + j] << " ";
        }
        fout << std::endl;
    }
    fout.close();
}

static double EIConfidence(const std::vector<int> &list)
{
    const int Num = list.size();
    const double MaxEI = -log(1.0 / Num);
    double Sumnlogn = 0.0;
    double Sumn = 0.0;
    for (int i = 0; i<Num; ++i) {
        int n = list.at(i);
        if (n>0) {
            Sumnlogn += n*log(n);
            Sumn += n;
        }
    }
    double EI = (Sumn*log(Sumn) - Sumnlogn) / (Sumn*1.0);
    return EI / MaxEI;
}
static double EIConfidence(const std::vector<double> &list)
{// list MUST be in order!
    const int BinNum = 10;
    const int N = list.size();
    const double min = list.at(0);
    const double max = list.at(N - 1);
    const double BinStep = abs(max - min) / (BinNum - 1);
    const double BeginVal = min < max ? min - BinStep*0.5 : max - BinStep*0.5;
    std::vector<int> NList;
    for (int i = 0; i<BinNum; i++)
        NList.push_back(0);
    for (int i = 0; i<N; ++i) {
        double v = list.at(i);
        int n = (v - BeginVal) / BinStep;
        NList[n]++;
    }
    double conf = EIConfidence(NList);
    NList.clear();
    return conf;
}


/////////////////////////////////
// ------- For Plane -------
/////////////////////////////////

const vcg::Plane3f _NON_PLANE = vcg::Plane3f(0.0, vcg::Point3f(0.0, 0.0, 0.0));

// [1] Detect Planes
enum FixedAxis {
    FixedAxis_X = 0,
    FixedAxis_Y = 1,
    FixedAxis_Z = 2
};
int HoughPlaneOne(
    vcg::Plane3f &Plane,
    const FixedAxis fix,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const double _intercept, const double _a, const double _s);
int HoughPlane(
    vcg::Plane3f &Plane,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const double _intercept, const double _a, const double _s);

// [2] Surface Points Verification
int AttachToPlane(
    std::vector<int> &PtIdx,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const vcg::Plane3f &plane,
    const double _TDis, const double _TAng);

// [3] Coplanar Separation
void PicMaxRegion(
    const std::vector<vcg::Point3f> &PointList,
    std::vector<int> &Index,
    const double _TDis);

// [4] LS Fit
double FinePlane(
    const std::vector<vcg::Point3f> &pointList,
    const std::vector<int> &planeVerList,
    vcg::Plane3f &plane);

// [5] Get Minimum-Bounding-Rectangle
bool PatchDimensionOne(
    std::vector<std::pair<double, int>> &proList,
    int &idxBegin, int &idxEnd, double r = 1.414);
double PatchDimensions(
    const std::vector<vcg::Point3f> &pts,
    const vcg::Point3f &nx,
    const vcg::Point3f &ny,
    vcg::Point3f &Pt_O, vcg::Point3f &Pt_Dx, vcg::Point3f &Pt_Dy,
    double *ex = 0, double *ey = 0);
ObjRect* ExtractMBR(
    CMeshO &mesh,
    const vcg::Plane3f &Plane,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<int> &IndexList,
    const std::vector<int> &PlaneVerList);
ObjCircle* CircleCheck(
    const ObjRect *rect,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<int> &PlaneVerList);

int DetectHTPlanes(
    std::vector<vcg::Plane3f> &Planes,
    const std::vector<vcg::Point3f> &pointList,
    const std::vector<vcg::Point3f> &normList,
    const double _intercept, const double _a, const double _s,
    const double TDis, const double TAng,
    const double TNPtsRatio = 0.0, const int TNPtsHard = 0,
    const int ExpPlaneNum = 0,
    std::vector<double> *errors = 0,
    std::vector<ObjPatch*> *pPlanes = 0,
    CMeshO *pMesh = 0, std::vector<int> *pIndexList = 0);

int ExtractPatches(
    CMeshO &mesh,
    std::vector<ObjPatch*> &patches,
    const std::vector<int> &indexList,
    const std::vector<vcg::Point3f> &pointList,    
    const int planeNum, const int *labels);


/////////////////////////////////
// ------- For Cube -------
/////////////////////////////////

inline bool IsParallel(const vcg::Point3f &L1, const vcg::Point3f &L2, const double AngTD);
inline bool IsPerpendicular(const vcg::Point3f &L1, const vcg::Point3f &L2, const double AngTD);
double PlaneIoU(const ObjRect *P1, const ObjRect *P2);
void RemovePlanes(
    std::vector<ObjRect*> &planes,
    const std::vector<ObjRect*> &planesTobeRemoved,
    const bool memRelease = false);

// [1] Infer Cube Faces
enum PlaneRelation {
    PlaneRelation_NoRetation  = 0x00,
    /*****************************************
     *           |               | 
     *           |      [U]p     |
     * ---------[O] -----[X]-----+---------
     *           | / / / / / / / |
     *  [L]eft  [Y] / /Rect / / /| [R]ight
     *           | / / / / / / / |
     * ----------+---------------+---------
     *           |    [B]ottom   |
     *           |               |
     *****************************************/
    PlaneRelation_Adjacency_U = 0x01,
    PlaneRelation_Adjacency_B = 0x03,
    PlaneRelation_Adjacency_L = 0x05,    
    PlaneRelation_Adjacency_R = 0x07,
    PlaneRelation_AtOppo      = 0x08,
};
bool IsOppoFaces(
    const ObjRect *P1, const ObjRect *P2,
    const double TAng, const double TIoU);
bool IsAdjacencyFaces(
    const ObjRect *P1, const ObjRect *P2,
    const double TRDis, const double TAng,
    PlaneRelation *adjType = 0);
PlaneRelation EstPlaneRelation(
    const ObjRect *P1, const ObjRect *P2,
    const double TRDis, const double TAng, const double TIoU);


bool BuildBox(
    ObjRect* cubePlane[6], ObjRect* Rect,
    const double TRDis, const double TAng, const double TIoU);
std::vector<ObjRect*> CubeFaceInferringOne(
    const std::vector<ObjRect*> &Rects,
    const double TRDis, const double TAng, const double TIoU);
int CubeFaceInferring(
    std::vector< std::vector<ObjRect*> > &CubeFaces,
    std::vector<ObjRect*> &Rects,
    const double TRDis, const double TAng, const double TIoU,
    const bool remove = true);


// [2] Estimate Cube
void RobustOrientation(
    const std::vector<ObjRect*> &CubePlanes,
    vcg::Point3f &NX, vcg::Point3f &NY, vcg::Point3f &NZ,
    const double TAng);

void RobustDimention(
    const std::vector<ObjRect*> &CubePlanes,
    const vcg::Point3f &NX, const vcg::Point3f &NY, const vcg::Point3f &NZ,
    vcg::Point3f &OP, vcg::Point3f &SIZE, vcg::Point3f &WEIGHTS,
    const double TAng);

void CubeMeasure(
    const std::vector<ObjRect*> &CubePlanes,
    ObjCube *cube, const double TAng);

// [3] Attach Planes To Cube
bool MergeToCube(
    CMeshO &mesh,
    const ObjRect *rect, const ObjCube *Cube,
    std::vector<ObjRect*> &planeSplit,   
    const double TAng, const double TDis);
int AttachToCube(
    CMeshO &mesh,
    std::vector<ObjRect*> &rects, 
    std::vector< std::vector<ObjRect*>> &CubeFaces,
    const std::vector<ObjCube*> &cubes,   
    const double TAng, const double TDis);

/////////////////////////////////
// ------- For Cylinder -------
/////////////////////////////////

// [1] Detect Symmetric Axis
bool DetectSymAxis(
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,    
    vcg::Point3f &PO, vcg::Point3f &NN,
    const double _a);

// [2] Attach Points to Cylinder and Estimate Radius
struct CylinderPointInfo {
    int index;
    double r;
    double l;
    CylinderPointInfo(const int _idx = 0, const double _r = 0.0, const double _l = 0.0)
        : index(_idx), r(_r), l(_l) {};
};
double EstRadius(const std::vector<double> &RList);
void AttachToCylinder(
    std::vector<int> &PtOnCylinder,
    ObjCylinder &cylinder,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,    
    const double _a, const double TR);

// [3] Others (Mainly Used in MCF-GCO)
double SignedDistanceCylinderPoint(
    const ObjCylinder &cyl, const vcg::Point3f &p);
double AngCylinderPoint(
    const ObjCylinder &cyl, const vcg::Point3f &p, const vcg::Point3f &n);
bool CylinderInlier(
	const ObjCylinder &cyl,
	const vcg::Point3f &p, const double TDis,
	const bool containInner = false);
bool CylinderInlier(
	const ObjCylinder &cyl,
	const vcg::Point3f &p, const double TDis,
	const vcg::Point3f &n, const double TAng, const bool containInner = false);
int CylinderInliers(
    const ObjCylinder &cyl,
    const std::vector<vcg::Point3f> &pointList,
    const std::vector<vcg::Point3f> &normList,
    const double TDis, const double TAng,
	const bool containInner = false,
    std::vector<int> *inlierIdx = 0);
bool CloseCylinders(
    const ObjCylinder &cyl1,
    const ObjCylinder &cyl2,
    const double TRDis = 0.1, const double TAng = 2);
ObjCylinder *EstCylinderTwoPoint(
    const vcg::Point3f p1, const vcg::Point3f n1,
    const vcg::Point3f p2, const vcg::Point3f n2,
    const double TDisDeviation = (1.0-0.6)/(1.0+0.6),
    const int TAngRequired = 15);
ObjCylinder *FineCylinder(
    const std::vector<vcg::Point3f> &pointList,
    const std::vector<int> &cylVerList,
    double &err);

double DetectCylinderRansac(
    const std::vector<vcg::Point3f> &pointList, 
    const std::vector<vcg::Point3f> &normList, 
    std::vector<ObjCylinder*> &cylCandidates,
    const double TDis, const double TAng,
    const int maxN = 1, const double inlierRatio = 0.1,
    const vcg::Box3f *constriantBox = 0);
double FineCylinderLength(
    CMeshO &mesh,
    ObjCylinder &cyl,
    std::vector<int> &inlierIdx);
int AttachToCylinder(
    CMeshO &mesh,
    std::vector<ObjCylinder*> &cylinders,
	const double TDis, const double TAng,
	const double TInlierRatio = 0.0);

////////////////////////////////////////////////
// ------- Multi-Model Fitting with GCO -------
////////////////////////////////////////////////
// MPFGCO : Multi-Plane Fitting with GCO
// MCFGCO : Multi-Cylinder Fitting with GCO
// MMFGCO : Multi-Model Fitting with GCO
// GCO    : Graph Cut Optimization
struct MPFGCOCost {
    int numLabel;
    int numSite;

    int *dataCost;
    int *smoothCost;
    int labelCost;

    MPFGCOCost() : numLabel(0), numSite(0), dataCost(0), smoothCost(0), labelCost(0) {};
    void memRelease() {
        numLabel = 0;
        numSite = 0;

        if (!dataCost)
            delete[] dataCost;
        dataCost = 0;
        if (!smoothCost)
            delete[] smoothCost;
        smoothCost = 0;
        labelCost = 0;
    }
};
struct MPFGCONeighbors {
    int numNeighbor;
    int numSite;

    int *neighborsCounts;
    int **neighborsIndexes;
    int **neighborsWeights;

    MPFGCONeighbors() : numNeighbor(0), numSite(0), neighborsCounts(0), neighborsIndexes(0), neighborsWeights(0) {};
    void memRelease() {
        numNeighbor = 0;
        numSite = 0;
        if (!neighborsCounts)
            delete[] neighborsCounts;
        neighborsCounts = 0;

        if (!neighborsIndexes) {
            if (!neighborsIndexes[0])
                delete[] neighborsIndexes[0];
            delete[] neighborsIndexes;
        }
        neighborsIndexes = 0;

        if (!neighborsWeights) {
            if (!neighborsWeights[0])
                delete[] neighborsWeights[0];
            delete[] neighborsWeights;
        }
        neighborsWeights = 0;;
    }
};

MPFGCONeighbors MPFGCOParseNeighbors(
    CMeshO &mesh,
    const std::vector<int> &ptIndex,
    const int lambda, const double delta,
    const int numNeighbors = 0,
    const double unit_a = 1.0);

// [ -- !!! IMPLEMENTS ARE TOO SIMILAR !!! -- ]
MPFGCOCost MPFGCOGeneratCost(
    const std::vector<ObjCylinder*> &cylinders,
    const std::vector<vcg::Point3f> &points,
    const std::vector<vcg::Point3f> &norms,
    const double unit_a = 1.0,
    const int cost_noise = 0,
    const int cost_label = 0);
MPFGCOCost MPFGCOGeneratCost(
    const std::vector<vcg::Plane3f> &planes,
    const std::vector<vcg::Point3f> &points,
    const std::vector<vcg::Point3f> &norms,
    const double unit_a = 1.0,
    const int cost_noise = 0,
    const int cost_label = 0);

// [ -- !!! IMPLEMENTS ARE TOO SIMILAR !!! -- ]
std::vector<double> GCOReEstimat(
    std::vector<vcg::Plane3f> &planes,
    const std::vector<vcg::Point3f> &pointList,
    const int *labels, const unsigned int TInlier = 0);
std::vector<double> GCOReEstimat(
    std::vector<ObjCylinder*> &planes,
    const std::vector<vcg::Point3f> &pointList,
    const int *labels, const unsigned int TInlier = 0);
#endif // !_POINT_CLOUD_FIT_UTIL_H_FILE_
