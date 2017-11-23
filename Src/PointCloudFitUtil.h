#ifndef _POINT_CLOUD_FIT_UTIL_H_FILE_
#include "PointCloudFit.h"
#include "utility/flog.h"
#ifdef _USE_OPENMP_
#include <omp.h>
#endif
#include <QTime>

#ifndef _DRTrans
#define _DRTrans
#define _D2R 0.017453292
#define _R2D 57.29577951
#define D2R(deg) (deg*_D2R)
#define R2D(rad) (rad*_R2D)
#endif

#define _ReportOut_

template<class T>
static void reportMat(T* mat, const int row, const int col, const char *file)
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
static double EIConfidence(const std::vector<std::pair<double, int>> &list)
{
    const int BinNum = 10;
    const int N = list.size();
    const double min = list.at(0).first;
    const double max = list.at(N - 1).first;
    const double BinStep = (max - min) / (BinNum - 1);
    const double BeginVal = min - BinStep*0.5;
    std::vector<int> NList;
    for (int i = 0; i<BinNum; i++)
        NList.push_back(0);
    for (int i = 0; i<N; ++i) {
        double v = list.at(i).first;
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
    const double _intercept, const double _a, const double _s);
int HoughPlane(
    vcg::Point4f &Plane,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const double _intercept, const double _a, const double _s);

// [2] Surface Points Verification
int AttachToPlane(
    std::vector<int> &PtOnPlane,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const vcg::Point4f &plane,
    const double _TDis, const double _TAng);

// [3] Coplanar Separation
void PicMaxRegion(
    const std::vector<vcg::Point3f> &PointList,
    std::vector<int> &index,
    const double _TDis);

// [4] LS Fit
double FineFit(
    const std::vector<vcg::Point3f> &pointList,
    const std::vector<int> &planeVerList,
    vcg::Point4f &plane);

// [5] Get Minimum-Bounding-Rectangle
double PatchDimensions(
    const std::vector<vcg::Point3f> & pts,
    const vcg::Point3f &nx,
    const vcg::Point3f &ny,
    vcg::Point3f &Pt_O, vcg::Point3f &Pt_Dx, vcg::Point3f &Pt_Dy,
    double &ex, double &ey);
void ExtractMBR(
    CMeshO &mesh,
    ObjPlane &APlane,
    const vcg::Point4f &Plane,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<int> &IndexList,
    const std::vector<int> &PlaneVerList);

std::vector<vcg::Point4f> DetectHTPlanes(
    const std::vector<vcg::Point3f> &pointList,
    const std::vector<vcg::Point3f> &normList,
    const double _intercept, const double _a, const double _s,
    const double TDis, const double TAng,
    const double TNPtsRatio = 0.0, const int TNPtsHard = 0,
    const int ExpPlaneNum = 0,
    std::vector<double> *errors = 0,
    std::vector<ObjPlane*> *pPlanes = 0,
    CMeshO *pMesh = 0, std::vector<int> *pIndexList = 0);

std::vector<ObjPlane*> ExtractPlanes(
    CMeshO &mesh,
    const std::vector<int> &indexList,
    const std::vector<vcg::Point3f> &pointList,    
    const int planeNUm, const int *labels);

// Plane Fitting with GCO
// MPFGCO : Multi-Plane Fitting with GCO
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
MPFGCOCost MPFGCOGeneratCost(
    const std::vector<vcg::Point4f> &planes,
    const std::vector<vcg::Point3f> &points, 
    const std::vector<vcg::Point3f> &norms,
    const double unit_a = 1.0,
    const int cost_noise = 0,
    const int cost_label = 0);
MPFGCONeighbors MPFGCOParseNeighbors(
    CMeshO &mesh,
    const std::vector<int> &ptIndex,
    const int lambda, const double delta,
    const int numNeighbors = 0,
    const double unit_a = 1.0);
std::vector<double> GCOReEstimat(
    std::vector<vcg::Point4f> &planes,
    const std::vector<vcg::Point3f> &pointList,
    const int *labels);

/////////////////////////////////
// ------- For Cube -------
/////////////////////////////////
void RemovePlanes(
    std::vector<ObjPlane*> &planes,
    const std::vector<ObjPlane*> &planesTobeRemoved,
    const bool memRelease = false);
bool IsParallel(const vcg::Point3f &L1, const vcg::Point3f &L2, const double AngTD);
bool IsPerpendicular(const vcg::Point3f &L1, const vcg::Point3f &L2, const double AngTD);

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
double PlaneIoU(const ObjPlane *P1, const ObjPlane *P2);
bool IsOppoFaces(
    const ObjPlane *P1, const ObjPlane *P2,
    const double TAng, const double TIoU);
bool IsAdjacencyFaces(
    const ObjPlane *P1, const ObjPlane *P2,
    const double TRDis, const double TAng,
    PlaneRelation *adjType = 0);
PlaneRelation EstPlaneRelation(
    const ObjPlane *P1, const ObjPlane *P2,
    const double TRDis, const double TAng, const double TIoU);

bool BuildBox(
    ObjPlane* cubePlane[6], ObjPlane* plane,
    const double TRDis, const double TAng, const double TIoU);
std::vector<ObjPlane*> CubeFaceInferringOne(
    const std::vector<ObjPlane*> &planes,
    const double TRDis, const double TAng, const double TIoU);

int CubeFaceInferring(
    std::vector< std::vector<ObjPlane*> > &cubefaces,
    std::vector<ObjPlane*> &planes,
    const double TRDis, const double TAng, const double TIoU,
    const bool remove = true);


// [2] Estimate Cube
void RobustOrientation(
    const std::vector<ObjPlane*> &CubePlanes,
    vcg::Point3f &NX, vcg::Point3f &NY, vcg::Point3f &NZ,
    const double TAng);

void RobustDimention(
    const std::vector<ObjPlane*> &CubePlanes,
    vcg::Point3f &NX, vcg::Point3f &NY, vcg::Point3f &NZ,
    vcg::Point3f &OP, vcg::Point3f &SIZE, vcg::Point3f &WEIGHTS,
    const double TAng);

ObjCube* CubeMeasure(const std::vector<ObjPlane*> &CubePlanes, const double TAng);

// [3] Attach Planes To Cube
bool BelongToCube(
    const ObjPlane *plane,
    const ObjCube *Cube,
    const double TAng);
int AttachToCube(
    std::vector<ObjPlane*> &planes, 
    std::vector< std::vector<ObjPlane*>> &CubeFaces,
    const std::vector<ObjCube*> &cubes,   
    const double TAng,
    const bool remove = true);

/////////////////////////////////
// ------- For Cylinder -------
/////////////////////////////////

// [1] Detect Symmetric Axis
bool DetectSymAxis(
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,
    const double _a,
    vcg::Point3f &PO, vcg::Point3f &NN);

// [2] Attach Points to Cylinder and Estimate Radius
struct CylinderPointInfo {
    int index;
    double r;
    double l;
    CylinderPointInfo(const int _idx = 0, const double _r = 0.0, const double _l = 0.0)
        : index(_idx), r(_r), l(_l) {};
};
double EstRadius(const std::vector<double> &RList, double &weight);
void AttachToCylinder(
    std::vector<int> &PtOnCylinder,
    ObjCylinder &cylinder,
    const std::vector<vcg::Point3f> &PointList,
    const std::vector<vcg::Point3f> &NormList,    
    const double _a, const double TR);

#endif // !_POINT_CLOUD_FIT_UTIL_H_FILE_
