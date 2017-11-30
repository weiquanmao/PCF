#ifndef _POINT_CLOUD_FIT_H_FILE_
#define _POINT_CLOUD_FIT_H_FILE_

#if 1
#define _RECON_DATA_ 1
#else
#define _SYN_DATA_ 1
#endif

#include "MeshDoc.h"
#include "GeometryObject.h"
#include "utility/flog.h"
#ifdef _USE_OPENMP_
#include <omp.h>
#endif

#include <QTime>

void doFailure();

#define PtAttri_GeoType "PerVerAttri_PtType"
typedef int PtType;
#define MaxModelNum 256
enum _PtType {
	Pt_Undefined  = 0x0000, // 00000
	Pt_Noise      = 0x0100, // 00001
	Pt_OnPlane    = 0x0200, // 00010
    Pt_OnCube     = 0x0600, // 00110
    Pt_OnCylinder = 0x0800, // 01000
    Pt_OnCone     = 0x1800  // 11000
};
int _ResetObjCode(_PtType type);
int _GetObjCode(_PtType type);

class PCFit
{
public:
	PCFit(const int nThread = 0);
	PCFit(const char *plyFilePath, const int nThread = 0);
	~PCFit();

	void printLogo();
	void printParams();

	void setThreadNum(const int nThread);
	void initMParams(const char *iniFile = 0);

	bool loadPly(const char *plyFilePath);
	bool savePly(const char *plyFilePath);
	void clearD();
	void inverseD();
	int deletePts(const int maskID);
	int keepPts(const int maskID);
	int recolorPts(const int maskID, const unsigned char r, const unsigned char g, const unsigned char b, const unsigned char a = 255);
	int deletePtType(const _PtType type);
	int keepPtType(const _PtType type);
	int recolorPtType(const _PtType type, const unsigned char r, const unsigned char g, const unsigned char b, const unsigned char a = 255);
	int resetType(const _PtType type);

	void autoColor();
    enum ProType {
        DoNothing                = 0x00,
        //---------------------
        OneStep_RemoveOutlier    = 0x01,
        OneStep_DetectCylinder   = 0x02,
        OneStep_DetectPlane      = 0x04,
        OneStep_DetectCube       = 0x08,
        //---------------------
        SpSteps_OnlyPlane        = 0x04,  // Without Remove Outlier Step
        SpSteps_ToPlane          = 0x06,  // Without Remove Outlier Step
        SpSteps_ToCube           = 0x0E,  // Without Remove Outlier Step
        Steps_OnlyPlane          = 0x05,
        Steps_ToPlane            = 0x07,
        Steps_ToCube             = 0x0F,
        Steps_All                = 0x0F,
        //---------------------
        
    };
	bool GEOFit(ProType proType = Steps_All, bool keepAttribute = true);
	ObjSet *getGEOObjSet() { return m_GEOObjSet; }
public:
	// All Thresholds
    // 基础参数-[基本不用改]
	int Threshold_NPts;                      // 整体点云数小于[Threshold_NPts]不进行处理
	// double RefA_Ratio;                       // 单位长度比例

    // 去噪参数-[基本不用改]
	int DeNoise_MaxIteration;                // 去噪迭代阈值, 默认100(<=0)
	int DeNoise_KNNNeighbors;                // 去噪时KNN最近邻个数(KNN)
    int DeNoise_GrowNeighbors;               // 去噪时KNN最近邻个数(区域生长)
	double DeNoise_DisRatioOfOutlier;        // 去噪距离比例阈值

    // 平面+圆柱检测参数
	double Precision_HT;                     // HT系数
    int Threshold_MaxModelNumPre;            // 圆检测前预提取最大平面数
    int Threshold_MaxModelNum;               // 最大模型个数
	double Threshold_NPtsPlane;              // 平面点个数比例阈值
    double Threshold_NPtsCylinder;           // 圆柱点个数阈值系数
	double Threshold_DisToSurface;           // 平面+圆柱检测距离阈值系数
	double Threshold_AngToSurface;           // 平面+圆柱检测角度阈值
    
    // 立方体推断参数
	double Threshold_PRAng;                  // 平面关系角度阈值
	double Threshold_PRDis;                  // 平面关系距离阈值系数
    double Threshold_PRIoU;                  // 平面关系交并比阈值
private:
    MeshDocument m_meshDoc;
    ObjSet *m_GEOObjSet;
    double m_refa;                           // 长度参考单位
private:

    vcg::Point3f GetPointList(
        std::vector<int> &indexList,
        std::vector<vcg::Point3f> &pointList,
        std::vector<vcg::Point3f> &normList,
        const bool moveToCenter = false);

	// Preproc
	int DeNoiseKNN();
	int DeNoiseRegGrw();
	bool PCADimensions(
        std::vector<vcg::Point3f> &PDirections, 
        vcg::Point3f &PSize);
	double Roughness();

    
	// Detect Plane
	std::vector<ObjPatch*> DetectPlanesHT(const int expPlaneNum);
	
	// Detect Cude
    std::vector<ObjCube*> DetectCubeFromPlanes(std::vector<ObjPatch*> &patches);
	
	// Detect Cylinder
	ObjCylinder* DetectCylinderSymAxis();

    // Detect By MMF-GCO
    // [ -- !!! IMPLEMENTS ARE TOO SIMILAR !!! -- ]
    std::vector<ObjPatch*> DetectPlanesGCO(const int expPlaneNum, const int iteration = -1);
    std::vector<ObjCylinder*> DetectCylinderGCO(const int expCylinderNum, const int iteration = -1);
};

#endif // !_POINT_CLOUD_FIT_H_FILE_
