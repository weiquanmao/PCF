#include "PointCloudFit.h"

#include <QSettings> 

#include <wrap/io_trimesh/io_mask.h>

const unsigned int  nColorChannel = 3;
const unsigned char Color_Gray[nColorChannel]           = { 100, 100, 100 };
const unsigned char Color_Solid[nColorChannel]          = { 255,   0,   0 };
const unsigned char Color_Solid_Cube[nColorChannel]     = { 255,   0,   0 };
const unsigned char Color_Solid_Cylinder[nColorChannel] = {   0,   0, 255 };
const unsigned char Color_Noise[nColorChannel]          = {   0, 255,   0 };
const unsigned char Color_Plane[10][nColorChannel] =
{
	{ 215,  55, 101 },
	{ 255, 245,  55 },
	{  97,  54, 130 },
	{ 194,  85,  38 },
	{   0, 192, 210 },
	{   0, 178, 111 },
	{ 250, 202,  87 },
	{ 255,  95,  61 },
	{ 128,  88, 189 },
	{ 149,  85,  66 }
};

void doFailure() {
	system("pause");
	exit(1);
}
int _code_plane_ = 0;
int _code_cube_ = 0;
int _code_cylinder_ = 0;
int _code_cone_ = 0;
int _ResetObjCode(_PtType type)
{
    int _old_ = -1;
    switch (type)
    {
    case Pt_OnPlane:
        _old_ = _code_plane_;
        _code_plane_ = 0;
        break;
    case Pt_OnCube:
        _old_ = _code_cube_;
        _code_cube_ = 0;
        break;
    case Pt_OnCylinder:
        _old_ = _code_cylinder_;
        _code_cylinder_ = 0;
        break;
    case Pt_OnCone:
        _old_ = _code_cone_;
        _code_cone_ = 0;
        break;
    default:
        break;
    }
    return _old_;
}
int _GetObjCode(_PtType type)
{
    int _code_ = -1;
    switch (type)
    {
    case Pt_OnPlane:
        _code_ =  _code_plane_++;
        break;
    case Pt_OnCube:
        _code_ =  _code_cube_++;
        break;
    case Pt_OnCylinder:
        _code_ =  _code_cylinder_++;
        break;
    case Pt_OnCone:
        _code_ = _code_cone_++;
        break;
    default:
        break;
    }
	if (_code_ > MaxModelNum)
		flog("[=_GetObjCode=]: [Warning] Model Number Exceed the Maximum %d.\n", MaxModelNum);

    return type+_code_;
}

PCFit::PCFit(const int nThread)
	: m_GEOObjSet(0)
	, m_refa(1.0)
{
	printLogo();
	initMParams();
	setThreadNum(nThread);
}
PCFit::PCFit(const char *plyFilePath, const int nThread)
	: m_GEOObjSet(0)
	, m_refa(1.0)
{
	printLogo();
	initMParams();
	setThreadNum(nThread);
	loadPly(plyFilePath);
}
PCFit::~PCFit()
{
	if (!m_GEOObjSet) {
		delete m_GEOObjSet;
        m_GEOObjSet = 0;
	}
}


void PCFit::printLogo()
{
	flog(
		"\n\n"
		"    ______  ______         ______  __  _______   |  PCFit - Point Cloud Fit   \n"
		"   /  _  / / ____/  ___   / ____/ / / /__  __/   |  Version 1.2.0             \n"
		"  / ____/ / /___   /__/  / ___/  / /    / /      |  November 2017 @ IPC.BUAA  \n"
		" /_/     /_____/        /_/     /_/    /_/       |  By WeiQM                  \n"
		"                                                 |  Email: weiqm@buaa.edu.cn  \n"
		"\n"
		"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
		"\n"
	);
}
void PCFit::printParams()
{
	flog(
		" +------------------------------------------------------+\n"
		" | PARAMETERS IN USE:                                   |\n"
		" |------------------------------------------------------|\n"
		" | [Threshold_NPts              ]:       < %000008d >   |\n"
        // " | [RefA_Ratio                  ]:       < %0008.3f >   |\n"
		" | [DeNoise_MaxIteration        ]:       < %000008d >   |\n"
		" | [DeNoise_KNNNeighbors        ]:       < %000008d >   |\n"
        " | [DeNoise_GrowNeighbors       ]:       < %000008d >   |\n"
        " | [DeNoise_DisRatioOfOutlier   ]:       < %0008.3f >   |\n"
		" | [Precision_HT                ]:       < %0008.3f >   |\n"
        " | [Threshold_MaxModelNumPre    ]:       < %000008d >   |\n"
        " | [Threshold_MaxModelNum       ]:       < %000008d >   |\n"
		" | [Threshold_NPtsPlane         ]:       < %0008.3f >   |\n"
        " | [Threshold_NPtsCylinder      ]:       < %0008.3f >   |\n"
		" | [Threshold_DisToSurface      ]:       < %0008.3f >   |\n"
		" | [Threshold_AngToSurface      ]:       < %0008.3f >   |\n"
		" | [Threshold_PRAng             ]:       < %0008.3f >   |\n"
		" | [Threshold_PRDis             ]:       < %0008.3f >   |\n"
        " | [Threshold_PRIoU             ]:       < %0008.3f >   |\n"		
		" +------------------------------------------------------+\n",
		Threshold_NPts, /*RefA_Ratio,*/
		DeNoise_MaxIteration, DeNoise_KNNNeighbors, DeNoise_GrowNeighbors, DeNoise_DisRatioOfOutlier,
		Precision_HT, Threshold_MaxModelNumPre, Threshold_MaxModelNum,
        Threshold_NPtsPlane, Threshold_NPtsCylinder,
		Threshold_DisToSurface, Threshold_AngToSurface, 
		Threshold_PRAng, Threshold_PRDis, Threshold_PRIoU);
}


void PCFit::setThreadNum(const int nThread)
{
#ifdef _USE_OPENMP_
	int nThreadMax = omp_get_max_threads();
	int nThreadUsed = (nThread <= 0 || nThread > nThreadMax) ? nThreadMax : nThread;
	omp_set_num_threads(nThreadUsed);
	flog(
		" +------------------------------------------------------+\n"
		" |          USE OPENMP: %02d THREAD(s) ARE USED.          |\n"
		" +------------------------------------------------------+\n",
		nThreadUsed);
#endif // _USE_OPENMP_
}
void PCFit::initMParams(const char *iniFile)
{
	Threshold_NPts    = 100;
    // RefA_Ratio        = 0.01;

	DeNoise_MaxIteration      = 0;
	DeNoise_KNNNeighbors      = 50;
    DeNoise_GrowNeighbors     = 20;
    DeNoise_DisRatioOfOutlier = 0.1;

    Precision_HT             = 0.01;
    Threshold_MaxModelNumPre = 10;
    Threshold_MaxModelNum    = 30;
    Threshold_NPtsPlane      = 0.05;
    Threshold_NPtsCylinder   = 0.1;
    Threshold_DisToSurface   = 2.0;
    Threshold_AngToSurface   = 20.0;

#if  _RECON_DATA_
    Threshold_AngToSurface = 30.0;
#elif _SYN_DATA_
    Threshold_AngToSurface = 15.0;
#endif // _SYN_DATA_
    
	Threshold_PRAng = 15.0;
	Threshold_PRDis = 0.50;
    Threshold_PRIoU = 0.60;

	if (!iniFile && QFileInfo(iniFile).exists()) {
		QSettings conf(iniFile, QSettings::IniFormat);
		QStringList keys = conf.allKeys();
		if (keys.contains("Threshold_NPts"))
			Threshold_NPts = conf.value("Threshold_NPts").toInt();
        //if (keys.contains("RefA_Ratio"))
        //    RefA_Ratio = conf.value("RefA_Ratio").toInt();

		if (keys.contains("DeNoise_MaxIteration"))
			DeNoise_MaxIteration = conf.value("DeNoise_MaxIteration").toInt();
		if (keys.contains("DeNoise_KNNNeighbors"))
			DeNoise_KNNNeighbors = conf.value("DeNoise_KNNNeighbors").toInt();
        if (keys.contains("DeNoise_GrowNeighbors"))
            DeNoise_GrowNeighbors = conf.value("DeNoise_GrowNeighbors").toInt();
		if (keys.contains("DeNoise_DisRatioOfOutlier"))
			DeNoise_DisRatioOfOutlier = conf.value("DeNoise_DisRatioOfOutlier").toDouble();
	
		if (keys.contains("Precision_HT"))
			Precision_HT = conf.value("Precision_HT").toDouble();
        if (keys.contains("Threshold_MaxModelNumPre"))
            Threshold_MaxModelNumPre = conf.value("Threshold_MaxModelNumPre").toInt();
        if (keys.contains("Threshold_MaxModelNum"))
            Threshold_MaxModelNum = conf.value("Threshold_MaxModelNum").toInt();

		if (keys.contains("Threshold_NPtsPlane"))
			Threshold_NPtsPlane = conf.value("Threshold_NPtsPlane").toDouble();
        if (keys.contains("Threshold_NPtsCylinder"))
            Threshold_NPtsCylinder = conf.value("Threshold_NPtsCylinder").toDouble();

		if (keys.contains("Threshold_DisToSurface"))
			Threshold_DisToSurface = conf.value("Threshold_DisToSurface").toDouble();
		if (keys.contains("Threshold_AngToSurface"))
			Threshold_AngToSurface = conf.value("Threshold_AngToSurface").toDouble();
        
		if (keys.contains("Threshold_PRAng"))
			Threshold_PRAng = conf.value("Threshold_PRAng").toDouble();
		if (keys.contains("Threshold_PRDis"))
			Threshold_PRDis = conf.value("Threshold_PRDis").toDouble();
        if (keys.contains("Threshold_PRIoU"))
            Threshold_PRIoU = conf.value("Threshold_PRIoU").toDouble();		
	}
	
}


bool PCFit::loadPly(const char * plyFilePath)
{
	QTime time;
	time.start();
	flog("\n\n[=LoadPly=]: -->> Loading Points from File: %s <<--  \n", plyFilePath);	
    //----[[
	bool bOpen = false;
	if (plyFilePath) {		
		bOpen = m_meshDoc.loadMesh(plyFilePath);
		if (bOpen && m_GEOObjSet != 0) {
			delete m_GEOObjSet;
            m_GEOObjSet = 0;
		}
	}
	//----]]
	if (bOpen)
		flog("[=LoadPly=]: Done, %d point(s) were loaded in %.4f seconds.\n", m_meshDoc.mesh->cm.vn, time.elapsed() / 1000.0);
	else {
		flog("[=LoadPly=]: [Error] Failed to open file %s for reading.\n", plyFilePath);
		doFailure();
	}
	return bOpen;
}
bool PCFit::savePly(const char *plyFilePath)
{
	QTime time;
	time.start();
	flog("\n\n[=SavePly=]: -->> Save Points as Ply File: %s <<--\n", plyFilePath);	
    //----[[
	bool bSave = m_meshDoc.saveMesh(plyFilePath, false);
	//----]]
	if (bSave)
		flog("[=SavePly=]: Done, %d point(s) were saved in %.4f seconds.\n", m_meshDoc.mesh->cm.vn, time.elapsed() / 1000.0);
	else {
		flog("[=SavePly=]: [Error] Failed to open file %s for writing.\n", plyFilePath);
		doFailure();
	}
	return bSave;
}

void PCFit::clearD()
{
	m_meshDoc.clearD();
}
void PCFit::inverseD()
{
	m_meshDoc.inverseD();
}

typedef bool(*_PtrFunc_)(const int a, const int b); 
typedef _PtrFunc_ _OP_;
inline bool _Comp_Equ_(const int a, const int b) { return a == b; }
inline bool _Comp_Con_(const int a, const int b) { return (a & b) == b; }
int _deletePts(CMeshO &mesh, _OP_ op, const int mask)
{
	int cnt = 0;
	if (vcg::tri::HasPerVertexAttribute(mesh, PtAttri_GeoType))
	{
		CMeshO::PerVertexAttributeHandle<PtType> type_hi =
			vcg::tri::Allocator<CMeshO>::GetPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
		CMeshO::VertexIterator vi;
		for (vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi)
			if (op(type_hi[vi],mask)) {
				cnt++;
				if (!vi->IsD())
					vcg::tri::Allocator<CMeshO>::DeleteVertex(mesh, *vi);
			}
	}
	return cnt;
}
int _keepPts(CMeshO &mesh, _OP_ op, const int mask)
{
	int cnt = 0;
	if (vcg::tri::HasPerVertexAttribute(mesh, PtAttri_GeoType))
	{
		CMeshO::PerVertexAttributeHandle<PtType> type_hi =
			vcg::tri::Allocator<CMeshO>::GetPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
		CMeshO::VertexIterator vi;
		for (vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi) {
			if (op(type_hi[vi],mask)) {
				cnt++;
				vi->ClearD();
			}
			else {
				vi->SetD();
			}
		}
	}
	mesh.vn = cnt;
	return cnt;
}
int _recolorPts(CMeshO &mesh, _OP_ op, const int mask, const unsigned char r, const unsigned char g, const unsigned char b, const unsigned char a)
{
	int cnt = 0;
	if (vcg::tri::HasPerVertexAttribute(mesh, PtAttri_GeoType))
	{
		CMeshO::PerVertexAttributeHandle<PtType> type_hi =
			vcg::tri::Allocator<CMeshO>::GetPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
		CMeshO::VertexIterator vi;
		for (vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi)
		{
			if (op(mask,type_hi[vi])) {
				vi->C().X() = r;
				vi->C().Y() = g;
				vi->C().Z() = b;
				vi->C().W() = a;
				++cnt;
			}
		}
	}
	return cnt;
}

int PCFit::deletePts(const int maskID)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
	_OP_ op = _Comp_Equ_;
	return _deletePts(mesh, op, maskID);
}
int PCFit::keepPts(const int maskID)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
	_OP_ op = _Comp_Equ_;
	return _keepPts(mesh, op, maskID);
}
int PCFit::recolorPts(const int maskID, const unsigned char r, const unsigned char g, const unsigned char b, const unsigned char a)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
	_OP_ op = _Comp_Equ_;
	return _recolorPts(mesh, op, maskID, r, g, b, a);
}
int PCFit::deletePtType(const _PtType type)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
	_OP_ op = _Comp_Con_;
	return _deletePts(mesh, op, type);
}
int PCFit::keepPtType(const _PtType type)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
	_OP_ op = _Comp_Con_;
	return _keepPts(mesh, op, type);
}
int PCFit::recolorPtType(const _PtType type, const unsigned char r, const unsigned char g, const unsigned char b, const unsigned char a)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
	_OP_ op = _Comp_Con_;
	return _recolorPts(mesh, op, type, r, g, b, a);
}
int PCFit::resetType(const _PtType type)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
	int cnt = 0;
	if (vcg::tri::HasPerVertexAttribute(mesh, PtAttri_GeoType))
	{
		CMeshO::PerVertexAttributeHandle<PtType> type_hi =
			vcg::tri::Allocator<CMeshO>::GetPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
		CMeshO::VertexIterator vi;
		for (vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi)
			if (!vi->IsD() && _Comp_Con_(type_hi[vi], type)) {
				type_hi[vi] = Pt_Undefined;
				cnt++;
			}
	}
	return cnt;
}

void PCFit::autoColor()
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
	if (vcg::tri::HasPerVertexAttribute(mesh, PtAttri_GeoType)) {
		CMeshO::PerVertexAttributeHandle<PtType> type_hi = 
			vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
		unsigned char PlaneColor[3];
		for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi)
		{
            int code = type_hi[vi];
			if (code == Pt_Undefined)
				continue;

			const unsigned char* pColor = Color_Gray;
			if (code == Pt_Noise)
				pColor = Color_Noise;
			else if (_Comp_Con_(code,Pt_OnCube))
				pColor = Color_Solid_Cube;
			else if (_Comp_Con_(code,Pt_OnCylinder))
				pColor = Color_Solid_Cylinder;
			else if (_Comp_Con_(code,Pt_OnPlane)) {
				PlaneColor[0] = Color_Plane[(code - Pt_OnPlane) % 10][0];
				PlaneColor[1] = Color_Plane[(code - Pt_OnPlane) % 10][1];
				PlaneColor[2] = Color_Plane[(code - Pt_OnPlane) % 10][2];
				pColor = PlaneColor;
			}
			for (int i = 0; i<nColorChannel; i++)
				(*vi).C().V(i) = pColor[i];
		}

	}
}

bool PCFit::GEOFit(ProType proType, bool keepAttribute)
{
    if (proType == DoNothing)
        return true;

	if (m_meshDoc.mesh == 0 || m_meshDoc.svn() < Threshold_NPts)
		return false;

	QTime time;
	CMeshO &mesh = m_meshDoc.mesh->cm;
    m_GEOObjSet = new ObjSet();

	// [+] Add Attribute
	bool bAttriAdded = false;
	{
		time.restart();
		flog("\n\n[=InitPtAttri=]: -->> Initialize Satellite Point Attribute <<--\n");		
        //----[[
		CMeshO::PerVertexAttributeHandle<PtType> type_hi;
		if (!vcg::tri::HasPerVertexAttribute(mesh, PtAttri_GeoType))
		{// Add Attribute
			flog("    >> Add Attri ... \n");
			type_hi = vcg::tri::Allocator<CMeshO>::GetPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
			bAttriAdded = true;
		}
		else
		{// Get Attribute
			flog("    >> Get Attri ... \n");
			type_hi = vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, PtAttri_GeoType);
			bAttriAdded = false;
		}
		flog("    >> Set Attri ... \n");
		for (int i = 0; i < mesh.VN(); i++)
		{// Set Attribute
			type_hi[i] = Pt_Undefined;
		}
		//-----]]
		flog("[=InitPtAttri=]: Done, %d point(s) were set in %.4f seconds.\n", mesh.vn, time.elapsed() / 1000.0);
	}


	// [1] Remove Outliers
    if ((proType & OneStep_RemoveOutlier) != 0)
    {
        time.restart();
        flog("\n\n[=DeNoise=]: -->> Remove Outliers <<--  \n");
        //----[[
#if 0 // DeNoise by KNN
        int nNoise = DeNoiseKNN();
#else // DeNoise by Region Grow
        int nNoise = DeNoiseRegGrw();
#endif
        //----]]
        flog("[=DeNoise=]: Done, %d outliers were removed in %.4f seconds.\n", nNoise, time.elapsed() / 1000.0);
    }
   
    // [2] Get Dimension Reference Unit
    {
		// [1.2] 
		time.restart();
		flog("\n\n[=RefSize=]: -->> Dimension Reference Unit Ana. <<--  \n");	
        //----[[
		vcg::Point3f PSize;
		std::vector<vcg::Point3f> PDirection;
        PCADimensions(PDirection, PSize);
        double PCASize_x1 = PSize.V(2)*0.01;
        double PCASize_x3 = PSize.V(2)*0.03;
        double RoughSize = Roughness();
        m_refa = RoughSize < PCASize_x1 ? PCASize_x1 : (RoughSize > PCASize_x3 ? PCASize_x3 : RoughSize);

		//----]]
		flog("[=RefSize=]: Done in %.4f seconds.\n", time.elapsed() / 1000.0);
	}
    
    // [3] Detect Cylinders
    std::vector<ObjCylinder*> objCylinder;
    std::vector<ObjPatch*> prePlanes;
    if ((proType & OneStep_DetectCylinder) != 0) 
    {
        // -- 3.1 Pre Planes Detect       
        {
            flog("\n\n[=PrePlaneFit_HT=]: -->> Try to Remove %d Planes by Hough Translation <<--  \n", Threshold_MaxModelNumPre);
            time.restart();
            //-------------------------------
            prePlanes = DetectPlanesHT(Threshold_MaxModelNumPre);
            //-------------------------------
            flog("[=PlaneFit_HT=]: Done, %d plane(s) were removed in %.4f seconds.\n", prePlanes.size(), time.elapsed() / 1000.0);
        }

        // -- 3.2 Detect Cylinder
        
        if (m_meshDoc.mesh->hasDataMask(vcg::tri::io::Mask::IOM_VERTNORMAL)) {
#if 0 // Detect Planes by Symmetric Axis Detection with Hough Transform
            flog("\n\n[=DetectCylinder_SA=]: -->> Try to Detect Cylinder by Symmetric Axis Detection <<--  \n");
            time.restart();
            //-------------------------------
            ObjCylinder *oneCly = DetectCylinderSymAxis();
            //-------------------------------
            if (oneCly != 0)
                objCylinder.push_back(oneCly);
            flog("[=DetectCylinder_SA=]: Done, %d cylinder is detected in %.4f seconds.\n", objCylinder.size(), time.elapsed() / 1000.0);
#else // Detect Planes by Multi-Model Fitting with GCO
            flog("\n\n[=DetectMCF_GCO=]: -->> Try to Detect Cylinder by Energy-based Multi-Model Fitting with GCO ... <<--  \n");
            time.restart();
            //-------------------------------
            objCylinder = DetectCylinderGCO(Threshold_MaxModelNum);
            //-------------------------------
            flog("[=DetectMCF_GCO=]: Done, %d cylinder is detected in %.4f seconds.\n", objCylinder.size(), time.elapsed() / 1000.0);

#endif
        }
        else {
            flog("\n\n[=DetectCylinder=]: [ ): ] Normals are needed to detected cylinder. \n");
        }
        // -- 3.3 Remove Pre-Plane
        resetType(Pt_OnPlane);
        
    }
    for (int i = 0; i < objCylinder.size(); ++i)
        m_GEOObjSet->m_SolidList.push_back(objCylinder.at(i));
    m_GEOObjSet->m_PlaneList.swap(prePlanes); // It May Be Useful

	// [4] Detect All Planes
    std::vector<ObjPatch*> planes;
    if ((proType & OneStep_DetectPlane) != 0)
	{
#if 1 // Detect Planes by Hough Transform
		flog("\n\n[=PlaneFit_HT=]: -->> Try to Detect %d Planes by Hough Translation <<--  \n", Threshold_MaxModelNum);
		time.restart();
		//-------------------------------
		planes = DetectPlanesHT(Threshold_MaxModelNum);
		//-------------------------------
		flog("[=PlaneFit_HT=]: Done, %d plane(s) were detected in %.4f seconds.\n", planes.size(), time.elapsed() / 1000.0);
#else // Detect Planes by Multi-Model Fitting with GCO
        flog("\n\n[=PlaneFit_MPFGCO=]: -->> Try to Fit by Multi-Planes Fitting in GCO with %d Expected Initial Planes <<--  \n", Threshold_MaxModelNum);
        time.restart();
        //-------------------------------
        planes = DetectPlanesGCO(Threshold_MaxModelNum, 10);
        //-------------------------------
        flog("[=PlaneFit_MPFGCO=]: Done, %d plane(s) were detected in %.4f seconds.\n", planes.size(), time.elapsed() / 1000.0);      
#endif
    }
   
	// [5] Detect Cube from Planes
	std::vector<ObjCube*> cubes;
    if ((proType & OneStep_DetectCube) != 0)
    {
        if (planes.size() >= 2) {
            flog("\n\n[=DetectCube=]: -->> Try to Detect Cube from %d Planes <<--  \n", planes.size());
            time.restart();
            //-------------------------------
            cubes = DetectCubeFromPlanes(planes);
            //-------------------------------
            flog("[=DetectCube=]: Done, %d cube(s) is detected in %.4f seconds, the other %d plane(s) are left.\n", cubes.size(), time.elapsed() / 1000.0, planes.size());
        }
    }
    for (int i = 0; i < cubes.size(); ++i)
        m_GEOObjSet->m_SolidList.push_back(cubes.at(i));
    
	// [6] Set Planes
    if ((proType & OneStep_DetectPlane) != 0)
	{
		flog("\n\n[=PlanesCheck=]: -->> %d plane(s) are left. << -- \n", planes.size());
        for (int i = 0; i < m_GEOObjSet->m_PlaneList.size(); ++i)
            delete m_GEOObjSet->m_PlaneList.at(i);
		m_GEOObjSet->m_PlaneList.swap(planes);
	}

	
	// [-] Remove Added Attribute
	if (bAttriAdded && !keepAttribute) {
		flog("\n\n[=CleanPtAttri=]: -->> Delete Point Sate Attribute <<--\n");
		vcg::tri::Allocator<CMeshO>::DeletePerVertexAttribute(mesh, PtAttri_GeoType);
	}
	
	return true;
}
