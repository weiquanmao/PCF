#include "PointCloudFit.h"

#include <QSettings> 

#include <wrap/io_trimesh/io_mask.h>

const unsigned int  nColorChannel = 3;
const unsigned char Color_Gray[nColorChannel] = { 100, 100, 100 };
const unsigned char Color_Solid[nColorChannel] = { 255,   0,   0 };
const unsigned char Color_Solid_Cube[nColorChannel] = { 255,   0,   0 };
const unsigned char Color_Solid_Cylinder[nColorChannel] = { 255,   0,   0 };
const unsigned char Color_Noise[nColorChannel] = { 0, 255,   0 };
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


PCFit::PCFit(const int nThread)
	: m_GEOObjSet(0)
	, m_refa(1.0)
{
	printLogo();
	initMParams();
	setThreadNum(nThread);
}
PCFit::PCFit(const char *PlyFilePath, const int nThread)
	: m_GEOObjSet(0)
	, m_refa(1.0)
{
	printLogo();
	initMParams();
	setThreadNum(nThread);
	loadPly(PlyFilePath);
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
		"   /  _  / / ____/  ___   / ____/ / / /__  __/   |  Ver. 0.4.0                \n"
		"  / ____/ / /___   /__/  / ___/  / /    / /      |  November 2017 @ IPC.BUAA \n"
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
		//" | [RefA_Ratio                  ]:       < %0008.3f >   |\n"
		" | [PlaneNum_Expected           ]:       < %000008d >   |\n"
		" | [DeNoise_MaxIteration        ]:       < %000008d >   |\n"
		" | [DeNoise_KNNNeighbors        ]:       < %000008d >   |\n"
        " | [DeNoise_GrowNeighbors       ]:       < %000008d >   |\n"
        " | [DeNoise_DisRatioOfOutlier   ]:       < %0008.3f >   |\n"
		" | [Precision_HT                ]:       < %0008.3f >   |\n"
		" | [Threshold_NPtsPlane         ]:       < %0008.3f >   |\n"
		" | [Threshold_DisToPlane        ]:       < %0008.3f >   |\n"
		" | [Threshold_AngToPlane        ]:       < %0008.3f >   |\n"		
		" | [Threshold_PRAng             ]:       < %0008.3f >   |\n"
		" | [Threshold_PRDis             ]:       < %0008.3f >   |\n"
		" | [Threshold_NPtsCyl           ]:       < %0008.3f >   |\n"
		" +------------------------------------------------------+\n",
		Threshold_NPts, /*RefA_Ratio,*/ PlaneNum_Expected,
		DeNoise_MaxIteration, DeNoise_KNNNeighbors, DeNoise_GrowNeighbors, DeNoise_DisRatioOfOutlier,
		Precision_HT, Threshold_NPtsPlane,
		Threshold_DisToPlane, Threshold_AngToPlane,
		Threshold_PRAng, Threshold_PRDis,
		Threshold_NPtsCyl);
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
	Threshold_NPts = 100;
	//RefA_Ratio = 0.01;
	PlaneNum_Expected = 20;

	DeNoise_MaxIteration = 0;
	DeNoise_KNNNeighbors = 50;
    DeNoise_GrowNeighbors = 20;
    DeNoise_DisRatioOfOutlier = 0.1;

	Precision_HT = 0.01;
	Threshold_NPtsPlane = 0.04;
	Threshold_DisToPlane = 4.0;
	Threshold_AngToPlane = 30.0;
	
	Threshold_PRAng = 15.0;
	Threshold_PRDis = 1.0/6.0;

	Threshold_NPtsCyl = 0.2;

	if (!iniFile && QFileInfo(iniFile).exists()) {
		QSettings conf(iniFile, QSettings::IniFormat);
		QStringList keys = conf.allKeys();
		if (keys.contains("Threshold_NPts"))
			Threshold_NPts = conf.value("Threshold_NPts").toInt();
		//if (keys.contains("RefA_Ratio"))
		//	RefA_Ratio = conf.value("RefA_Ratio").toDouble();
		if (keys.contains("PlaneNum_Expected"))
			PlaneNum_Expected = conf.value("PlaneNum_Expected").toInt();

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
		if (keys.contains("Threshold_NPtsPlane"))
			Threshold_NPtsPlane = conf.value("Threshold_NPtsPlane").toDouble();
		if (keys.contains("Threshold_DisToPlane"))
			Threshold_DisToPlane = conf.value("Threshold_DisToPlane").toDouble();
		if (keys.contains("Threshold_AngToPlane"))
			Threshold_AngToPlane = conf.value("Threshold_AngToPlane").toDouble();

		if (keys.contains("Threshold_PRAng"))
			Threshold_PRAng = conf.value("Threshold_PRAng").toDouble();
		if (keys.contains("Threshold_PRDis"))
			Threshold_PRDis = conf.value("Threshold_PRDis").toDouble();

		if (keys.contains("Threshold_NPtsCyl"))
			Threshold_NPtsCyl = conf.value("Threshold_NPtsCyl").toDouble();
		
	}
	
}


bool PCFit::loadPly(const char * PlyFilePath)
{
	QTime time;
	time.start();
	flog("\n\n[=LoadPly=]: -->> Loading Points from File: %s <<--  \n", PlyFilePath);	
    //----[[
	bool bOpen = false;
	if (PlyFilePath) {		
		bOpen = m_meshDoc.loadMesh(PlyFilePath);
		if (bOpen && m_GEOObjSet != 0) {
			delete m_GEOObjSet;
            m_GEOObjSet = 0;
		}
	}
	//----]]
	if (bOpen)
		flog("[=LoadPly=]: Done, %d point(s) were loaded in %.4f seconds.\n", m_meshDoc.mesh->cm.vn, time.elapsed() / 1000.0);
	else {
		flog("[=LoadPly=]: [Error] Failed to open file %s for reading.\n", PlyFilePath);
		doFailure();
	}
	return bOpen;
}
bool PCFit::savePly(const char *PlyFilePath)
{
	QTime time;
	time.start();
	flog("\n\n[=SavePly=]: -->> Svae Points as Ply File: %s <<--\n", PlyFilePath);	
    //----[[
	bool bSave = m_meshDoc.saveMesh(PlyFilePath, false);
	//----]]
	if (bSave)
		flog("[=SavePly=]: Done, %d point(s) were saved in %.4f seconds.\n", m_meshDoc.mesh->cm.vn, time.elapsed() / 1000.0);
	else {
		flog("[=SavePly=]: [Error] Failed to open file %s for writing.\n", PlyFilePath);
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

int PCFit::deletePts(const int mask)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
	int cnt = 0;
	if (vcg::tri::HasPerVertexAttribute(mesh, _MyPtAttri))
	{
		CMeshO::PerVertexAttributeHandle<PtType> type_hi =
			vcg::tri::Allocator<CMeshO>::GetPerVertexAttribute<PtType>(mesh, _MyPtAttri);
		CMeshO::VertexIterator vi;
		for (vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi)
			if ((type_hi[vi] & mask) != 0) {
				cnt++;
				if (!vi->IsD())
					vcg::tri::Allocator<CMeshO>::DeleteVertex(mesh, *vi);
			}
	}
	return cnt;
}
int PCFit::keepPts(const int mask)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
	int cnt = 0;
	if (vcg::tri::HasPerVertexAttribute(mesh, _MyPtAttri))
	{
		CMeshO::PerVertexAttributeHandle<PtType> type_hi =
			vcg::tri::Allocator<CMeshO>::GetPerVertexAttribute<PtType>(mesh, _MyPtAttri);
		CMeshO::VertexIterator vi;
		for (vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi) {
			if ((type_hi[vi] & mask) != 0) {
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

int PCFit::recolorPts(const int mask, const unsigned char R, const unsigned char G, const unsigned char B, const unsigned char A)
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
	int cnt = 0;
	if (vcg::tri::HasPerVertexAttribute(mesh, _MyPtAttri))
	{
		CMeshO::PerVertexAttributeHandle<PtType> type_hi =
			vcg::tri::Allocator<CMeshO>::GetPerVertexAttribute<PtType>(mesh, _MyPtAttri);
		CMeshO::VertexIterator vi;
		for (vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi)
		{
			if ( (mask == 0 && type_hi[vi] == 0) ||
                (mask & type_hi[vi]) != 0) {
				vi->C().X() = R;
				vi->C().Y() = G;
				vi->C().Z() = B;
				vi->C().W() = A;
				++cnt;
			}
		}
	}
	return cnt;
}

void PCFit::autoColor()
{
	CMeshO &mesh = m_meshDoc.mesh->cm;
	if (vcg::tri::HasPerVertexAttribute(mesh, _MyPtAttri)) {
		CMeshO::PerVertexAttributeHandle<PtType> type_hi = 
			vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, _MyPtAttri);
		unsigned char PlaneColor[3];
		for (CMeshO::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi)
		{
			if (type_hi[vi] == Pt_Undefined || vi->IsD())
				continue;

			const unsigned char* pColor = Color_Gray;
			if (type_hi[vi] == Pt_OnCube)
				pColor = Color_Solid_Cube;
			else if (type_hi[vi] == Pt_OnCylinder)
				pColor = Color_Solid_Cylinder;
			else if (type_hi[vi] == Pt_Noise)
				pColor = Color_Noise;
			else if (type_hi[vi] >= Pt_OnPlane) {
				PlaneColor[0] = Color_Plane[type_hi[vi] % 10][0];
				PlaneColor[1] = Color_Plane[type_hi[vi] % 10][1];
				PlaneColor[2] = Color_Plane[type_hi[vi] % 10][2];
				pColor = PlaneColor;
			}
			for (int i = 0; i<nColorChannel; i++)
				(*vi).C().V(i) = pColor[i];
		}

	}
}

bool PCFit::Fit_Sate(bool keepAttribute)
{
	if (m_meshDoc.mesh == 0 || m_meshDoc.svn() < Threshold_NPts)
		return false;

	QTime time;
	CMeshO &mesh = m_meshDoc.mesh->cm;
    m_GEOObjSet = new ObjSet();

	// [A+] Add Attribute
	bool bAttriAdded = false;
	{
		time.restart();
		flog("\n\n[=InitPtAttri=]: -->> Initialize Satellite Point Attribute <<--\n");		
        //----[[
		CMeshO::PerVertexAttributeHandle<PtType> type_hi;
		if (!vcg::tri::HasPerVertexAttribute(mesh, _MyPtAttri))
		{// Add It
			flog("    >> Add Attri ... \n");
			type_hi = vcg::tri::Allocator<CMeshO>::GetPerVertexAttribute<PtType>(mesh, _MyPtAttri);
			bAttriAdded = true;
		}
		else
		{// Get It
			flog("    >> Get Attri ... \n");
			type_hi = vcg::tri::Allocator<CMeshO>::FindPerVertexAttribute<PtType>(mesh, _MyPtAttri);
			bAttriAdded = false;
		}
		flog("    >> Set Attri ... \n");
		for (int i = 0; i < mesh.VN(); i++)
		{// Set It
			type_hi[i] = Pt_Undefined;
		}
		//-----]]
		flog("[=InitPtAttri=]: Done, %d point(s) were setted in %.4f seconds.\n", mesh.vn, time.elapsed() / 1000.0);
	}



	// [1] A Little Preprocessing
	{
		// [1.1] Remove Outliers
		time.restart();
		flog("\n\n[=DeNoise=]: -->> Remove Outliers <<--  \n");	
        //----[[
#if 0 // DeNoise by KNN
		int nNoise = DeNoiseKNN();
#else // DeNoise by Region Grow
		int nNoise = DeNoiseRegGrw();
#endif
		//----]]
		flog("[=DeNoise=]: Done, %d outlier(s ) were removed in %.4f seconds.\n", nNoise, time.elapsed()/1000.0);


		// [1.2] Get Dimension Reference Unit
		time.restart();
		flog("\n\n[=RefSize=]: -->> Dimension Reference Unit Ana. <<--  \n");	
        //----[[
		vcg::Point3f PSize;
		std::vector<vcg::Point3f> PDirection;

		// PCADimensionAna(PSize, PDirection, true);
		// m_refa = PSize.V(2)*RefA_Ratio;
		m_refa = RoughnessAna(true);

		//----]]
		flog("[=RefSize=]: Done in %.4f seconds.\n", time.elapsed() / 1000.0);
	}



	// -- 2. Find All Planes
	std::vector<ObjPlane*> planes;
	{
#if 1
		flog("\n\n[=PlaneFit_HT=]: -->> Try to Fit %d Planes by Hough Translation <<--  \n", PlaneNum_Expected);
		time.restart();
		//-------------------------------
		planes = DetectPlanesHT(PlaneNum_Expected);
		//-------------------------------
		flog("[=PlaneFit_HT=]: Done, %d plane(s) were detected in %.4f seconds.\n", planes.size(), time.elapsed() / 1000.0);
#else
        flog("\n\n[=PlaneFit_MPFGCO=]: -->> Try to Fit by Multi-Planes Fitting in GCO with %d Expected Initial Planes <<--  \n", PlaneNum_Expected);
        time.restart();
        //-------------------------------
        planes = DetectPlanesGCO(PlaneNum_Expected, 10);
        //-------------------------------
        flog("[=PlaneFit_MPFGCO=]: Done, %d plane(s) were detected in %.4f seconds.\n", planes.size(), time.elapsed() / 1000.0);      
#endif
    }

	// -- 3. Judge Cube Planes
	std::vector<ObjCube*> cubes;
	if (planes.size() >= 2) {
		flog("\n\n[=DetectCube=]: -->> Try to Detect Cube from %d Planes <<--  \n", planes.size());
		time.restart();
		//-------------------------------
        cubes = DetectCubeFromPlanes(planes);
		//-------------------------------
		flog("[=DetectCube=]: Done, %d cube(s) is detected in %.4f seconds, the other %d plane(s) are left.\n", cubes.size(), time.elapsed() / 1000.0, planes.size());
	}
    for (int i = 0; i < cubes.size(); ++i)
        m_GEOObjSet->m_SolidList.push_back(cubes.at(i));



	// -- 4. Find Cyclinde
    ObjCylinder *objCylinder = 0;
	if ( m_meshDoc.mesh->hasDataMask(vcg::tri::io::Mask::IOM_VERTNORMAL) ) {
#if 1
		flog("\n\n[=DetectCylinder_SA=]: -->> Try to Detect Cylinder by Symmetric Axis <<--  \n");
		time.restart();
		//-------------------------------
        objCylinder = DetectCylinderSymAxis();
		//-------------------------------
		if (objCylinder != 0)
			flog("[=DetectCylinder_SA=]: Done, cylinder is detected in %.4f seconds.\n", time.elapsed() / 1000.0);
		else
			flog("[=DetectCylinder_SA=]: No cylinder is detected. Elapsed %.4f seconds.\n", time.elapsed() / 1000.0);
#else
#endif
    }
	else {
		flog("\n\n[=DetectCylinder=]: [ ): ] Normals are needed to detected cylinder. \n");
	}
    if (objCylinder != 0)
        m_GEOObjSet->m_SolidList.push_back(objCylinder);


	// -- 5. Set SailbordList
	{
		flog("\n\n[=PlanesCheck=]: %d plane(s) are left. \n", planes.size());
		m_GEOObjSet->m_PlaneList.swap(planes);
	}


	// -- *Remove Added Attribute
	if (bAttriAdded && !keepAttribute) {
		flog("\n\n[=CleanPtAttri=]: -->> Delete Point Sate Attribute <<--\n");
		vcg::tri::Allocator<CMeshO>::DeletePerVertexAttribute(mesh, _MyPtAttri);
	}
	
	return true;
}