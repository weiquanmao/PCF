#include "MeshDoc.h"
#include <wrap/io_trimesh/export.h>
#include <wrap/io_trimesh/import.h>


MeshModel::MeshModel(QString fullFileName, QString labelName)
	: cm()
{
  Clear();
  if( !fullFileName.isEmpty() )   
	  _fullPathFileName=fullFileName;
  if( !labelName.isEmpty() )     
	  _label=labelName;
}
MeshModel::~MeshModel(void)
{
	Clear();
}
void MeshModel::Clear()
{
	_currentMask = vcg::tri::io::Mask::IOM_NONE;
	cm.Clear();
	CMeshO::VertContainer tempVert;
	CMeshO::FaceContainer tempFace;
	CMeshO::EdgeContainer tempEdge;
	cm.vert.swap(tempVert);
	cm.face.swap(tempFace);
	cm.edge.swap(tempEdge);
}

MeshDocument::MeshDocument()
{
    mesh=0;
    busy=false;
}
MeshDocument::~MeshDocument()
{
	delMesh();
}

bool MeshDocument::loadMesh(QString Path, QString label)
{
	if(Path.isEmpty())
		return false;

	QFileInfo fi(Path);
	QString extension = fi.suffix();
	if (! (extension.toUpper() == tr("PLY") ||
		   extension.toUpper() == tr("STL") ||
		   extension.toUpper() == tr("OBJ")))
		return false;
	
	if(!fi.exists() || !fi.isReadable())
        return false;


	std::string filename = Path.toLocal8Bit().data();
	int mask = 0;
	if (label.isEmpty())
		label = fi.baseName();
	MeshModel *newMesh = new MeshModel(Path, label);
	int ret = vcg::tri::io::Importer<CMeshO>::Open(newMesh->cm, filename.c_str(), mask);
	if ( ret == 0 || (ret<6 && extension.toUpper() == tr("OBJ")) )
	{
		newMesh->Enable(mask);
		vcg::tri::UpdateBounding<CMeshO>::Box(newMesh->cm);
		delMesh();
		mesh = newMesh;
		return true;
	}
	else {
        printf(vcg::tri::io::Exporter<CMeshO>::ErrorMsg(ret));
		delete newMesh;
		return false;	
	}
}
bool MeshDocument::saveMesh(QString Path, bool bSaveAll)
{
	if (bSaveAll)
		clearD();

	std::string filename = Path.toLocal8Bit().data();
	int mask = mesh->mask();
	int result = vcg::tri::io::Exporter<CMeshO>::Save(mesh->cm,filename.c_str(),mask);
    if (result != 0) {
        printf(vcg::tri::io::Exporter<CMeshO>::ErrorMsg(result));
        return false;
    }
    else
		return true;
}
void MeshDocument::delMesh()
{
	if (mesh != 0) {
		delete mesh;
		mesh = 0;
	}
}

int MeshDocument::clearD()
{
	int cnt = 0;
	for (CMeshO::VertexIterator vi = mesh->cm.vert.begin(); vi != mesh->cm.vert.end(); ++vi) {
		vi->ClearD();
		++cnt;
	}
	mesh->cm.vn = cnt;
	return cnt;
}
int MeshDocument::inverseD()
{
	int cnt = 0;
	for (CMeshO::VertexIterator vi = mesh->cm.vert.begin(); vi != mesh->cm.vert.end(); ++vi) {
		if (vi->IsD()) {
			vi->ClearD();
			++cnt;
		}
		else
			vi->SetD();
	}
	mesh->cm.vn = cnt;
	return cnt;
}