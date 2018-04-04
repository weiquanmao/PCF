#ifndef __MESHMODEL_H
#define __MESHMODEL_H


#include <QFileInfo>
#include <vcg/complex/complex.h>

class CVertexO;
class CFaceO;
class CEdgeO;

class CUsedTypesO : public vcg::UsedTypes <
	vcg::Use<CVertexO>::AsVertexType,
	vcg::Use<CEdgeO  >::AsEdgeType,
	vcg::Use<CFaceO  >::AsFaceType 
>{};
class CVertexO  : public vcg::Vertex< 
	CUsedTypesO,
	vcg::vertex::Coord3f,           /* 12b */
	vcg::vertex::Normal3f,          /* 12b */	
	vcg::vertex::Color4b,           /*  4b */
	vcg::vertex::BitFlags           /*  4b */
>{};
class CFaceO : public vcg::Face<  
	CUsedTypesO,	
	vcg::face::VertexRef,            /*12b */
	vcg::face::BitFlags              /* 4b */
> {};
class CEdgeO : public vcg::Edge<
	CUsedTypesO
> {};
class CMeshO : public vcg::tri::TriMesh< 
	std::vector<CVertexO>,
	std::vector<CFaceO>,
	std::vector<CEdgeO>
> {};


class MeshModel
{
public:
	CMeshO cm;

private:
	int _currentMask;
	QString _fullPathFileName;
	QString _label;

public:
	MeshModel(QString fullFileName, QString labelName);
	~MeshModel(void);

	void Enable(int mask) { _currentMask |= mask; }
	bool hasDataMask(int maskToBeTested) { return ((_currentMask & maskToBeTested) != 0); }
	int mask() {return _currentMask;}

	void Clear();
	QString label() const
	{
		if(_label.isEmpty()) 
			return shortName(); 
		else
			return _label;
	};
	QString shortName() const { return QFileInfo(_fullPathFileName).fileName(); };
	QString fullName() const { return _fullPathFileName; };
};

class MeshDocument : public QObject
{
public:
	MeshModel *mesh;

private:
    bool busy;

public:

	MeshDocument();
	~MeshDocument();

	bool loadMesh(QString Path, QString Label = "");
	bool saveMesh(QString Path, bool bSaveAll = true);
	void delMesh();

	bool isBusy() { return busy;}
	void setBusy(bool _busy) {	busy=_busy;}

	int clearD();
	int inverseD();

	int svn() 
	{
		if (mesh != 0) return mesh->cm.vn;
		else return 0;
	}
	int sfn() 
	{
		if (mesh != 0) return mesh->cm.fn;
		else return 0;
	}
	int sen()
	{
		if (mesh != 0) return mesh->cm.en;
		else return 0;
	}
	vcg::Box3f bbox()
	{
		vcg::Box3f FullBBox;
		if (mesh != 0)
			FullBBox.Add(mesh->cm.bbox);
		return FullBBox;
	}
};

#endif // __MESHMODEL_H_
