#ifndef _GEOMETRY_OBJECT_H_FILE_
#define _GEOMETRY_OBJECT_H_FILE_

#include <vector>
#include <vcg/space/deprecated_point3.h>

class ObjPlane
{
public:
	float height() const { return m_dY.Norm(); }
	float width() const { return m_dX.Norm(); }
	float thickness() const { return m_dZ.Norm(); }

	vcg::Point3f m_pO;
	vcg::Point3f m_dX;
	vcg::Point3f m_dY;
	vcg::Point3f m_dZ;
	vcg::Point3f m_sizeConfidence;
	vcg::Point3f m_N;
	double m_varN;
	const int m_PlaneIndex;
    ObjPlane(int index)
		: m_PlaneIndex(index) {
		m_N.SetZero();
		m_pO.SetZero();
		m_dX.SetZero();
		m_dY.SetZero();
		m_dZ.SetZero();
		m_sizeConfidence.SetZero();
		m_varN = 0.0;
	}
};
class ObjSolid
{
public:
	enum SolidType {
		Solid_Cube = 1,
		Solid_Cylinder = 2,
        Solid_Cone = 3
	};
	virtual SolidType type() const = 0;
};
class ObjCube : public ObjSolid
{
public:
	vcg::Point3f m_pO;
	vcg::Point3f m_dX;
	vcg::Point3f m_dY;
	vcg::Point3f m_dZ;
	vcg::Point3f m_sizeConfidence;
	double DimX() const { return m_dX.Norm(); }
	double DimY() const { return m_dY.Norm(); }
	double DimZ() const { return m_dZ.Norm(); }
    ObjCube()
		: ObjSolid() {
		m_sizeConfidence.SetZero();
		m_pO.SetZero();
		m_dX.SetZero();
		m_dY.SetZero();
		m_dZ.SetZero();
	}
    SolidType type() const { return Solid_Cube; }
};
class ObjCylinder : public ObjSolid
{
public:
	double m_radius;
	double m_length;
	double m_radius_weight;
	double m_length_weight;
	vcg::Point3f m_pO;
	vcg::Point3f m_N;
    ObjCylinder()
		: ObjSolid() {
		m_radius = 0.0;
		m_length = 0.0;
		m_N.SetZero();
		m_pO.SetZero();
	}
    SolidType type() const { return Solid_Cylinder; }
};
class ObjSet {
public:
	std::vector<ObjPlane*> m_PlaneList;
    std::vector<ObjSolid*> m_SolidList;
    ObjSet() {}
	~ObjSet() {
		clean();
	}

    bool containSolid(ObjSolid *solid) {
        if (solid == 0)
            return false;
        for (int i = 0; i<m_SolidList.size(); i++) {
            if (m_SolidList.at(i) == solid)
                return true;
        }
        return false;
    }
	bool containBord(ObjPlane *plane) {
		if (plane == 0)
			return false;
		for (int i = 0; i<m_PlaneList.size(); i++) {
			if (m_PlaneList.at(i) == plane)
				return true;
		}
		return false;
	}
    void delSolid(ObjSolid *solid) {
        std::vector<ObjSolid*>::iterator it;
        for (it = m_SolidList.begin(); it != m_SolidList.end(); it++) {
            if (*it == solid) {
                m_SolidList.erase(it);
                break;
            }
        }

    }
	void delPlane(ObjPlane *plane) {
		std::vector<ObjPlane*>::iterator it;
		for (it = m_PlaneList.begin(); it != m_PlaneList.end(); it++) {
			if (*it == plane) {
                m_PlaneList.erase(it);
				break;
			}
		}

	}
	void cleanSolid() {
        for (int i = 0; i<m_SolidList.size(); i++)
            delete m_SolidList.at(i);
        m_SolidList.clear();
	}
	void cleanPlanes() {
		for (int i = 0; i<m_PlaneList.size(); i++)
			delete m_PlaneList.at(i);
        m_PlaneList.clear();
	}
	void clean() {
        cleanPlanes();
        cleanSolid();
	}
};

void SaveObjSet(ObjSet *objSet, const char *file);
ObjSet* LoadObjSet(const char *file);

#endif // !_GEOMETRY_OBJECT_H_FILE_
