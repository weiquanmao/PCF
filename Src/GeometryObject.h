#ifndef _GEOMETRY_OBJECT_H_FILE_
#define _GEOMETRY_OBJECT_H_FILE_

#include <vector>
#include <vcg/space/deprecated_point3.h>

enum GeoObjType {
    GeoObj_Patch = 0x100,
    Patch_Rectangle = 0x111,
    Patch_Circle = 0x121,
    Patch_Arbitary = 0x1FF,

    GeoObj_Solid = 0x200,
    Solid_Cube = 0x211,
    Solid_Cylinder = 0x212,
    Solid_Cone = 0x221
};

class GeoObj
{
public:
    virtual GeoObjType type() const = 0;
    vcg::Point3f m_O;
    const int m_index;
    GeoObj(int index) : m_index(index) { m_O.SetZero(); }
};

class ObjPatch : public GeoObj
{
public:
    virtual GeoObjType type() const { return GeoObj_Patch; }

    vcg::Point3f m_N;
    double m_varN;
    ObjPatch(int index)
        : GeoObj(index) {
        m_N.SetZero();
        m_varN = 0.0;
    }
};
class ObjRect : public ObjPatch
{
public:
    GeoObjType type() const { return Patch_Rectangle; }

    float height() const { return m_AY.Norm(); }
    float width() const { return m_AX.Norm(); }

    vcg::Point3f m_AX;
    vcg::Point3f m_AY;
    double  m_EIConfX;
    double  m_EIConfY;

    ObjRect(int index)
        : ObjPatch(index) {
        m_AX.SetZero();
        m_AY.SetZero();
        m_EIConfX = 0.0;
        m_EIConfY = 0.0;
    }
};
class ObjCircle : public ObjPatch
{
public:
    GeoObjType type() const { return Patch_Circle; }
    double  m_radius;
    ObjCircle(int index) : ObjPatch(index) { m_radius = 0.0; }
};

class ObjSolid : public GeoObj
{
public:
    virtual GeoObjType type() const { return GeoObj_Solid; }
    ObjSolid(int index) : GeoObj(index) { }
};
class ObjCube : public ObjSolid
{
public:
    GeoObjType type() const { return Solid_Cube; }

    double DimX() const { return m_AX.Norm(); }
    double DimY() const { return m_AY.Norm(); }
    double DimZ() const { return m_AZ.Norm(); }

    vcg::Point3f m_AX;
    vcg::Point3f m_AY;
    vcg::Point3f m_AZ;
    double  m_EIConfX;
    double  m_EIConfY;
    double  m_EIConfZ;
    ObjCube(int index)
        : ObjSolid(index) {
        m_EIConfX = 0.0;
        m_EIConfY = 0.0;
        m_EIConfZ = 0.0;
        m_AX.SetZero();
        m_AY.SetZero();
        m_AZ.SetZero();
    }
};
class ObjCylinder : public ObjSolid
{
public:
    GeoObjType type() const { return Solid_Cylinder; }

    double m_radius;
    double m_length;
    vcg::Point3f m_N;
    ObjCylinder(int index)
        : ObjSolid(index) {
        m_N.SetZero();
        m_radius = 0.0;
        m_length = 0.0;
    }

};
class ObjSet {
public:
    std::vector<ObjPatch*> m_PlaneList;
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
    bool containPlane(ObjPatch *plane) {
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
    void delPlane(ObjPatch *plane) {
        std::vector<ObjPatch*>::iterator it;
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