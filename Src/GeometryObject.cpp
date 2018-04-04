#include "GeometryObject.h"
#include <fstream>
#include <iomanip>
#include <string>

#define GeometryObjectHead \
"# --------  COMMENTS --------\n" \
"# GeometryObject 2.1         \n" \
"# by PCF                     \n" \
"# Comment by start with '#'. \n" \
"#                            \n" \
"# Type code for Solid:       \n" \
"#   +------------------------\n" \
"#   | [0] : Have no solid    \n" \
"#   | [1] : Cube             \n" \
"#   | [2] : Cylinder         \n" \
"#   | [3] : Cone             \n" \
"#   +------------------------\n" \
"# Type code for Plane:       \n" \
"#   +------------------------\n" \
"#   | [0] : Have no plane    \n" \
"#   | [1] : Rectangle        \n" \
"#   | [2] : Circle           \n" \
"#   +------------------------\n" \
"#                            \n" \
"# Record in following order: \n" \
"#   +------------------------\n" \
"#   | 1 : Solids             \n" \
"#   | 2 : Planes             \n" \
"#   +------------------------\n" \
"#                            \n" \
"# Record Format              \n" \
"#   +------------------------\n" \
"#   | <Number of Solid(s)>   \n" \
"#   | --[For Cuboid]--       \n" \
"#   | <Solid Type := 1> <ID> \n" \
"#   | <O.X> <O.Y> <O.Z>      \n" \
"#   | <X.X> <X.Y> <X.Z> <EX> \n" \
"#   | <Y.X> <Y.Y> <Y.Z> <EY> \n" \
"#   | <Z.X> <Z.Y> <Z.Z> <EZ> \n" \
"#   | --[For Cylinder]--     \n" \
"#   | <Solid Type := 2> <ID> \n" \
"#   | <O.X> <O.Y> <O.Z>      \n" \
"#   | <N.X> <N.Y> <N.Z>      \n" \
"#   | <Radius> <Length>      \n" \
"#   | ...                    \n" \
"#   |                        \n" \
"#   | <Number of Plane(s)>   \n" \
"#   | --[For Rectangle]--    \n" \
"#   | <Plane Type := 1> <ID> \n" \
"#   | <O.X> <O.Y> <O.Z>      \n" \
"#   | <N.X> <N.Y> <N.Z> <Var>\n" \
"#   | <X.X> <X.Y> <X.Z> <EX> \n" \
"#   | <Y.X> <Y.Y> <Y.Z> <EY> \n" \
"#   | --[For Circle]--       \n" \
"#   | <Plane Type := 2> <ID> \n" \
"#   | <O.X> <O.Y> <O.Z> <R>  \n" \
"#   | <N.X> <N.Y> <N.Z> <Var>\n" \
"#   | ...                    \n" \
"#   +------------------------\n" \
"# ----- END OF COMMENTS -----\n" \


void SaveObjSet(ObjSet *objSet, const char *file)
{
    if (objSet != 0 && file != 0) {
        std::ofstream out(file, std::ios::out);
        out << std::fixed << std::setprecision(10);
        out << GeometryObjectHead << '\n';

        out << objSet->m_SolidList.size() << '\n';
        for (int i = 0; i<objSet->m_SolidList.size(); i++) {
            ObjSolid *solid = objSet->m_SolidList.at(i);
            if (solid->type() == Solid_Cube) { // Cube
                ObjCube *cube = (ObjCube*)solid;
                out << "1" << '\t' << cube->m_index << "\n";
                out << cube->m_O.X() << '\t' << cube->m_O.Y() << '\t' << cube->m_O.Z() << '\n'
                    << cube->m_AX.X() << '\t' << cube->m_AX.Y() << '\t' << cube->m_AX.Z() << '\t' << cube->m_EIConfX << '\n'
                    << cube->m_AY.X() << '\t' << cube->m_AY.Y() << '\t' << cube->m_AY.Z() << '\t' << cube->m_EIConfY << '\n'
                    << cube->m_AZ.X() << '\t' << cube->m_AZ.Y() << '\t' << cube->m_AZ.Z() << '\t' << cube->m_EIConfZ << '\n';
            }
            if (solid->type() == Solid_Cylinder) { // Cylinder
                ObjCylinder *cyl = (ObjCylinder*)solid;
                out << "2" << '\t' << cyl->m_index << "\n";
                out << cyl->m_O.X() << '\t' << cyl->m_O.Y() << '\t' << cyl->m_O.Z() << '\n'
                    << cyl->m_N.X() << '\t' << cyl->m_N.Y() << '\t' << cyl->m_N.Z() << '\n'
                    << cyl->m_radius << '\t' << cyl->m_length << '\n';
            }
            if (solid->type() == Solid_Cone) { // Cone
                out << "3" << "\n";
            }
        }


        out << objSet->m_PlaneList.size() << '\n';
        for (int i = 0; i<objSet->m_PlaneList.size(); i++) {

            ObjPatch *patch = objSet->m_PlaneList.at(i);
            if (patch->type() == Patch_Rectangle) { // Rectangle
                ObjRect *rect = (ObjRect*)patch;
                out << "1" << '\t' << rect->m_index << "\n";
                out << rect->m_O.X() << '\t' << rect->m_O.Y() << '\t' << rect->m_O.Z() << '\n'
                    << rect->m_N.X() << '\t' << rect->m_N.Y() << '\t' << rect->m_N.Z() << '\t' << rect->m_varN << '\n'
                    << rect->m_AX.X() << '\t' << rect->m_AX.Y() << '\t' << rect->m_AX.Z() << '\t' << rect->m_EIConfX << '\n'
                    << rect->m_AY.X() << '\t' << rect->m_AY.Y() << '\t' << rect->m_AY.Z() << '\t' << rect->m_EIConfY << '\n';
            }
            if (patch->type() == Patch_Circle) { // Circle
                ObjCircle *cir = (ObjCircle*)patch;
                out << "2" << '\t' << cir->m_index << "\n";
                out << cir->m_O.X() << '\t' << cir->m_O.Y() << '\t' << cir->m_O.Z() << '\t' << cir->m_radius << '\n'
                    << cir->m_N.X() << '\t' << cir->m_N.Y() << '\t' << cir->m_N.Z() << '\t' << cir->m_varN << '\n';
            }
            if (patch->type() == Patch_Arbitary) { // Arbitary
                out << "3" << "\n";
            }

        }
        out.close();
    }
}
ObjSet* LoadObjSet(const char *file)
{
    if (file != 0) {
        std::ifstream in(file, std::ios::in);
        if (!in.is_open())
            return 0;
        ObjSet *objSet = new ObjSet();

        std::string taken;
        while (in.peek() == '#') {
            std::getline(in, taken);
            if (in.peek() == '\n')
                in.get();
        }
        int nSolid;
        in >> nSolid;
        for (int i = 0; i<nSolid; i++) {
            int solidType, nCode;
            in >> solidType >> nCode;
            if (solidType == 1) { // Cube
                ObjCube *cube = new ObjCube(nCode);
                in >> cube->m_O.X() >> cube->m_O.Y() >> cube->m_O.Z()
                    >> cube->m_AX.X() >> cube->m_AX.Y() >> cube->m_AX.Z() >> cube->m_EIConfX
                    >> cube->m_AY.X() >> cube->m_AY.Y() >> cube->m_AY.Z() >> cube->m_EIConfY
                    >> cube->m_AZ.X() >> cube->m_AZ.Y() >> cube->m_AZ.Z() >> cube->m_EIConfZ;
                objSet->m_SolidList.push_back(cube);
            }
            if (solidType == 2) { // Sylinder
                ObjCylinder *cyl = new ObjCylinder(nCode);
                in >> cyl->m_O.X() >> cyl->m_O.Y() >> cyl->m_O.Z()
                    >> cyl->m_N.X() >> cyl->m_N.Y() >> cyl->m_N.Z()
                    >> cyl->m_radius >> cyl->m_length;
                objSet->m_SolidList.push_back(cyl);
            }
            if (solidType == 3) { // Cone
                                  //
            }
        }


        int nPatch;
        in >> nPatch;
        for (int i = 0; i<nPatch; i++) { // Patch
            int patchType, nCode;
            in >> patchType >> nCode;
            if (patchType == 1) { // Rectangle
                ObjRect *rect = new ObjRect(nCode);
                in >> rect->m_O.X() >> rect->m_O.Y() >> rect->m_O.Z()
                    >> rect->m_N.X() >> rect->m_N.Y() >> rect->m_N.Z() >> rect->m_varN
                    >> rect->m_AX.X() >> rect->m_AX.Y() >> rect->m_AX.Z() >> rect->m_EIConfX
                    >> rect->m_AY.X() >> rect->m_AY.Y() >> rect->m_AY.Z() >> rect->m_EIConfY;

                objSet->m_PlaneList.push_back(rect);
            }
            if (patchType == 2) { // Circle
                ObjCircle *cir = new ObjCircle(nCode);
                in >> cir->m_O.X() >> cir->m_O.Y() >> cir->m_O.Z() >> cir->m_radius
                    >> cir->m_N.X() >> cir->m_N.Y() >> cir->m_N.Z() >> cir->m_varN;

                objSet->m_PlaneList.push_back(cir);
            }
            if (patchType == 3) { // Arbitary
                                  //
            }
        }
        in.close();
        return objSet;
    }
    return 0;
}
