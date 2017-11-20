#include "GeometryObject.h"
#include <fstream>
#include <iomanip>
#include <string>

void SaveObjSet(ObjSet *objSet, const char *file)
{
	if ( objSet != 0 && 
        (!objSet->m_SolidList.empty() || !objSet->m_PlaneList.empty()) &&
		file != 0) {
		std::ofstream out(file, std::ios::out);
		out << std::fixed << std::setprecision(10);
        out << "#GeometryObject 1.0" << '\n'
            << "#by PCF" << '\n'
            << "#Sold" << '\n'
            << "# [0] : Have no solid" << '\n'
            << "# [1] : Cube" << '\n' 
            << "# [2] : Cylinder" << '\n'
            << "# [3] : Cone" << '\n'
			<< "#SBs" << '\n'
			<< '\n';

        out << objSet->m_SolidList.size() << '\n';
        for (int i = 0; i<objSet->m_SolidList.size(); i++) {
            ObjSolid *solid = objSet->m_SolidList.at(i);
            if (solid->type() == ObjSolid::Solid_Cube) {
                ObjCube *cube = (ObjCube*)solid;
                out << "1" << "\n";
                out << cube->m_pO.X() << '\t' << cube->m_pO.Y() << '\t' << cube->m_pO.Z() << '\n'
                    << cube->m_dX.X() << '\t' << cube->m_dX.Y() << '\t' << cube->m_dX.Z() << cube->m_sizeConfidence.X() << '\n'
                    << cube->m_dY.X() << '\t' << cube->m_dY.Y() << '\t' << cube->m_dY.Z() << cube->m_sizeConfidence.Y() << '\n'
                    << cube->m_dZ.X() << '\t' << cube->m_dZ.Y() << '\t' << cube->m_dZ.Z() << cube->m_sizeConfidence.Z() << '\n';
            }
            if (solid->type() == ObjSolid::Solid_Cylinder) {
                ObjCylinder *cyl = (ObjCylinder*)solid;
                out << "2" << "\n";
                out << cyl->m_pO.X() << '\t' << cyl->m_pO.Y() << '\t' << cyl->m_pO.Z() << '\n'
                    << cyl->m_N.X() << '\t' << cyl->m_N.Y() << '\t' << cyl->m_N.Z() << '\n'
                    << cyl->m_length << '\t' << cyl->m_length_weight << '\n'
                    << cyl->m_radius << '\t' << cyl->m_radius_weight << '\n';
            }
            if (solid->type() == ObjSolid::Solid_Cone) {
                //ObjCone *cone = (ObjCone*)solid;
                out << "3" << "\n";
                //out << cyl->m_pO.X() << '\t' << cyl->m_pO.Y() << '\t' << cyl->m_pO.Z() << '\n'
                //    << cyl->m_N.X() << '\t' << cyl->m_N.Y() << '\t' << cyl->m_N.Z() << '\n'
                //    << cyl->m_length << '\t' << cyl->m_length_weight << '\n'
                //    << cyl->m_radius << '\t' << cyl->m_radius_weight << '\n';
            }
        }

		
		out << objSet->m_PlaneList.size() << '\n';
		for (int i = 0; i<objSet->m_PlaneList.size(); i++) {
			ObjPlane *ple = objSet->m_PlaneList.at(i);
			out << ple->m_PlaneIndex << '\t' << ple->m_varN << '\n';
			out << ple->m_pO.X() << '\t' << ple->m_pO.Y() << '\t' << ple->m_pO.Z() << '\n'
				<< ple->m_N.X() << '\t' << ple->m_N.Y() << '\t' << ple->m_N.Z() << '\n'
				<< ple->m_dX.X() << '\t' << ple->m_dX.Y() << '\t' << ple->m_dX.Z() << ple->m_sizeConfidence.X() << '\n'
				<< ple->m_dY.X() << '\t' << ple->m_dY.Y() << '\t' << ple->m_dY.Z() << ple->m_sizeConfidence.Y() << '\n'
				<< 0.0 << '\t' << 0.0 << '\t' << 0.0 << 0.0 << '\n';
		}
		out.close();
	}
}
ObjSet* LoadObjSet(const char *file)
{
	if (file) {
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
            int solidType;
            in >> solidType;
            if (solidType == 1) { // Cube
                ObjCube *cube = new ObjCube();
                in >> cube->m_pO.X() >> cube->m_pO.Y() >> cube->m_pO.Z()
                    >> cube->m_dX.X() >> cube->m_dX.Y() >> cube->m_dX.Z() >> cube->m_sizeConfidence.X()
                    >> cube->m_dY.X() >> cube->m_dY.Y() >> cube->m_dY.Z() >> cube->m_sizeConfidence.Y()
                    >> cube->m_dZ.X() >> cube->m_dZ.Y() >> cube->m_dZ.Z() >> cube->m_sizeConfidence.Z();
                objSet->m_SolidList.push_back(cube);
            }
            if (solidType == 2) { // Sylinder
                ObjCylinder *cyl = new ObjCylinder();
                in >> cyl->m_pO.X() >> cyl->m_pO.Y() >> cyl->m_pO.Z()
                    >> cyl->m_N.X() >> cyl->m_N.Y() >> cyl->m_N.Z()
                    >> cyl->m_length >> cyl->m_length_weight
                    >> cyl->m_radius >> cyl->m_radius_weight;
                objSet->m_SolidList.push_back(cyl);
            }
            if (solidType == 3) { // Cone
                //
            }
        }

		
		int nSB;
		in >> nSB;
		for (int i = 0; i<nSB; i++) {
			int nCode;
			double nVar;
			in >> nCode >> nVar;
			ObjPlane *ple = new ObjPlane(nCode);
            ple->m_varN = nVar;

			in >> ple->m_pO.X() >> ple->m_pO.Y() >> ple->m_pO.Z()
				>> ple->m_N.X() >> ple->m_N.Y() >> ple->m_N.Z()
				>> ple->m_dX.X() >> ple->m_dX.Y() >> ple->m_dX.Z() >> ple->m_sizeConfidence.X()
				>> ple->m_dY.X() >> ple->m_dY.Y() >> ple->m_dY.Z() >> ple->m_sizeConfidence.Y()
				>> ple->m_dZ.X() >> ple->m_dZ.Y() >> ple->m_dZ.Z() >> ple->m_sizeConfidence.Z();

            objSet->m_PlaneList.push_back(ple);
		}
		in.close();
		return objSet;
	}
	return 0;
}
