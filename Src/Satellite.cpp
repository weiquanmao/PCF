#include "Satellite.h"
#include <fstream>
#include <iomanip>
#include <string>

void SaveSatellite(Satellite *sate, const char *file)
{
	if (sate != 0 && (sate->m_Mainbody != 0 || !sate->m_SailboradList.empty()) &&
		!file) {
		std::ofstream out(file, std::ios::out);
		out << std::fixed << std::setprecision(10);
		out << "#SateStruct 1.0" << '\n'
			<< "#by LEORecon" << '\n'
			<< "#Mainbody 1:Cube 2:Cylinder" << '\n'
			<< "#SBs" << '\n'
			<< '\n';

		if (sate->m_Mainbody == 0)
			out << "0" << '\n';
		else {
			if (sate->m_Mainbody->type() == Mainbody::MainBody_Cube) {
				MainBodyCube *cube = (MainBodyCube*)sate->m_Mainbody;
				out << "1" << ' ' << "1" << "\n";
				out << cube->m_pO.X() << '\t' << cube->m_pO.Y() << '\t' << cube->m_pO.Z() << '\n'
					<< cube->m_dX.X() << '\t' << cube->m_dX.Y() << '\t' << cube->m_dX.Z() << cube->m_sizeConfidence.X() << '\n'
					<< cube->m_dY.X() << '\t' << cube->m_dY.Y() << '\t' << cube->m_dY.Z() << cube->m_sizeConfidence.Y() << '\n'
					<< cube->m_dZ.X() << '\t' << cube->m_dZ.Y() << '\t' << cube->m_dZ.Z() << cube->m_sizeConfidence.Z() << '\n';
			}
			if (sate->m_Mainbody->type() == Mainbody::MainBody_Cylinder) {
				MainBodyCylinder *cyl = (MainBodyCylinder*)sate->m_Mainbody;
				out << "1" << ' ' << "2" << "\n";
				out << cyl->m_pO.X() << '\t' << cyl->m_pO.Y() << '\t' << cyl->m_pO.Z() << '\n'
					<< cyl->m_N.X() << '\t' << cyl->m_N.Y() << '\t' << cyl->m_N.Z() << '\n'
					<< cyl->m_length << '\t' << cyl->m_length_weight << '\n'
					<< cyl->m_radius << '\t' << cyl->m_radius_weight << '\n';
			}
		}

		out << sate->m_SailboradList.size() << '\n';
		for (int i = 0; i<sate->m_SailboradList.size(); i++) {
			Sailboard *sb = sate->m_SailboradList.at(i);
			out << sb->m_PlaneIndex << '\t' << sb->m_varN << '\n';
			out << sb->m_pO.X() << '\t' << sb->m_pO.Y() << '\t' << sb->m_pO.Z() << '\n'
				<< sb->m_N.X() << '\t' << sb->m_N.Y() << '\t' << sb->m_N.Z() << '\n'
				<< sb->m_dX.X() << '\t' << sb->m_dX.Y() << '\t' << sb->m_dX.Z() << sb->m_sizeConfidence.X() << '\n'
				<< sb->m_dY.X() << '\t' << sb->m_dY.Y() << '\t' << sb->m_dY.Z() << sb->m_sizeConfidence.Y() << '\n'
				<< 0.0 << '\t' << 0.0 << '\t' << 0.0 << 0.0 << '\n';
		}
		out.close();
	}
}
Satellite* LoadSatellite(const char *file)
{
	if (!file) {
		std::ifstream in(file, std::ios::in);
		if (!in.is_open())
			return 0;
		Satellite *sate = new Satellite();

		std::string taken;
		while (in.peek() == '#') {
			std::getline(in, taken);
			if (in.peek() == '\n')
				in.get();
		}
		int nMain, typeMain;
		in >> nMain;
		if (nMain == 1)
		{
			in >> typeMain;
			if (typeMain == 1) {
				MainBodyCube *cube = new MainBodyCube();
				in >> cube->m_pO.X() >> cube->m_pO.Y() >> cube->m_pO.Z()
					>> cube->m_dX.X() >> cube->m_dX.Y() >> cube->m_dX.Z() >> cube->m_sizeConfidence.X()
					>> cube->m_dY.X() >> cube->m_dY.Y() >> cube->m_dY.Z() >> cube->m_sizeConfidence.Y()
					>> cube->m_dZ.X() >> cube->m_dZ.Y() >> cube->m_dZ.Z() >> cube->m_sizeConfidence.Z();
				sate->m_Mainbody = cube;
			}
			if (typeMain == 2) {
				MainBodyCylinder *cyl = new MainBodyCylinder();
				in >> cyl->m_pO.X() >> cyl->m_pO.Y() >> cyl->m_pO.Z()
					>> cyl->m_N.X() >> cyl->m_N.Y() >> cyl->m_N.Z()
					>> cyl->m_length >> cyl->m_length_weight
					>> cyl->m_radius >> cyl->m_radius_weight;
				sate->m_Mainbody = cyl;
			}
		}

		int nSB;
		in >> nSB;
		for (int i = 0; i<nSB; i++) {
			int nCode;
			double nVar;
			in >> nCode >> nVar;
			Sailboard *sb = new Sailboard(nCode);
			sb->m_varN = nVar;

			in >> sb->m_pO.X() >> sb->m_pO.Y() >> sb->m_pO.Z()
				>> sb->m_N.X() >> sb->m_N.Y() >> sb->m_N.Z()
				>> sb->m_dX.X() >> sb->m_dX.Y() >> sb->m_dX.Z() >> sb->m_sizeConfidence.X()
				>> sb->m_dY.X() >> sb->m_dY.Y() >> sb->m_dY.Z() >> sb->m_sizeConfidence.Y()
				>> sb->m_dZ.X() >> sb->m_dZ.Y() >> sb->m_dZ.Z() >> sb->m_sizeConfidence.Z();

			sate->m_SailboradList.push_back(sb);
		}
		in.close();
		return sate;
	}
	return 0;
}
