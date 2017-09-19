#ifndef _SATELLITE_H_FILE_
#define _SATELLITE_H_FILE_

#include <vector>
#include <vcg/space/deprecated_point3.h>

class Sailboard
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
	Sailboard(int index)
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
class Mainbody
{
public:
	enum MainBodyType {
		MainBody_Cube,
		MainBody_Cylinder
	};
	virtual MainBodyType type() const = 0;
};
class MainBodyCube : public Mainbody
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
	MainBodyCube()
		: Mainbody() {
		m_sizeConfidence.SetZero();
		m_pO.SetZero();
		m_dX.SetZero();
		m_dY.SetZero();
		m_dZ.SetZero();
	}
	MainBodyType type() const { return MainBody_Cube; }
};
class MainBodyCylinder : public Mainbody
{
public:
	double m_radius;
	double m_length;
	double m_radius_weight;
	double m_length_weight;
	vcg::Point3f m_pO;
	vcg::Point3f m_N;
	MainBodyCylinder()
		: Mainbody() {
		m_radius = 0.0;
		m_length = 0.0;
		m_N.SetZero();
		m_pO.SetZero();
	}
	MainBodyType type() const { return MainBody_Cylinder; }
};
class Satellite {
public:
	std::vector<Sailboard*> m_SailboradList;
	Mainbody *m_Mainbody;
	Satellite() :m_Mainbody(0) {}
	~Satellite() {
		clean();
	}

	void setMainbody(Mainbody *mainBody) {
		cleanMainbody();
		m_Mainbody = mainBody;
	}
	bool containBord(Sailboard *board) {
		if (board == 0)
			return false;
		for (int i = 0; i<m_SailboradList.size(); i++) {
			if (m_SailboradList.at(i) == board)
				return true;
		}
		return false;
	}
	void addSailbord(Sailboard *board) {
		if (board != 0 && !containBord(board))
			m_SailboradList.push_back(board);
	}
	void delSailbord(Sailboard *board) {
		std::vector<Sailboard*>::iterator it;
		for (it = m_SailboradList.begin(); it != m_SailboradList.end(); it++) {
			if (*it == board) {
				m_SailboradList.erase(it);
				break;
			}
		}

	}
	void cleanMainbody() {
		if (m_Mainbody != 0)
			delete m_Mainbody;
	}
	void cleanSailbord() {
		for (int i = 0; i<m_SailboradList.size(); i++)
			delete m_SailboradList.at(i);
		m_SailboradList.clear();
	}
	void clean() {
		cleanMainbody();
		cleanSailbord();
	}
};

void SaveSatellite(Satellite *sate, const char *file);
Satellite* LoadSatellite(const char *file);

#endif // !_SATELLITE_H_FILE_
