#include "PointCloudFit.h"
#include <windows.h>
#include <iostream>
#include <string>
#include <io.h>

using std::cout;
using std::endl;
using std::string;
#pragma comment(lib,"winmm.lib")

int main(int argc, char *argv[])
{
	const string ResultDir = "../../Result";
	const string ModelDir = "../../3DModels";
	const string szPath = ModelDir + "/*.ply";

	intptr_t hFile = 0;
	struct _finddata_t fileinfo;	
	if ((hFile = _findfirst(szPath.c_str(), &fileinfo)) != -1)
	{
		do
		{
			if ((fileinfo.attrib & _A_ARCH))
			{
				const string FileName = fileinfo.name;
				const string BaseName = FileName.substr(0, FileName.rfind('.'));
				const string FilePath = ModelDir + "/" + FileName;			
				const string FileNameOpt = ResultDir + "/" + BaseName + "_Opt.ply";
				const string FileNameRej = ResultDir + "/" + BaseName + "_Rej.ply";
				const string SateStruct = ResultDir + "/" + BaseName + ".struct";

				cout
					<< "===============================\n"
					<< "[ ----- Point Cloud Fit ----- ]\n"
					<< "[+] " << FilePath << endl
					<< "[-] " << FileNameOpt << endl
					<< "[-] " << FileNameRej << endl
					<< "[-] " << SateStruct << endl
					<< "===============================\n";

				DWORD t1, t2;
				t1 = timeGetTime();
				cout << "\n--------------------------------------------------------------------------------------------\n";
				PCFit Fit;
				Fit.printParams();
				Fit.loadPly(FilePath.c_str());
				Fit.Fit_Sate();
				Fit.clearD();
				Fit.recolorPts(Pt_Noise, 255, 0, 0);
				Fit.savePly(FileNameOpt.c_str());
				Fit.keepPts(Pt_Noise);
				Fit.savePly(FileNameRej.c_str());
				SaveSatellite(Fit.getSate(), SateStruct.c_str());
				cout << "\n--------------------------------------------------------------------------------------------\n";
				t2 = timeGetTime();
				cout << "\n\n"
					<< "=============================\n"
					<< " # Elapsed Time: " << (t2 - t1)*1.0 / 1000 << " Seconds.\n"
					<< "=============================\n"
					<< "\n\n";	
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}

	system("pause");
	return 0;
}
