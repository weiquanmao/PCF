#include "PointCloudFit.h"
#include "utility/flog.h"
#include "utility/melody.h"
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
#if 0 // Reconstructed Models
	const string InputDir = "../../../Data/PCModelRecon";
	const string OutputDir = "../../../Data/ResultRecon";
#else // Synthesized Models
    const string InputDir = "../../../Data/PCModelSyn/50K";
    const string OutputDir = "../../../Data/ResultSyn";
#endif
	const string szPath = InputDir + "/*.ply";

	intptr_t hFile = 0;
	struct _finddata_t fileinfo;	
	if ((hFile = _findfirst(szPath.c_str(), &fileinfo)) != -1)
	{
		do
		{
			if ((fileinfo.attrib & _A_ARCH))
			{
				const string FileName = fileinfo.name;
                SetConsoleTitleA(FileName.c_str());
#if 0 // Test Only One
				if (FileName != "dsp.ply")
					continue;
#endif
				const string BaseName = FileName.substr(0, FileName.rfind('.'));
				const string FilePath = InputDir + "/" + FileName;
				const string FileNameOpt = OutputDir + "/" + BaseName + "_Opt.ply";
				const string FileNameRej = OutputDir + "/" + BaseName + "_Rej.ply";
				const string SateStruct = OutputDir + "/" + BaseName + ".struct";

                clog
					<< "===============================\n"
					<< "[ ----- Point Cloud Fit ----- ]\n"
					<< "[+] " << FilePath << '\n'
					<< "[-] " << FileNameOpt << '\n'
					<< "[-] " << FileNameRej << '\n'
					<< "[-] " << SateStruct << '\n'
					<< "===============================\n";
                
				DWORD t1, t2;
				t1 = timeGetTime();
                clog << "\n--------------------------------------------------------------------------------------------\n";
				PCFit Fit;
				Fit.printParams();
				Fit.loadPly(FilePath.c_str());
				Fit.Fit_Sate();
				Fit.clearD();
				Fit.recolorPts(Pt_Noise, 255, 0, 0);
				Fit.savePly(FileNameOpt.c_str());
				Fit.autoColor();
				Fit.savePly(FileNameRej.c_str());
				SaveObjSet(Fit.getGEOObjSet(), SateStruct.c_str());
                clog << "\n--------------------------------------------------------------------------------------------\n";
				t2 = timeGetTime();
                clog << "\n\n"
					<< "=============================\n"
					<< " # Elapsed Time: " << (t2 - t1)*1.0 / 1000 << " Seconds.\n"
					<< "=============================\n"
					<< "\n\n";	
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
    
    // MelodyPlay_Notice();
    SetConsoleTitleA("-- Finished --");
	system("pause");
	return 0;
}
