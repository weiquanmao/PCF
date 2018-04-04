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

const std::string ReconInputDir = "../TestData/Ply/Recon";
const std::string ReconOutputDir = "../TestData/Ply/Recon";

const std::string SynInputDir = "../TestData/Ply/Syn/50K_00U_00D";
const std::string SynOutputDir = "../TestData/Ply/Syn/50K_00U_00D";

PCFit::ProType proType = PCFit::Steps_All;

int ProFolder(const string &InputDir, const string &OutputDir, const double &TAng = -1)
{
    DWORD T1, T2;
    T1 = timeGetTime();

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
                const string BaseName = FileName.substr(0, FileName.rfind('.'));
                const string FilePath = InputDir + "/" + FileName;
                const string FileNameOpt = OutputDir + "/" + BaseName + "_Opt.ply";
                const string FileNameRej = OutputDir + "/" + BaseName + "_Rej.ply";
                const string SateStruct = OutputDir + "/" + BaseName + ".struct";

                SetConsoleTitleA(FilePath.c_str());

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
                if (TAng > 0)
                    Fit.Threshold_AngToSurface = TAng;

                Fit.printParams();
                Fit.loadPly(FilePath.c_str());
                Fit.GEOFit(proType);
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
    T2 = timeGetTime();
    clog << "\n\n"
        << "===========================================\n"
        << " # Total Elapsed Time: " << (T2 - T1)*1.0 / 1000 << " Seconds.\n"
        << "===========================================\n"
        << "\n\n";

    return 0;
}

int main(int argc, char *argv[])
{

    string outLog;

    outLog = SynOutputDir + "/Log.txt";
    setLogFile(outLog.c_str());
    ProFolder(SynInputDir, SynOutputDir, 15.0);

    outLog = ReconOutputDir + "/Log.txt";
    setLogFile(outLog.c_str());
    ProFolder(ReconInputDir, ReconOutputDir, 30.0);

    MelodyPlay_Notice();
    system("pause");
    return 0;

}
  