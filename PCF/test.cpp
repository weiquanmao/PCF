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

const std::string ReconInputDir = "../../../Data/PCModelRecon";
const std::string ReconOutputDir = "../../../Data/ResultRecon";

const std::string SynInputDir = "../../../Data/PCModelSyn";
const std::string SynOutputDir = "../../../Data/ResultSyn";

enum ProMe { ProRecon, ProSynAll, ProSynOne };


const ProMe proMe = ProRecon;
PCFit::ProType proType = PCFit::Steps_All;
const bool bOnlySpecialTar = false;

// For Special Target
const char SpecialTar[] = "spot.ply";
// For One Syn
const int num_one = 50000;
const int ndis_one = 0;
const int nang_one = 5;
// For All Syn
const int _K = 5;
const int _M = 5;
const int _N = 5;
unsigned int num_all[_K] = { 50000, 20000, 10000, 5000, 2000 };
unsigned int ndis_all[_M] = { 0, 1, 2, 4, 8 };
unsigned int nang_all[_N] = { 0, 5, 10, 15, 30 };




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
                if (bOnlySpecialTar && FileName != SpecialTar)
                    continue;
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
void testSynOne()
{
    DWORD T1, T2;
    T1 = timeGetTime();

    char _folderName[1024];
    if (num_one > 1000 && (num_one % 1000) == 0)
        sprintf(_folderName, "%dK_%02dU_%02dD", num_one / 1000, ndis_one, nang_one);
    else
        sprintf(_folderName, "%d_%02dU_%02dD", num_one / 1000, ndis_one, nang_one);

    const string folderName = _folderName;
    const string inFolder = SynInputDir + "/" + folderName;
    const string outFolder = SynOutputDir + "/" + folderName;
    const string outLog = SynOutputDir + "/" + folderName + "/Log.txt";
    cout
        << endl
        << "===============================================================================\n"
        << "===============================================================================\n"
        << "[<] " << inFolder << "\n"
        << "[>] " << outFolder << "\n"
        << "[>] " << outLog << "\n"
        << "===============================================================================\n"
        << "===============================================================================\n"
        << endl;
    CreateDirectoryA(outFolder.c_str(), NULL);
    double TAng = 15;
    if (nang_one >= 10)
        TAng = 30;
    setLogFile(outLog.c_str());
    ProFolder(inFolder, outFolder, TAng);

    cout << endl
        << "===============================================================================\n"
        << "===============================================================================\n"
        << " # Total Elapsed Time: " << (T2 - T1)*1.0 / 1000 << " Seconds.\n"
        << "===============================================================================\n"
        << "===============================================================================\n"
        << endl;

}
void testSynAll()
{
    DWORD T1, T2;
    T1 = timeGetTime();

    const unsigned nnn = _K*_M*_N;
    unsigned int _num, _ndis, _nang;

    for (int idx = 0; idx < nnn; ++idx) {
        unsigned k = idx / (_M*_N);
        unsigned i = (idx % (_M*_N)) / _N;
        unsigned j = idx % _N;
        unsigned _num = num_all[k];
        unsigned _ndis = ndis_all[i];
        unsigned _nang = nang_all[j];

        char _folderName[1024];
        if (_num > 1000 && (_num % 1000) == 0)
            sprintf(_folderName, "%dK_%02dU_%02dD", _num / 1000, _ndis, _nang);
        else
            sprintf(_folderName, "%d_%02dU_%02dD", _num / 1000, _ndis, _nang);

        const string folderName = _folderName;
        const string inFolder = SynInputDir + "/" + folderName;
        const string outFolder = SynOutputDir + "/" + folderName;
        const string outLog = SynOutputDir + "/" + folderName + "/Log.txt";
        cout
            << endl
            << "===============================================================================\n"
            << "===============================================================================\n"
            << "[<] " << inFolder << "\n"
            << "[>] " << outFolder << "\n"
            << "[>] " << outLog << "\n"
            << "===============================================================================\n"
            << "===============================================================================\n"
            << endl;
        CreateDirectoryA(outFolder.c_str(), NULL);
        double TAng = 15;
        if (_nang >= 10)
            TAng = 30;
        setLogFile(outLog.c_str());
        ProFolder(inFolder, outFolder, TAng);
    }

    cout << endl
        << "===============================================================================\n"
        << "===============================================================================\n"
        << " # Total Elapsed Time: " << (T2 - T1)*1.0 / 1000 << " Seconds.\n"
        << "===============================================================================\n"
        << "===============================================================================\n"
        << endl;
}
int main(int argc, char *argv[])
{
    switch (proMe)
    {
    case ProRecon: {
        const string outLog = ReconOutputDir + "/" + "/Log.txt";
        setLogFile(outLog.c_str());
        ProFolder(ReconInputDir, ReconOutputDir, 30.0);
        break;
    }
    case ProSynOne:
        testSynOne();
        break;
    case ProSynAll:
        testSynAll();
        break;
    default:
        break;
    }
    SetConsoleTitleA("-- Finished --");
    MelodyPlay_Notice();
    system("pause");
    return 0;
}

