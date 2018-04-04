#include "struct2ply.h"

#include <string>
#include <io.h>

using namespace std;

//const string InDir = "../TestData/Ply/Recon";
//const string OutDir = "../TestData/Ply/Recon";

const string InDir = "../TestData/Ply/Syn/50K_00U_00D";
const string OutDir = "../TestData/Ply/Syn/50K_00U_00D";

int main(int argc, char *argv[])
{
    
    const string szPath = InDir + "/*.struct";
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
                const string inFile = InDir + "/" + FileName;
                const string plyFile = OutDir + "/" + BaseName + "_Opt.ply";
                const string outFile = OutDir + "/" + BaseName + "_sample.ply";
                //----
                int K = PlyPts(plyFile);
                struct2ply(inFile, K, outFile);
                //----
            }
        } while (_findnext(hFile, &fileinfo) == 0);
    }
    system("pause");
    return 0;
}
