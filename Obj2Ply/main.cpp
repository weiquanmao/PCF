#include <stdio.h>
#include <omp.h>

#include "obj2ply.h"

const std::string InputDir = "../../../Data/[MeshModelPruned]";
const std::string OutputDir = "../../../Data/SynData";

int main()
{
	const int _K = 5;
	const int _M = 5;
	const int _N = 5;
	unsigned int num[_K] = { 50000, 20000, 10000, 5000, 2000 };
	int ndis[_M] = { 0, 1, 2, 4, 8 };
	int nang[_N] = { 0, 5, 10, 15, 30 };

	//const int _K = 5;
	//const int _M = 5;
	//const int _N = 1;
	//unsigned int num[_K] = { 50000, 20000, 10000, 5000, 2000 };
	//int ndis[_M] = { 0, 1, 2, 4, 8 };
	//int nang[_N] = { 0 };

	const unsigned nnn = _K*_M*_N;
	unsigned int _num, _ndis, _nang;

	int nThreadMax = omp_get_max_threads();
	omp_set_num_threads(nThreadMax);
#pragma omp parallel for
	for (int idx = 0; idx < nnn; ++idx) {
		unsigned k = idx / (_M*_N);
		unsigned i = (idx % (_M*_N)) / _N;
		unsigned j = idx % _N;
		unsigned _num = num[k];
		unsigned _ndis = ndis[i];
		unsigned _nang = nang[j];
		obj2ply(InputDir, OutputDir, _num, _ndis, _nang);
	}
    system("pause");
    return 0;
}
