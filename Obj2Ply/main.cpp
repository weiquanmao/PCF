#include <stdio.h>
#include <omp.h>

#include "obj2ply.h"

const std::string InputDir = "../TestData/Obj";
const std::string OutputDir = "../TestData/Ply/Syn";

int main()
{
#if 1
    obj2ply(InputDir, OutputDir, 50000);
#else
	const int _K = 5;
	const int _M = 4;
	const int _N = 4;
	unsigned int num[_K] = { 50000, 20000, 10000, 5000, 2000 };
	int ndis[_M] = { 0, 1, 2, 4 };
	int nang[_N] = { 0, 5, 10, 15 };

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
#endif
    system("pause");
    return 0;
}
