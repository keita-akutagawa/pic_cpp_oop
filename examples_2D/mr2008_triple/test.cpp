#include <omp.h>
#include <vector>


int main()
{
    int size = 100000000;
    std::vector<double> a(size, 0.0);
    std::vector<double> b(size, 10.0);
    std::vector<double> c(size, 110.0);

    #pragma omp parallel for
    for (int i = 0; i < size; i++) {
        a[i] = b[i] + c[i];
    }

}


