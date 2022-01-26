
#include <stdio.h>
#include "blast.hpp"

int main() {
    using namespace blast;

    printf("Hello World!\n");

    Mat3 m;
    m(0, 0) = 1;
    m(1, 0) = 2;
    m(2, 0) = 3;

    print(m);

    printf("\nPress enter to exit.\n");
    getchar();
    return 0;
}
