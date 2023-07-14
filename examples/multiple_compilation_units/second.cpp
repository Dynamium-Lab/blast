#include "second.h"

int something_that_uses_blast(){
    blast::Array a(3);
    a[0] = 1;
    a[1] = 2;
    a[2] = 3;
    return (int)blast::sum(a);
}
