#define PROFILER 1
#include "blast.h"
#include "blast_utilities.hpp"
#include <iostream>


int main() {
    blast::begin_profile();
    using blast::print;
    using std::cout;
    using std::endl;
    {
        blast_time_function;
        {
            blast_time_block("Sleep");
            Sleep(1000);
        }
    }
    blast::end_profile();
    return 0;
}

blast_profiler_end_compilation_unit;
