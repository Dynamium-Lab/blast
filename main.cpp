#define PROFILER 1
#include "blast.h"
#include "blast_utilities.hpp"
#include <chrono>
#include <iostream>
#include <thread>

int main() {
    blast::begin_profile();

    using namespace std::chrono_literals;
    using blast::print;
    using std::cout;
    using std::endl;

    {
        blast_time_function;
        {
            blast_time_block("Sleep");
            std::this_thread::sleep_for(1s);
        }
    }

    blast::end_profile();
    return 0;
}

blast_profiler_end_compilation_unit;
