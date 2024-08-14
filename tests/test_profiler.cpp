#include "blast.h"
#include <chrono>
#include <iostream>
#include <thread>

int main() {
    using namespace std::chrono_literals;

    blast::begin_profile();
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
