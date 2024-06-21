#define PROFILER 1
#include "blast.h"
#include "blast_utilities.h"
#include <iostream>


int main() {
    blast::BeginProfile();
    using blast::print;
    using std::cout;
    using std::endl;
    {
        TimeBlock("Main");
        {
            TimeBlock("Sleep");
            Sleep(1000);
        }
    }
    blast::EndAndPrintProfile();
    return 0;
}
ProfilerEndOfCompilationUnit;
