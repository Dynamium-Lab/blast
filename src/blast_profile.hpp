#pragma once

#if _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <intrin.h>
#include <Windows.h>
#undef NOMINMAX
#undef WIN32_LEAN_AND_MEAN
#endif

namespace blast {
#define ArrayCount(Array) (sizeof(Array)/sizeof((Array)[0]))


#if _WIN32

static u64 GetOSTimerFreq(void) {
    LARGE_INTEGER Freq;
    QueryPerformanceFrequency(&Freq);
    return Freq.QuadPart;
}

static u64 ReadOSTimer(void) {
    LARGE_INTEGER Value;
    QueryPerformanceCounter(&Value);
    return Value.QuadPart;
}

#else

#include <x86intrin.h>
#include <sys/time.h>

static u64 GetOSTimerFreq(void) {
    return 1000000;
}

static u64 ReadOSTimer(void) {
    // NOTE(casey): The "struct" keyword is not necessary here when compiling in C++,
    // but just in case anyone is using this file from C, I include it.
    struct timeval Value;
    gettimeofday(&Value, 0);

    u64 Result = GetOSTimerFreq()*(u64)Value.tv_sec + (u64)Value.tv_usec;
    return Result;
}

#endif

inline u64 ReadCPUTimer(void) {
    return __rdtsc();
}

static u64 EstimateCPUTimerFreq(void) {
    u64 MillisecondsToWait = 100;
    u64 OSFreq = GetOSTimerFreq();

    u64 CPUStart = ReadCPUTimer();
    u64 OSStart = ReadOSTimer();
    u64 OSEnd = 0;
    u64 OSElapsed = 0;
    u64 OSWaitTime = OSFreq * MillisecondsToWait / 1000;
    while(OSElapsed < OSWaitTime) {
        OSEnd = ReadOSTimer();
        OSElapsed = OSEnd - OSStart;
    }

    u64 CPUEnd = ReadCPUTimer();
    u64 CPUElapsed = CPUEnd - CPUStart;

    u64 CPUFreq = 0;
    if(OSElapsed) {
        CPUFreq = OSFreq * CPUElapsed / OSElapsed;
    }

    return CPUFreq;
}

#ifndef PROFILER
#define PROFILER 0
#endif

#ifndef READ_BLOCK_TIMER
#define READ_BLOCK_TIMER ReadCPUTimer
#endif

#if PROFILER

struct profile_anchor {
    u64 TSCElapsedExclusive;
    u64 TSCElapsedInclusive;
    u64 HitCount;
    u64 ProcessedByteCount;
    char const *Label;
};
static profile_anchor GlobalProfilerAnchors[4096];
static u32 GlobalProfilerParent;

struct profile_block {
    profile_block(char const *Label_, u32 AnchorIndex_, u64 ByteCount) {
        ParentIndex = GlobalProfilerParent;

        AnchorIndex = AnchorIndex_;
        Label = Label_;

        profile_anchor *Anchor = GlobalProfilerAnchors + AnchorIndex;
        OldTSCElapsedInclusive = Anchor->TSCElapsedInclusive;
        Anchor->ProcessedByteCount += ByteCount;

        GlobalProfilerParent = AnchorIndex;
        StartTSC = READ_BLOCK_TIMER();
    }

    ~profile_block(void) {
        u64 Elapsed = READ_BLOCK_TIMER() - StartTSC;
        GlobalProfilerParent = ParentIndex;

        profile_anchor *Parent = GlobalProfilerAnchors + ParentIndex;
        profile_anchor *Anchor = GlobalProfilerAnchors + AnchorIndex;

        Parent->TSCElapsedExclusive -= Elapsed;
        Anchor->TSCElapsedExclusive += Elapsed;
        Anchor->TSCElapsedInclusive = OldTSCElapsedInclusive + Elapsed;
        ++Anchor->HitCount;

        /* NOTE(casey): This write happens every time solely because there is no
           straightforward way in C++ to have the same ease-of-use. In a better programming
           language, it would be simple to have the anchor points gathered and labeled at compile
           time, and this repetative write would be eliminated. */
        Anchor->Label = Label;
    }

    char const *Label;
    u64 OldTSCElapsedInclusive;
    u64 StartTSC;
    u32 ParentIndex;
    u32 AnchorIndex;
};

#define NameConcat2(A, B) A##B
#define NameConcat(A, B) NameConcat2(A, B)
#define TimeBandwidth(Name, bytes) blast::profile_block NameConcat(Block, __LINE__)(Name, __COUNTER__ + 1, bytes);
#define ProfilerEndOfCompilationUnit static_assert(__COUNTER__ < ArrayCount(blast::GlobalProfilerAnchors), "Number of profile points exceeds size of profiler::Anchors array")

static void PrintTimeElapsed(u64 TotalTSCElapsed, u64 TimerFreq, profile_anchor *Anchor) {
    double Percent = 100.0 * ((double)Anchor->TSCElapsedExclusive / (double)TotalTSCElapsed);
    printf("  %s[%llu]: %llu (%.2f%%", Anchor->Label, Anchor->HitCount, Anchor->TSCElapsedExclusive, Percent);
    if(Anchor->TSCElapsedInclusive != Anchor->TSCElapsedExclusive) {
        double PercentWithChildren = 100.0 * ((double)Anchor->TSCElapsedInclusive / (double)TotalTSCElapsed);
        printf(", %.2f%% w/children", PercentWithChildren);
    }
    printf(")");

    if(Anchor->ProcessedByteCount) {
        double Megabyte = 1024.0f*1024.0f;
        double Gigabyte = Megabyte*1024.0f;

        double Seconds = (double)Anchor->TSCElapsedInclusive / (double)TimerFreq;
        double BytesPerSecond = (double)Anchor->ProcessedByteCount / Seconds;
        double Megabytes = (double)Anchor->ProcessedByteCount / (double)Megabyte;
        double GigabytesPerSecond = BytesPerSecond / Gigabyte;

        printf("  %.3fmb at %.2fgb/s", Megabytes, GigabytesPerSecond);
    }

    printf("\n");
}

static void PrintAnchorData(u64 TotalCPUElapsed, u64 TimerFreq) {
    for(u32 AnchorIndex = 0; AnchorIndex < ArrayCount(GlobalProfilerAnchors); ++AnchorIndex) {
        profile_anchor *Anchor = GlobalProfilerAnchors + AnchorIndex;
        if(Anchor->TSCElapsedInclusive) {
            PrintTimeElapsed(TotalCPUElapsed, TimerFreq, Anchor);
        }
    }
}

#else

#define TimeBandwidth(...)
#define PrintAnchorData(...)
#define ProfilerEndOfCompilationUnit

#endif

struct profiler {
    u64 StartTSC;
    u64 EndTSC;
};
static profiler GlobalProfiler;

#define TimeBlock(Name) TimeBandwidth(Name, 0)
#define TimeFunction TimeBlock(__func__)

static u64 EstimateBlockTimerFreq(void) {
    (void)&EstimateCPUTimerFreq; // NOTE(casey): This has to be voided here to prevent compilers from warning us that it is not used

    u64 MillisecondsToWait = 100;
    u64 OSFreq = GetOSTimerFreq();

    u64 BlockStart = READ_BLOCK_TIMER();
    u64 OSStart = ReadOSTimer();
    u64 OSEnd = 0;
    u64 OSElapsed = 0;
    u64 OSWaitTime = OSFreq * MillisecondsToWait / 1000;
    while(OSElapsed < OSWaitTime) {
        OSEnd = ReadOSTimer();
        OSElapsed = OSEnd - OSStart;
    }

    u64 BlockEnd = READ_BLOCK_TIMER();
    u64 BlockElapsed = BlockEnd - BlockStart;

    u64 BlockFreq = 0;
    if(OSElapsed) {
        BlockFreq = OSFreq * BlockElapsed / OSElapsed;
    }

    return BlockFreq;
}

static void BeginProfile(void) {
    GlobalProfiler.StartTSC = READ_BLOCK_TIMER();
}

static void EndAndPrintProfile() {
    GlobalProfiler.EndTSC = READ_BLOCK_TIMER();
    u64 TimerFreq = EstimateBlockTimerFreq();

    u64 TotalTSCElapsed = GlobalProfiler.EndTSC - GlobalProfiler.StartTSC;

    if(TimerFreq) {
        printf("\nTotal time: %0.4fms (timer freq %llu)\n", 1000.0 * (double)TotalTSCElapsed / (double)TimerFreq, TimerFreq);
    }

    PrintAnchorData(TotalTSCElapsed, TimerFreq);
}

} // namespace blast
