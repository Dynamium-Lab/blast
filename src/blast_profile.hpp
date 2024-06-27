#pragma once

#if _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <intrin.h>
#include <Windows.h>
#undef NOMINMAX
#undef WIN32_LEAN_AND_MEAN
#else
#include <x86intrin.h>
#include <sys/time.h>
#endif

namespace blast {
#define array_count(array) (sizeof(array)/sizeof((array)[0]))


#if _WIN32
static u64 get_os_freq(void) {
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    return freq.QuadPart;
}
static u64 read_os_timer(void) {
    LARGE_INTEGER time;
    QueryPerformanceCounter(&time);
    return time.QuadPart;
}

#else
static u64 get_os_freq(void) {
    return 1000000;
}
static u64 read_os_timer(void) {
    timeval time;
    gettimeofday(&time, 0);
    u64 result = get_os_freq()*(u64)time.tv_sec + (u64)time.tv_usec;
    return result;
}
#endif

inline u64 read_cpu_timer(void) {
    return __rdtsc();
}

#ifndef PROFILER
#define PROFILER 0
#endif

#ifndef READ_BLOCK_TIMER
#define READ_BLOCK_TIMER read_cpu_timer
#endif

#if PROFILER

struct ProfileAnchor {
    u64 TSC_elapsed_exclusive;
    u64 TSC_elapsed_inclusive;
    u64 hit_count;
    u64 byte_count;
    char const *label;
};
static ProfileAnchor g_profiler_anchors[4096];
static u32 g_profiler_parent;

struct ProfileBlock {
    ProfileBlock(char const *label_, u32 anchor_id_, u64 byte_count) {
        parent_id = g_profiler_parent;

        anchor_id = anchor_id_;
        label = label_;

        ProfileAnchor *anchor = g_profiler_anchors + anchor_id;
        old_TSC_elapsed_inclusive = anchor->TSC_elapsed_inclusive;
        anchor->byte_count += byte_count;

        g_profiler_parent = anchor_id;
        start_TSC = READ_BLOCK_TIMER();
    }

    ~ProfileBlock(void) {
        u64 elapsed = READ_BLOCK_TIMER() - start_TSC;
        g_profiler_parent = parent_id;

        ProfileAnchor *parent = g_profiler_anchors + parent_id;
        ProfileAnchor *anchor = g_profiler_anchors + anchor_id;

        parent->TSC_elapsed_exclusive -= elapsed;
        anchor->TSC_elapsed_exclusive += elapsed;
        anchor->TSC_elapsed_inclusive = old_TSC_elapsed_inclusive + elapsed;
        anchor->hit_count++;
        anchor->label = label;
    }

    char const *label;
    u64 old_TSC_elapsed_inclusive;
    u64 start_TSC;
    u32 parent_id;
    u32 anchor_id;
};

#define NameConcat2(A, B) A##B
#define NameConcat(A, B) NameConcat2(A, B)
#define blast_time_bandwidth(Name, bytes) blast::ProfileBlock NameConcat(Block, __LINE__)(Name, __COUNTER__ + 1, bytes);
#define blast_profiler_end_compilation_unit static_assert(__COUNTER__ < array_count(blast::g_profiler_anchors), "Number of profile points exceeds size of profiler::anchors array")

static void print_time_elapsed(u64 total_elapsed, u64 freq, ProfileAnchor *anchor) {
    double percent = 100.0 * ((double)anchor->TSC_elapsed_exclusive / (double)total_elapsed);
    printf("  %s[%I64u]: %I64u (%.2f%%", anchor->label, anchor->hit_count, anchor->TSC_elapsed_exclusive, percent);
    if(anchor->TSC_elapsed_inclusive != anchor->TSC_elapsed_exclusive) {
        double percent_w_children = 100.0 * ((double)anchor->TSC_elapsed_inclusive / (double)total_elapsed);
        printf(", %.2f%% w/children", percent_w_children);
    }
    printf(")");

    if(anchor->byte_count) {
        double MB = 1024.0*1024.0;
        double GB = MB*1024.0;

        double seconds = (double)anchor->TSC_elapsed_inclusive / (double)freq;
        double bytes_per_second = (double)anchor->byte_count / seconds;
        double megabytes = (double)anchor->byte_count / MB;
        double GB_per_second = bytes_per_second / GB;

        printf("  %.3fmb at %.2fgb/s", megabytes, GB_per_second);
    }

    printf("\n");
}

static void print_anchor(u64 total_CPU_elapsed, u64 freq) {
    for(u32 anchor_id = 0; anchor_id < array_count(g_profiler_anchors); ++anchor_id) {
        ProfileAnchor *anchor = g_profiler_anchors + anchor_id;
        if(anchor->TSC_elapsed_inclusive) {
            print_time_elapsed(total_CPU_elapsed, freq, anchor);
        }
    }
}

#else

#define blast_time_bandwidth(...)
#define print_anchor(...)
#define blast_profiler_end_compilation_unit

#endif

struct Profiler {
    u64 start_TSC;
    u64 end_TSC;
};
static Profiler g_profiler;

#define blast_time_block(Name) blast_time_bandwidth(Name, 0)
#define blast_time_function blast_time_block(__func__)

static u64 estimate_block_freq(void) {
    u64 ms_to_wait = 100;
    u64 os_freq = get_os_freq();

    u64 block_start = READ_BLOCK_TIMER();
    u64 os_start = read_os_timer();
    u64 os_end = 0;
    u64 os_elapsed = 0;
    u64 os_wait_time = os_freq * ms_to_wait / 1000;
    while(os_elapsed < os_wait_time) {
        os_end = read_os_timer();
        os_elapsed = os_end - os_start;
    }

    u64 block_end = READ_BLOCK_TIMER();
    u64 block_elapsed = block_end - block_start;

    u64 block_freq = 0;
    if(os_elapsed) {
        block_freq = os_freq * block_elapsed / os_elapsed;
    }

    return block_freq;
}

static void begin_profile(void) {
    g_profiler.start_TSC = READ_BLOCK_TIMER();
}

static void end_profile() {
    g_profiler.end_TSC = READ_BLOCK_TIMER();
    u64 freq = estimate_block_freq();

    u64 total_elapsed = g_profiler.end_TSC - g_profiler.start_TSC;

    if(freq) {
        printf("\nTotal time: %0.4fms (timer freq %I64u)\n", 1000.0 * (double)total_elapsed / (double)freq, freq);
    }

    print_anchor(total_elapsed, freq);
}

} // namespace blast
