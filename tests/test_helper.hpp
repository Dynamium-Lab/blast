#pragma once

#include <blast>

// Assert is active under the same condition as in blast_error.hpp.
#if defined(BLAST_DEBUG) || defined(_DEBUG)
#define BLAST_TEST_ASSERTS_ON 1
#else
#define BLAST_TEST_ASSERTS_ON 0
#endif

// A failing Assert aborts the process, so it is checked in a forked child (POSIX only).
#if BLAST_TEST_ASSERTS_ON && !defined(_WIN32)
#define BLAST_TEST_DEATH_TESTS 1
#include <csignal>
#include <sys/wait.h>
#include <unistd.h>
#else
#define BLAST_TEST_DEATH_TESTS 0
#endif

namespace blast::test {

#if BLAST_TEST_DEATH_TESTS
/** @brief True when f aborts (a failing Assert), run in a child process. */
template<typename F>
bool aborts(F&& f) {
  const pid_t pid = fork();
  if (pid == 0) {
    (void) !freopen("/dev/null", "w", stderr); // hide the expected assertion message
    f();
    _exit(0);
  }
  int status = 0;
  waitpid(pid, &status, 0);
  return WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT;
}
#endif

// Absolute tolerance for "exact" numeric checks in tests, scaled to the active
// precision of blast::real. Float carries ~7 significant digits (machine eps
// ~1.2e-7), so the 1e-9 absolute bound that is appropriate for double is simply
// unattainable in single precision. The quantities checked here are O(1)..O(1e2),
// so 1e-4 comfortably absorbs float rounding while still catching real errors.
#if BLAST_USE_DOUBLES
constexpr real abs_tol = 1e-9;
#else
constexpr real abs_tol = (real) 1e-4;
#endif

// Looser tolerance for inherently approximate checks: finite-difference vs.
// analytical derivatives and accumulated basis-function sums. Float finite
// differencing with a tiny dt loses most of its significant digits to
// cancellation, so it needs a far looser bound than double. Double keeps the
// historical 1e-5 (== BLAST_EPSILON) so those checks stay as strict as before.
#if BLAST_USE_DOUBLES
constexpr real approx_tol = 1e-5;
#else
constexpr real approx_tol = (real) 1e-1;
#endif

} // namespace blast::test
