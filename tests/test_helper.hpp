#pragma once

#include <blast>

namespace blast::test {

// Absolute tolerance for "exact" numeric checks in tests, scaled to the active
// precision of blast::real. Float carries ~7 significant digits (machine eps
// ~1.2e-7), so the 1e-9 absolute bound that is appropriate for double is simply
// unattainable in single precision. The quantities checked here are O(1)..O(1e2),
// so 1e-4 comfortably absorbs float rounding while still catching real errors.
#if BLAST_USE_DOUBLES
constexpr real abs_tol = 1e-9;
#else
constexpr real abs_tol = 1e-4;
#endif

// Looser tolerance for inherently approximate checks: finite-difference vs.
// analytical derivatives and accumulated basis-function sums. Float finite
// differencing with a tiny dt loses most of its significant digits to
// cancellation, so it needs a far looser bound than double. Double keeps the
// historical 1e-5 (== BLAST_EPSILON) so those checks stay as strict as before.
#if BLAST_USE_DOUBLES
constexpr real approx_tol = 1e-5;
#else
constexpr real approx_tol = 1e-1;
#endif

} // namespace blast::test
