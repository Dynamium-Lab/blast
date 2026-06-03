# C++ API reference

This reference is generated from the C++ headers via Doxygen + Breathe. Documentation
quality tracks the Doxygen comments in the source — types with `/** @brief ... */` blocks
render richest. The full namespace is `blast::`.

```{note}
Header conversion to Doxygen comments is ongoing. `blast_manipulator.hpp` is fully
documented; other modules currently expose plain `//` comments and will improve as
headers are annotated.
```

## Manipulator

:::{doxygenclass} blast::Manipulator
:members:
:::

:::{doxygenstruct} blast::ManipulatorLimits
:members:
:::

:::{doxygenstruct} blast::ManipulatorKinematics
:members:
:::

:::{doxygenstruct} blast::ManipulatorDynamics
:members:
:::

:::{doxygenstruct} blast::ManipulatorCapsules
:members:
:::

:::{doxygenstruct} blast::Tool
:members:
:::

## World & collision primitives

:::{doxygenstruct} blast::World
:members:
:::

:::{doxygenstruct} blast::Box
:members:
:::

:::{doxygenstruct} blast::Sphere
:members:
:::

:::{doxygenstruct} blast::Capsule
:members:
:::

:::{doxygenstruct} blast::DynamicBox
:members:
:::

:::{doxygenstruct} blast::DynamicSphere
:members:
:::

:::{doxygenstruct} blast::DynamicCapsule
:members:
:::

:::{doxygenstruct} blast::DynamicDoor
:members:
:::

## Optimization

:::{doxygenstruct} blast::Optimization
:members:
:::

:::{doxygenstruct} blast::Objective
:members:
:::

:::{doxygenstruct} blast::Guess
:members:
:::

:::{doxygenstruct} blast::ConstraintSelection
:members:
:::

:::{doxygenstruct} blast::Result
:members:
:::

:::{doxygenenum} blast::OptimizationMethod
:::

:::{doxygenenum} blast::ConstraintType
:::

## Task & trajectory

:::{doxygenstruct} blast::Task
:members:
:::

:::{doxygenstruct} blast::Bspline
:members:
:::

```{note}
Free functions such as `blast::optimize` live in `blast_optimization.hpp`. Once they
carry Doxygen comments, add `{doxygenfunction}` directives here (a disambiguating
signature is required for overloaded functions).
```
