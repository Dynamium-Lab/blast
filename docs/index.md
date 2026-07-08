# Blast

> 1. Make it work
> 2. Make it faster
> 3. Repeat step 2

**Blast** is a fast C++ trajectory optimization library for robot manipulators. Given a
robot model and a motion task (start pose, goal pose, and constraints), Blast finds a
smooth, time-optimal, collision-free B-spline trajectory that satisfies the following limits

- joint position
- joint velocity
- joint acceleration
- joint torque
- tool speed

It is **header-only** from the consumer's perspective (`#include <blast>`) and ships with
first-class **Python bindings** (`import blast`).

Blast is the software accompanying the ICRA 2026 paper *"Accelerating Trajectory
Optimization by Exploiting B-Spline Gradient Structure."* See [Citing Blast](about/citing.md).

```{toctree}
:caption: Getting started
:maxdepth: 2

getting-started/installation
getting-started/quickstart-cpp
getting-started/quickstart-python
```

```{toctree}
:caption: Guides
:maxdepth: 2

guides/concepts
guides/constraints
guides/objectives-and-guesses
```

```{toctree}
:caption: API reference
:maxdepth: 2

api/cpp/index
api/python/index
```

```{toctree}
:caption: About
:maxdepth: 1

about/citing
about/contributing
about/changelog
```
