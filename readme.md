# Blast
A motion optimization header library that generates smooth trajectories that respect kinematic and dynamic constraints. This library does not use the GPU so it should build easily with most modern compilers.

# Build
To build a project with Blast, make sure the header files are in the path and include `blast.hpp`. All other files will be included in the correct order.

```c++
#include "blast.hpp"
```

# Examples
A set of examples are provided in the `examples` folder. 