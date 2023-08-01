# Blast wishlist 2023
## Must have
- No dependancy optimization
    - SQP implementation
- Fast collision distance detection
    - 4 basic primitives
- Tested dynamic models for supported robots
    - Link6
    - Gen3 7dof
    - UR5e
- User friendly API
    - Cleanup Manipulator
        - Define required API
    - Dissociate trajectory data from functions
- More robust Error handling
- Basic example
- Documentation for external API

## Nice to have
- Hide trajectory generation details and setup
- Hide optimization details and setup
- Cleanup
- GPU
    - Fast constraints computation
    - PSO
    - Gray wolf algorithm
- Collision detection
    - Collision for general shapes
    - Example scene
    - Broad-narrow phase
    - Incremental detection for speed
- Exposed objective function and constraint API
- Add math functionality
- Easy interface with MoveIt
- Generic manipulator
- Read scene from file

## Definitely will not have but would be awesome
- Automatic scene construction from 3D scan

## Super no!
- Super slow
- Super buggy with misleading or no error messages
- Lying documentation
- Random number generator to refuse a request
- Misleading function names (compute_trajectory gives fart noise)