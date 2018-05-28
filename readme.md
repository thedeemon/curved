# Inside Curved Spaces
Visualizing curved 3D spaces.

We take a sphere or torus or some other curved surface and consider it a 2D space with curvature, then add one more dimension, orthogonal to those two, to get something like S²xR, a curved 3D space. Then render it by casting rays following geodesics of this space.

See:

https://www.youtube.com/watch?v=s_PNYf4qVKc

![Screenshot](https://bitbucket.org/thedeemon/curved/downloads/curved-screenshot-sphere.jpg)

Controls:

WASD to move and turn, IJKL to rotate the 3D view of the surface, O/U to zoom in/out. 

You can change R, radius of the sphere/torus/whatever (in pixels), higher values of R make everything look less curved. And you can change vision range - how far the rays should go, relative to the eye-screen distance.

"Dynamic step" is an option to make longer steps during ray walking when away from camera, it makes rendering much faster but a bit less accurate.

### More screenshots:

#### Torus
Its outer part has positive curvature similar to a sphere or ellipsoid, but the inner part looks quite different, the curvature is negative there. Worth a look before going to wormholes.

![torus](https://bitbucket.org/thedeemon/curved/downloads/curved-screenshot-torus.jpg)

#### Black hole
We take 2D section of Schwarzschild metric and find a 2D surface which, being embedded in ordinary 3D Euclidean space, has the same metric. This is pretty much what you usually see in popular clips showing spacetime curvature around black hole.

![black hole](https://bitbucket.org/thedeemon/curved/downloads/curved-screenshot-blackhole.jpg)

#### Wormhole 
A similar geometry, continuted to the other size. Two infinite 2D planes connected at one place via such wormhole.

![wormhole](https://bitbucket.org/thedeemon/curved/downloads/curved-screenshot-wormhole.jpg)

## How it's done
To describe a geometry we just write the surface equation: how 2D coordinates (u,v) are mapped into 3D coordinates (x,y,z). For instance, in case of sphere it's
```
Expr[] sphereEq() {
    auto R = new Var("R");
    return [mul(R, mul(new Cos("u"), new Cos("v"))),
            mul(R, new Sin("v")),
            mul(R, mul(new Sin("u"), new Cos("v")))  ];
}
```
The equation is given as AST,  a symbolic expression, similar to what one writes in SymPy. We have a little symbolic computation module that can do arithmetics, dot products and partial derivatives. So the main program takes this one equation in symbolic form and automatically computes basis vectors, metric tensor, its inverse, then Christoffel symbols, simplify the expressions, and in the end generate source code that is immediately used to automatically create optimized routines for computing all this stuff for different surfaces with different metrics. All these symbolic computations happen at compile time, thanks to D's compile time function execution feature. Then it's just a matter of numeric integration to walk along geodesics.

## Download and run
You can find Linux and Windows binaries in the Download section.

To run Linux version you need libsdl2 installed (not necessarily dev version).

If you want to build from source you need

Dub (package manager and build tool) and preferrably LDC2 (compiler)

https://github.com/ldc-developers/ldc/releases

(archive with LDC already includes Dub)

Build command on 64-bit Linux:

dub -b release --compiler=ldc2

Build command on Windows:

dub -b release --compiler=ldmd2 --arch=x86_64

## License 
MIT
