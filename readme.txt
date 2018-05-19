To run Linux version you need libsdl2 installed (not necessarily dev version).

To build from source you need
Dub (package manager and build tool) and preferrably LDC2 (compiler)
https://github.com/ldc-developers/ldc/releases
(archive with LDC already includes Dub)

Build command on 64-bit Linux:
dub -b release --compiler=ldc2

Build command on Windows:
dub -b release --compiler=ldmd2 --arch=x86_64



