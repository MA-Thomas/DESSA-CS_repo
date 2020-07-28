To run this code:

1. Save the directory in Go/src
  (e.g. in Windows   C:\Go\src\DESSA_CS__FiniteParticleRepresentation_numerical)

2. Navigate to this directory in the terminal

3. The command is (Windows or Linux):
  go run numerical_integrations_concurrency.go ExpErfc_TaylorApprox.go

Note: It is highly recommended to use this directory structure.
Golang will create Go/src during installation, and it wants to
see all code within Go/src. Go also expects to see a main package and a function
called main. Both files in the command in (3) belong to the main package.
The file numerical_integrations_concurrency.go contains the main function.
The file ExpErfc_TaylorApprox.go contains an auxiliary function used in the
computations.

If you would like to suppress outputs being written to files, comment out
from lines 151-177.

If you would like to add functionality, I suggest adding functions to new files
(also package main). In that case, to run the main function
which calls any new functions, simply add the file containing the new functions
to go run ...
