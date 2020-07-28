To run this code:

1. Save the directory in Go/src
  (e.g. C:\Go\src\DESSA_CS__PointParticleRepresentation_numerical)

2. Navigate to this directory in the terminal

3. The command is:
  go run numerical_integrations_concurrency.go

Note: It is highly recommended to use this directory structure.
Golang will create Go/src during installation, and it wants to
see all code within Go/src. Go also expects to see a main package and a function
called main. The file in the command in (3) belongs to the main package.
The file numerical_integrations_concurrency.go contains the main function.


If you would like to suppress outputs being written to files, comment out
from lines 125-152.

If you would like to add functionality, I suggest adding functions to new files
(also package main). In that case, to run the main function 
which calls any new functions, simply add the file containing the new functions
to go run ...
