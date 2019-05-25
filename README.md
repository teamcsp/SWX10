# SWX10

:skull:  
X10 Implementation of the [Smith-Waterman Algorithm](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)

Contains Sequential and Parallel implementation
- For parallel implementation, refer to this [paper](https://github.com/khongyew/SWX10/blob/master/paper/Parallel%20Implementation%20of%20the%20Smith-Waterman%20Algorithm.pdf)
 
Requires X10 to compile and execute:
- Read http://x10-lang.org/releases/x10-release-262.html
- X10 can be compiled to Java or native (C++) executable

Possible to use X10 with Eclipse IDE:
- Install x10 extensions onto plain Eclipse using **update site**
- See more at http://x10-lang.org/documentation/x10dt-installation.html (X10DT Installation via Eclipse Update Manager)

## Compile and run
Quick notes  
Java backend:
- compile using `x10c Match.x10`
- run using `x10 Match.x10`

C++ backend:
- compile using `x10c++ -o Match Match.x10`
- run using `./Match`
