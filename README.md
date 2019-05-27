# SWX10

:skull:  

X10 Implementation of the [Smith-Waterman Algorithm](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)

Source code contains Sequential and Parallel implementation
- For information on the parallel algorithm, refer to this [paper](https://github.com/khongyew/SWX10/blob/master/paper/Parallel%20Implementation%20of%20the%20Smith-Waterman%20Algorithm.pdf)
 
Requires X10 to compile and execute:
- Read http://x10-lang.org/releases/x10-release-262.html
- X10 can be compiled to Java or native (C++) executable

Possible to use X10 with Eclipse IDE:
- Install x10 extensions onto plain Eclipse using **update site**
- See more at http://x10-lang.org/documentation/x10dt-installation.html (X10DT Installation via Eclipse Update Manager)

## Prerequisites
1. use a linux environment (e.g. Ubuntu 18.04)
2. download the X10 binaries from http://x10-lang.org/releases/x10-release-262.html
3. extract the contents and put them into `/opt/x10`
4. put the following lines at the end of `~/.bashrc`
```
# Settings for X10
export PATH=$PATH:/opt/x10/bin
export LD_LIBRARY_PATH=/opt/x10/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/x10/stdlib/lib:$LD_LIBRARY_PATH
export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
```
5. run `source ~/.bashrc` or reboot the pc

## Compile and run
Quick notes  

Java backend:
- compile using `x10c Match.x10`
- run using `x10 Match.x10`

C++ backend:
- compile using `x10c++ -o Match Match.x10`
- run using `./Match`
