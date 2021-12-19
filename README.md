Purpose: This is an example repository of how to use CMake, C++, and google tests to get code coverage of your code base. 

Included: 
-As an example, the source code included uses algorithms being taught in Computation Geometry at Washington University in St. Louis in the Fall 2021 semester.
-CMake is used to make static and dynamic libraries as well as executables and unit tests.
-Github workflow is used for continuous integration and continuous deployment (CI-CD) (building, testing, coverage analysis, delivering)
-Code coverage is automatically calculated in the Ubuntu Debug job.

Output (On every push)
-The source code will be built for {Windows, MacOs, Ubuntu} for {Release, Debug}
-Unit tests (google tests) will be run in each configuration.
  -Test result summary will be uploaded as an artifact for Ubuntu Release. 
  -This repository does not use OS specific items, so only one test results was uploaded, but this could be tinkered with if need be.
-The executables, tests, libraies and associated library headers will be installed with cmake
  -They are uploaded as an artifact for each Release configuration.
  -Debug configuration has extra instrumentation linked in that I did not necessarily want to upload, but this could be tinkered with.
-Code coverage will be computed only on the Ubuntu Debug configuration. 
  -Debug configurations are linked with extra instrumentation that is not in the release installation.
  -lcov and genhtml is used to generate the code coverage
  -I have no plans to do code coverage for other OS's and compilers as this project does not have any lines of code that are OS or compiler dependent. 

