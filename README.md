# **Introduction**

This software is a back-end (JNI/C++11 based) part of the `SCOTS2SR` (<https://github.com/ivan-zapreev/SCOTS2SR>) tool for generating functional representations of `SCOTSv2.0` (<https://gitlab.lrz.de/matthias/SCOTSv0.2>) BDD controllers.

# **Dependencies**

This project is dependent on:

1. GNU Scientific Library - `GSL` (<https://www.gnu.org/software/gsl/>)
2. Tools for BDD controller determinization - `SCOTS2C` (<https://github.com/ivan-zapreev/SCOTS2C>)
3. JNI project for `SCOTS2SR ` - `SCOTS2JNI` (<https://github.com/ivan-zapreev/SCOTS2JNI>)

# **Required tools**

In order to build the project one requires to have:

1. Netbeans version 8.2 or later in its version containing: Java/JDK and C++
2. Java version 1.8 or later
3. C++ version 11 or later
3. GNU Compiler Collection (GCC/G++) 5.5.0 or later
4. XCode SDK with c++ version 4.2.1 or later (on Max OS X)

# **Build instructions**

Before the project can be build `SCOTS2C` and `SCOTS2JNI` are to be downloaded and build in the folders next the the folder containing this project. The directory structure is assumed to be as follows:

```
$ ls -al
drwxr-xr-x  10 user  staff   320 May 15 10:22 .
drwxr-xr-x   8 user  staff   256 Feb 20 08:41 ..
drwxr-xr-x  13 user  staff   416 May 14 14:30 SCOTS2C
drwxr-xr-x  13 user  staff   416 May 15 14:04 SCOTS2DLL
drwxr-xr-x   8 user  staff   256 May  7 12:12 SCOTS2JNI
```
Where `SCOTS2DLL` is storing this project. Further one needs to:

1. Build and install `GSL` version 2.4 or later from <ftp://ftp.gnu.org/gnu/gsl/>. The assumed installation folders are `/opt/local/` for Max OS X and `/user/local/` for Linux/Ubuntu
2. Build `SCOTS2C` following instructions in <https://github.com/ivan-zapreev/SCOTS2C>
3. Build `SCOTS2JNI` following instructions in <https://github.com/ivan-zapreev/SCOTS2JNI>

Further one requires to

1. Open the `SCOTS2DLL ` project in Netbeans
2. Choose proper build configuration, depending on the platform: `Release_Mac` or `Release_Linux`
3. Open the Project properties dialog
4. Go to the `Build` and then `C++ Compiler` options
5. Open the `Include Directories` dialog
6. Update the JDK `include` and `include/<platform>` directories to point to the JDK folders at your file system. This is needed for inclusion of the JNI related header files.
7. For the used platform, the step above is to be done for both `Release` and `Debug` configurations.
8. Close the project properties dialog
3. Run project `build` from the Netbeans IDE

# **Resulting binary**

Building the project results in a dynamic library generated into `./dist`. Depending on the platform this will be either `libSCOTS2DLL.dylib` (on Max OS X) or `libSCOTS2DLL.so` (on Linux/Ubuntu). The library is to be loaded from `SCOTS2SR`'s Java interface on its start-up.
