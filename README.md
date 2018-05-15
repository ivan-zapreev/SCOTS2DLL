1. Make sure the Java JNI includes are correct
2. Make sure the CUDD installation inlude is correct
3. The Projects are located in the same folder so relative paths work
4. GSL (https://www.gnu.org/software/gsl/) >= v2.4 (ftp://ftp.gnu.org/gnu/gsl/) is installed
5. Choose proper confogiration in to run Mac/Linux
6. Make sure the java paths are correct in the project properties.
7. Do not need to resolve the project issues if the proper configuration is chosen (?)
8. Static linking GSL on linux, somehow the dynamic version fails to load when JNI library is loaded, however things go well with dynamic version of CUDD (?????? why ?????)