# SiMaNo: Simulated Annealing Modeler

## Introduction
SiMaNo (SIMulated Anealing Modeler) is a C++-software to generate runtime performance models for programs. It is based on principles of [Extra-P](https://www.scalasca.org/software/extra-p/download.html) and was presented in [BUR18]. Besides models for programs with polynomial runtime behavior, meaning that the runtime increases polynomially with the size of an input parameter *p*, SiMaNo is also able to model exponential runtime behavior where the runtime grows exponentially in *p*.

The model generation process is based on an SMP-parallelized simulated annealing procedure.

## Compilation
#### Prerequistes
* Microsoft Visual Studio 2017 (or higher) or gcc/g++ > 5.5 (-std=c++-17 supported)
* [Eigen Library](http://eigen.tuxfamily.org/index.php?title=Main_Page)
* Optional: [NAG Library](https://www.nag.com/content/nag-library)
*
### Windows
SiMaNo has successfully been tested with Visual Studio 2017/2019, Eigen 3.7 and NAG C-library 24/25/26.

Open the Visual Studio and Build the Release (x64) configuration. Debug (x64) is also supported. The path to Eigen must be integrated into Windows' $(INCLUDE) system variable or has to be added to the addional include directorys setting in the VS project properties.

If you prefer to employ the NAG library instead of Eigen you have to:
* Add the path to the NAG headers to $(INCLUDE) or in the VS project settings
* Add the path to the NAG libraries to $(LIB) or in the VS project settings
* Add *nagc_nag_MT.lib* and *legacy_stdio_definitions.lib* as additional dependencies to the linker in the VS project settings
* Define the preprocessor macro `USE_NAG`in the VS project settings

### Linux
SiMaNo has successfully been tested with GCC 5.5, 6.3, 7.1, 8.3, 9.0 and 10.1 with Eigen 3.7.

Make sure that the Eigen-headers are included in $CPATH environment variable or addpend the pass to the `INC` variable in the `Makefile`. Then type `make` while being in the Makefile directory.


## Usage

## References
<table style="border:0px">
<tr>
    <td valign="top"><a name="ref-MACH20"></a>[BUR18]</td>
    <td>Burger, Michael ; Bischof, Christian ; Calotoiu, Alexandru ; Wunderer, Thomas ; Wolf, Felix:
       <a href=https://apps.fz-juelich.de/jsc-pubsystem/pub-webpages/general/get_attach.php?pubid=1861>Exploring the Performance Envelope of the LLL Algorithm.</a></td>
</tr>
<tr>
    <td valign="top"><a name="ref-MACH20"></a>[BUR20]</td>
    <td>Burger, Michael ; Bischof, Christian ; Nguyen, Giang Nam :
       <a href=https://apps.fz-juelich.de/jsc-pubsystem/pub-webpages/general/get_attach.php?pubid=1861>Developing Models for the Runtime of Programs
With Exponential Runtime Behavior.</a></td>
</tr>
</table>
