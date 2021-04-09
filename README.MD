# SimAnMo: Simulated Annealing Modeler

## Introduction
SimAnMo (SIMulated Anealing Modeler) is a C++-software to generate runtime performance models for programs. It is based on principles of [Extra-P](https://www.scalasca.org/software/extra-p/download.html) and was presented in [BUR18]. Besides models for programs with polynomial runtime behavior, meaning that the runtime increases polynomially with the size of an input parameter *p*, SiMaNo is also able to model exponential runtime behavior where the runtime grows exponentially in *p*. The models are exported as LaTeX documents.

The model generation process is based on an SMP-parallelized simulated annealing procedure. Training data is required in files of the following format:

    (	p_1	;	val_1	)
    (	p_2	;	val_2	)
    (	p_3	;	val_3	)
    ...
    (	p_l	;	val_l	)

The first value p_i in the brackets is the value of the input parameter p for the respective mesaurement. val_i gives the runtime in our case.
You can also add additional *k* measurements for evaluation. Then, you have to add them in the following way after the *l* training point with *ADDPOINTS* as mandatory indicator for the software:

    #ADDPOINTS
    (	p_(l+1)	;	val_(l+1)	)
    ....
    (	p_(l+k)	;	val_(l+k)	)

Comments can be added by using a *#*, e.g.,

    (	p_1	;	val_1	)
    #(	p_2	;	val_2	)
    (	p_3	;	val_3	)

which will ignore the second point for the model generation.
Example files are given in the *inputs* directory.

## Compilation
### Prerequistes
* Microsoft Visual Studio 2017 (or higher) or gcc/g++ 6.0 (or higher, -std=c++-17 must be supported)
* [Eigen Library](http://eigen.tuxfamily.org/index.php?title=Main_Page)
* Optional: [NAG Library](https://www.nag.com/content/nag-library)

### Windows
SiMaNo has successfully been tested with Visual Studio 2017/2019, Eigen 3.3 and NAG C-library 24/25/26.

Open the Visual Studio and build the Release (x64) configuration. Debug (x64) is also supported. The path to Eigen must be integrated into Windows' %INCLUDE% system variable or has to be added to the addional include directories in the VS project settings.

The Release (x64) binaries are generated in the main project directory in the subfolder *bin*. Debug configuration is built in the directory set in the VS project properties.

If you prefer to employ the NAG library instead of Eigen you have to:
* Add the path to the NAG headers to %INCLUDE% or in the VS project settings
* Add the path to the NAG libraries to %LIB% or in the VS project settings
* Add *nagc_nag_MT.lib* and *legacy_stdio_definitions.lib* as additional dependencies to the linker in the VS project settings
* Define the preprocessor macro `USE_NAG` in the VS project settings

### Linux
SiMaNo has successfully been tested with GCC  8.3, 8.4, 9.0 and 10.1 with Eigen 3.3.

Make sure that the Eigen-headers are included in $CPATH environment variable or addpend the path to the `INC` variable in the `Makefile`. Then type `make` while being in the Makefile directory. The binary is generated in the main project directory in the subfolder *bin*.


## Usage
First, we describe all calling parameter and second, the code changes regquired to change the model type and the cost metric which is minimized related to the training data.

### Execution and Parameter
The application is highly configurable. The main usage (depending on the OS) is:

`(./) SiMaNo(.exe) [options]  --inputfile PATH_AND_NAME_OF_INPUTFILE --gl --outpath PATH_WHERE_TO_PLACE_OUTPUT --texfile NAME_OF_TEX_AND_PDF_FILES`

There are many command line parameter which influence the behavior. These are:
#### General and Annealing Process
* `--help` / -h Print the help
* `--inputfile / -i + PATH_TO_INPUT_FILE/FILNAME` MANDATORY. The path and name to/of the input file
* `--outfile / -o + PATH_TO_OUT_FILE` The path to the output files (default=folder of executable)
* `--texfile / -t + NAME_OF_FILES`  MANDATORY. The name of output tex and pdf files
* `--number_of_threads / --nt + INT` How many threads anneal in parallel (default=1)
* `--number_of_trials / --tr + INT`  How many repetitions of the annealing process (default=1)
* `--ann_target_temp / --att + FLOAT` Target temperature for annealing (default=1e-9)
* `--ann_cooling_rate / --acr + FLOAT` Temperature degredation per iteration (default=0.99)
* `--ann_steps / --as + INT` How many steps are performed per temperature (default=15)
* `--ann_steps_wo_mod / --awm + INT` Heuristic: Stop if no improvement after this number of steps (default=200000)
* `--ann_steps_backtrack / --abt + INT`Heuristic: Backtrack per thread after this number of steps (default=20000)

#### Model and Cost Calculation Configuration
* `--costcalc_type / --cct + STRING` Decides whether to use RSS (STRING="rsscost") or raRSD (STRING="rarsdcost") as cost metric to minimize (default="rsscost")

#### LaTeX Printing
* `--genlatex / --gl` Activate the LaTeX report generation (default=false)
* `--openpdf / --op` Generated pdf file is automatically opened with pdfxchange at --pathtopdfxchange (default=false)
* `--pathtopdfxchange` If --openpdf is set, pdf file is automatically opened with pdfxchange at this path
* `--logy` The y-axis in the prediction graph is scaled logarithmically
* `--print_confidence / --pc` Print the confidence interval in the predictiion (default=false)
* `--confidence_interval / --ci + FLOAT` Set size of confidence interval when printing it (default=0.0)
* `--print_cost_details / --pcd` Print details of cost development during annealing (default=false)

#### Parameter Depending on Model Type Employed
##### Polynomial-logarithmic (pol-log) Model Configuration [Extra-P Standard Model]
*`--max_log_range / --melog + FLOAT` Maximum exponent for lorarithms (default=4.00)
*`--max_pol_range / --mepol + FLOAT` Maximum exponent for polynoms (default=6.00)

##### Linear-logarithmic (lin-log) Model Configuration [Linear Regression]
*`--create_lin_log / --ll` Create a lin-log model if set (default=false)
*`--base_lin_log / --bll + INT` Basis that is used for lin-log-model (default=2)

##### Exponential (exp) Model Configuration [Extended Extra-P Model]
* `--max_exp_range / --meexp + INT` Maximum coefficient in exponent (default=4.00)

##### Exponential-polynomial (exp-pol) Model Configuration [Introduced in [BUR20]]
* `--exp_pol_min_coeff / --epmic + FLOAT` Minimum coefficient in the exponent of exp-pol models (default=0.01)
* `-exp_pol_max_coeff / --epmac + FLOAT` Maximum coefficient in the exponent of exp-pol models (default=2.0)
* `--exp_pol_min_exp / --epmie + FLOAT` Minimum exponent in the exponent of exp-pol models (default=0.5)
* `--exp_pol_max_exp / --epmae + FLOAT` Maximum exponent in the exponent of exp-pol models (default=3.0)

##### Factorial (fac) Model Configuration [Introduced in [BUR21]]
* `--fac_pol_min_coeff / --fpmic + FLOAT` Minimum coefficient in the exponent of the polynomial part of the fac models (default=-0.5)
* `--fac_pol_max_coeff / --fpmac + FLOAT` Maximum coefficient in the exponent of the polynomial part of the fac models (default=3.0)
* `--fac_log_min_coeff  / --flmic + FLOAT` Minimum coefficient in the exponent of the logarithmic part of the fac models (default=1e-3)
* `--fac_log_max_coeff  / --flmac + FLOAT` Maximum coefficient in the exponent of the logarithmic part of the fac models (default=1.5)

### Different Models and Metrics
#### Changing the Models
SimAnMo supports the four types presented in the parameter section. If you want to create a model of that shape for your training data you have to edit the *SimulatedAnnealingExtraP.cpp*. There, go to the end of the `main`-function, which is the end of the file, and change the template parameter for the call of `annealingManager`. You can use:
* `ExtraPSolution`: Extra-P model as employed in [BUR18]
* `ExponentialSolution`: Extended Extra-P model
* `ExponentialPolynomSolution`: The exponential-polynomial model presented in [BUR20]

Lin-log models can be created automatically by using the `--create_lin_log / --ll` calling parameter.

#### Canging the Metric
SimAnMo supports two different metrics in its current version *RSS* and the relative average Residual Sum of Differences *raRSD* as defined in [Bur20]. To exchange the metric employed, you have to edit the function `annealingManager` in *SimulatedAnnealingExtraP.cpp*. There, search the call of `doAnnealing` an change the second template parameter. Use `nnrRSSCostCalculator` for *raRSD* and *SSCostCalculator* for *RSS*.


## Reproducing Paper Results
To reproduce the results presented in [BUR20], execute the respective script (*bat* for Windows, *sh* for Linux) in the *script* folder. It will generate the pdf-files in *outputs* directory (it overwrites the example files shipped with the code). These scripts can also be employed for a general test of functionality.

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
       Developing Models for the Runtime of Programs
With Exponential Runtime Behavior.</td>
</tr>
<tr>
    <td valign="top"><a name="ref-MACH20"></a>[BUR21]</td>
    <td>Burger, Michael ; Bischof, Christian ; Nguyen, Giang Nam :
       SimAnMo: A Parallelized Runtime Model Generator.</td>
</tr>
</table>

## Contact
michael.burger@tu-darmstadt.de
giang_nam.nguyen@tu-darmstadt.de
