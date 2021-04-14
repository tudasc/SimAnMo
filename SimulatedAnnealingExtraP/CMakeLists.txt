
add_library(SimAnMoLib AlglibLinearRegressionFinder.cpp
CostCalculator.cpp
EigenParameterEstimator.cpp
ExponentialPolynomSolution.cpp
ExponentialSolution.cpp
ExtraPSolution.cpp
FactorialSolution.cpp
GeneralParameterEstimator.cpp
LatexPrinter.cpp
LibraryInterface.cpp
LinearRegressionFinder.cpp
LinearSolution.cpp
MeasurementDB.cpp
MeasurementDBReader.cpp
nnrRSSCostCalculator.cpp
ParameterEstimator.cpp
QualityLogger.cpp
raRSDParameterEstimator.cpp
RSSCostCalculator.cpp
SimAnMoTools.cpp
SimulatedAnnealingExtraP.cpp
Solution.cpp
SolutionModifier.cpp
StartSolutionFinder.cpp
TemperatureInitializer.cpp)

add_subdirectory(alglib)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(SimAnMoLib PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(SimAnMoLib PUBLIC 
	AlgLib	
)

target_include_directories(SimAnMoLib PUBLIC
   "${PROJECT_BINARY_DIR}"
	alglib
#                           ${PROJECT_SOURCE_DIR}/SimulatedAnnealingExtraP
                           )
