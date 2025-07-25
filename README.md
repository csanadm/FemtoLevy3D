# Correlation function calculation with a 3D Lévy source and Coulomb FSI

 Bose-Einstein correlation function integral calculation based on a 3D, not necessarily spherical Lévy source, incorporating the Coulomb final-state interaction. Most of the integrals can be performed analytically, only very well behaving integrals remain, similarly to <a href="https://github.com/csanadm/CoulCorrLevyIntegral/">https://github.com/csanadm/CoulCorrLevyIntegral/</a>.

## Description
This package contains a calculation for quantum-statistical correlation functions, including Coulomb-correction, based on analytic results and final numerical integration. For the calculation to work, the `boost` library is needed, although the numerical integral can be programmed by the user as well. One plotting code and the fit test requires `ROOT` (latter also `Minuit2`), and several `python` libraries are required for the `python` plotters as well.

## File content

### Basics
- [**README.md**](README.md): This README file
- [**Makefile**](Makefile): Using `make all`, it will create all executables (which will have the suffix `.exe`), and requires the existence of a `deps` directory (to store the dependencies)

### Libraries
- [**Levy3D_CoulCalc.cpp**](Levy3D_CoulCalc.cpp): The main calculator class, containing the formulas and the final integral (via the Gauss-Kronrod method of `boost`)
- [**Levy3D_CoulCalc.h**](Levy3D_CoulCalc.h): Header file for the `Levy3D_CoulCalc` class
- [**HypCalculator.cpp**](HypCalculator.cpp): The functions in this class calculate the hypergeometric function 2F1.
- [**HypCalculator.h**](HypCalculator.h): Header file for the `HypCalculator` class
- [**functions.cpp**](functions.cpp): Auxiliary functions, such as the Gamma function
- [**functions.h**](functions.h): Header file for `functions.cpp`
- [**basics.h**](basics.h): Basic constants, needed for all calculations
- [**my_includes.h**](my_includes.h): A few general ROOT and other includes

### Testing of the libraries
- [**coulcorrtest.cc**](coulcorrtest.cc): An example code for testing the library
- [**coulcorrtest.py**](coulcorrtest.py): A `python` plotter for the test result
- [**coulcorrcompare.cc**](coulcorrcompare.cc): Comparing the new 3D calculation to an approximation where the Coulomb correction is calculated based on a spherically symmetric source
- [**coulcorrcompare.py**](coulcorrcompare.py): A `python` plotter for comparing the full and the approximative result, from `coulcorrcompare.cc`
- [**coulcorrcompare.betatest.py**](coulcorrcompare.betatest.py): A `python` plotter for the $\beta_{T}$-dependence of the relative error, calculated via `coulcorrcompare.cc`
- [**coulcorrcompare.2D.cc**](coulcorrcompare.2D.cc): Comparing the new 3D calculation to the approximation, a 2D version
- [**coulcorrcompare.2D.plot.cc**](coulcorrcompare.2D.plot.cc): A `ROOT` plotter for the results of `coulcorrcompare.2D.cc`


## Running the codes
- A `deps` directory should be created, this is where `make` stores the dependencies.
- A `make all` command creates all executables (they will have the `.exe` suffix), or a specific `make <codename>.exe` command produces just the given executable (here `<codename>` can be for example `coulcorrtest` or `coulcorrtestcompare`).
- One shall run `./coulcorrtest.exe > coulcorrtest.out` and then `python coulcorrtest.py` to produce a test plot (can edit `coulcorrtest.cc` to make the plot for different parameters).
- The code `coulcorrtestcompare.cc` should also be run as  `./coulcorrtestcompare.exe > coulcorrtestcompare.out`, which then can be plotted as `python coulcorrtestcompare.py`.

## Example results

### An example result
<img width="600" alt="coulcorrtest" src="https://github.com/user-attachments/assets/b1ad9c6c-1c84-4adc-8dc0-bcabbaf783ba" />

### Comparision between the full calculation and an approximation
<img width="600" alt="coulcorrcompare" src="https://github.com/user-attachments/assets/152acf79-b01a-45ac-9ba4-e40082c92e24" />

### Difference between the full calculation and an approximation
<img width="400" alt="coulcorrcompare betatest" src="https://github.com/user-attachments/assets/1611abd3-e52a-45cd-a19f-04c858229eed" />

### 2D differences
<img width="300" alt="coulcorrcompare 2D Cdiff qOqL" src="https://github.com/user-attachments/assets/8c23fb02-1834-4091-87f6-7927f03c8947" />
<img width="300" alt="coulcorrcompare 2D Cdiff qTqL" src="https://github.com/user-attachments/assets/74764955-847a-4394-b6ba-aba27b0f1e59" />
<img width="300" alt="coulcorrcompare 2D Cdiff qSqL" src="https://github.com/user-attachments/assets/86b274e0-432f-4e63-b4f1-c6879d0223a3" />
<img width="300" alt="coulcorrcompare 2D Cdiff qOqS" src="https://github.com/user-attachments/assets/90c65a73-d4f6-418e-941a-0b8eca558096" />
