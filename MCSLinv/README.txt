MCSLinv: Inverse Monte Carlo algorithm 
Version 1.0
July 2017

National Institute of Science and Technology
Gaithersburg, Maryland

Written by:	
Richelle Streater, Colorado School of Mines, NIST SURF Student
Anne-Michelle Lieberson, Sherwood High School, NIST SHIP Student
Zachary Levine, National Institute of Standards and Technology
	
With help from Adam Pintar
Experimental data provided by Catherine Cooksey and Paul Lemaillet 


I. Overview

A. Program details

Language: C++ (2011 or 2014)/OpenMP, Mathematica for Post-Processing
Contact details: zlevine@nist.gov
Companion program: MCMLpar: A parallelized Monte Carlo light propagation simulation. MCSLinv includes MCMLpar, but MCMLpar is also provided alone [1, 2].

B. Summary:

MCSLinv solves for the attenuation and scattering coefficients of a single-layer material. Its development involved the implementation of a parallel version of MCML [3], a popular Monte  Carlo light propagation model. The program uses importance sampling to evaluate the parallelized MCML, or MCMLpar, over a grid of possible parameter values. It scores each combination of parameters by calculating the profile log-likelihood, using experimental data from an input file. MCSLinv finds the optimal parameters by narrowing the search region several times, and finally by fitting a paraboloid to the data. 


II. Directions for use:

A. Variables and units:

The only unit is an inverse length. While we prefer to enter the thickness in millimeters and get the scattering parameters in inverse millimeters, if the user chooses, e.g., to enter the length unit in centimeters, the scattering parameters will be output in inverse centimeters.

n- index of refraction
g- Henyey-Greenstein anisotropy parameter
t- thickness (mm)
mua- inverse interaction length for absorption
mus- inverse interaction length for scattering
mut- inverse interaction length for either absorption or scattering (total)
etaa- mua/mut, ratio of absorption events to total events

B. Input files:

1. input.txt: This file contains all user-defined parameter inputs. The first three inputs, the n, g, and t values, MUST be updated by the user to be consistent with the parameters from the experimental or numerical input data that the inverse-solving function will fit. The next six parameters could be updated if the etaa and mut values are known to be out of the current range, or if the function is returning the following error: "Error: Maximum out of bounds (from discMax.cpp)". The starting number of particles and number of search iterations should be updated if the inverse-solving function is not converging properly. The number of processors should be set to however many threads the computer can handle. The radius of the detector should match with the real physical detector radius if comparing to an experiment, or the input detector radius if comparing to numerical forward data. To neglect the finite radius effect, set the radius input to any large number, such as 10000000. Note: MCML neglects the finite radius effect, so if MCML is providing the forward data, set the radius input to a large number [3].

2. exp.txt: This file contains experimental ARS curves. 

i.   Experimental input: Input intensity should be normalized to a measurement with no sample present. Input should be averaged over polarizations and each Angle-Resolved Scattering (ARS) value should be adjusted for solid angle by a factor of A/R^2, where A is the detector aperture area and R is the detector radius. Data should be in central-angle form, in which the detector is centered around round angle numbers like 0 degrees, 5 degrees, 10 degrees, etc., and each ARS value is collected around that angle. For instance, the angle should range from -2.5 to 2.5 degrees, then from 2.5 to 7.5 degrees, etc. See (ii) if the data is not in this form.

ii.   Numerical input: The output data from MCMLpar will be in the correct form to input to MCSLinv [2]. MCMLpar was written by the same team as MCSLinv, and MCMLpar can be downloaded from the same location as MCSLinv. The output data from Wang and Jacques’ MCML will not be in central-angle form [3]. To convert to central-angle form, collect ARS data for twice as many points as desired. For instance, if a 36-point ARS curve is desired (5 degree angle ranges), collect 72 points of ARS data, so that data is stored between 0 and 2.5 degrees, then 2.5 and 5 degrees, etc. Double the ARS count between 0 and 2.5, then add together consecutive counts in pairs to form the rest of the list. For instance, add together the contents from 2.5 to 5 degrees and from 5 to 7.5 degrees and store the sum as the second element. Then add together the contents from 7.5 to 10 degrees and from 10 to 12.5 degrees and store this sum as the third element. This should result in a total of 36 list elements. The fixARS.cpp function in the source code also will do this.

3. expExample.txt and inputExample.txt: These are two example inputs. With these inputs, the program should run for only 1-5 minutes and should output a mus value around 2.75 mm and a mua value around 0.0085 mm.

C. Output files:

1. MCSLoutput.csv: This file will include the time elapsed, the maximum likelihood estimate for mua and mus, and the paraboloid parameters for Mathematica to input.

2. MCSLoutput.nb: This Mathematica notebook will use MCSLoutput.csv. The notebook will plot the maximum likelihood estimate for mua and mus as a point, with an elliptical confidence interval (confidence level can be changed by the user).

3. MCSLoutputExample.csv: This is an example output file. The program should output something similar to this, given the input files inputExample.txt and expExample.txt.

D. External packages:

1. SPRNG 2.0b: This is the best version of SPRNG to download for Windows, since a Windows package is available from NAADSM [4]. SPRNG 2.0b is also available for Linux and Mac [5]. SPRNG 2.0b is the default; no additional flags are necessary. libsprng.a or some other SPRNG library must be linked.

2. SPRNG 5.0: Linux and Mac users may wish to download the latest version of SPRNG. The flag “-D SPRNGFIVE” must be used if the user has SPRNG 5.0 instead of SPRNG 2.0b.

3. OpenMP: The code uses OpenMP for parallelization. OpenMP can generally be downloaded with MinGW. The OpenMP library must be linked and the flag -fopenmp must be included.

4. Eigen: Eigen is used in fitting the paraboloid to the data [6]. The search directory may need to be updated to include the Eigen folder.

5. MinGW libraries: The following libraries must be in the same directory as the .exe for it to work: libgomp-1.dll, libstdc++-6.dll, and libgcc_s_dw2-1.dll. These are provided with the .exe, and can be found at: http://www.mingw.org/


III. Program contents

A. Directories:

1. analyticTests: Contains Mathematica notebook with analytical single-scatter tests.
2. bin: Contains executable program
3. dataIn: Contains user input files
4. dataOut: Contains output files and a Mathematica notebook to graph results.
5. obj: Contains .o files

IV. Design

A. Design choices:

Style: Variables and functions use lowerCamelCase. Constants use UPPER_CASE. Class definitions use UpperCamelCase.

main.cpp: Coordinates etaa and mut were introduced, where mut = mua + mus (mua is a layer’s attenuation coefficient and mus is a layer’s scattering coefficient), and etaa = mua/mut. The advantage here is that etaa samples from a Bernoulli distribution, which requires less computation than a sample from an exponential distribution. Also, these two distributions are independent, which means that instead of constantly weighting on a grid of possible parameters with Importance Sampling, the program only weights two vectors of possible parameters until the particle escapes the medium.

main.cpp: Particle, layer, and weight classes were introduced. The layer class is convenient because it stores all reference values for the layer that a particle is in. Although the inverse part of the program only currently works for single-layer material, MCMLpar works for multi-layer material, so the program creates a vector of layers at the beginning. To introduce multi-layers, changes to the forward program would not be necessary, but some changes to the weight class and the updateInterval function would. The particle class includes a layer as a data member because the particle constantly needs information about the layer that it exists in. This might create unboxing slow-down. The weight class is convenient because it groups all weight elements together: a scalar weight, to account for the probability of the particle attenuating in the material, and mut and etaa weight vectors, to provide importance sampling weights for these two parameters. The particle class also includes a weight as a member. It would be possible to unroll all weight data members and functions into the particle class.

main.cpp: Vectors were used for almost all sets of values. This design was used because vectors can be resized more easily than arrays, but it is possible to change some vectors in the program into arrays.

main.cpp: The median value of the etaa and mut vectors is the reference point for Importance Sampling. It is possible to speed up the program by using a reference value that runs faster, but this could result in lower precision in the ARS curve.

main.cpp: Some functions were initialized in main, within the parallel for loop, to avoid race conditions. There may be other ways of accomplishing this.

contour.cpp: Instead of simply returning the likelihood function and allowing the user to find the contours and the maximum likelihood estimates, the program does it automatically. This allows the program to return a very small, readable output file and lets the user quickly receive an answer and evaluate the Mathematica notebook. However, it requires the user to download Eigen.

detect.cpp: Instead of just collecting output angles, this function measures where the photon would intersect with a detector given a finite detector radius. To turn off this feature, the user only needs to set the detector radius to be very large.

findRegion.cpp: The program currently updates search region by comparing each log-likelihood value to the 99.9999% confidence level chi-square value. This is the current design because it is the simplest, but other test designs included fitting a paraboloid to the likelihood surface and solving for the maximum and minimum values for the elliptical contour at the same chi square value, as well as fitting a paraboloid to the likelihood surface and updating the search region by comparing the paraboloid values to the chi-square value. These methods were all functional and yielded the same results.

findRegion.cpp: The shrinking of the confidence region was limited by a factor of two. In other words, if the likelihood surface had a contour with a bounding etaa or mut region that was less than half of the previous one, the program would increase these bounds evenly about the discrete maximum. This was implemented after the program started “losing” the real etaa and mut values by choosing a search region too small to include them. After the 2x limit was implemented, the function found the true parameters every time, given a valid experimental input.

fixARS.cpp: In most experiments, researchers set the detector at evenly spaced central angles (for example, 0, 5, 10, …). The program corrects for this by collecting data for twice as many angle slices and adding the slices up in the proper way to imitate central angle data. The program also switches the first and last elements in the ARS vector, the second and second-to last elements, etc. This resulted in the transmission angles ranging from 90-180, while the reflection angles ranged from 0-90. The only purpose of this was to create more intuitive ARS curves, based on the idea that most humans read from left to right and are most comfortable with light travelling from left to right.

propagate.cpp: When implementing Importance Sampling, the program must update weights every time it samples from etaa and mut distributions. It turns out that the mut PDF is a piecewise function, since the particle can hit boundaries. For every distance value sampled that would not move the particle beyond the boundary, the mut PDF is the exponential function of a typical exponential distribution. For every distance sampled that would move the particle beyond the boundary, the mut PDF is a single point that holds the value of the exponential CDF, evaluated from the boundary distance to infinity. Because of this, a particle must update its vector of mut weights with a different function, depending on whether it hits a boundary or not.

roulette.cpp: This function uses the same design and values as the roulette function in MCML, including where and when it is called.

scattFunction.cpp: Here, the program takes advantage of the fact that it gains a degree of freedom from the problem’s phi-independence. The math operations in this function are equivalent to building an orthogonal triad, where the particle’s direction vector is the first column in the triad, and multiplying the triad by the Euler rotation matrix. The phi independence allows for freedom in choosing one of three sets of math operations. By finding the index of the smallest element in the particle’s direction vector and using this index to decide which of the three sets of math operations to use, the program avoids roundoff error.

searchRegion.cpp: The choice of using the chi-square value at 99.9999% confidence was based on a desire to narrow down the search region while also preserving enough of the likelihood surface to fit a paraboloid that yields a good confidence region contour. The program’s functionality did not change very much, however, when this value was varied.

specularR.cpp: The presence of this function stays consistent with the MCML convention. While it does not make a big difference to the forward Monte Carlo simulation, if the boundary function was called first instead of specular, it would be harder to compare to MCML and some strange behavior might pop up, such as the inverse search region converging at a different number of particles for different indices of refraction.

updateInterval.cpp: updateInterval and several other functions return Booleans. This allows the entire program to stop when an error occurs, such as mismatched vector sizes or a lack of a discrete maximum within the edges of the search region. 


V. Troubleshooting:

A. Possible program errors: 

1. "Error: saddle point (from checkEigenVals.cpp)" or "Error: minimum (from checkEigenVals.cpp)": This error means that the paraboloidal fit to the likelihood surface does not have a maximum. If everything else looks good up to this point, there is probably something wrong with Eigen.

2. "Error: Maximum out of bounds (from discMax.cpp)": This error comes up most often when the inverse-solving function cannot achieve the input ARS curve, no matter what parameters it uses. The most likely cause is a mismatch between input parameters and the true parameters of the experimental or numerical input data. Check to make sure that the t, n, and g values are correct and are all in the correct units. Also, make sure that the input data is in central-angle form (see section I.A.2 of this file for more details). This error could also be the result of a bad starting search region: for instance, the true mua is 0.05 and the starting mua range is between 0.001 and 0.01. Try expanding the starting search region. Or, the error could be caused by the program “losing” the true mua, mus values. Try increasing the starting number of particles.

3. "Error: File did not open (in setParameters.cpp)." or "File reading error (in setParameters.cpp)." or "Error: File did not open (in fileToVec.cpp).": An input file is probably incorrectly named or missing, or for some reason the program cannot find it. There could also be a problem with the format of the input file.

4. "Error: x actual out of y wanted streams spawned!": This means that the parallelization did not work. It is possible that something went wrong with SPRNG, or the user must change numProc in the input file.

5. "Error: could not decrease search region (from findRegion.cpp) ": The program could not narrow the search region down. This could be due to a starting search region that does not include the correct value: for instance, the true mua is 0.05 and the starting mua range is between 0.001 and 0.01.

6. Compiler not recognizing nullptr or initialization lists: User must use a compiler that is compatible with C++ 2011 or C++ 2014.


VI. Works Cited:

[1]    Z.H. Levine, R.H. Streater, A.-M. R. Lieberson, A.L. Pintar, C.C. Cooksey, and P. Lemaillet, “Algorithm for Rapid Determination of Optical Scattering Parameters”, unpublished.

[2]    R.H Streater, A-M. R. Lieberson, A.L. Pintar, Z.H. Levine, “MCMLpar and MCSLinv: Parallel MCML and an Inverse Monte Carlo Algorithm to Calculate Optical Scattering Parameters”, unpublished

[3]    L. Wang, S.L. Jacques, L. Zheng, “MCML-Monte Carlo modeling of light transport in multi layered tissues”, Computer Methods and Programs in Biomedicine 47, 1995, pp. 131-146.

[4]    NAADSM Development Team, “SPRNG for Microsoft Windows and NAADSM”, http://www.naadsm.org/opensource/sprng, 2013, accessed 28 July 2017.

[5]    M. Mascagni, A. Srinivasan, “Algorithm 806: SPRNG: a scalable library for pseudorandom number generation”, ACM Transactions on Mathematica Software 26, 2000, pp. 436-461; http://sprng.org, accessed 28 July 2017.

[6]    G. Guennebaud, B. Jacob, “Eigen v3”, http://eigen.tuxfamily.org, 2010, accessed 28 July 2017.
