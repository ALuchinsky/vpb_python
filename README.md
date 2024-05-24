Python realization for Vectorized Persistence Block vectorization

The goal is to implement TDAvec library in python using R as a source. The following functions needs to be translated:
* `computeVPB`:     A Vector Summary of the Persistence Block
* `computeECC`:     A Vector Summary of the Euler Characteristic Curve
* `computeNL`:      A Vector Summary of the Normalized Life Curve
* `computePES`:     A Vector Summary of the Persistent Entropy Summary Function
* `computePI`:      A Vector Summary of the Persistence Surface
* `computePL`:      A Vector Summary of the Persistence Landscape Function
* `computePS`:      A Vector Summary of the Persistence Silhouette Function
* `computeVAB`:     A Vector Summary of the Betti Curve

For each of the functions I should:
* translate the code from R to python in the python notebook
* create some results in R
* Make sure that python results agree with R
* translate the code into pyx format
* make sure that results agree with R


File structure:
* R/: all R files live here
  * `TDA_vec.Rmd`: creates files to compare with using TDAvec R library
