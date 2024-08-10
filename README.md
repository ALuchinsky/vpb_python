# Python realization for Vectorized Persistence Block vectorization



The goal is to implement TDAvec library in python using R as a source (see GitHub repository https://github.com/alexey-luchinsky/TDAvec)

## TODO

* Create remote repository
* Find other python libraries for vectorization
* Incorporate my into ML workflow
* Start the JOSS paper

The following functions were translated from R and checked against R results:
* `computeVPB`:     A Vector Summary of the Persistence Block
* `computePL`:      A Vector Summary of the Persistence Landscape Function
* `computePS`:      A Vector Summary of the Persistence Silhouette
* `computeNL`:      A Vector Summary of the Normalized Life Curve
* `computeVAB`:     A Vector Summary of the Betti Curve
* `computeECC`:     A Vector Summary of the Euler Characteristic Curve
* `computePES`:     A Vector Summary of the Persistent Entropy Summary Function
* `computePI`:      A Vector Summary of the Persistence Surface



## File structure:
* python/: python files live here
    * TDAvec.pyx: Cython version of the package
    * TDAvec_pyx_test.ipynb: Jupyter notebook, comparison with R results
    * setup.py: Cython make file. See TDAvec_pyx_test.ipynb for instructions
* R/: R files live here
    * TDA_vec.Rmd: extracting test results for different functions