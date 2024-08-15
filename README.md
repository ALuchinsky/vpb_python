# Python realization for Vectorized Persistence Block vectorization

The goal is to implement TDAvec library in python using R as a source (see GitHub repository https://github.com/alexey-luchinsky/TDAvec)

## Other python libraries:

* skit-tda: https://docs.scikit-tda.org/en/latest/libraries.html
* Giotto-TDA: https://giotto-ai.github.io/gtda-docs/0.5.1/index.html
* Geometry Understanding in Higher Dimensions: https://github.com/GUDHI, https://gudhi.inria.fr/introduction/
* Dionysus 2: https://www.mrzv.org/software/dionysus2/




## TODO

* Remove outliers
* Rewrite paper
* Send to Umar

Done tasks

* Create remote repository: done Aug 10, 2024, 09:55
* Start the JOSS paper template: Aug 10, 2024, 10:05
* Find other python libraries for vectorization: Aug 10, 2024, 16:39
* Incorporate my into ML workflow
* Redo simulations: done Aug 15, 2024, 07:38
* Correct table: done Aug 15, 2024, 07:38
* Add dim0+dim1 set of predictors: done Aug 15, 2024, 08:50


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
    * `TDAvec.pyx`: Cython version of the package
    * `TDAvec_pyx_test.ipynb`: Jupyter notebook, comparison with R results
    * `setup.py`: Cython make file. See TDAvec_pyx_test.ipynb for instructions
* R/: R files live here
    * `TDA_vec.Rmd`: extracting test results for different functions