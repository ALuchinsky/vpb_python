# Python realization for Vectorized Persistence Block vectorization

The goal is to implement TDAvec library in python using R as a source (see GitHub repository https://github.com/alexey-luchinsky/TDAvec)

## How to run:

1. create python enviroment and install all packages from requirements.txt
2. run `python setup.py build_ext --inplace` from the `python/` directory to compile pyx file
3. launch `R/TDA_vec.Rmd` to get R results
4. run through `test_TDAvectorizer.ipynb` notebook to see python results and save figures
5. launch docker
6. Run `source make_pdf.sh` from `paper` directory to create pdf

## File structure:
* python/: python files live here
    * `TDAvec.pyx`: Cython version of the package. 
    * `TDAvectorizer.py`: TDAvectorizer class
    * `test_TDAvectorizer.ipynb`: Jupyter notebook, testing TDAvectorizer
    * `setup.py`: Cython make file. See TDAvec_pyx_test.ipynb for instructions
    * `TDAvec_pyx_test.ipynb`: Jupyter notebook, comparison with R results
* R/: R files live here
    * `TDA_vec.Rmd`: extracting test results for different functions
* paper/: paper files live here

## TODO

* FDA in R
* Quick FDA in Python
* FDA: compare to R
* add FDA
* reread/rewrite paper

### Done tasks

* DONE: Aug 10, 2024, 09:55: Create remote repository: 
* DONE: Aug 10, 2024, 10:05: Start the JOSS paper template
* DONE: Aug 10, 2024, 16:39: Find other python libraries for vectorization
* DONE: Aug 11, 2024, 16:39: Incorporate my into ML workflow
* DONE: Aug 15, 2024, 07:38: Redo simulations: 
* DONE: Aug 15, 2024, 07:38: Correct table: 
* DONE: Aug 15, 2024, 08:50: Add dim0+dim1 set of predictors: 
* DONE: Aug 15, 2024, 11:54: Rewrite paper: 
* DONE: Aug 15, 2024, 13:00: Remove outliers: 
* DONE: Aug 15, 2024, 16:53: Send to Umar:
* DONE: Aug 16, 2024, 07:56: FDA in R 


## MISC

### Other python libraries:

* skit-tda: https://docs.scikit-tda.org/en/latest/libraries.html
* Giotto-TDA: https://giotto-ai.github.io/gtda-docs/0.5.1/index.html
* Geometry Understanding in Higher Dimensions: https://github.com/GUDHI, https://gudhi.inria.fr/introduction/
* Dionysus 2: https://www.mrzv.org/software/dionysus2/



The following functions were translated from R and checked against R results:
* `computeVPB`:     A Vector Summary of the Persistence Block
* `computePL`:      A Vector Summary of the Persistence Landscape Function
* `computePS`:      A Vector Summary of the Persistence Silhouette
* `computeNL`:      A Vector Summary of the Normalized Life Curve
* `computeVAB`:     A Vector Summary of the Betti Curve
* `computeECC`:     A Vector Summary of the Euler Characteristic Curve
* `computePES`:     A Vector Summary of the Persistent Entropy Summary Function
* `computePI`:      A Vector Summary of the Persistence Surface



