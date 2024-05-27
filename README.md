Python realization for Vectorized Persistence Block vectorization


The goal is to implement TDAvec library in python using R as a source (see GitHub repository https://github.com/alexey-luchinsky/TDAvec)


The following functions needs to be translated:
* `computeVPB`:     A Vector Summary of the Persistence Block

    done

* `computePL`:      A Vector Summary of the Persistence Landscape Function

    done

* `computePS`:      A Vector Summary of the Persistence Silhouette 

    done

* `computeNL`:      A Vector Summary of the Normalized Life Curve

    done

* `computeVAB`:     A Vector Summary of the Betti Curve

    done

* `computeECC`:     A Vector Summary of the Euler Characteristic Curve

    done

* `computePES`:     A Vector Summary of the Persistent Entropy Summary Function

    done

* `computePI`:      A Vector Summary of the Persistence Surface

    done

For each of the functions I should:
* translate the code from R to python in the python notebook
* create some results in R
* Make sure that python results agree with R
* translate the code into pyx format
* make sure that results agree with R

Finally, I should collect all files into one library

## File structure:
* python/: python files live here
    * TDAvec.pyx: Cython version of the package
    * TDAvec_pyx_test.ipynb: Jupyter notebook, comparison with R results
    * setup.py: Cython make file. See TDAvec_pyx_test.ipynb for instructions
* R/: R files live here
    * TDA_vec.Rmd: extracting test results for different functions