Python realization for Vectorized Persistence Block vectorization

The goal is to implement TDAvec library in python using R as a source. The following functions needs to be translated:
* `computeVPB`:     A Vector Summary of the Persistence Block

    done

* `computePL`:      A Vector Summary of the Persistence Landscape Function

    done

* `computeECC`:     A Vector Summary of the Euler Characteristic Curve
* `computeNL`:      A Vector Summary of the Normalized Life Curve
* `computePES`:     A Vector Summary of the Persistent Entropy Summary Function
* `computePI`:      A Vector Summary of the Persistence Surface
* `computePS`:      A Vector Summary of the Persistence Silhouette Function
* `computeVAB`:     A Vector Summary of the Betti Curve

For each of the functions I should:
* translate the code from R to python in the python notebook
* create some results in R
* Make sure that python results agree with R
* translate the code into pyx format
* make sure that results agree with R

Finally, I should collect all files into one library

## File structure:
* `src/`: all source files live here
    * `computeVPB/`: folder related to computeVPB function:
        * `computeVPB.R`: R file to produce sample files to compare with
        * `computeVPB.ipynb`: python notebook that implements computeVPB function in python
        * `computeVPB.pyx`: cython version of computeVPB function
        * `computeVPB_test_pyx.ipynb`: python notebook that compares cython and R results
