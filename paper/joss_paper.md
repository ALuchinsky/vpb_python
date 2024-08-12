---
title: "Python Implemetation of VPB package"
tags:
authors:
affiliations:
date: Aug 10, 2024
---

# 1) Summary

aaa

# 2) Statement of Need

In the table below you can see list of some existing vectorization methods.

Method | this | Giotto-TDA | GUDHI | Persim | TDA | TDAvec
------- | ---- | --------- | ----- | ------ |  -- | ------
BC     | V     |  V        |  V    |        |     |  V
EEC     | V    |           |       |        |     |  V
NLC    | V    |            |       |        |     |  V
PES    | V    |            | V     |        |     |  V
PS     | V    |  V         | V     |        |  V  |  V
PL     | V    |  V         | V     |   V    |  V  |  V
PI     | V    |  V         | V     |   V    |     |  V
VPB    | V    |            |       |        |     |  V

See last sections for references.

# 3) Software Details

As it was noticed in the previous section, lots of different vectorization methods of the can be found in the literature. For a given persistence diagram
$$
PD = \{(b_i, d_i\}_{i=1}^N
$$
we can consider such quantities as

1) **Betti Curves**, when each value of dimension $d$ corresponds to function
$$
\beta_d(t) = \sum_{i=1}^N I_{[b_i, d_i)}(t),
$$
where $I_{[b,d)}(t)$ stands for the indicator function, which is equal to unity on the region $b\le t<d$ and vanishes outside.
In the following we will refer to this vectorization as **BC**.

2) **Euler Characteristic Curve**, which is a linear combination of the Betti Curves
$$
\chi(t) = \sum_{d} (-1)^d \beta_d(t)
$$
In the following it will be referred to as **EEC**.

3) **Normalized Line Curve**, where for each dimension $d$ we have
$$
s_d(t) = \sum_{i=1}^N \frac{d_i - b_i}{L} I_{[b_i, d_i]}(t),
$$
where
$$
L = \sum_{i=1}^N (d_i - b_i)
$$
In the following it will be referred to as **NLC**.

4) **Persistence Entropy Summary** function
$$
S_d(t) = \sum_{i=1}^N \frac{d_i-b_i}{L} \log_2\left(\frac{d_i - b_i}{L}\right) I_{[b_i, d_i)}(t)
$$
In the following it will be referred to as **PES**.

5) **Persistence Silhouette** function
$$
\phi_p(t) = \frac{1}{\sum_i|d_i-b_i|^p} \sum_{i=1}^N |d_i - b_i|^p \Lambda_i(t),
$$
where $p$ is a hyper-parameter of the model and a block triangle function $\Lambda$ is defined as
$$
\Lambda_i(t) = \begin{cases}
t - b_i \qquad & b_i \le t \le (b_i + d_i)/2,\\
d_i - t & (b_i + d_i)/2 \le t \le d_i,\\
0 & \textrm{otherwise}
\end{cases}
$$
In the following it will be referred to as **PS**.


5) **Persistence Landscape** function (PL in the following)
$$
\lambda_k(t) = \mathrm{arg}\max_{1\le i\le N} \Lambda_i(t)
$$

6) **Persistence Image** function (PI in the following)
$$
\rho(x, y) = \sum_{i=1}^N f(b_i, p_i) \phi_{(b_i, p_i)}(x,y),
$$
where
$$
\phi_{(b_i,d_i)}(x,y) = \frac{1}{\sqrt{2\pi\sigma^2}}\exp{-\frac{(x-b_i)^2 + (y-p_i)^2}{2\sigma^2}}
$$
is a Gauss distribution centered at the point $(b_i, p_i=d_i-b_i)$ and
$$
f(b, p) = w(p) = \begin{cases}
0 \qquad & p \le 0\\
p/p_{max} & 0 \le p \le p_{max},\\
1 & p \ge p_{max}.
\end{cases}
$$

7) **Vectorized Persistence Block** (VPB in the following)
$$
V(x, y) = \sum_{i=1}^N I_{E(b_i, p_i)}(x, y),
$$
where the indicator function is different from zero on the rectangle
$$
E(b_i, p_i) = \left[b_i - \frac{\lambda_i}{2}; b_i + \frac{\lambda_i}{2}\right] \times  \left[p_i - \frac{\lambda_i}{2}; p_i + \frac{\lambda_i}{2}\right]
$$

All defined above functions are now elements of some vector spaces and can be used in theoretical statistical analysis. In practical calculations, however, it is useful to digitize them and consider the values on some discrete 1-dimensional or 2-dimensional grids.



References:

Good overview of the existing libraries can be found here: https://github.com/FatemehTarashi/awesome-tda?tab=readme-ov-file#python

Below are some extractions:
* Python:
    * Giotto-TDA: https://giotto-ai.github.io/gtda-docs/0.5.1/index.html
    * GUDHI: https://gudhi.inria.fr/
    * Persim: https://persim.scikit-tda.org/en/latest/
* R
    * TDA: https://cran.r-project.org/web/packages/TDA/
    * TDAstats: https://github.com/rrrlw/TDAstats
    * TDAvec: https://cran.r-project.org/web/packages/TDAvec/index.html

# 4) Program Workflow

# 5) Conclusion

# 6) Acknowledgements

# 7) References


    