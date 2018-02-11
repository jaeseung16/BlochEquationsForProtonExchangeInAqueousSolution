---
title: Bloch equations for proton exchange reactions in an aqueous solution
layout: page
---

# Bloch equations for proton exchange reactions in an aqueous solution

The MATLAB codes to produce Figures 1, 2, and 4 in Ref. [1].

[1] Lee J-S, Regatte RR, Jerschow A. Bloch equations for proton exchange reactions in an aqueous solution. Concepts Magn Reson Part A. 2017;e21397. https://doi.org/10.1002/cmr.a.21397

#### The following are some corrections for the errors in Ref. [1].
1. In the 4th paragraph of Section 6, '10<sup>pK<sub>a</sub></sup> - 15.7' should be '10<sup>pK<sub>a</sub> - 15.7</sup>'.
2. In Eq. (38), HPO<sub>4</sub><sup>2-</sup> should not be boldfaced.
3. In Section 9, there are two citations of Figure 5D. They should be Figure 5F.

$$
\begin{align*}
  & \phi(x,y) = \phi \left(\sum_{i=1}^n x_ie_i, \sum_{j=1}^n y_je_j \right)
  = \sum_{i=1}^n \sum_{j=1}^n x_i y_j \phi(e_i, e_j) = \\
  & (x_1, \ldots, x_n) \left( \begin{array}{ccc}
      \phi(e_1, e_1) & \cdots & \phi(e_1, e_n) \\
      \vdots & \ddots & \vdots \\
      \phi(e_n, e_1) & \cdots & \phi(e_n, e_n)
    \end{array} \right)
  \left( \begin{array}{c}
      y_1 \\
      \vdots \\
      y_n
    \end{array} \right)
\end{align*}
$$
