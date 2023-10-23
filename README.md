# Efficient_LocationRobust_TDOALoc_MPR
The signal emitted by an acoustic source may be propagating in an environment in which the speed is not known, such as in solid or ocean. Localization of such a source through observing the signal by a number of sensors requires joint estimation with the propagation speed.  This work applies the nullspace projection approach to the pseudo-linear formulation for the localization problem to obtain a closed-form solution, which is refined by error-compensation to reach the final estimation.  In contrast to the methods from the literature that are either suboptimal or computationally demanding, the proposed method is both statistically and computationally efficient, and is shown analytically to achieve the Cramer-Rao Lower Bound accuracy.

>[1] Y. Sun, K. C. Ho, Y. Yang, and L. Chen, "An asymptotically optimal estimator for source location and propagation speed by TDOA," *IEEE Signal Process. Lett.*, vol. 30, pp. 1037-1041, Aug. 2023.

**If you use any of the following codes in your research, please cite the paper as a reference in your publication. Thank you!**

## An asymptotically optimal estimator for source location and propagation speed by TDOA (10/23/2023)

### <u>Reference</u>
>Y. Sun, K. C. Ho, Y. Yang, and L. Chen, "An asymptotically optimal estimator for source location and propagation speed by TDOA," *IEEE Signal Process. Lett.*, vol. 30, pp. 1037-1041, Aug. 2023.

### Code List:
- Closed-Form Solution: TDOALoc_Proj_UPS
- Cramer Rao Bound: CRLB_TDOA_UPS
- Analytical Covariance: COV_Proj_UPS
- Example: main
