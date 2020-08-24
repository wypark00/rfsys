# Codes to solve a linear rational expectations model

## gensysx

*gensysx* solves a linear rational expectations model using the ordinary Schur decomposition when the model is invertible (non-singular) and the QZ (generalized Schur) decomposition when the model is not invertible (singular). 

**Note that the order of outputs has changed from that of original gensys by Chris Sims.**

Date: May 20, 2020.

### Update (August 1, 2020)

*gensysx* now returns the intermediate outputs if `zy=false` for efficient evaluation of the likelihood of a linear DSGE model. If `zy=false`, a real form of the Schur decomposition or the QZ decomposition is used and the solution only for the stable block is returned. If you want *gensysx* to return the solution for the whole vector of the endogenous variables, set `zy=true`.

### Download

- Right click *gensysx.m* above and choose "Save the link as..." or click the green "Code" button on the top to download everything.
- [An example of system reduction described in Lee and Park (2020b)](/efflkh/gensysx)

## rfsys

*rfsys* solves a reduced-form linear rational expectations model using the ordinary Schur decomposition. The code is optimized to take advantage of the ordinary Schur decomposition.

Date: May 20, 2020.

### Update (August 1, 2020)

*rfsys* now returns the intermediate outputs if `zy=false` for efficient evaluation of the likelihood of a linear DSGE model. If `zy=false`, the solution only for the stable block is returned. If you want *rfsys* to return the solution for the whole vector of the endogenous variables, set `zy=true`.

### Download

- Right click *rfsys.m* above and choose "Save the link as..." click the green "Code" button on the top to download everything.
- [An example of system reduction described in Lee and Park (2020b)](/efflkh/rfsys)

## References

- Lee and Park, 2020a, "[Solving Reduced-form Linear Rational Expectations Models](https://drive.google.com/file/d/1cRdCQWVO3J1u7F06hJ0WMrWdZh_T6gQT/view?usp=sharing)," Working paper.
- Lee and Park, 2020b, "An Efficient Likelihood Evaluation of Dynamic Stochastic General Equilibrium Models," Working paper.
- Sims, 2002, "[Solving Linear Rational Expectations Models](https://doi.org/10.1023/A:1020517101123)," *Computational Economics*, 20(1-2), 1-20.
