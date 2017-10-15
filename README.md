# The MOSOO algorithm for Multi-Objective Black-Box Optimization

This repository hosts the code for the multi-objective simultaneous optimistic optimisation algorithm. Paper can be found [here](https://arxiv.org/pdf/1612.08412.pdf).


# Contents

The algorithm is implemented in `MATLAB` and there is also a `C`-implementation (though not thoroughly tested). Scripts to generate the theoretical bounds and figures are also presented.


To run a demo, in `MATLAB`, `cd` to the `mosoo-demo` directory and execute the following
```
>>run_MOSOO_Demo
```


# Citation

If you write a scientific paper describing research that made use of this code, please consider citing the following paper:
```
@article{ALDUJAILI2018159,
title = "Multi-Objective Simultaneous Optimistic Optimization",
journal = "Information Sciences",
volume = "424",
number = "Supplement C",
pages = "159 - 174",
year = "2018",
issn = "0020-0255",
doi = "https://doi.org/10.1016/j.ins.2017.09.066",
url = "http://www.sciencedirect.com/science/article/pii/S0020025517309854",
author = "Abdullah Al-Dujaili and S. Suresh",
keywords = "Multi-objective optimization",
keywords = "Optimistic methods",
keywords = "Multi-armed bandits",
keywords = "Simultaneous optimistic optimization",
keywords = "Finite-time analysis",
keywords = "Asymptotic analysis"
}
```

