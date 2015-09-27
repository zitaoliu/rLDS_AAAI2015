# Matlab code for learning regularized linear dynamical systems

This code contains the learning procedures for regularized linear dynamical systems (rLDS) with two different regularizations (nuclear norm and group lasso).

This code is built on top of the [Kalman Filter Matlab package](http://www.cs.ubc.ca/~murphyk/Software/Kalman/kalman.html) from [Kevin P Murphy](http://www.cs.ubc.ca/~murphyk/).

The group lasso SOCP solver is developed by [Francis Bach](http://www.di.ens.fr/~fbach/), which can be download from [http://www.di.ens.fr/~fbach/grouplasso/](http://www.di.ens.fr/~fbach/grouplasso/).

Thanks very much to Kevin and Francis for sharing the code.

More mathematical details about this code, please check out our paper:

@inproceedings{liu2015regularized,
  title={A Regularized Linear Dynamical System Framework for Multivariate Time Series Analysis},
  author={Liu, Zitao and Hauskrecht, Milos},
  booktitle={Proceedings of the Twenty-Ninth AAAI Conference on Artificial Intelligence. AAAI Conference on Artificial Intelligence},
  volume={2015},
  pages={1798 - 1804},
  year={2015}
}
