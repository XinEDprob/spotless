### Description
This repo is a fork of https://github.com/anirudhamajumdar/spotless. We added the support of power cone optimization with Mosek (https://www.mosek.com/) for paper ```A new characterization of symmetric H^+-tensors and M-tensors```.


### Installation instructions

1. git clone the repo with 
```
https://github.com/XinEDprob/spotless.git
``` 
2. intall the spotless tool with ```spot_install.m```
3. install Tensor Toolbox for MATLAB v3.1 for MATLAB from https://www.tensortoolbox.org/index.html
4. download ```allcomb(varargin)``` from Matlab file exchange: https://www.mathworks.com/matlabcentral/fileexchange/10064-allcomb-varargin. Unzip the downloaded file and put the unzipped in folder spotless.
5. install Mosek 9.3 from https://www.mosek.com
6. (optional) if you want to compare our methods with the one proposed in TenEig (https://epubs.siam.org/doi/abs/10.1137/15M1010725), you also need to download and install TenEig from https://users.math.msu.edu/users/chenlipi/teneig.html

### Running experiments

The implementation of our methods for Table 1-3 and Example 4-7 are contained in ```tensor_paper/experiments.m```

For implementation of Example 8 
- ```tensor_paper/Fan_product_M_tensors.m``` is for the computation of H-eigenvalues for our methods. 
- ```lb_H_eigenvalue_Fan_Prodcut.m``` is to compute the lower bound with the H-eigenvalues from ```tensor_paper/Fan_product_M_tensors.m```.
- ```other_method_H_eigen_M_tensor.m``` is for the computation of H-eigenvalues for other methods in the comparison.
- ```other_method_H_eigen_Fan.m``` is to compute the lower bound with the H-eigenvalues from ```other_method_H_eigen_M_tensor.m```

The experiments with TenEig (https://users.math.msu.edu/users/chenlipi/teneig.html) (for part of Table 2) are contained in ```tensor_paper/homotopy_H_eig.m```

