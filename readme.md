Spitz - Optimal edge weight reduction with targeted spectral radius
======================
This code implements the Spitz method as published in XXXXX [link]

Check `example.m` for an example usage

The methods expects three variables: 
* Tthe adjacency matrix `A` 
* The spectral radius `beta`
* An option structure `options` (optional)

First import two folders
```matlab
addpath('projections/'); 
addpath('lbfgs/');
``` 

### The edges problem:
To solve the problem using the dykstra method:
```matlab
X_dykstra = near_bounded_sparse(A, beta, options);
``` 
To solve the problem using the interior points method:
```matlab
X_intp = near_bounded_interior_p(A, beta, options);
```
To solve the problem using the dynamical importance method:
```matlab
X_dynamical_importance = binary_deletion_dynamical_importance(A, beta, options);
``` 

### The vertices problem:
For the vertices problem you also need to create a "vertices influence matrix":
```matlab
P = createVerticesInfluenceMatrix(A, 'equal_weight');
``` 
To solve the problem using the dykstra method:
```matlab
X_dykstra_v = near_bounded_sparse_vertices(A, beta, P, options);
``` 
To solve the problem using the interior points method:
```matlab
X_intp_v = near_bounded_interior_p_vertices(A, beta, P, options);
``` 
To solve the problem using the dynamical importance method:
```matlab
X_dynamical_importance_v = binary_deletion_dynamical_importance_vertices(A, beta, options);
``` 



To reproduce the experiment you can download the networks from koblenz network data:
Go to the `data` folder and run the script `download_networks.m`. 
It will create a subfolder under data called net_mat with the networks.

The `projections` folder contains projection function that are used by the methods.

The `data` folder contains data for the real world networks.

The `lbfgs` folder contains functions to run lbfgs with box constraints provided by Stephen Becker and Peter Carbonetto
http://www.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb--l-bfgs-b--mex-wrapper

