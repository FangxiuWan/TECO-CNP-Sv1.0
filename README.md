# TECO-CNP Sv1.0

This respository contains the code that accompanies the article 
TECO-CNP Sv1.0: A coupled carbon-nitrogen-phosphorus model with data assimilation for subtropical forests
Wan et al., 2025 [DOI]. 

Please refer to the README file for the steps to follow before running the model.

# Check compiler
This code runs with gfortran compile, e.g.
```console
@TECO-CNP~: gfortran --version
GNU Fortran (GCC) 4.8.5 20150623 (Red Hat 4.8.5-44)
Copyright (C) 2015 Free Software Foundation, Inc.
GNU Fortran comes with NO WARRANTY, to the extent permitted by law.
You may redistribute copies of GNU Fortran
under the terms of the GNU General Public License.
For more information about these matters, see the file named COPYING
```

# Install the Dependencies
To run the TECO-CNP model, you need to install the following packages:

- **LAPACK**: A library for solving complex linear algebra problems, providing efficient routines for matrix computations.
- **BLAS**: Basic Linear Algebra Subprograms, a foundational library that provides basic vector and matrix operations.
- **LAPACKE**: A C interface to LAPACK, which allows easier integration of LAPACK functionality into C/C++ programs.

# Create compiled Fortran modules and run model
To generate the compiled Fortran modules and run model, please execute the runtecocnp.sh script in your terminal using the following command
During the peer review stage, a one-year forcing dataset has been made available to demonstrate the model's functionality.

```console
@TECO-CNP~: ./runtecocnp.sh
```
# State of Input files 
The example files for driving model is provided, only the format is referable and the data can obtain by corresponding the authors.
Note: Example files containing forcing data and flux data are provided for reference. Please replace these with your own data as needed.

