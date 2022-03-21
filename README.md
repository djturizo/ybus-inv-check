# Admittance Matrix Invertibility Checking for Balanced Power Systems
This repository provides a proof-of-concept MATLAB R2012b code for certifying the invertibility of a balanced power system admittance matrix in time complexity linear with respect to the system size. The repository includes test cases from the OPF benchmark of [PGLib](https://github.com/power-grid-lib/pglib-opf).

The program is based on the results developed in the following article:

D. Turizo and D. Molzahn, "Invertibility Conditions for the Admittance Matrices of Balanced Power Systems," *IEEE Trans. Power Syst.*, submitted. See the [arXiv version](https://arxiv.org/abs/2012.04087).

The code description assumes that you read the article, as it makes multiple references to the theorems developed there.

The program does not compute the rank of the admittance matrix when there are no shunts. In order to do so, it is necessary to directly compute the rank of the generalized incidence matrix. Such computation can be done in linear time by means of a thorough implementation using graph theory, but such level of complexity is beyond the scope of this proof-of-concept program. 

If you ever find this repository useful, please cite the previous article.


# Code Description
Below you can find the description of the core functions of the program. The repository also includes additional auxiliary functions that are required by the core functions to work properly.

## `check_inv.m`
Check invertibility of admittance matrix

### Syntax

    flag = check_inv(mpc, tol);
    
This function attempts to apply Theorem 1 of the paper to certify the invertibility of the admittance matrix of a power system with `N` nodes and `L` lines. `mpc` is a struct with all the power system information, in [MATPOWER](https://github.com/MATPOWER/matpower) format. `tol` is a tolerance used for numerical comparations (mostly to determine wheter a given values is zero or not). `check_inv` returns an integer `flag`, which takes one of the following values:

* SUCCESS VALUES:
  * `1`: According to Theorem 1 the admittance matrix is invertible.
* CORNER-CASE VALUES:
  * `0`: There are no shunts, so according to Theorem 1 the rank of the admittance matrix equals the rank of the generalized incidence matrix. The matrix may be singular.
* FAILURE VALUES:
  * `-1`: Assumption 1 is violated: after reduction some series branches are zero or non-finite.
  * `-2`: Assumption 2 is violated: some taps ratio are zero.
  * `-4`: Passivity assumption of Theorem 2 is violated: after reduction some admittances were found to have negative conductance.
  * `-5`: Connectivity assumption of Theorem 1 is violated: the network graph is not connected.
  * `-6`: Some reactive component of the network does not have any of the required topologies of Theorem 3, the algorithm cannot proceed.
  * `-7`: Some reactive component of the network satisfies the topology of condition 1) of Theorem 3, but none of the conditions are satisfied anyway, so the algorithm cannot proceed.

A failure value indicates that the invertibility of the admittance matrix cannot be asserted. A sucess value indicates that the function sucessfully certified that the matrix is invertible. This function runs in `O(N+L)` time.


## `test_pglib.m`
Script for running the program on PGLib test cases

### Syntax

    test_pglib;
    
This script searches for all test cases stored in the `\pglib` folder, and runs `check_inv` on each of them. The test cases must be stored in the `\pglib` folder as `.m` functions. The functions must require no inputs, and must return as a single output an `mpc` structure in [MATPOWER](https://github.com/MATPOWER/matpower) format. The script uses a tolerance of `1e-12` for all calls of `check_inv`. After running the function over all test cases, the script stores the results in the file `results.mat`. The file stores the following variables (`C` is the number of test cases found):

* `cases`: `C` x `1` cell array. The i-th cell is a string with the filename of the i-th test case.
* `flags`: `C` x `1` vector. The i-th entry is the `flag` output of `check_inv` for the i-th test case.
* `times`: `C` x `1` vector. The i-th entry is the execution time of `check_inv` for the i-th test case.
* `nodes`: `C` x `1` vector. The i-th entry is the number of nodes of the i-th test case.
* `lines`: `C` x `1` vector. The i-th entry is the number of lines of the i-th test case.

After storing the results in `results.mat`, the script calls the function `print_results`.


## `print_results.m`
Print results of PGLib test cases in  console

### Syntax

    print_results();
    print_results(ignore_flags);
    
This function loads the results stored in `results.mat` and prints them in console. Optional argument `ignore_flags` is a vector whose entries are `flag` values of the function `check_inv` that will NOT be printed. This means that test cases with one of the `flag` values of `ignore_flags` will not be printed. When no input is provided no flags are ignored, and all test cases are printed.


### Sample output

    >> print_results();
    ________________________________________________________________________________
    |                                |       |       |Thm.   |Is YN |     |        |
    |Test case                       |N      |L      |holds? |inv.? |Flag |Time[s] |
    +--------------------------------+-------+-------+-------+------+-----+--------+
    |pglib_opf_case10000_goc         |10000  |13193  |Yes    |Yes   |1    |2.48    |
    |pglib_opf_case10480_goc         |10480  |18559  |Yes    |Yes   |1    |2.76    |
    |pglib_opf_case118_ieee          |118    |186    |Yes    |Yes   |1    |0.04    |
    |pglib_opf_case1354_pegase       |1354   |1991   |Yes    |Yes   |1    |0.25    |
    |pglib_opf_case13659_pegase      |13659  |20467  |No     |-     |-4   |1.18    |
    |pglib_opf_case14_ieee           |14     |20     |No     |-     |-6   |0.00    |
    +--------------------------------+-------+-------+-------+------+-----+--------+
    
    >>




