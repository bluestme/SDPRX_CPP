# SDPRX_CPP
SDPRX is a statistical method for cross-population prediction of complex traits. It integrates GWAS summary statistics and LD matrices from two populations (EUR and non-EUR) to compuate polygenic risk scores. The original SDPRX is implemented with Python at https://github.com/eldronzhou/SDPRX and SDPRX_CPP is the C++ version of SDPRX, which has been proven to significantly improve the computation efficiency compared with the original Python version.

# Install

The SDPRX can be downloaded via
```
git clone https://github.com/bluestme/SDPRX_CPP.git
```
The SDPRX_CPP has been tested under both linux x86 and macOS with Apple Silicon environments.
## x86 Linux 
If you plan to run SDPR_CPP on a linux system with a modern intel processor, you may use the precompiled binary `SDPRX`. Please make sure that dynamic libraries `gsl/lib/libgsl.so` and `MKL/lib/libmkl_rt.so` are not changed, otherwise SDPRX is unable to load the libraries. If you are not able to run the MKL library, you can use openBLAS instead as described in the section below. If you encounter errors reporting cannot load these two dynamic libraries, you may want to try export `LD_LIBRARY_PATH=$LD_LIBRARY_PATH:gsl/lib/:MKL/lib/`.

If you want to compile SDPRX from the source for best performance, you need a C/C++ compiler like g++ (tested under version 4.8.5), GSL (version 2.60) and MKL library (version 2017). For convenience, we redistribute the compiled GSL and MKL building on our Intel Xeon Processors in the repo. To install, type `make`. If this version does not work, please report the error to the issue page. If the issue is related to GSL, you may want to download the source code of GSL and compile it yourself. For details about downloading and installing GSL, please refer to [this page](https://www.gnu.org/software/gsl/) and [this page](https://www.gnu.org/software/gsl/doc/html/usage.html#compiling-and-linking). A tutorial of installing gsl is given on [this page](https://coral.ise.lehigh.edu/jild13/2016/07/11/hello/). If you have other versions of MKL library, please refer to [this manual](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html) for linking advice. On AMD CPU, we recommend using OpenBLAS, please see [here](https://github.com/eldronzhou/SDPR/issues/3) for instructions.

## Apple Silicon
Installation of gsl library is necessary if SDPRX_CPP is run on the new MacBook with Apple Silicon, which should be installed via `brew`
```
brew install gsl
```
Then you may end the version implemented for Apple Silicon with Arm chip architechture.
```
cd ./SDPRX_CPP-Apple
make
```
to install the package.

# Quickstart
SDPRX_CPP can be run from the command line. To see the full list of options, please type
```
./SDPRX -h
```
SDPRX provides two functions: 

(1) Estimating and paritioning the reference LD matrix out of 2 population with their shared individuals

(2) perform MCMC to estimate the posterior effect sizes for each SNP of 2 populations. 


# Make the reference LD

# Running SDPRX_CPP
Important options for running mcmc are:
- N1 (required): Sample size of the EUR summary statistics.
- N2 (required): Sample size of the non-EUR summary statistics.
- ss1 (required): Path to the EUR summary statistics.
- ss2 (required): Path to the non-EUR summary statistics.
- out1 (required) path to the output file for population 1 containing estimated effect sizes, end with .txt recommended
- out2 (required) path to the output file for population 2 containing estimated effect sizes, end with .txt recommended
- chr (required) chromsome to work on. Currently support 1-22. Recommend to run in parallel.
- load_ld (required): Path to the referecence LD directory.
- valid (required): Path to the bim file for the testing dataset, including the .bim suffix.
- rho (required): Trans-ethnic genetic correlation output by PopCorn between 0 and 1. Default is 0.8.
- force_shared (required): Whether to force sharing of effect sizes between populations. Default is True.
- n_threads (optional): number of threads to use. Default is 1.

An example to run the mcmc is
```
../SDPRX -mcmc -ref_dir test -ss1 test/SDPRX_EUR.txt -ss2 test/SDPRX_EAS.txt -valid test/Ukb_imp_v2.bim -N1 885541 -N2 116404 -out1 ./output1.txt -out2 ./output2.txt -chr 21
```
using the test file provided

## Summary Statistics
The summary statistics should at least contain following columns with the same name (order of the column is not important).
```
SNP	A1	A2	BETA	P
rs737657        A       G       -2.044  0.0409
rs7086391       T       C       -2.257  0.024
rs1983865       T       C       3.652   0.00026
...
```
where SNP is the marker name, A1 is the effect allele, A2 is the alternative allele, BETA is the regression coefficient for quantitative traits or log odds ratio for binary traits, and P is the p value.

# References
Zhou G, Chen T, Zhao H. SDPRX: A statistical method for cross-population prediction of complex traits. Am J Hum Genet. 2023 Jan 5;110(1):13-22. doi: 10.1016/j.ajhg.2022.11.007.
