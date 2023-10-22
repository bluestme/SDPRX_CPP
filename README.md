# SDPRX_CPP
SDPRX is a statistical method for cross-population prediction of complex traits. It integrates GWAS summary statistics and LD matrices from two populations (EUR and non-EUR) to compuate polygenic risk scores. The original SDPRX is implemented with Python at https://github.com/eldronzhou/SDPRX and SDPRX_CPP is the C++ version of SDPRX, which has been proven to significantlty improved the computation efficiency compared with the original Python version.

# Install

The SDPRX can be downloaded via
```
git clone https://github.com/bluestme/SDPRX_CPP.git
```
The SDPRX_CPP has been tested under both linux x86 and macOS arm environments. It works under both environments with different makefiles.
## x86 Linux 
If you plan to run SDPR_CPP on a linux system with a modern intel processor, you may use the precompiled binary `SDPRX`. Please make sure that dynamic libraries `gsl/lib/libgsl.so` and `MKL/lib/libmkl_rt.so` are not changed, otherwise SDPRX is unable to load the libraries. If you are not able to run the MKL library, you can use openBLAS instead as described in the section below. If you encounter errors reporting cannot load these two dynamic libraries, you may want to try export `LD_LIBRARY_PATH=$LD_LIBRARY_PATH:gsl/lib/:MKL/lib/`.

## Apple Silicon and arm architecture chips
Installation of gsl library is necessary if SDPRX_CPP is run on the new MacBook with Apple Silicon, which should be installed via `brew`
```
brew install gsl
```

# Quickstart
SDPRX_CPP can be run from the command line. To see the full list of options, please type
```
./SDPRX
```
SDPRX provides two functions: 
(1) Estimating and paritioning the reference LD matrix out of 2 population 
(2) perform MCMC to estimate the posterior effect sizes for each SNP of 2 populations. We provide an example usage for the test dataset:


# Make the reference LD

# Running SDPRX_CPP
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
