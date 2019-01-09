===================================================================================================================================================================================================================================================
Before installing the package CInLPN (Causal Inference in a Network of Latent Processes), please make sure you have the newest R version or R version >= 3.4.3 and make sure you have the followed packages or software installed on your computer:
===================================================================================================================================================================================================================================================
1- "Rcpp": this package is useful to build R packages that contain C++ code.
2- "Rtools": if you are on windows system, this R additional software is useful to build R packages.
3- "devtools": this package can be used to install packages from github repository
4- "marqLevAlgParallel": this package is used for the optimization. It can be installed from github repository as follows : 
   devtools::install_github("VivianePhilipps/marqLevAlgParallel")


==========================
Installing CInLPN package
==========================
The installation is done from github repository with the following command:
devtools::install_github("Bachirtadde/CInLPN")

=====
Note
=====
When installing CInLPN package, the dependencies packages RcppArmadillo, spline2, and MASS are also installed if not installed yet.
