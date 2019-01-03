===================================================================================================================================================================================================================================================
Before installing the package CInLPN (Causal Inference in a Network of Latent Processes), please make sure you have the newest R version or R version >= 3.4.3 and make sure you have the followed packages or software installed on your computer:
===================================================================================================================================================================================================================================================
1- "Rcpp": this package is usefull to build R package that contains C++ code.
2- "Rtools": if you are on windows system, this R additionnal software is usefull to build R package.
3- "devtools": this package can be used to install package from github repository
4- "marqLevAlgParallel": this package is used for optimization task. It can be installed from github repository as follows : 
   devtools::install_github("VivianePhilipps/marqLevAlgParallel")


==========================
Installing CInLPN package
==========================
The installation is done from github repository with the followed command:
devtools::install_github("Bachirtadde/CInLPN")

=====
Note
=====
By installing CInLPN package, the dependencies packages RcppArmadillo, spline2, and MASS are also installed if not installed yet.