*The main analysis file is:

Analysis_Surprised.do

note:
*Figure 2 has dependencies on files generated in another directory.
* In particular, One section of the file file depends on having
 
 biasGVT-0-0.dta

*as these have player specific bias-adjustments
*The file was created by the .ipynb Jupyter file in 


*This file is meant to be opened, and run block-by-block
*to check the results & comments of the paper.


*other dependencies include
gilovichshooters.dta
genregressionvars_v2.do
3ptcode.do

*Note, the "old" folder contains the previous approach to generate bias-adjustments, i.e. using simulation.
BIAS2.do
BIAS3.do


