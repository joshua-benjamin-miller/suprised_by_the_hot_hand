This folder consists of the files necessary to reproduce Figure 1 in Section 2

The numbers have been cross-validated by simulation, and by enumeration.

How to reproduce the figure.
1. GraphProp.m builds the figure, it requires
	a. mainfigure.m 
	b. Ecprop.m
2. Ecprop.m requires that the data points PropFine.mat, it is the second, outside expected value.
3. PropFine.mat is produced with ProcessDataProp.m
4. ProcessDataProp.m requires the two data files PropData and PropDataFine.
5. PropData and PropDataFine are produced using ForServer.m
6. ForServer.m runs the function WriteProp.m twice to generate data
7. WriteProp.m calls the conditional expected value function Eprop.m many times to build the data.
8. Eprop.m is the basic formula for the conditional expected value.