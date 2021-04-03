
%Version 1.0
%Date: 19-June-2015


%The code ForServer.m generated 2 mat files
%we need to bring the data together
%Begin with the following data sets in the directory,
%and end with 'PropFine.mat'
data1 = load('PropData.mat');
data2 = load('PropDataFine.mat');

nvaluesfine = 2:1:30;
nvaluescoarse = 40:10:100;
kvalues = 2:5;
DATA=cell(length(kvalues),length(nvaluesfine)+length(nvaluescoarse));

for k=1:length(kvalues)
    for n=1:length(nvaluesfine)
    DATA{k,n}=data2.DATA{k,n}
    end
end


for k=1:length(kvalues)
    for n=1:length(nvaluescoarse)
        n1=n+length(nvaluesfine);
    DATA{k,n1}=data1.DATA{k,n+3} %start at n=40
    end
end


nvalues=[2:1:30 40:10:100];

nindex = @(n) find(nvalues==n);
kindex = @(k) find(kvalues==k);

clear data1 data2 k ans n n1
filename = 'PropFine.mat';

save(filename)

