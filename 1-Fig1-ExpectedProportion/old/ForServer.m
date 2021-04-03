%Version 1.0
%Date: 19-June-2015


%This file needs to be run on the server as it takes a couple hours
%It runs WriteProp.m twice
%When complete we have two files:  (1) 'PropData.mat', (2) 'PropDataFine.mat'
%To run it use the UNIX command
%nohup \matlab -nojvm -nodisplay < ForServer.m >& ForServer.log & 
addpath('C:\Users\Miller\Dropbox\josh\work\projects\Adam\HotHand\Data\statTheory\posted')
rehash path
WriteProp(0)
WriteProp(1)