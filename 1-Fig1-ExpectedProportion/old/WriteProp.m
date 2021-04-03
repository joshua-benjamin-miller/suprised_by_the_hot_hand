function WriteProp(fine)
%Version 1.0
%Date: 19-June-2015


%Description: This builds the a dataset used later for graphing purposes
%
%For n<=30 we use every integer
%
% disp('----------------------------')
% disp('----------------------------')
  
% 
%   rehash path
 %Matlab instructions 
 %nohup \matlab -nojvm -nodisplay < WriteProp.m >& WriteProp.log & 
%jobs 13087

% play withn vlaues and k values, filename = 'PropDataFine.mat';
addpath('C:\Users\Miller\Dropbox\josh\work\projects\Adam\HotHand\Data\statTheory\posted')
%%%%%
if fine == 0
    nvalues = 10:10:100;
elseif fine == 1
    nvalues = 2:1:30;
else
    fprintf('ERROR: Please enter fine= 0 or 1\n')
    return
end
kvalues = 2:5;
%nindex = @(n) [n==5]+2*[n==10];
nindex = @(n) find(nvalues==n);
kindex = @(k) k-1;
DATA=cell(length(kvalues),length(nvalues));



 formatSpec = '****** nvalues = 3:1:30; ******\n';
        fprintf(formatSpec)
   formatSpec = '****** kvalues = 2:5; ******\n';
   fprintf(formatSpec)
for n=nvalues
    for k=kvalues
        formatSpec = '****** n = %4.0f, k = %4.0f ******\n';
        fprintf(formatSpec,n,k)
        DATAkn=[];
        for n1=0:n
            formatSpec = 'n1 = %4.0f\n';
            fprintf(formatSpec,n1)
            tic
            if n1>=k
            [prop nCounted nUncounted]=Eprop( n1,n,k );
            else
                prop=0;
                nCounted=0;
                nUncounted=nchoosek(n,n1);

            end
            time=toc;
            DATAkn=[DATAkn; [k n n1 prop nCounted nUncounted time]];
        end
        DATA{kindex(k),nindex(n)}=DATAkn;
    end
end

clear diff DATAkn k kvalues n n1 nCounted nUncounted nvalues time
if fine == 0
    filename = 'PropData.mat';
elseif fine == 1
    filename = 'PropDataFine.mat';
else
    fprintf('Enter fine= 0 or 1')
    exit
end

save(filename)

end




% 
% clear
% load DiffPropData.mat
% p=.80;
% k=3;
% n=10;
% lookup=DATA{kindex(k),nindex(n)};
% n1 =lookup(:,3);
% diff=lookup(:,4);
% n=n*ones(size(n1));
% nCounted=lookup(:,5);
% nUncounted=lookup(:,6);
% 
% %Average over the counted sequences, and subt
% total=sum(nCounted.*(p.^n1.*(1-p).^(n-n1)).*diff)/(1-sum(nUncounted.*(p.^n1.*(1-p).^(n-n1))))
% 
% %Equivalent to the same divided by their total probability.
% total=sum(nCounted.*(p.^n1.*(1-p).^(n-n1)).*diff)/(sum(nCounted.*(p.^n1.*(1-p).^(n-n1))))
