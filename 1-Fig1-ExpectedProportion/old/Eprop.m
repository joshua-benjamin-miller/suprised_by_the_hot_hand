function [f nC nU]= Eprop( n1,n,k )

% %Version 2.0 
%Date: 19-June-2015


%Description: This code computes the result of the formula reported in Appendix E,
%Theorem 5. The formula computes the exectaed value of the proportion of
%successes on trials that are immediately preceded by k consecutive successes,
%conditional on (i) there being n_1 successes in the sequence, and
%(ii) there being at least one instance of a trial immediately precede by 
%k consecutive successes.

%The  output:
%   1. f: the expected proportion given n1 hits in n trials.  It is the arithmetic average across all nchoosek(n,n1) sequences
%that are measurable, so we need to subtract out those sequences which are not measureable.
%   2. nU: the number of sequences in which the proportion cannot be
%   calculated because the denominator is zero.  This is typically a small
%   fraction of the sequences and can be checked, e.g. for n1=25,n=50,k=3:
%       [f A B]=Eprop( n1,n,k );
%       B/nchoosek(n,n1)=0.0024
%
% NEW TO Version 2.0: improved efficiency by combining cases and using
% binomial identities.

% With this function for each n we can build a look up table with rows n1 and columns k% 
% we have average proportions,

%Execution Notes:
% 1. With floating point arithmatic, nchoosek and like functions cease to be exact at high values
%   they remain relatively close, as can be seen when we sum the
%   probabilities and we are within 10^(-16) of the true probabilities.
%2. %for k=5 and larger the code pretty slow, perhaps the indexing isn't
%efficient?  We do have combinitorial explosion for large k, which was the
%reason this approach was taken to begin with.


%%%%%%%%%%%%%

%For testing purposes
% disp('----------------------------')
% disp('----------------------------')
%  addpath('C:\Users\Miller\Dropbox\josh\work\projects\Adam\HotHand\Data\statTheory')
%  rehash path
   % n=13;
 %   n1=7;
  % k=3;
n0=n-n1;


if n1>n
    error('n1 must be less than n');
    return,
end

if k>n1
    error('k must be less than n1');
    return,
end

if k<2
    error('k must be greater than or equal to 2')
    return
end
cumseqnotcounted=0;
cumprop=0;



%disp(' ****** Case 1 ****** ')
% Count the sequences for which the proportion is not
% defined
%%This has been cross checked several times to make sure it sums accurately.
% Use Verson 1.0 code pieces below to cross-check

%disp(' ****** Case 1i ****** ')
%Case 1i: No streaks of k or more:
tic
nseq1i=0;
for r1=1:min(n1,n0+1)
    ulim=min( floor( (n1-r1)/(k-1) ), r1);
    nmult=0;
    for t1=0:ulim
        t2=n1-r1-t1*(k-1);
        nmult=nmult+ (-1)^t1*nchoosek(r1,t1)*nchoosek(r1-1+t2,r1-1);
    end
    nseq1i=nseq1i+nmult*nchoosek(n0+1,r1);
    
end

%disp(' ****** Case 1ii ****** ')
%Case 1ii
%Here we count the number of sequences that end in a run of k, with all
%other runs less than k length.
nseq1ii=0;
if n1>k
    for r1=2:min(n1-k+1,n0+1)
        ulim=min( floor( (n1-k-(r1-1))/(k-1) ), r1-1);
        nmult=0;
        for t1=0:ulim
            t2=n1-k-(r1-1)-t1*(k-1);
            nmult=nmult+ (-1)^t1*nchoosek(r1-1,t1)*nchoosek(r1-2+t2,r1-2);
        end
        nseq1ii=nseq1ii+nmult*nchoosek(n0,r1-1);
    end
elseif n1==k
    nseq1ii=1; %only one way to end with k 1s in a row
else
    error('Error in control flow')
    return
end
cumseqnotcounted=nseq1i+nseq1ii;


%disp(' ****** Case 2 ****** ')
% pk=(fk-s1k)/(fk-1)
%Case 2: End in run of 1s of k+ && s1k>=1 (excluding the case of s1k=1 and
%rk=1)
%more)

%We have multidimensional indices, and the number of dimensions depend on k, lets put them into a single dimension so
%we need only one loop (This could be done more efficiently with nested loops that depend on each other).
UpperLimitsFixedLength1runs= floor(n1*ones(1,k-1)./[1:k-1]); %The higher number of runs for each fixed run length
IndexCountFixedLength1runs=ones(1,k-1)+UpperLimitsFixedLength1runs; %include the first zero term
UpperLimitFixedLength1runsLinearSum=prod(IndexCountFixedLength1runs);
runcell = cell( 1, k-1 );
for i=1:UpperLimitFixedLength1runsLinearSum
    [ runcell{:} ] = ind2sub(IndexCountFixedLength1runs,i);  %go from linear index to multi-dimensional index
    FixedLength1runs=cell2mat(runcell)-ones(1,k-1); %Instead of 1:floor(n1/j)+1 it is 0:floor(n1/j), [r11 r12.. r1k-1]
    nOnesInFixedLength1runs=sum( ( 1:length(FixedLength1runs) ).*FixedLength1runs );
   
    if nOnesInFixedLength1runs<n1-k %we could reduce the number of loops before by taking away k, but lets just add control flow, don't lose much time.
        for s1k=1:floor( (n1-nOnesInFixedLength1runs)/k )  %s1k must start at 1
             fstat = freqstat(FixedLength1runs,s1k,n1);
             prop1 = (fstat-s1k)/(fstat-1);  %prop when ends with run of k or more
             prop0 = (fstat-s1k)/(fstat);  %prop there is at least 1 zero in the final k positions
             avprop=  ( s1k/( n0+1 ) )*prop1+( (n0+1-s1k)/( n0+1 ) )*prop0;
             cumprop = cumprop + NseqsRunsOfOneType( FixedLength1runs,s1k,n1,n )*avprop;
        end
    end
end


%Return Values
%Note for large n1,n0, floating point arithmetic means 

nU=cumseqnotcounted;
nC=nchoosek(n,n1)-nU;
f=cumprop/nC;
% 
% % TESTING
% % Below output for testing, to make sure things add correctly
% % Confirm the partitition of the r11,r12,s13 sequences, that total count, divided by
% % nchoosek(n,n1), sums to 1.
% 
% 
% disp('-----------------------------')
% disp('************  Final Count ****************')
% disp('-----------------------------')
% 
% disp('Counted Sequences:')
% Tseq
% 
% disp('UnCounted Sequences:')
% cumseqnotcounted
% 
% disp('Total Sequences:')
% mycount=Tseq+cumseqnotcounted
% 
% disp('True Count')
% disp('Combination Function')
% truecount=nchoosek(n,n1)
% 
% 
% disp('Absolute Difference')
% diff=abs(truecount-mycount)
% 
% disp('Relative Difference: This is the key check') 
% pdiff=abs(1-mycount/truecount)
% 
% disp('-----------------------------')
% disp('************  Expected Proportion ****************')
% 
% disp('Actual Fraction 1s out of n')
% n1/n
% 
% 
% disp('Expected proportion')
% cumprop/Tseq


end

% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %                   VERSION 1.0
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [f nC nU]= Eprop( n1,n,k )
% 
% %Date: 24-May-2015
% %Author: Prof. Joshua B. Miller
% 
% %Description: This code computes
% %   1. f: the expected proportion given n1 hits in n trials.  It is the arithmetic average across all nchoosek(n,n1) sequences
% %that are measurable, so we need to subtract out those sequences which are not measureable.
% %   2. nU: the number of sequences in which the proportion cannot be
% %   calculated because the denominator is zero.  This is typically a small
% %   fraction of the sequences and can be checked, e.g. for n1=25,n=50,k=3:
% %       [A B]=Eprop( n1,n,k );
% %       B/nchoosek(n,n1)=0.0024
% 
% % With this function for each n we can build a look up table with rows n1 and columns k% 
% % we have average proportions, 
% 
% %Execution Notes:
% % 1. With floating point arithmatic, nchoosek and like functions cease to be exact at high values
% %   they remain relatively close, as can be seen when we sum the
% %   probabilities and we are within 10^(-16) of the true probabilities.
% %2. %for k=5 and larger the code pretty slow, perhaps the indexing isn't
% %efficient?  We do have combinitorial explosion for large k, which was the
% %reason this approach was taken to begin with.
% 
% 
% %%%%%%%%%%%%%
% 
% disp('----------------------------')
% disp('----------------------------')
%  addpath('C:\Users\Miller\Dropbox\josh\work\projects\Adam\HotHand\Data\statTheory')
%  rehash path
%    % n=13;
%  %   n1=7;
%   % k=3;
% n0=n-n1;
% 
% 
% if n1>=n
%     error('n1 must be less than n');
%     return,
% end
% 
% if k>n1
%     error('k must be less than n1');
%     return,
% end
% 
% cumseqnotcounted=0;
% cumprop=0;
% Tseq=0;
% nseq=0;
% 
% 
% 
% %Case 1: pk= not observable
% disp(' ****** Case 1i ****** ')
% % pk= not observable
% 
% %************  proportion not measureable************
% 
% %We have multidimensional indices, and the number of dimensions depend on k, lets put them into a single dimension so
% %we need only one loop.
% % note r1j=0,1,2,...,floor(n1/j) for j=1,2,...,k-1
% UpperLimitsFixedLength1runs= floor(n1*ones(1,k-1)./[1:k-1]); %The higher number of runs for each fixed run length
% IndexCountFixedLength1runs=ones(1,k-1)+UpperLimitsFixedLength1runs; %include the first zero term ( matrix consisting k-1 nonneg integer, each integer representing the number of different values the runs can take)
% UpperLimitFixedLength1runsLinearSum=prod(IndexCountFixedLength1runs);
% runcell = cell( 1, k-1 );
% for i=1:UpperLimitFixedLength1runsLinearSum
%     [ runcell{:} ] = ind2sub(IndexCountFixedLength1runs,i); %go from linear index to multi-dimensional index
%     FixedLength1runs=cell2mat(runcell)-ones(1,k-1); %Instead of 1:floor(n1/j)+1 it is 0:floor(n1/j)
%     
%     nOnesInFixedLength1runs=sum( ( 1:length(FixedLength1runs) ).*FixedLength1runs );
%         
%     
%     if nOnesInFixedLength1runs == n1   %
%         s1k=0;  %it cannot be greater, so no need to loop over values of s1k
%         nseq=nseq+NseqsRunsOfOneType( FixedLength1runs,s1k,n1,n );
%     end
%      %[i1 i2] = [i1 i2]
% end
% 
% 
% 
% 
% Tseq=Tseq+nseq;
% 
% nseq=0;
% 
% disp(' ****** Case 1ii ****** ')
% % pk= not observable
% %Case 1ii: End in run of 1s of k+ && s1k=1 &&  nOnesInFixedLength1runs=n1-k (i.e. r1k=1)
% %************  proportion not measureable************
% 
% %We have multidimensional indices, and the number of dimensions depend on k, lets put them into a single dimension so
% %we need only one loop.
% UpperLimitsFixedLength1runs= floor(n1*ones(1,k-1)./[1:k-1]); %The higher number of runs for each fixed run length
% IndexCountFixedLength1runs=ones(1,k-1)+UpperLimitsFixedLength1runs; %include the first zero term
% UpperLimitFixedLength1runsLinearSum=prod(IndexCountFixedLength1runs);
% runcell = cell( 1, k-1 );
% for i=1:UpperLimitFixedLength1runsLinearSum
%     [ runcell{:} ] = ind2sub(IndexCountFixedLength1runs,i); %go from linear index to multi-dimensional index
%     FixedLength1runs=cell2mat(runcell)-ones(1,k-1); %Instead of 1:floor(n1/j)+1 it is 0:floor(n1/j)
%     nOnesInFixedLength1runs=sum( ( 1:length(FixedLength1runs) ).*FixedLength1runs );
%         
%     
% 
%     if nOnesInFixedLength1runs == n1 - k  %
%        s1k=1; %all the remaining k ones must form a run of length k, so no need to loop over values of s1k
%         %r
%         nseq=nseq+( s1k/( n0+1 ) )*NseqsRunsOfOneType( FixedLength1runs,s1k,n1,n );
%     end
%      %[i1 i2] = [i1 i2]
% end
% 
% 
% 
% 
% 
% Tseq=Tseq+nseq;
% cumseqnotcounted=Tseq;
% Tseq=0;
% nseq=0;
% 
% 
% 
% disp(' ****** Case 2 ****** ')
% % pk=(fk-s1k)/(fk-1)
% %Case 2: End in run of 1s of k+ && s1k>=1 (excluding the case of s1k=1 and
% %rk=1)
% %more)
% 
% %We have multidimensional indices, and the number of dimensions depend on k, lets put them into a single dimension so
% %we need only one loop.
% UpperLimitsFixedLength1runs= floor(n1*ones(1,k-1)./[1:k-1]); %The higher number of runs for each fixed run length
% IndexCountFixedLength1runs=ones(1,k-1)+UpperLimitsFixedLength1runs; %include the first zero term
% UpperLimitFixedLength1runsLinearSum=prod(IndexCountFixedLength1runs);
% runcell = cell( 1, k-1 );
% for i=1:UpperLimitFixedLength1runsLinearSum
%     [ runcell{:} ] = ind2sub(IndexCountFixedLength1runs,i);  %go from linear index to multi-dimensional index
%     FixedLength1runs=cell2mat(runcell)-ones(1,k-1); %Instead of 1:floor(n1/j)+1 it is 0:floor(n1/j), [r11 r12.. r1k-1]
%     nOnesInFixedLength1runs=sum( ( 1:length(FixedLength1runs) ).*FixedLength1runs );
%         
%     for s1k=1:floor( (n1-nOnesInFixedLength1runs)/k )
%         %Don't nee
%          if nOnesInFixedLength1runs < n1 - k  %
%              nseq=nseq+( s1k/( n0+1 ) )*NseqsRunsOfOneType( FixedLength1runs,s1k,n1,n );  %Same as 2ii   
%              fstat = freqstat(FixedLength1runs,s1k,n1);
%              mstat = (fstat-s1k)/(fstat-1);  %because ends with run of k or more
%              cumprop = cumprop + mstat* ( s1k/( n0+1 ) )*NseqsRunsOfOneType( FixedLength1runs,s1k,n1,n );
%          end
%      end
% end
% 
% 
% 
% Tseq=Tseq+nseq;
% 
% nseq=0;
% 
% % OLD CASE ii & iii: used to distinguish between 
% % disp(' ****** Case 2ii ****** ')
% % %Case 2ii:  End in run of 1s of k+ && s1k=1 &&  nOnesInFixedLength1runs<n1-k (i.e. r1k=0)
% % %more)
% % 
% % %We have multidimensional indices, and the number of dimensions depend on k, lets put them into a single dimension so
% % %we need only one loop.
% % UpperLimitsFixedLength1runs= floor(n1*ones(1,k-1)./[1:k-1]); %The higher number of runs for each fixed run length
% % IndexCountFixedLength1runs=ones(1,k-1)+UpperLimitsFixedLength1runs; %include the first zero term
% % UpperLimitFixedLength1runsLinearSum=prod(IndexCountFixedLength1runs);
% % runcell = cell( 1, k-1 );
% % for i=1:UpperLimitFixedLength1runsLinearSum
% %     [ runcell{:} ] = ind2sub(IndexCountFixedLength1runs,i);
% %     FixedLength1runs=cell2mat(runcell)-ones(1,k-1); %Instead of 1:floor(n1/j)+1 it is 0:floor(n1/j)
% %     nOnesInFixedLength1runs=sum( ( 1:length(FixedLength1runs) ).*FixedLength1runs );
% %         
% %     s1k=1;
% % 
% %     if nOnesInFixedLength1runs < n1 - k  %
% %         nseq=nseq+( s1k/( n0+1 ) )*NseqsRunsOfOneType( FixedLength1runs,s1k,n1,n );
% %         fstat = freqstat(FixedLength1runs,s1k,n1);
% %         mstat = (fstat-s1k)/(fstat-1); %because ends with run of k or more
% %         cumprop = cumprop + mstat* ( s1k/( n0+1 ) )*NseqsRunsOfOneType( FixedLength1runs,s1k,n1,n );
% %     end
% % end
% % 
% % 
% % 
% % Tseq=Tseq+nseq;
% % 
% % nseq=0;
% % 
% % disp(' ****** Case 2iii ****** ')
% % %Case 2iii: End in run of 1s of k+ && s1k>=2 s
% % %more)
% % 
% % %Actually, can combine with case ii!! check it. 
% % 
% % 
% % %We have multidimensional indices, and the number of dimensions depend on k, lets put them into a single dimension so
% % %we need only one loop.
% % UpperLimitsFixedLength1runs= floor(n1*ones(1,k-1)./[1:k-1]); %The higher number of runs for each fixed run length
% % IndexCountFixedLength1runs=ones(1,k-1)+UpperLimitsFixedLength1runs; %include the first zero term
% % UpperLimitFixedLength1runsLinearSum=prod(IndexCountFixedLength1runs);
% % runcell = cell( 1, k-1 );
% % for i=1:UpperLimitFixedLength1runsLinearSum
% %     [ runcell{:} ] = ind2sub(IndexCountFixedLength1runs,i);
% %     FixedLength1runs=cell2mat(runcell)-ones(1,k-1); %Instead of 1:floor(n1/j)+1 it is 0:floor(n1/j)
% %     nOnesInFixedLength1runs=sum( ( 1:length(FixedLength1runs) ).*FixedLength1runs );
% %         
% %     for s1k=2:floor( (n1-nOnesInFixedLength1runs)/k )
% %         %Don't nee
% %          nseq=nseq+( s1k/( n0+1 ) )*NseqsRunsOfOneType( FixedLength1runs,s1k,n1,n );  %Same as 2ii   
% %          fstat = freqstat(FixedLength1runs,s1k,n1);
% %          mstat = (fstat-s1k)/(fstat-1);  %because ends with run of k or more
% %          cumprop = cumprop + mstat* ( s1k/( n0+1 ) )*NseqsRunsOfOneType( FixedLength1runs,s1k,n1,n );
% %      end
% % end
% 
% 
% disp(' ****** Case 3i ****** ')
% %Case 3i:  End in run of 1s of less than k && s1k>=1 
% %less)
% 
% %We have multidimensional indices, and the number of dimensions depend on k, lets put them into a single dimension so
% %we need only one loop.
% UpperLimitsFixedLength1runs= floor(n1*ones(1,k-1)./[1:k-1]); %The higher number of runs for each fixed run length
% IndexCountFixedLength1runs=ones(1,k-1)+UpperLimitsFixedLength1runs; %include the first zero term
% UpperLimitFixedLength1runsLinearSum=prod(IndexCountFixedLength1runs);
% runcell = cell( 1, k-1 );
% for i=1:UpperLimitFixedLength1runsLinearSum
%     [ runcell{:} ] = ind2sub(IndexCountFixedLength1runs,i);
%     FixedLength1runs=cell2mat(runcell)-ones(1,k-1); %Instead of 1:floor(n1/j)+1 it is 0:floor(n1/j)
%     nOnesInFixedLength1runs=sum( ( 1:length(FixedLength1runs) ).*FixedLength1runs );
%     
%      
% 
%    for s1k=1:floor( (n1-nOnesInFixedLength1runs)/k )
%      r1=sum([FixedLength1runs s1k]); %total runs
%      nseq=nseq+( (r1-s1k)/( n0+1 ) )*NseqsRunsOfOneType( FixedLength1runs,s1k,n1,n );
% 
% 
%      fstat = freqstat(FixedLength1runs,s1k,n1);
%      mstat = (fstat-s1k)/fstat;  %because ends with run of k-1 or less
%      cumprop = cumprop + mstat*( (r1-s1k)/( n0+1 ) )*NseqsRunsOfOneType( FixedLength1runs,s1k,n1,n );
%    end
% 
% end
% 
% 
% 
% Tseq=Tseq+nseq;
% 
% nseq=0;
% 
% disp(' ****** Case 3ii ****** ')
% %Case 3ii: End in run of 0s  && s1k>=1 
% 
% %We have multidimensional indices, and the number of dimensions depend on k, lets put them into a single dimension so
% %we need only one loop.
% UpperLimitsFixedLength1runs= floor(n1*ones(1,k-1)./[1:k-1]); %The higher number of runs for each fixed run length
% IndexCountFixedLength1runs=ones(1,k-1)+UpperLimitsFixedLength1runs; %include the first zero term
% UpperLimitFixedLength1runsLinearSum=prod(IndexCountFixedLength1runs);
% runcell = cell( 1, k-1 );
% for i=1:UpperLimitFixedLength1runsLinearSum
%     [ runcell{:} ] = ind2sub(IndexCountFixedLength1runs,i);
%     FixedLength1runs=cell2mat(runcell)-ones(1,k-1); %Instead of 1:floor(n1/j)+1 it is 0:floor(n1/j)
%     nOnesInFixedLength1runs=sum( ( 1:length(FixedLength1runs) ).*FixedLength1runs );
% 
%     for s1k=1:floor( (n1-nOnesInFixedLength1runs)/k )
%         r1=sum([FixedLength1runs s1k]); %total runs
%         if r1<=n0
%           %  disp('adding the term')
%            % runs
%             nseq=nseq+( (n0-r1+1)/( n0+1 ) )*NseqsRunsOfOneType( FixedLength1runs,s1k,n1,n );
% 
%             fstat = freqstat(FixedLength1runs,s1k,n1);
%              mstat = (fstat-s1k)/fstat;  %because ends with run of k-1 or less
%              cumprop = cumprop + mstat*( (n0-r1+1)/( n0+1 ) )*NseqsRunsOfOneType( FixedLength1runs,s1k,n1,n );
%         end
%     end    
% 
% end
% 
% 
% 
% 
% Tseq=Tseq+nseq;
% 
% %Return Values
% f=cumprop/Tseq;
% nU=cumseqnotcounted;
% nC=Tseq;
% 
% % 
% % % TESTING
% % % Below output for testing, to make sure things add correctly
% % % Confirm the partitition of the r11,r12,s13 sequences, that total count, divided by
% % % nchoosek(n,n1), sums to 1.
% % 
% % 
% % disp('-----------------------------')
% % disp('************  Final Count ****************')
% % disp('-----------------------------')
% % 
% % disp('Counted Sequences:')
% % Tseq
% % 
% % disp('UnCounted Sequences:')
% % cumseqnotcounted
% % 
% % disp('Total Sequences:')
% % mycount=Tseq+cumseqnotcounted
% % 
% % disp('True Count')
% % disp('Combination Function')
% % truecount=nchoosek(n,n1)
% % 
% % 
% % disp('Absolute Difference')
% % diff=abs(truecount-mycount)
% % 
% % disp('Relative Difference: This is the key check') 
% % pdiff=abs(1-mycount/truecount)
% % 
% % disp('-----------------------------')
% % disp('************  Expected Proportion ****************')
% % 
% % disp('Actual Fraction 1s out of n')
% % n1/n
% % 
% % 
% % disp('Expected proportion')
% % cumprop/Tseq
% 
% 
% 
% 
% 
% 
% 
% end
