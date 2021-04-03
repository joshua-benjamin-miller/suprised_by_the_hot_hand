*
*Author: Joshua B. Miller
*Send comments/questions to: joshua.benjamin.miller@gmail.com




**1.  Define the program
* When you run this program e.g. gensequences, nshots(10) nstreak(3)
* It enumerates every possible sequence of 10 shots, and calculates the difference
* in hit rate after hitting 3+ and missing 3+
cap program drop gensequences
program define gensequences
	syntax [, nshots(integer 5) nstreak(integer 2) freqstats(string)]
	quietly{
	if `nshots'>15{
	  	noi disp ""
		noi disp " *********************************"
		noi disp " ****** ERROR: TOO MANY SHOTS, THIS WILL TAKE TOO LONG!!
		noi disp " *********************************"
		exit
	
	}
	clear
	
	*local nshots=6
	*local nstreak=2
	local Lseq=`nshots'
	local Lstreak=`nstreak'
	*local Lseq=10
	*local Lstreak=2
	local Nseq=2^`Lseq'						// # of possibilities of sequences of length `n'
	set obs `Nseq'							// Set # of observations = # of sequences
	gen decimal=_n-1						// Decimal number
	gen str`Lseq' binary=""					// generate a string variable of length Lseq
	forvalues i=1/`Nseq'{
		di "hi"
		local x=decimal[`i']
		inbase 2 `x'						// convert integer to based 2
		replace binary=r(base) in `i'		// pull base 2
	}
	*deci decimal, f(10) t(2) gen(binary)	// Use stata module "deci" to convert decimal to binary
	*format %`Lseq'.0f binary
	*gen str`Lseq' seq0=string(binary)		// The the string version on the binary
	gen gap=`Lseq'-strlen(binary)				// The stings will not all have length Lseq, find the gap
	gen seq= gap*"0" + binary					// Construction the sequence

	

	drop decimal binary gap
	*the sequence strings are ready now

	forvalues i=1/`Lseq'{					//parse string and convert to numeric
		gen var`i'=real(substr(seq,`i',1))


	}
	*

	disp _N
	gen seqid=_n
	reshape long var, i(seqid) j(shot)  //put each sequence in long form so we can calculate easier
	rename var make
	gen streaktype=.
	*gen cirstreaktype=.


	*Calculate hit rate after Hitting `Lstreak'+ and after missing `Lstreak'+
	if `Lstreak' == 1{
		by seqid:replace streaktype=make[_n-1] if _n>1
		*by seqid:replace cirstreaktype=streaktype
	}
	else if `Lstreak' == 2{
		by seqid:replace streaktype=make[_n-1] if _n>2&make[_n-1]==make[_n-2]
	}
	else if `Lstreak' == 3{
		by seqid:replace streaktype=make[_n-1] if _n>3&make[_n-1]==make[_n-2] &make[_n-2]==make[_n-3]
	}
	else if `Lstreak' == 4{
		by seqid:replace streaktype=make[_n-1] if _n>4&make[_n-1]==make[_n-2] &make[_n-2]==make[_n-3] &make[_n-3]==make[_n-4]
	}
	else if `Lstreak' == 5{
		by seqid:replace streaktype=make[_n-1] if _n>5&make[_n-1]==make[_n-2] &make[_n-2]==make[_n-3] &make[_n-3]==make[_n-4] &make[_n-4]==make[_n-5]
	}
	else{
		noi disp "Lstreak =`Lstreak'"
		noi disp ""
		noi disp " *********************************"
		noi disp " ****** ERROR: Can only look at streaks of 1, 2 or 3 "
		noi disp " *********************************"
		exit
	}
	
	by seqid: gen endS=(make[_n]!=make[_n-1])&make[_n-1]==1 if _n>1
	by seqid: replace endS = make==1 if _n==_N
	by seqid: egen nSruns=total(endS)
	
	by seqid: gen alternation=(make[_n]!=make[_n-1]) if _n>1
	by seqid: egen nruns=total(alternation)
	replace nruns=nruns+1
	drop alternation endS
	
	gen dHstreak=.
	gen dMstreak=.
	gen dNonstreak=.
		replace dHstreak=(streaktype==1)
		replace dMstreak=(streaktype==0)
		replace dNonstreak=(streaktype==.)
	by seqid: egen nHstreak=total(dHstreak==1)
	by seqid: egen nCstreak=total(dMstreak==1)
	by seqid: egen nNonstreak=total(dNonstreak==1)
	by seqid: egen nHitHstreak=total(make==1&dHstreak==1)
	by seqid: egen nHitCstreak=total(make==1&dMstreak==1)
	by seqid: egen nHitNonstreak=total(make==1&dNonstreak==1)
	gen fgHstreak=nHitHstreak/nHstreak
	gen fgCstreak=nHitCstreak/nCstreak
	gen fgNotHstreak=(nHitCstreak+nHitNonstreak)/( nCstreak + nNonstreak )
	gen fgNotCstreak=(nHitHstreak+nHitNonstreak)/( nHstreak + nNonstreak )
	if "`freqstats'" == "yes"{
			by seqid:gen S1=make[_n]
			by seqid:gen S2=make[_n] if _n>1&make[_n]==make[_n-1]
			by seqid:gen S3=make[_n] if _n>2&make[_n]==make[_n-1] &make[_n-1]==make[_n-2]
			by seqid:gen S4=make[_n] if _n>3&make[_n]==make[_n-1] &make[_n-1]==make[_n-2] &make[_n-2]==make[_n-3]
			forvalues i=1/4{
				by seqid: egen nHf`i'=total(S`i'==1)
				by seqid: egen nCf`i'=total(S`i'==0)
				drop S`i'
			}
			
	}

	gen cons=.
	gen slope=.
	levelsof seqid,local(seqids)
	foreach seqid of local seqids{
		reg make dHstreak if seqid==`seqid' 
		replace cons=el(r(table),1,2) if  seqid==`seqid'
		replace slope=el(r(table),1,1) if  seqid==`seqid'
	}
	
	drop streaktype dHstreak dMstreak dNonstreak
	by seqid: egen nhits=total(make==1)
	
	** The predictor
	gen cum=.
	gen guess=1 if shot==1
	gen predrate=.
	gen pay=.
	
	by seqid:replace cum = (guess==make) if shot==1
	by seqid:replace predrate=cum  if shot==1
	by seqid:replace pay=predrate-.5 if shot ==1
	
	sort seqid shot

	levelsof shot,local(shots)
	foreach shot of local shots{
	    disp 
		if `shot'>1 {
			by seqid:replace guess=1 if pay[`shot'-1]<=0 & shot==`shot'
			by seqid:replace cum=cum[`shot'-1]+(guess==make) if shot==`shot'
			by seqid:replace predrate=cum/_n*(guess==1)+predrate[_n-1]*(guess==.) if shot==`shot'
			by seqid:replace pay=predrate-.5 if shot==`shot'
		}
	
	}
	by seqid: replace predrate=predrate[_N]
	by seqid:gen predictorgain=pay[_N]
	drop cum guess pay
	

	reshape wide make, i(seqid) j(shot)  // Return to short form


	sort nhits seqid
	replace seqid=_n

	*order nhits fgHstreak fgCstreak
	gen diff=fgHstreak-fgCstreak
	order nhits seq diff fgHstreak fgNotHstreak fgCstreak fgNotCstreak
	drop make1-make`Lseq'
	sort nhits diff
	*list
	}
end
*

capture:program drop pHH
    program define pHH, rclass
		syntax [, nshots(integer 5) ]
		tempname sum term weight pHH
		scalar `weight' = 0
		scalar `sum' = 0
		forvalues k=1/`nshots'{
			 if `k' == 1{
				scalar `weight' = comb(`nshots',`k') - 1
			 
			 }
			 else {
				scalar `weight' = comb(`nshots',`k')
			 
			 }
			
			 scalar `term' = `weight'*(`k'-1)/(`nshots'-1)
			
			 scalar `sum' = `sum' + `term'
			 
			 * noi di "k = " `k'
			* noi di "weight = " `weight'
			 *	 noi di "term = " `term'
		}
		*noi di "sum = " `sum'
		scalar `pHH' = `sum'/(2^`nshots'-2) // for 1 back, 2 sequences cannot be counted
		return scalar pHH = `pHH'
	end
	
capture:program drop probabilities
    program define probabilities, rclass
		syntax [, nshots(integer 5) ]
		tempname sumHH sumHT sumdiff termHH termHT termdiff weightHH weightHT weightdiff pHHk pHTk pHH pHT diff
		
		scalar `sumHH' = 0
		scalar `sumHT' = 0
		scalar `sumdiff' = 0
		forvalues k=0/`nshots'{
			 if `k' == 0{
			 	scalar `weightHH' = 0
				scalar `weightHT' = comb(`nshots',`k')
			}
			 else if `k' == 1{
				scalar `weightHH' = comb(`nshots',`k') - 1
				scalar `weightHT' = comb(`nshots',`k')
			 }
			 else if `k' == `nshots' -1 {
				scalar `weightHH' = comb(`nshots',`k') 
				scalar `weightHT' = comb(`nshots',`k') - 1
			 }
			 else if `k' == `nshots' {
				scalar `weightHH' = comb(`nshots',`k') 
				scalar `weightHT' = 0
			 }
			 else {
				scalar `weightHH' = comb(`nshots',`k')
				scalar `weightHT' = `weightHH'
			 
			 }
			 scalar `weightdiff' = min(`weightHH',`weightHT') // works only because they differ by 1, won't work for 2+ and 3+
		
		*Compute probabilities and terms:
			if `weightHH' != 0{
				scalar `pHHk' = (`k'-1)/(`nshots'-1)
				scalar `termHH' = `weightHH'*`pHHk'
			}
			else{
				scalar `termHH' = 0
			}
			
			if `weightHT' != 0{
				scalar `pHTk' = (`k')/(`nshots'-1)
				scalar `termHT' = `weightHT'*`pHTk'
			}
			else{
				scalar `termHT' = 0
			}
			
			scalar `sumHH' = `sumHH' + `termHH'
			scalar `sumHT' = `sumHT' + `termHT'
			*The weight is tricky because pHHk and pHTk are computed differently
			*Gotta pull out certain sequences and recalculate
			*first modify pHHk with one fewer sequence for k=N-1
			if `k' == 1{
				scalar `weightHH' = comb(`nshots',`k') - 1
				scalar `weightHT' = comb(`nshots',`k')
			 }
			 else if `k' == `nshots' -1 {
				scalar `pHHk' = ( (`nshots'-2)*(`k'-2)/(`nshots'-2)+1 )/(`nshots' - 1)
				scalar `weightHT' = comb(`nshots',`k') - 1
			 }
			
			/*
			if `weightdiff' != 0{
				scalar `termdiff' = `weightdiff'*(`pHHk'- `pHTk')
			
			}
			else{
				scalar `termdiff' = 0
			
			}
			*/
			 
			 * scalar `sumdiff' = `sumdiff' + `termdiff'
			 * noi di "k = " `k'
			* noi di "weight = " `weight'
			 *	 noi di "term = " `term'
		}
		*noi di "sum = " `sum'
		scalar `pHH' = `sumHH'/(2^`nshots'-2)
		scalar `pHT' = `sumHT'/(2^`nshots'-2)
		return scalar pHH = `pHH'
		return scalar pHT= `pHT'
	end

capture:program drop expected_diff
    program define expected_diff, rclass
	 syntax [, nshots(integer 5) streak_length(integer 2) p(real .5)  ]
	
	gensequences, nshots(`nshots') nstreak(`streak_length')
	gen prob=`p'^(nhits)*(1-`p')^(`n'-nhits)
		drop if diff==.
	gen wdiff=prob*diff
	egen tprob=total(prob)
	egen twdiff=total(wdiff)
	gen avdiff=twdiff/tprob
	return scalar e_diff= avdiff[1]
end


* COPY-PASTE of one of the simulations codes i had lying around

cap program drop subcurry
program define subcurry, rclass
	syntax [, nshots(integer 100) fg(real .5) fgcurry(real .10) ncurry(integer 20)]
clear
	set obs `nshots'
	gen rand1=runiform()
	gen rand2=runiform()
	gen period=_n
	gen make=(rand1<=`fg')
	
	*gen hot=((make[_n-1]==1&make[_n-2]==1&make[_n-3]==1) & rand2<`randomhot' )
	replace make =(rand2<=`fgcurry') if period>= `nshots'/2-`ncurry'/2 & period<`nshots'/2+`ncurry'/2

		
	gen dMade1ormore=(make[_n-1]==1) if _n>1
	gen dMissed1ormore=(make[_n-1]==0) if _n>1
	
	
	gen dMade2ormore=(make[_n-1]==1&make[_n-2]==1) if _n>2
	gen dMissed2ormore=(make[_n-1]==0&make[_n-2]==0) if _n>2
	
	gen dMade3ormore=(make[_n-1]==1&make[_n-2]==1&make[_n-3]==1) if _n>3
	gen dMissed3ormore=(make[_n-1]==0&make[_n-2]==0&make[_n-3]==0) if _n>3
	
	summ make if dMade3ormore==1
			return scalar fg3H=r(mean)
			summ make if dMade2ormore==1
			return scalar fg2H=r(mean)
			summ make if dMade1ormore==1
			return scalar fg1H=r(mean)
			
			summ make if dMissed3ormore==1
			return scalar fg3M=r(mean)
			summ make if dMissed2ormore==1
			return scalar fg2M=r(mean)
			summ make if dMissed1ormore==1
			return scalar fg1M=r(mean)
			
	
end
subcurry, fg(.50) fgcurry(.50) ncurry(20)

* 40, 60 for 20 of 100 shots, you will find nothing
* 35, 45 for 20 of 100 shots, you will find -7
* 35, 45 for 20 of 200 shots, you will find -4
* 35, 45 for 40 of 200 shots, you will find -2
* 40, 60 for 20 of 100 shots, you will find nothing

simulate  fg3M=r(fg3M) fg2M=r(fg2M) fg1M=r(fg1M) ///
			fg3H=r(fg3H) fg2H=r(fg2H) fg1H=r(fg1H) ///
					, reps(1000) : subcurry, nshots(100) fg(.4) fgcurry(.6) ncurry(40)
		
exit
*************************************************************

************BELOW HERE DOES't EXECUTE UNLESS SELECTED INTERACTIVELY


************************************************************
***INDEPENDENLTY VERIFY MATLAB CODE
**Test Diff


expected_diff, nshots(5) streak_length(2) p(.5)

**Test prop
		local n=10
		local p=.25
		local k=4
		gensequences, nshots(`n') nstreak(`k')
		drop if fgHstreak==.
		gen prob=`p'^(nhits)*(1-`p')^(`n'-nhits)
		gen wprop=prob*fgHstreak
		egen tprob=total(prob)
		egen twprop=total(wprop)
		gen avprop=twprop/tprob

	local n=5
	local p=.5
	local k=3
	gensequences, nshots(`n') nstreak(`k')

	collapse fgHstreak (count) C=fgHstreak (count) N=seqid, by(nhits)
		

** 2. Explore

*Below we see with a session of 10 shots, there are 2^10=1024 possible performance outcomes
*Of these only 148 will allow us to compare Hit3+ with Miss3+

gensequences, nshots(8) nstreak(1)
count
keep if diff!=.
count
list

local nshots = 6
gensequences, nshots(`nshots') nstreak(3)
gen pmk=(nhits-3)/(`nshots'-3) if nhits>=3
by nhits: summ fgHstreak pmk


local nshots = 15
local lstreak =2
gensequences, nshots(`nshots') nstreak(`lstreak')
gen palternation=nHitCstreak/nCstreak
summ palternation

local nshots = 10
gensequences, nshots(`nshots') nstreak(5)
gen palternation=nHitCstreak/nCstreak
summ palternation


*Show the alternation rate is still bad if heads and tails are not distinguished
gensequences, nshots(10) nstreak(3)
gen Nalt=nHitCstreak+(nHstreak-nHitHstreak)
gen altrate=Nalt/(nHstreak+nCstreak)
summ altrate

***CALCULATE

*Lets enumerate all the sequences
*What we learned:
* H4/H3- M4/M3 is very far away from the diff in proportions because of the end point issue

*In addition, we verify the formula for 1 back
local nshots = 10
local lstreak =1
gensequences, nshots(`nshots') nstreak(`lstreak') freqstats("no")

gen pHC=fgCstreak
gen pHH=fgHstreak
cap drop pHCe pHHe
*gen pHCe = (1-nCf2/(nCf1))
*gen pHHe = nHf2/(nHf1)
*gen diff1=pHH-pHC

*order nhits seq diff pHH pHHe pHC pHCe

*Lets get the actual expected 
sum pHH pHC diff

pHH, nshots(`nshots')
return list

probabilities, nshots(`nshots')
return list

gen pmk=(nhits-3)/(`nshots'-3) if nhits>=3
by nhits: summ fgHstreak pmk

*This is the expected value of the statistic based on a 50-50 bernoulli shooter
summ diff

*But it will be true even if they aren't Benoulli, to see this
*Lets group "sessions" (sequence) by # of hits
*Within each group each sequence is equally likely regardless of the bernoulli prior
*So we can average across sequences in each group with the same # of hits, 
*for each group we see these averages are always negative, 
*so even if they have difference probabilities
*there will always be a bias.  The quantity of the bias will depeend on the prior, of course.
by nhits: summ diff

exit
*We can play around with the other stats, if we want
by nhits: summ fgHstreak fgCstreak

by nhits: summ fgHstreak if fgCstreak==.

by nhits: summ fgCstreak if fgHstreak==.

summ fgHstreak if fgCstreak==.

summ fgCstreak if fgHstreak==.

summ fgHstreak if fgCstreak!=.

summ fgCstreak if fgHstreak!=.
summ diff


summ fgHstreak

summ fgCstreak


*Calculate expected bias under various base rate assumptions
cap drop rate prob tprob probwdiff
gen rate=.9
gen prob=rate^nhits*(1-rate)^(6-nhits)
egen tprob=total(prob)
gen probwdiff=(prob/tprob)*diff
summ probwdiff
