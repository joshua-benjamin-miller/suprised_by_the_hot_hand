* Aggregate & Inidvidual Conditional probability tests


* Version 2.0

capture:program drop cptest
    program define cptest, rclass
	syntax [, nstreak(integer 3) vsmiss(integer 0) dropcondition(string)]
		preserve
			sort sid session period // Permute may misorder periods
			tempvar hotstate coldstate
			
			if `nstreak' == 2{
				by sid session:gen `hotstate'=(make[_n-1]==1&make[_n-2]==1) if _n>2
				if `vsmiss' == 1{
					by sid session:gen `coldstate' =(make[_n-1]==0&make[_n-2]==0) if _n>2
				}
			}
			else{		
				by sid session:gen `hotstate'=(make[_n-1]==1&make[_n-2]==1&make[_n-3]==1) if _n>3
				if `vsmiss' == 1{
					by sid session:gen `coldstate'=(make[_n-1]==0&make[_n-2]==0&make[_n-3]==0) if _n>3		
				}
			}
			
			if `vsmiss' == 1{
			    drop if `hotstate'==0 & `coldstate'==0
				drop if `coldstate'==.
			}
			
		
			drop if `hotstate'==.
			
			cap drop if `dropcondition'
		
			collapse (mean) av=make (count) n=make, by(sid `hotstate')
			gen fg_h=.
			gen fg_c=.
			gen n_h=. 	
			gen n_c =.

			sort sid `hotstate'
			by sid:replace fg_c=av[1]
			by sid:replace fg_h=av[2]
			by sid:replace n_c=n[1]
			by sid:replace n_h=n[2]
			drop if `hotstate'==0  //gets rid of duplicates & guys who don't have observations when hotstate== 1
			drop `hotstate' av n
			*signrank fg_c=fg_h
			
			gen diff=fg_h-fg_c
			*
			ttest diff ==0
	
			return scalar mean=r(mu_1)
			return scalar tstat=r(t)
			restore
		end

	
capture:program drop momtest
    program define momtest, rclass
	syntax [, nstreak(integer 3) vsmiss(integer 0)]
		inspect sid
			if r(N_unique)!=1{
				di "ERROR: only 1 sid allowed"
				exit
			}
		preserve
			
			sort sid session period // Permute may misorder periods
			tempvar hotstate coldstate
			
			if `nstreak' == 2{
				by sid session:gen `hotstate'=(make[_n-1]==1&make[_n-2]==1) if _n>2
				if `vsmiss' == 1{
					by sid session:gen `coldstate' =(make[_n-1]==0&make[_n-2]==0) if _n>2
				}
			}
			else{		
				by sid session:gen `hotstate'=(make[_n-1]==1&make[_n-2]==1&make[_n-3]==1) if _n>3
				if `vsmiss' == 1{
					by sid session:gen `coldstate'=(make[_n-1]==0&make[_n-2]==0&make[_n-3]==0) if _n>3		
				}
			}
			
			if `vsmiss' == 1{
			    drop if `hotstate'==0 & `coldstate'==0
				drop if `coldstate'==.
			}
			prtest make,by(`hotstate')
			return list
			local diff=r(P_2)-r(P_1)
			di `diff'
			return scalar fghot=r(P_2)
			return scalar diff=r(P_2)-r(P_1)
			return scalar tstat=r(z)
		end

		
	

*** SHOW THE BIAS IN 3+

*This code is copied from the power analysis
cap program drop genshots
program define genshots
	syntax [, nshots(integer 900) nsids(integer 10) fgpercent(real .5) ]
clear
	noi display "local obs = `nshots'*`nsids'"
	local obs = `nshots'*`nsids'
	*noi display "HI"
	set obs `obs'
	*noi display "HI"
	gen rand1=runiform()
	scalar fg=`fgpercent'
	*noi display "HI"
	gen sid=floor((_n-1)/`nshots')+1
	by sid,sort: gen period=_n
	gen make=(rand1<=scalar(fg))

	gen session=1
	sort sid period
	
	di "GENERATE STREAK 6 VARIABLES"
	
	by sid:gen dMade4ormore=(make[_n-1]==1&make[_n-2]==1&make[_n-3]==1&make[_n-4]==1) if _n>4
	by sid:gen dMissed4ormore=(make[_n-1]==0&make[_n-2]==0&make[_n-3]==0&make[_n-4]==0) if _n>4
	
	by sid:gen dMade3ormore=(make[_n-1]==1&make[_n-2]==1&make[_n-3]==1) if _n>3
	by sid:gen dMissed3ormore=(make[_n-1]==0&make[_n-2]==0&make[_n-3]==0) if _n>3
	
	by sid:gen dMade2ormore=(make[_n-1]==1&make[_n-2]==1) if _n>2
	by sid:gen dMissed2ormore=(make[_n-1]==0&make[_n-2]==0) if _n>2
	
end

*This code takes the 3pt data, expands it into `expand' copies
* The generates random shots bassed on a player's overall hit rate
cap program drop genshots3pt
program define genshots3pt
	syntax, [expand(integer 1) dropcutoff(integer 0) ]
	noi display "HI"
	clear
	use $allstar
	*local expand=3
	replace sid = sid-1000 
	by sid,sort:egen fg=mean(make)
	drop make
	count
	local nobs=r(N)
	
	summ sid
	local nsids=r(max)
	
	gen id=_n
	expand `expand', generate(xcopy)
	sort id xcopy
	by id:gen replica=_n-1
	
	
	
	replace sid=sid+replica*`nsids'
	sort sid session period
	drop id xcopy
	noi display "HI"
	drop if nshots<`dropcutoff'
	gen rand1=runiform()

	gen make=(rand1<=fg)

	sort sid session period
	by sid session:gen dMade3ormore=(make[_n-1]==1&make[_n-2]==1&make[_n-3]==1) if _n>3
	by sid session:gen dMissed3ormore=(make[_n-1]==0&make[_n-2]==0&make[_n-3]==0) if _n>3
	by sid session:gen dMade2ormore=(make[_n-1]==1&make[_n-2]==1) if _n>2
	by sid session:gen dMissed2ormore=(make[_n-1]==0&make[_n-2]==0) if _n>2
	
end

*Program computes the different in hit rates, to get at the bias
cap program drop testbias
program define testbias
	syntax [, twoplus(integer 0) alldata(integer 0)  streak(integer 3) ]
	preserve
	qui{
		gen treatment=.
		local k=`streak'
		
		if `twoplus'==1{
			local k=2
		} // end if
		if `alldata'==0{ //Madek+ vs. Missed k+
			if `k' == 2{ 
				replace treatment=dMade2ormore
					drop if dMade2ormore==. | dMissed2ormore ==.
					drop if dMissed2ormore==0 & dMade2ormore==0
			}
			else if `k'==3{
				replace treatment=dMade3ormore
					drop if dMade3ormore==. | dMissed3ormore ==.
					drop if dMissed3ormore==0 & dMade3ormore==0
			}
			else if `k'==4 {
				replace treatment=dMade4ormore
					drop if dMade4ormore==. | dMissed4ormore ==.
					drop if dMissed4ormore==0 & dMade4ormore==0
			}
			else{
					di "ERROR: k defined improperly = `k'"
			}
		}
		else if `alldata' == 1{ //Madek+ vs. other
			if `k' == 2{ 
				replace treatment=dMade2ormore
					drop if dMade2ormore==. 

			}
			else if `k'==3{
				replace treatment=dMade3ormore
					drop if dMade3ormore==. 

			}
			else if `k'==4{
				replace treatment=dMade4ormore
					drop if dMade4ormore==. 
	
			}
			else{
					di "ERROR: k defined improperly = `k'"
			}	
		}
		else if `alldata' == 2{ // Missedk+ vs. other
			if `k' == 2{
				replace treatment=dMissed2ormore
					drop if dMissed2ormore==.
				}
			else if `k' == 3{
				replace treatment=dMissed3ormore
					drop if dMissed3ormore==. 
			}
			else if `k' == 4{
				replace treatment=dMissed4ormore
					drop if dMissed4ormore==. 
			}
			else{
				di "ERROR: k defined improperly = `k'"
			}
		}
		else{
			di "ERROR: alldata defined improperly"
			di "ERROR: alldata  = `alldata'"
			exit
		}
		
			collapse (mean) av=make (count) n=make, by(sid treatment)
			gen fg_T1=.
			gen fg_T0=.
			gen n_T1=.
			gen n_T0 =.

			sort sid treatment
			by sid:replace fg_T0=av[1]
			by sid:replace fg_T1=av[2]
			by sid:replace n_T0=n[1]
			by sid:replace n_T1=n[2]
			drop if n_T0==.|n_T1==. // drop if we can't measure both
			drop if treatment==0 //drop the duplicates
			drop treatment av n
			
			gen ehot = fg_T1-fg_T0
	}
		qui: count if (fg_T1>fg_T0)
		*di r(N)/_N
		summ ehot

		*hist ehot
		
		
		*ttest fg_c=fg_h
		*signrank fg_c=fg_h
		
		
	restore	
end





**1.  Define the program
cap program drop gensequences
program define gensequences
	syntax [, nshots(integer 5) nstreak(integer 2)]
	quietly{
	if `nshots'>15{
	  	noi disp ""
		noi disp " *********************************"
		noi disp " ****** ERROR: TOO MANY SHOTS, THIS WILL TAKE TOO LONG!!
		noi disp " *********************************"
		exit
	
	}
	clear
	local Lseq=`nshots'
	local Lstreak=`nstreak'
	*local Lseq=10
	*local Lstreak=2
	local Nseq=2^`Lseq'						// # of possibilities of sequences of length `n'
	set obs `Nseq'							// Set # of observations = # of sequences
	gen decimal=_n-1						// Decimal number
	gen str`Lseq' binary=""
	forvalues i=1/`Nseq'{
		di "hi"
		local x=decimal[`i']
		inbase 2 `x'						// convert to binary
		replace binary=r(base) in `i'		//record binary value
	}
	*deci decimal, f(10) t(2) gen(binary)	// Use stata module "deci" to convert decimal to binary
	*format %`Lseq'.0f binary
	*gen str`Lseq' seq0=string(binary)		// The the string version on the binary
	gen gap=`Lseq'-strlen(binary)				// The stings will not all have length Lseq, find the gap
	gen seq= gap*"0" + binary					// Construction the sequence

	

	drop decimal binary gap

	forvalues i=1/`Lseq'{					//parse string and convert to numeric
		gen var`i'=real(substr(seq,`i',1))


	}
	*

	disp _N
	gen seqid=_n
	reshape long var, i(seqid) j(shot)  //put each sequence in long form so we can calculate easier
	rename var make
	gen streakoutcome=.

	*Calculate hit rate after Hitting `Lstreak'+ and after missing `Lstreak'+
	if `Lstreak' == 1{
		by seqid:replace streakoutcome=make[_n-1] if _n>1
	}
	else if `Lstreak' == 2{
		by seqid:replace streakoutcome=make[_n-1] if _n>2&make[_n-1]==make[_n-2]
	}
	else if `Lstreak' == 3{
		by seqid:replace streakoutcome=make[_n-1] if _n>3&make[_n-1]==make[_n-2] &make[_n-2]==make[_n-3]
	}
	else{
		noi disp ""
		noi disp " *********************************"
		noi disp " ****** ERROR: Can only look at streaks of 1, 2 or 3 "
		noi disp " *********************************"
		exit
	}
	by seqid: egen nHstreak=total(streakoutcome==1)
	by seqid: egen nCstreak=total(streakoutcome==0)
	by seqid: egen nHitHstreak=total(make==1&streakoutcome==1)
	by seqid: egen nHitCstreak=total(make==1&streakoutcome==0)
	gen fgHstreak=nHitHstreak/nHstreak
	gen fgCstreak=nHitCstreak/nCstreak
	
	
	drop streakoutcome
	by seqid: egen nhits=total(make==1)

	reshape wide make, i(seqid) j(shot)  // Return to short form


	sort nhits seqid
	replace seqid=_n

	*order nhits fgHstreak fgCstreak
	gen diff=fgHstreak-fgCstreak
	order nhits seq diff fgHstreak fgCstreak
	drop make1-make`Lseq'
	sort nhits diff
	*list
	}
end
