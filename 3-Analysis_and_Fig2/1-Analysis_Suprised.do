*Date: 18-Aug-2017

*Analysis for "Surprised by the Hot Hand Fallacy"

*This file is annated and intended to be run by highlights relevant sections
*and then Ctrl+D to run it, not ran through from beginning to end.

*1. Gxplaining the original test, player by player
* Why it is biased	
	
	*First point out that 3 of GVTs stats, autocorrelation, runs, and slope of regression (diff in proportion) are all the same
	*We have a proof in the old paper, buts its simpler just to simulate
	*Note we don't have an _rc code below, which we should have (like in the testbias program further below)
	cap program drop genshots
	program define genshots, rclass
		syntax [, nshots(integer 900) fgpercent(real .5) ]
		clear
			noi display "local obs = `nshots'*`nsids'"
			local obs = `nshots'
			set obs `obs'
			gen rand1=runiform()
			scalar fg=`fgpercent'
			gen period=_n
			gen make=(rand1<=scalar(fg))

			gen dMade1ormore=(make[_n-1]==1) if _n>1
			gen dMissed1ormore=(make[_n-1]==0) if _n>1
			
			*Surprised Proves the following 2 are identical
			reg make dMade1ormore
			return scalar beta=_b[dMade1ormore]
			
			summ make if dMade1ormore==1
				local pH=r(mean)
			summ make if dMissed1ormore==1
				local pM=r(mean)
			return scalar diff=`pH'-`pM'
			
			*Cold shower proves the following 2 are identical (for permuations of n/2 hits out of n attempts)
			*Otherwise they are nearly perfectly correlated
			runtest make, threshold(.5)
			return scalar runs=r(n_runs)
			
			tsset period
			corrgram make, lags(1)
			return scalar autocorr=r(ac1)
			
				
			*Here is the old school way of measuring serial correlation
			gen x1=make if _n<_N
			gen x2=make[_n+1]
			corr x1 x2
			return scalar rho=r(rho)		
	end

	genshots, nshots(300) 
	return list
	
	simulate beta=r(beta) rho=r(rho) autocorr=r(autocorr) runs=r(runs) diff=r(diff), reps(1000): genshots, nshots(100) 
	
	*when simulation is done we can see the correlation.
	graph matrix autocorr runs diff

	*Step 1: Naive Correction (Pooled test)
	clear
	use gilovichshooters.dta
	do genregressionvars_v2.do  //generate the variables that are used in the regression
		genregressionvars
			replace dMade4ormore=0 if period<5	
			replace dMade3ormore=0 if period<4
			replace dMade2ormore=0 if period<3
			replace dMissed4ormore=0 if period<5
			replace dMissed3ormore=0 if period<4
			replace dMissed2ormore=0 if period<3
		*If control is perfect and they are 50% shooters:
		*there no reason to systematically underestimate the better shooters
		*but as long as some players have a better chance than others because of our error
		*there will be an aggregation bias
			prtest make if (dMade4ormore==1|dMissed4ormore==1),by(dMade4ormore)
			prtest make if (dMade3ormore==1|dMissed3ormore==1),by(dMade3ormore)
			ttest make if (dMade3ormore==1|dMissed3ormore==1),by(dMade3ormore)
			prtest make if (dMade2ormore==1|dMissed2ormore==1),by(dMade2ormore)
			
			prtest make ,by(dMade3ormore)
			prtest make, by(dMade2ormore)
			
			
		*right, but they didn't do a good job positioning them
		clear
		use gilovichshooters.dta
	
		levelsof sid, local(sidlist) // create text list of sid values
			matrix A=J(wordcount("`sidlist'"),3,.)
			*Fill the Matrix
			local it = 0
			
			foreach sid in `sidlist' { //iterate over each sid
				local ++it
				 prtest make=.5 if sid==`sid'
				 matrix A[`it',1]=`sid'
				 matrix A[`it',2]=r(P_1)
				  matrix A[`it',3]=normal(-abs(r(z)))
			}
		clear
		svmat A, names(col)
		count if c3<.05
		*12 out of 26 are different than .5, didn't do such a great job
				 

*Step 2: Bias table.  Please See AnalysisWriteUP4_v2 for how to creat latex tables:
* Plus player level test
		clear
			discard
			use gilovichshooters.dta
	
			do genregressionvars_v2.do  //generate the variables that are used in the regression
				genregressionvars	
		
		*2.1 Clean up
			*THIS IS AN IMPORTANT MODICICATION,
			replace dMade4ormore=0 if period<5
			replace dMade3ormore=0 if period<4
			replace dMade2ormore=0 if period<3
			replace dMissed4ormore=0 if period<5
			replace dMissed3ormore=0 if period<4
			replace dMissed2ormore=0 if period<3
			
				by sid, sort:egen nshots=total(make!=.)
				by sid, sort:egen nhits=total(make==1)
				gen fgp=nhits/nshots
				forvalues k=2/4{
					by sid, sort:egen nhit`k'=total(dMade`k'ormore==1)
					by sid, sort:egen nmiss`k'=total(dMissed`k'ormore==1)
					by sid, sort:egen nhit_hit`k'=total(dMade`k'ormore==1 & make==1)
					by sid, sort:egen nhit_miss`k'=total(dMissed`k'ormore==1& make==1)
				}
		
			*table name dMade4ormore, contents(mean make ) format(%3.2f)
				preserve
					drop if dMade4ormore==0&dMissed4ormore==0 //just compares make3+ vs miss3+
					collapse make nshots, by(sid dMade4ormore)
					sort sid dMade4ormore
					by sid:gen fgc=make[1]
					by sid:gen fgh=make[2]
					gen diff=fgh-fgc
					keep if dMade4ormore==1
					tabstat diff nshots, by(sid) format( %3.2f)
					merge 1:1 sid using biasGVT-4-0.dta
					gen cdiff=diff-bias
					
					mkmat diff cdiff, matrix(Dh4m4)
				restore

				preserve
					drop if dMade3ormore==0&dMissed3ormore==0 //just compares make3+ vs miss3+
					collapse make nshots, by(sid dMade3ormore)
					sort sid dMade3ormore
					by sid:gen fgc=make[1]
					by sid:gen fgh=make[2]
					gen diff=fgh-fgc
					keep if dMade3ormore==1
					tabstat diff nshots, by(sid) format( %3.2f)
					merge 1:1 sid using biasGVT-0-0.dta
					gen cdiff=diff-bias
					
					mkmat diff cdiff, matrix(Dh3m3)
				restore
				
				preserve
					drop if dMade2ormore==0&dMissed2ormore==0 //just compares make2+ vs miss3+
					collapse make nshots, by(sid dMade2ormore)
					sort sid dMade2ormore
					by sid:gen fgc=make[1]
					by sid:gen fgh=make[2]
					gen diff=fgh-fgc
					keep if dMade2ormore==1
					*tabstat diff nshots, by(sid) format( %3.2f)
					merge 1:1 sid using biasGVT-1-0.dta
					gen cdiff=diff-bias
					
					mkmat diff cdiff, matrix(Dh2m2)
				restore
				
				
				
				collapse (mean) make  nhit_hit3 nhit_miss3 nmiss3   nhit3 (sum) dMade2ormore dMissed2ormore dMade3ormore dMissed3ormore (mean) nshots, by(sid)
				svmat Dh4m4
				svmat Dh3m3
				svmat Dh2m2
				
				gen phit_hit3=nhit_hit3/nhit3
				gen phit_miss3=nhit_miss3/nmiss3
				
				tabstat make Dh4m41 Dh4m42 Dh3m31 Dh3m32 Dh2m21 Dh2m22 nshots , by(sid) format( %3.2f) save
				return list
	
		*But here are the results
		*Binomial test on medians works for 3+, but not for 2+ or 4+
		
		*ttest with 26 observations on means is significant on 3+, but not 2+ or 4+
				ttest Dh4m42=0
				ttest Dh3m32=0
				ttest Dh2m22=0
		*Let's be honest, we are being selective. But luckily the straightforward test helps us.
		
		
	*So lets do a player-by-player test, bias corrected
	*e.g
	*Let's make sure our tests ate sensible
	*for the t-test
			clear
			discard
			use gilovichshooters.dta
	
			do genregressionvars_v2.do  //generate the variables that are used in the regression
				genregressionvars	
		ttest make if sid==110&(dMade2ormore==1|dMissed2ormore==1),by(dMade2ormore)
		
		*Let's make sure stata is computing things as I would expect:
		
		*Let's check that the t-value is sensible
		disp (r(mu_1)-r(mu_2))/r(se) 
		*check
		
		*Let's check the one-sided left p-value
		disp t(r(df_t),r(t))
		*check
		
	
		*Let's check the one-sided right p-value 
		disp ttail(r(df_t),r(t))
		*check
		
		*Let's check the confidence interval on the difference:
		disp r(mu_1)-r(mu_2)+r(se)*invttail(r(df_t),.05/2)
		*check
		
		*Okay so we can compute the t value for the bias adjusted difference (here we change the sign of the difference
		local diff2=r(mu_2)-r(mu_1)- (Dh2m2[6,1]-Dh2m2[6,2])
		local t2=(`diff2')/r(se)  //backing out the transformed value, because r(t)=(r(mu_1)-r(mu_2))/sd
		disp "Upper CI = " `diff2'+r(se)*invttail(r(df_t),.05/2)
		disp "Lower CI = " `diff2'-r(se)*invttail(r(df_t),.05/2)
		
		
		disp  "Pvalues = " ttail(r(df_t) ,`t2') //the pvalue
		
		*This should replicate with the Wald (F-test) in the regression framework
		
		reg make dMade2ormore if sid==106&(dMade2ormore==1|dMissed2ormore==1)
		test dMade2ormore=0
		local tstat=sign(_b[dMade2ormore])*sqrt(r(F))
		disp `tstat'
		*Stata works! It gives the same t-stat, good.
		
		*And more directly we can perform a 1-sided t-test using the Wald stat
		test dMade2ormore=Dh2m2[6,1]-Dh2m2[6,2]
		local tstat=sign(_b[dMade2ormore])*sqrt(r(F))
		disp `tstat'
		disp  ttail(r(df_r) ,`tstat')
	
	*2.2 Now we can proceed to performing a player-by-player bias-adjusted ttest
		
		*load bias data for each player from monte-carlo simulation, the file has unfortunate naming
		clear
		use biasGVT-1-0.dta
		mkmat sid bias, matrix(bias2)
		clear
		use biasGVT-0-0.dta
		mkmat sid bias, matrix(bias3)
		clear
		use biasGVT-4-0.dta
		mkmat sid bias, matrix(bias4)
		clear
		clear
			discard
			use gilovichshooters.dta
	
			do genregressionvars_v2.do  //generate the variables that are used in the regression
				genregressionvars	
				
			*THIS IS AN IMPORTANT MODICICATION,
			replace dMade4ormore=0 if period<5
			replace dMade3ormore=0 if period<4
			replace dMade2ormore=0 if period<3
			replace dMissed4ormore=0 if period<5
			replace dMissed3ormore=0 if period<4
			replace dMissed2ormore=0 if period<3
			
				by sid, sort:egen nshots=total(make!=.)
				by sid, sort:egen nhits=total(make==1)
				gen fgp=nhits/nshots
				forvalues k=2/4{
					by sid, sort:egen nhit`k'=total(dMade`k'ormore==1)
					by sid, sort:egen nmiss`k'=total(dMissed`k'ormore==1)
					by sid, sort:egen nhit_hit`k'=total(dMade`k'ormore==1 & make==1)
					by sid, sort:egen nhit_miss`k'=total(dMissed`k'ormore==1& make==1)
				}
			
		levelsof sid, local(sidlist) // create text list of sid values
			local nstreaklengths = 3
			local nvariables = 4
			local ncols=1+ `nstreaklengths'*`nvariables'
			matrix A=J(wordcount("`sidlist'"),`ncols',.)
			local matrixcommand "matrix colnames A = sid "
			disp "`matrixcommand'"
			forvalues k=2(1)4{
				local matrixcommand "`matrixcommand' diff`k' bias`k' se`k' df`k'"
			}
			*disp "`matrixcommand'"
			`matrixcommand'
			*Fill the Matrix
			local it = 0
			
			foreach sid in `sidlist' { //iterate over each sid
			
			
				local ++it  // variable to control row of the matrix
				
				noi di "Sid # `sid'"
				noi di "Row # `it'"
			
				matrix A[`it',1] = `sid'
			
				preserve
					local it2=0
					forvalues j=2(1)4{		
						capture: ttest make if sid==`sid'&(dMade`j'ormore==1|dMissed`j'ormore==1),by(dMade`j'ormore)
							if _rc!=0{
								disp "*************** ERROR ************"
								disp "Sid =`sid'; k = `j'"
								local diff`j'=.
								local bias`j'=.
								local se`j' = .
								local df`j'=.
							}
							else{
								local diff`j'=r(mu_2)-r(mu_1)
								local bias`j'=bias`j'[`it',2]
								local se`j' = r(se)
								local df`j'=r(df_t)
								/*
								local diff`j'=r(mu_2)-r(mu_1)- bias`j'[`it',2]
								local t`j'=(`diff`j'')/r(se) //transform t value
								local p`j' = ttail(r(df_t) ,`t`j'')  //transform z value
								local ui`j'=`diff`j''+r(se)*invttail(r(df_t),.05/2)
								local li`j'=`diff`j''-r(se)*invttail(r(df_t),.05/2)
								
								
								if `t`j''==.{
									disp "*************** CALC ERROR ************"  //for debugging
									disp "Sid =`sid'; k = `j'"
									disp "r(mu_2) = " r(mu_2)
									disp "r(mu_1) = " r(mu_1)
									disp "r(t) = " r(t)
									disp "r(se) = " r(se)
									disp "r(N_1) = " r(N_1)
									disp "r(N_2) = " r(N_2)
									disp "bias`j'[`it',2] = " bias`j'[`it',2]
								
								}
								*/
							}
							
							local ++it2
							matrix A[`it',4*`it2'-2]= `diff`j''
							matrix A[`it',4*`it2'-1]=`bias`j''
							matrix A[`it',4*`it2']=`se`j''
							matrix A[`it',4*`it2'+1]=`df`j''
							
					}

									
				restore

			}
					

						clear
						svmat A, names(col)
				
			
			*Now we can look at the player-by-player tests
			forvalues k=2(1)4{
				cap:gen adiff`k'=diff`k'-bias`k'
				cap:gen at`k'=adiff`k'/se`k'
				cap:gen p`k'=ttail(df`k' ,at`k')
				count if p`k'<.05
				local nsig=r(N)
				count if p`k'<.
				disp "******* Streak k = " `k' " ************"
				disp "Number Shooters Eligible = " r(N)
				disp " Number significant = " `nsig'
				disp "p-value = " 1-binomial(r(N),`nsig'-1,.05)
				disp "*********************************"
			}
			
			*create upper and lower intervals, before and after bias shift
			forvalues k=2(1)4{
				cap:gen ui`k'=diff`k'+se`k'*invttail(df`k',.05/2)
				cap:gen li`k'=diff`k'-se`k'*invttail(df`k',.05/2)
				cap:gen aui`k'=adiff`k'+se`k'*invttail(df`k',.05/2)
				cap:gen ali`k'=adiff`k'-se`k'*invttail(df`k',.05/2)
				cap:gen aui_se`k'=adiff`k'+se`k'
				cap:gen ali_se`k'=adiff`k'-se`k'
			
			}
			
			*Of course the confidence intervals are rather wide, given the small samples:	
			local k=3
			sort diff`k'
			cap drop id`k'
			gen id`k'=_n
			
			summ diff`k'
			local mean_diff=r(mean)
			summ adiff`k'
			local mean_adiff=r(mean)
			
		
			twoway (rcap ali`k' aui`k' id`k', ) (scatter adiff`k' id`k',mlcolor(red) msymbol(circle) mfcolor(red)), ///
				ytitle(Difference) yscale(range(-.7 .7)) yline(0) yline(`mean_adiff', lpattern(dash)  lcolor(red)) ylabel(-.7(.1).7, angle(horizontal)) ///
				xtitle(Shooter) legend(off) scheme(s1mono) name(adiff`k',replace)
			twoway (rcap ali`k' aui`k' id`k', ) (scatter adiff`k' id`k',mlcolor(red) msymbol(circle) mfcolor(red)), ///
				ytitle(Difference) yscale(range(-1.1 1)) yline(0) yline(`mean_adiff', lpattern(dash)  lcolor(red)) ylabel(-1(.5)1, angle(horizontal)) ///
				xtitle(Shooter) legend(off) scheme(s1mono) name(adiff`k',replace)
				
			twoway (rcap li`k' ui`k' id`k', ) (scatter diff`k' id`k', msymbol(circle) mfcolor(gs9)), ///
				ytitle(Difference) yscale(range(-1.1 1)) yline(0) yline(`mean_diff', lpattern(dash)) ylabel(-1(.5)1,  angle(horizontal)) ///
				xtitle(Shooter) legend(off) scheme(s1mono) name(diff`k',replace)
				
						
			twoway (scatter diff`k' id`k',  msymbol(circle) mfcolor(gs9)) (rcap li`k' ui`k' id`k') (scatter adiff`k' id`k', mcolor(red) msymbol(circle)), ///
				ytitle(Difference) yscale(range(-1.1 1)) yline(0) yline(`mean_diff', lpattern(dash)) yline(`mean_adiff', lpattern(dash)  lcolor(red)) ylabel(-1(.5)1,  angle(horizontal)) ///
				xtitle(Shooter) legend(off) scheme(s1mono) name(both`k',replace)
		
				
				twoway (scatter diff`k' id`k') (rcap li`k' ui`k' id`k') (scatter adiff`k' id`k', mcolor(red) msymbol(circle)) (rcap ali`k' aui`k' id`k', lcolor(red)), ///
				ytitle(Difference) yscale(range(-1 1)) yline(0) yline(`mean_diff', lpattern(dash)) yline(`mean_adiff', lpattern(dash)  lcolor(red)) ylabel(-1(.5)1) ///
				xtitle(Shooter) legend(off) scheme(s1mono)
				
				local k=3
			sort diff`k'
			cap drop id`k'
			gen id`k'=_n
			
			summ diff`k'
			local mean_diff=r(mean)
			summ adiff`k'
			local mean_adiff=r(mean)
			
			replace ali`k' = ali`k'*100
			replace aui`k' = aui`k'*100
			replace	 ali_se`k' = ali_se`k'*100
			replace aui_se`k' = aui_se`k'*100
			replace adiff`k' = 100*adiff`k'
			
			local k=3
			
			twoway (scatter adiff`k' id`k') ///
				(rspike ali_se`k' aui_se`k' id`k', lwidth(medthick)) ///
				(rspike ali`k' aui`k' id`k', lwidth(vthin)), ///
				ytitle("Bias-corrected difference" "(percentage points)" ) yline(0, lwidth(thin)) ///
				ylabel(-100(50)100,  angle(horizontal)) ///
				xtitle(Shooter) ///
				legend(order(2 "+/- S.E." 3 "95% CI")) scheme(s1manual)
		

			
			
				
*Step 3: Permutation Test
	*For pooled test, may want to consider allowing certain data to be dropped.
	*There were some general purpose functions, but damn, they didn't seem to work
	*So I try to re-do the code from scratch here.
	
	*Rather than go into all the files that do the permutations and modify
	*Let's just repeat here and document
	
	*Step 1: analagous to whats in anPermute##__#
		clear
		set more off
		program drop _all
		*Instead of using separately defined programs, lets define a program right here
		*Instead of using individual_statistics_avsidfirst_v6.do
		* which defines all sorts of statistics (used in previous Permutation files)
		*lets just stick with the ones GVT used,
		*Define  program to compute:
		*1. difference, z value, whether it is a success
		*2. importantly, for the simulation, we need to handle errors
		*and the simulation needs a value that's not are not emptys
		
		*this program assumes the subject has been dropped
		capture:program drop streakdiff
			program define streakdiff, rclass
			syntax [, kstreak(integer 3)  ]
				*display "****PROGRAM*****"
				preserve
					gen nlast1=.
					by sid session,sort:replace nlast1=make[_n-1] if _n>1
					forvalues i=2(1)4{
						gen nlast`i'=.
						local k=`i'-1
						by sid,sort:replace nlast`i' = nlast`k' + make[_n-`i'] if _n>`i'
					} // end forvalues

				*create indicators for making i or more shots, or exactly i previous shots
					forvalues i=1(1)4{
						gen dMade`i'ormore=.
						by sid session,sort:replace dMade`i'ormore=(nlast`i'==`i') if _n>`i'
						
						gen dMissed`i'ormore=.
						by sid session,sort:replace dMissed`i'ormore=(nlast`i'==0) if _n>`i'
					}
					
					*record hit rate after hit 3+
					summ make if dMade`kstreak'ormore==1
					local nH=r(N)
					if `nH' >0 {
						local pH=r(mean)
					
					}
					*record hit rate after miss 3+
					summ make if dMissed`kstreak'ormore==1
					local nM=r(N)
					if `nM' > 0{
						local pM=r(mean)
					}
					if min(`nH',`nM')==0{
						return scalar success=0
						return scalar diff=-99
						return scalar pH=-99
						return scalar pM=-99
					 }
					 else{
						local diff= `pH'-`pM'
						return scalar success=1
						return scalar diff=`diff'
						return scalar pH=`pH'
						return scalar pM=`pM'
					 }
					*scalar drop pH pM //or else something weird happens in the permutation, I don't know why!
			restore
			*keep sid-exp //restore the data so we can regenerate the streak variables again
		end
		
		*Load in GVT data
		clear
		use gilovichshooters.dta
		drop if make==.	 //missing can create unexpected problems.
	

		
		*let's test to make sure return values make sense
		preserve
			keep if sid==211
			streakdiff, kstreak(3)
			return list
			
			do genregressionvars_v2.do  //generate the variables that are used in the regression
			genregressionvars
			capture:prtest make if dMade2ormore==1|dMissed2ormore==1,by(dMissed2ormore)
			display _rc
			return list
		restore



		local it=1
		timer clear	
		levelsof sid, local(sidlevels)	
		* Use excel =RANDBETWEEN(1,2147483647)
		set seed 855570539
		
		*matrix reps = J(12,1,5)
		*di el(reps,`it',1)
		*"el(reps,`it',1)"
		*matrix reps[2,1]=51500
		forvalues k=2(1)4{
			foreach sid in `sidlevels'{
							timer on `it'
							*local sid=212
							*local k =2
							preserve
								keep if sid==`sid' // it's import to drop, because we only want to permute on player's sequence
								noi: di "sid==`sid'"
								noi: di "streak==`k'"
								local reps=50000  //generate more than we need
								*if `sid' == 102 {
								*	local reps=51500
								*}
								streakdiff,kstreak(`k')
								return list
								if r(success)==1{
									permute make diff=r(diff) pH=r(pH) pM=r(pM) success=r(success) ,saving(shooters/Permute_GVT_diff`k'_Sid`sid'.dta) ///
										reps(`reps'): streakdiff, kstreak(`k') 
									}
							restore
							
							timer off `it'
							timer list
							local ++it

			} // end 	foreach sid in `sidlevels'
		} //end forvalues k=2(1)4
*

		
				
*Let's go through each dataset and
*	 1. cut the data to 10,000 each
*	 2. create standardized measure
		clear
		* Use excel =RANDBETWEEN(1,2147483647)
	
		local min=40000 //the number of permutations to use
			set seed 127605791
			
			matrix A=J(100,2,.) //to record and check the # of successful permutations
			local it=0
		
		forvalues k=2(1)4{
				clear
				use gilovichshooters.dta
			levelsof sid, local(sidlevels)
		
				*calculate just the subjects that are involved.
				levelsof sid, local(sidlevels)			
				foreach sid in `sidlevels' {
						preserve
						keep if sid==`sid'
						streakdiff,kstreak(`k')
						if r(success)==1{ // No error, i.e. the subject had data recorded for him/her
							local sidlevels`k' "`sidlevels`k'' `sid'"
						} //end  if _rc==0{
						else{
						
						} // end else
						restore
				} // end foreach sid in `sidlevels' 
		
		
			foreach sid in `sidlevels`k''{
					clear
					local ++it
					use shooters/Permute_GVT_diff`k'_Sid`sid'.dta
					display "**COUNT**"
					display "`sid' = " `sid'
					display "`k' = " `k'

					count if success==1
					display "count = " r(N)
					matrix A[`it',1]=`sid'
					matrix A[`it',2]=r(N)
					drop if success==0
					gen u1=runiform()
					gen u2=runiform() // in case there are ties with y
					sort u1 u2
					drop u1 u2
					keep if _n<=`min' // sid 204 (and a few other sids) doesn't have 10,000 permutations with 3+ (this will allow us to aggregate)	
					drop success
					*now standardize each stat based on the null mean and standard deviation
						ds 
						local varlist=r(varlist)
					foreach stat of varlist `varlist' {
						if ``stat'[permute]'>-99{						
							summarize `stat'
							gen s`stat' =(`stat'-r(mean))/r(sd) // standardize the realized statisitic from the permututed shots
							local x = (``stat'[permute]'-r(mean))/r(sd) // standardize the realized statisitic from the actual shots
						
							*Add meta-deta for the .dta file, the characteristics, now we have standardized values
							char s`stat'[permute] `x'
							char s`stat'[expression] "r(s`stat')"
							char s`stat'[coleq] "_"
							char s`stat'[colname] "s`stat'"
						}
						else{
							noi di "*** NOT INCLUDING Sid `sid'  Stat `stat'"
								gen s`stat' = 0 // standardize the realized statisitic from the permututed shots							
								char s`stat'[permute] .	//make it empty because it was in error anyway.
								char s`stat'[expression] "r(s`stat')"
								char s`stat'[coleq] "_"
								char s`stat'[colname] "s`stat'"		
						}
					}	// end 						foreach stat of varlist `varlist' {
					
					save temp/tempPermute_GVT_diff`k'_Sid`sid'.dta,replace

			
			} // end foreach sid in `sidlevels`k''{
		} //end 		forvalues k=2(1)4{

*

**Let's Aggregate the data sets
*And create a new permutation file		

		* Use excel =RANDBETWEEN(1,2147483647)
		set seed 1115185510
		local min=40000 //the number of permutations to use
		
		
		forvalues k=2(1)4{
				*which sids actually have esimated values
				clear
				use gilovichshooters.dta

				*Create a list of the subset of sids with data
				levelsof sid, local(sidlevels)			
				foreach sid in `sidlevels' {
						preserve
						keep if sid==`sid'
						streakdiff,kstreak(`k')
						if r(success)==1{ // No error, i.e. the subject had data recorded for him/her
							local sidlevels`k' "`sidlevels`k'' `sid'"
						} //end  if _rc==0{
						else{
						
						} // end else
						restore
				} // end foreach sid in `sidlevels' {
				
		
				clear
				use temp/tempPermute_GVT_diff`k'_Sid101.dta
				*Retrieve variable names
				ds
				local varlist=r(varlist)	//get variable names

				
				foreach stat of local varlist { // initialize sums
					local sum`stat'=0
				}
	
				
				*read in data, and sum the estimated values across subjects, for each variable
				local nsid=0 //count # of subjects for the average
				foreach sid in `sidlevels`k''{

					clear
					local ++nsid
					use temp/tempPermute_GVT_diff`k'_Sid`sid'.dta
					foreach stat of varlist `varlist' {
						local sum`stat' = `sum`stat'' + ``stat'[permute]' 
					}
	
				} // end foreach sid in `sidlevels`k''{
				
					
				*avereage the estimated values (will be meta-data later)
				foreach stat of local varlist  {
					local av`stat'=`sum`stat''/`nsid'
				}
				
				*Just below here we create a data set with `min' rows.  
				*we build each row interatively by setting the first `min' rows, doubling the size by apppending a shooter
				*then we add it together and drop the second half, then repeat
				local i=0	
				foreach sid in `sidlevels`k''{
					display "**Appending**"
					display "`sid' = " `sid'
					display "`k' = " `k'
					local ++i
					if 	`i'==1 { //the first sid's data set we read in
						clear
						use temp/tempPermute_GVT_diff`k'_Sid`sid'.dta
						gen block=1
					}
					else{
						append using temp/tempPermute_GVT_diff`k'_Sid`sid'.dta
						noi di "bring in: temp/tempPermute_GVT_diff`k'_Sid`sid'.dta" 
						replace block=2 if _n>`min'  
						*noi di "line 87"
						gen u1=`min'+runiform()  // make sure we only permute within sid strata 
						gen u2=runiform()
						replace u1=_n if block==1 // don't re-order the block 1
						sort u1 u2
						drop u1 u2
						foreach stat of varlist `varlist'  {
							replace `stat'=`stat' + `stat'[_n+`min']  //add previous sum, to the current sid in block 2
						} 
						keep if _n<=`min'  //drop block 2
					} // end else
				} //end 	foreach sid in `sidlevels`k''{	
				
				*Fix the Meta-Data
				foreach stat of varlist `varlist' {
			
					noi di "replace `stat'=`stat'/`nsid'"
					replace `stat'=`stat'/`nsid'
					char `stat'[permute] `av`stat''
				}
				drop block
				save temp/tempPermute_GVT_diff`k',replace

		
} //end forvalues k=2(1)4{

*
*pooled permutation tests
*They are not symmetric, and its a one-sided hypothesis
permute using "temp/tempPermute_GVT_diff2.dta", right
permute using "temp/tempPermute_GVT_diff3.dta", right
permute using "temp/tempPermute_GVT_diff4.dta", right


*Player by player permutation tests		
	local i=0
		for		matrix drop _all	
	values k=2(1)4{
			clear
			use gilovichshooters.dta
			*Create a list of the subset of sids with data
			levelsof sid, local(sidlevels)			
			foreach sid in `sidlevels' {
					preserve
					keep if sid==`sid'
					streakdiff,kstreak(`k')
					if r(success)==1{ // No error, i.e. the subject had data recorded for him/her
						local sidlevels`k' "`sidlevels`k'' `sid'"
					} //end  if _rc==0{
					else{
					
					} // end else
					restore
			} // end foreach sid in `sidlevels' {
			clear
		
			foreach sid in `sidlevels`k''{		
				local ++i
				display "*********** ITERATION # `i' *************"
				*Note: the histogram isn't symmetric, and it is a one-sided hypothesis
				permute diff using "temp/tempPermute_GVT_diff`k'_Sid`sid'.dta", right
				foreach  vec in `k' `sid' r(b) r(p) r(ci) { //8 potential columns
					*di "group=`group'"	
					matrix pvalues`k'_`sid' = (nullmat(pvalues`k'_`sid'), `vec'')  // nullmat(pvalues) allows us to start with an empty matrix.
				} // end foreach  vec in `k' `sid' r(b) r(p) r(ci) 
				if `i'==1{
					matrix pvalues=pvalues`k'_`sid' 
				}
				else {
					matrix pvalues=pvalues\pvalues`k'_`sid' 
				} // end if `i'==1{
				
			} // end 	foreach sid in `sidlevels`k''{	
		} // end 	forvalues k=2(1)4{
				
				matrix colnames pvalues = k sid t p l u //now we label the column names
					*Check if matrix is build correctly
					noi matrix list pvalues
					
					clear 
					svmat pvalues, names(col)
					
					count if k==2&p<.
					count if k==2&p<.05
					
					count if k==3&p<.
					count if k==3&p<.05
					
					count if k==4&p<.
					count if k==4&p<.05
					
*------------------------------

* A bias-correct CI check:
*With a test to make sure the intervals have the correct coverage properties
cap program drop testCI
program define testCI, rclass
	syntax [, fg(string) k(integer 3) dgp(string) ]
clear
	if "`dgp'"=="feedback"{
		genshots_posfeedback, dfg(.15)
	}
	else{
		genshots_GVT, fg(`fg')
	}
	genDiffCI, k(`k')
	return scalar ci_rad=r(ci_rad)
	return scalar diff=r(diff)
end

cap program drop genDiffCI
program define genDiffCI, rclass
	syntax [, k(integer 3) confidence(real .975) ]
	preserve
			*build command to make the streak variables
		*local k=3
		local make="make[_n-1]==1"
		local miss="make[_n-1]==0"
		forvalues i=2(1)`k'{
			local make="`make' & make[_n-`i']==1"
			local miss="`miss' & make[_n-`i']==0"
		}
		*disp "`make' "  "||| `miss'"
		cap gen dMade`k'ormore=`make' if period>`k'
		cap gen dMissed`k'ormore=`miss' if period>`k'
		
		scalar drop _all
		scalar tvar=0
		local nplayers=0
		local dropsids ""
		levelsof sid, local(sidlist) // create text list of sid values		
			foreach sid in `sidlist' { //iterate over each sid
				cap: ttest make if sid==`sid' &(dMade`k'ormore==1| dMissed`k'ormore==1),by(dMade`k'ormore)
				if _rc==0 & r(se)!=.{
						disp "diff= " -r(t)*r(se) ", std. error= " r(se) 
						local nplayers=`nplayers'+1 // computing the total # of players
						scalar tvar=scalar(tvar)+ r(se)^2	// computing the total variance		
				}
				else{
					drop if sid==`sid' // drop players with no variance
				}
			}
			disp "Total Variance " tvar
		scalar var=(1/`nplayers')^2*scalar(tvar)  //compute the variance of the average difference acorros players
		
		return scalar ci_rad=invnormal((`confidence'+1)/2)*(scalar(var))^(1/2)
		return scalar se=(scalar(var))^(1/2)
		drop if dMade`k'ormore==0& dMissed`k'ormore==0 | dMade`k'ormore==. | dMissed`k'ormore==.
		collapse make, by(sid dMade`k'ormore dMissed`k'ormore)
		collapse make, by(dMade`k'ormore dMissed`k'ormore)
		return scalar diff=make[2]-make[1]
		restore
end

	
testCI, fg(design) k(3)
return list
simulate diff=r(diff) ci_rad=r(ci_rad),reps(1000): testCI, fg(design) k(3)
simulate diff=r(diff) ,reps(1000): testCI, fg(design) k(3)

simulate diff=r(diff) ci_rad=r(ci_rad),reps(10000): testCI, dgp(feedback)

*Now do the adjustment for real
		*the mean is .1259292  from taking the average of adiff3 in the code above
		* or the average in Dh3m32
		*Let's generage the width of the CIs
		clear
		use gilovichshooters.dta
		genDiffCI, confidence(.95)
		return list
		*z-test
		
		*.0913594527
		
		*k=2 the mean is .0537492  by taking the average of adiff2
			clear
		use gilovichshooters.dta
		genDiffCI, k(2) confidence(.95)
		return list
		disp 1-normal(.0537492/r(se))
	
		*k=4 the mean is .102 by taking the average of adiff2
			clear
		use gilovichshooters.dta
		genDiffCI, k(4) confidence(.95)
		return list
		disp 1-normal(.102/r(se))
		*
		clear
		use corunashooters_all.dta
		drop if dDoubleSession==1|dPanel8==0
		genDiffCI, k(3) confidence(.95)
		return list
		
		clear
		use jagsdata.dta
		genDiffCI, k(3) confidence(.95)
		return list
		
		clear
		discard
		 use 3pt_v12.dta
		 keep if nshots>99
		 	genDiffCI, k(3) confidence(.95)
		return list

cap program drop genshots_posfeedback
program define genshots_posfeedback
	syntax [, nshots(integer 100) nsids(integer 26) fgn(real .5) dfg(real .10) ]
clear
	local obs=`nshots'*`nsids'
	set obs `obs'
	gen rand1=runiform()
	gen period=mod(_n-1,`nshots')+1
	gen sid=100+floor((_n-period)/`nshots')+1
	scalar fgn=`fgn'
	scalar dfg=`dfg'
	scalar fgh=scalar(fgn)+scalar(dfg)/2
	scalar fgc=scalar(fgn)-scalar(dfg)/2
	gen make=(rand1<=scalar(fgn)) if period<4
	replace make =(rand1<=scalar(fgc)&make[_n-1]==0&make[_n-2]==0&make[_n-3]==0| rand1<=scalar(fgn)& !((make[_n-1]==1&make[_n-2]==1&make[_n-3]==1)|(make[_n-1]==0&make[_n-2]==0&make[_n-3]==0)) |rand1<=scalar(fgh)&make[_n-1]==1&make[_n-2]==1&make[_n-3]==1)

end

* A bias-corrected pooled test, in a fixed-effects context

*First, lets see how biased a it is given the observed field goal percentages
cap program drop genshots_GVT
program define genshots_GVT
	syntax [, fg(string) nshots(integer 100) nsids(integer 26)]
	*This code regenerates GVT's shot outcomes, randomly 
	clear
	local obs=`nshots'*`nsids'
	set obs `obs'
	gen rand1=runiform()
	gen period=mod(_n-1,`nshots')+1
	gen sid=100+floor((_n-period)/`nshots')+1
	replace sid=200+sid-100-14 if sid>114
	
	*These shooters had empty observations, because they ended early
	drop if period>90 & sid==104
	drop if period>75& sid==107
	drop if period>50 & sid==108
	
	gen fg=.
	if "`fg'"=="design"{
		replace fg=.5
	}
	else if "`fg'"=="actual"{
		replace fg=.54 if sid==101
		replace fg=.35 if sid==102
		replace fg=.6 if sid==103
		replace fg=.4 if sid==104
		replace fg=.42 if sid==105
		replace fg=.57 if sid==106
		replace fg=.56 if sid==107
		replace fg=.5 if sid==108
		replace fg=.54 if sid==109
		replace fg=.6 if sid==110
		replace fg=.58 if sid==111
		replace fg=.44 if sid==112
		replace fg=.61 if sid==113
		replace fg=.59 if sid==114
		replace fg=.48 if sid==201
		replace fg=.34 if sid==202
		replace fg=.39 if sid==203
		replace fg=.32 if sid==204
		replace fg=.36 if sid==205
		replace fg=.46 if sid==206
		replace fg=.41 if sid==207
		replace fg=.53 if sid==208
		replace fg=.45 if sid==209
		replace fg=.46 if sid==210
		replace fg=.53 if sid==211
		replace fg=.25 if sid==212
	}
	else{
		display "ERROR: control flow issues"
		error
		exit
	}
	
	gen make=(rand1<=fg)
end

cap program drop anGVT
program define anGVT, rclass
	syntax [, fg(string) k(integer 3) command(string) ]
	*This code regenerates GVT's shot outcomes, randomly 
	*Then runs a fixed effect regression (biased), and returns the coefficient estimate 
	*hitting k or more in row
	clear
	genshots_GVT, fg(`fg')
	*Flexible, to add any command
	`command'
	local streak=3
	local make="make[_n-1]==1"
	local miss="make[_n-1]==0"
	forvalues i=2(1)`k'{
		local make="`make' & make[_n-`i']==1"
		local miss="`miss' & make[_n-`i']==0"
	}
	*disp "`make' "  "||| `miss'"
	gen dMade`k'ormore=`make' if period>`k'
	gen dMissed`k'ormore=`miss' if period>`k'
	
	reg make dMade`k'ormore i.sid if dMade`k'ormore==1| dMissed`k'ormore==1
	return scalar diff=el(r(table),1,1)
	/*
	drop if dMade3ormore==0& dMissed3ormore==0 | dMade3ormore==. | dMissed3ormore==.
	collapse make, by(dMade3ormore dMissed3ormore)
	return scalar diff=make[2]-make[1]
	*/
end

*reg make dMade3ormore i.sid if dMade3

*bias is =-.0401 for k=3; -.0211 for k=2; -.07647 for k==4
simulate diff=r(diff),reps(10000): anGVT, fg(actual) k(4)


*Calculate the bias
*For sid 109 only, k=3,4
simulate diff=r(diff),reps(10000) saving(MS_Bias_only109_k4.dta): anGVT, fg(design) k(4) command("keep if sid == 109") 
simulate diff=r(diff),reps(10000) saving(MS_Bias_without109_k3.dta): anGVT, fg(design) k(3) command("drop if sid == 109") 

*For GVT without sid 109, k=3,4
simulate diff=r(diff),reps(10000) saving(MS_Bias_without109_k3.dta): anGVT, fg(design) k(4) command("drop if sid == 109")

*bias is -=.042 for k=3; -.02 for k=2; -.0823 for k==4
simulate diff=r(diff),reps(10000): anGVT, fg(design) k(3) command("drop if sid == 109")
summarize, meanonly
local bias =  r(mean) 
clear
use gilovichshooters.dta
do genregressionvars_v2.do  //generate the variables that are used in the regression
	genregressionvars	
drop if sid ==109
xtset sid period 
xtreg make dMade3ormore if dMade3ormore==1| dMissed3ormore==1, fe
	test dMade3ormore=-.0429150
			local tstat=sign(_b[dMade3ormore])*sqrt(r(F))
			disp `tstat'
			disp  ttail(r(df_r) ,`tstat')	




simulate diff=r(diff),reps(10000): anGVT, fg(design) k(4)

clear
discard
use gilovichshooters.dta
do genregressionvars_v2.do  //generate the variables that are used in the regression
	genregressionvars		
	tsset sid period
			
global dummys=""
levelsof sid, local(sidlist) // create text list of sid values
foreach sid in `sidlist' { //iterate over each sid
	gen s`sid'=sid==`sid'
	if `sid' !=212{  // otherwise control variables are linearly dependent
	global dummys="$dummys s`sid'"
	}
}	
disp "$dummys"

*Pooled with k==3
reg make dMade3ormore if dMade3ormore==1| dMissed3ormore==1
reg make dMade3ormore i.sid if dMade3ormore==1| dMissed3ormore==1
reg make dMade3ormore $dummys if dMade3ormore==1| dMissed3ormore==1, robust
brl make  dMade3ormore  $dummys  if dMade3ormore==1| dMissed3ormore==1, cluster(sid)

*hypothesis test
	test dMade3ormore=-.04
			local tstat=sign(_b[dMade3ormore])*sqrt(r(F))
			disp `tstat'
			disp  ttail(r(df_r) ,`tstat')	
			
*Pooled with k==2
reg make dMade2ormore i.sid if dMade2ormore==1| dMissed2ormore==1
reg make dMade2ormore $dummys if dMade2ormore==1| dMissed2ormore==1, robust
brl make  dMade2ormore  $dummys  if dMade2ormore==1| dMissed2ormore==1, cluster(sid)

*hypothesis test
	test dMade2ormore=-.021
			local tstat=sign(_b[dMade2ormore])*sqrt(r(F))
			disp `tstat'
			disp  ttail(r(df_r) ,`tstat')	

*Pooled with k==4
reg make dMade4ormore i.sid if dMade4ormore==1| dMissed4ormore==1
reg make dMade4ormore $dummys if dMade4ormore==1| dMissed4ormore==1, robust
brl make  dMade4ormore  $dummys  if dMade4ormore==1| dMissed4ormore==1, cluster(sid)

*hypothesis test
	test dMade4ormore=-.076
			local tstat=sign(_b[dMade4ormore])*sqrt(r(F))
			disp `tstat'
			disp  ttail(r(df_r) ,`tstat')	
			
			
* Need to do the same for 3 point data.

			
clear
discard
 use 3pt_v12.dta
 keep if nshots>99
do genregressionvars_v2.do  //generate the variables that are used in the regression
	genregressionvars		

			
global dummys=""
levelsof sid, local(sidlist) // create text list of sid values
foreach sid in `sidlist' { //iterate over each sid
	gen s`sid'=sid==`sid'
	if `sid' !=101{  // otherwise control variables are linearly dependent
	global dummys="$dummys s`sid'"
	}
}	
disp "$dummys"

*Pooled with k==3
reg make dMade3ormore i.sid if dMade3ormore==1| dMissed3ormore==1, cluster(sid)
reg make dMade3ormore $dummys if dMade3ormore==1| dMissed3ormore==1, robust
brl make  dMade3ormore  $dummys  if dMade3ormore==1| dMissed3ormore==1, cluster(sid)
	 
			
			
*----------------------------- Analysis of beliefs

*Load regimeswitch further below
*this shows the long run correlation
	cap program drop predictm
	program define predictm, rclass
	genshots_regimeswitch, staynormal(.96) fg(.5) frachot(.15) dfg(.10)
	gen pmake=state==1
	corr pmake make
	return scalar rho=r(rho)
	end

	simulate rho=r(rho), reps(10000): predictm

	simulate rho=r(rho), reps(10000): predictm2

	cap program drop predictm2
	program define predictm2, rclass
	clear
	set obs 100
	gen state=_n<16
	gen rand1=runiform()
	gen make=(rand1<=.55)*(state==1)+(rand1<=.45)*(state==0)

	gen pmake=state==1
	corr pmake make
	return scalar rho=r(rho)
	end



clear
use "C:\Users\Miller\Dropbox\josh\work\projects\Adam\HotHand\Data\gilovichshooters.dta" 
gen betm=.
gen beth=.
by sid, sort: replace beth=1 if make[_n-1]==1&make[_n-2]==1 &make[_n-3]==1
by sid, sort: replace betm=0 if make[_n-1]==0&make[_n-2]==0 &make[_n-3]==0
by sid, sort:egen nbeth=total(beth==1)
by sid, sort:egen nbetm=total(betm==0)
by sid, sort:egen fgh=mean(make)
gen fgm=1-fgh

by sid, sort: gen winh=beth==make if beth<.
by sid, sort: gen winm=betm==make if betm<.

by sid, sort:egen mwinm=mean(winm)
by sid, sort:egen mwinh=mean(winh)
replace mwinm=0 if mwinm==.
replace mwinh=0 if mwinh==.

gen wsuccess=(mwinh-fgh)*nbeth + (mwinm-fgm)*nbetm
summ wsuccess if period==1

gen totalbets=nbeth+nbetm

keep if period==1
keep sid nbeth nbetm

mkmat sid nbeth nbetm, matrix(X)
matrix list X


*Put down the bets at random in the data set.
clear
use "C:\Users\Miller\Dropbox\josh\work\projects\Adam\HotHand\Data\gilovichshooters.dta" 
gen betm=.
gen beth=.
forvalues i=1(1)26{
	replace beth=1 if sid==X[`i',1] & period<=X[`i',2]
	replace betm=0 if sid==X[`i',1] & period>X[`i',2]& period<=X[`i',2]+X[`i',3]
}

gen winm=.
gen winh=.

*assuming beth and betm exist, this function calculates winnings 
capture:program drop winnings
			program define winnings, rclass

		gen winm=.
		gen winh=.
		by sid, sort: replace winh=beth==make if beth<.
		by sid, sort: replace winm=betm==make if betm<.

		by sid, sort:egen mwinm=mean(winm)
		by sid, sort:egen mwinh=mean(winh)
		replace mwinm=0 if mwinm==.
		replace mwinh=0 if mwinh==.

		gen wsuccess=(mwinh-fgh)*nbeth + (mwinm-fgm)*nbetm
		summ wsuccess if period==1

		return scalar profit=r(mean)*r(N)
		
		*drop the variables
		drop winh winm mwinm mwinh wsuccess
end


clear
use "C:\Users\Miller\Dropbox\josh\work\projects\Adam\HotHand\Data\gilovichshooters.dta" 
gen betm=.
gen beth=.
by sid, sort: replace beth=1 if make[_n-1]==1&make[_n-2]==1 &make[_n-3]==1
by sid, sort: replace betm=0 if make[_n-1]==0&make[_n-2]==0 &make[_n-3]==0
by sid, sort:egen nbeth=total(beth==1)
by sid, sort:egen nbetm=total(betm==0)
by sid, sort:egen fgh=mean(make)
gen fgm=1-fgh

permute make profit=r(profit), strata(sid) reps(1000) : winnings

*If they know a player's overall fg%
clear
use "C:\Users\Miller\Dropbox\josh\work\projects\Adam\HotHand\Data\gilovichshooters.dta" 
by sid, sort:egen fgh=mean(make)
gen fgm=1-fgh
gen betm=.
gen beth=.
gen dHit3=make[_n-1]==1&make[_n-2]==1 &make[_n-3]==1
gen dMiss3=make[_n-1]==0&make[_n-2]==0 &make[_n-3]==0
by sid, sort: replace beth=1 if dHit3==1
by sid, sort: replace betm=0 if dMiss3==1
by sid, sort:replace beth=1 if beth==. &  dMiss3==0& dHit3==0 & fgh>.5
by sid, sort:replace betm=0 if betm==. &  dMiss3==0& dHit3==0 & fgh<.5
*by sid, sort:replace betm=0 if betm==. &  dMiss3==0& dHit3==0 
by sid, sort:egen nbeth=total(beth==1)
by sid, sort:egen nbetm=total(betm==0)

 winnings
gen predict=.
replace predict=beth if beth!=.
replace predict=betm if betm!=.

gen rho=.
by sid, sort: gen backpredict=predict[_n+1]

gen backrho=.
forvalues i=1(1)26{
	corr predict make if sid==X[`i',1]
	replace rho=r(rho) if sid==X[`i',1]
	
	corr backpredict make if sid==X[`i',1]
	replace backrho=r(rho) if sid==X[`i',1]
}
permute make profit=r(profit), strata(sid) reps(1000) right: winnings


*------------ POWER

cap program drop genshots_3boost
program define genshots_3boost
	syntax [, nshots(integer 100) nsids(integer 26) fgpercent(real .5) dfgpercent(real .10) ]
clear
	local obs=`nshots'*`nsids'
	set obs `obs'
	gen rand1=runiform()
	scalar fg=`fgpercent'
	scalar dfg=`dfgpercent'
	scalar fghot=fg+dfg/2
	scalar fgcold=fg-dfg/2
	gen period=mod(_n,`nshots')
	gen sid=floor((_n-period)/`nshots')+1
	
	gen make=(rand1<=scalar(fg)) if mod(period,`nshots')<4
	replace make =(rand1<=scalar(fg)& !(make[_n-1]==make[_n-2]& make[_n-2]==make[_n-3]) |rand1<=scalar(fghot)&make[_n-1]==1&make[_n-2]==1&make[_n-3]==1| rand1<=scalar(fgcold)&make[_n-1]==0&make[_n-2]==0&make[_n-3]==0) if period>3

	gen dMade3ormore2=(make[_n-1]==1&make[_n-2]==1&make[_n-3]==1) if period>3
	gen dMissed3ormore2=(make[_n-1]==0&make[_n-2]==0&make[_n-3]==0) if period>3
end

genshots_3boost, nshots(100) nsids(26) fgpercent(.5) dfgpercent(.13)



cap program drop bintest
program define bintest, rclass
	syntax [,  biasadj(real .0794) ]
		*THIS IS A Simplied version of what is in PowerAnalysis.do
			replace dMade3ormore=0 if period<4
			replace dMissed3ormore=0 if period<4

		
		levelsof sid, local(sidlist) // create text list of sid values
			local ntests = 0
			local nreject =0
			
			foreach sid in `sidlist' { //iterate over each sid
			
			
				 // variable to control row of the matrix
				capture: ttest make if sid==`sid'&(dMade3ormore==1|dMissed3ormore==1),by(dMade3ormore)
							if _rc!=0{
					
							}
							else{
								local ++ntests
								local t_adj=(r(mu_2)-r(mu_1)+`biasadj')/r(se)
								local pvalue=ttail(r(df_t),`t_adj')
								local reject=cond(`pvalue'<=.05,1,0)
								local nreject=`nreject'+ `reject'
				
							}
			}
			disp "Number tests= `ntests', Number reject= `nreject'"
			
			local binpval= 1-binomial(`ntests',`nreject'-1,.05)
			return scalar pvalue=`binpval'
			return scalar reject=cond(`binpval'<=.05,1,0)
end


genshots_3boost, nshots(100) nsids(26) fgpercent(.5) dfgpercent(.13)

	
cap program drop power_experiment
program define power_experiment, rclass
syntax [, nshots(integer 100) nsids(integer 26) fgpercent(real .5) dfgpercent(real .10) biasadj(real .0794) ]

		genshots_3boost, nshots(`nshots') nsids(`nsids') fgpercent(`fgpercent') dfgpercent(`dfgpercent') 
		bintest, biasadj(`biasadj')
			return scalar pvalue=r(pvalue)
			return scalar reject=r(reject)
end

power_experiment, nshots(100) nsids(26) fgpercent(.5) dfgpercent(.13)  biasadj(.0794)
return list

set more off
simulate pvalue=r(pvalue) reject=r(reject), saving("C:\Users\Miller\Dropbox\josh\work\projects\Adam\HotHand\Data\statTheory\temp.dta") reps(10000) seed(130937061): power_experiment, nshots(100) nsids(26) fgpercent(.5) dfgpercent(.13)  biasadj(.0794)
			
			
* It has moderate power, .75
					
					
* Bias
cap program drop steph_curry
program define steph_curry, rclass
syntax [, nshots(integer 100) nsids(integer 1) fgpercent(real .5) dfgpercent(real .10)  ]

		genshots_3boost, nshots(`nshots') nsids(`nsids') fgpercent(`fgpercent') dfgpercent(`dfgpercent') 
		capture: ttest make if (dMade3ormore==1|dMissed3ormore==1),by(dMade3ormore)
			return scalar diff=r(mu_2)-r(mu_1)

end
steph_curry,  nshots(100) nsids(1) fgpercent(.5) dfgpercent(.10) 
simulate diff=r(diff), reps(2000): steph_curry,  nshots(100) nsids(1) fgpercent(.5) dfgpercent(.08) 


***************************
****** Showing the bias correction is conservative

* Note: for the regime shift the bias depends on the stationary distribution and 
*		the difference in probrability (transition probabilities in the hot state betwen .8 and .1 don't matter much)

* Note: for the 3 boost, 
****************************
clear
*Lets summarize what we have found:
use "C:\Users\Miller\Dropbox\josh\work\projects\Adam\HotHand\Data\BiasRobustSimulation.dta" 


*Notice that the bias in the regime shift model depends principally on the truediff, and the underlying fg% (fraction of time in hotstate 10-20 & transistion probabilities 80-100 p hot don't influence it much)
twoway (line bias fg), by(dgp dfg stayn frachot)


collapse diff bias, by(dgp fg dfg)
sort dgp dfg fg bias diff
order dgp dfg fg bias diff
export excel using "statTheory/BiasRobustSimulation.xlsx",replace // export to matlab


mata
	function probstayhot(probstaynormal,frachot)
	{
		p = 1- ( (1-frachot)/frachot )* (1-probstaynormal)
		return(p)
	}
end

mata: probstayhot(.99,.10)
mata: probstayhot(.98,.10)

mata: probstayhot(.99,.15)
mata: probstayhot(.97,.15)

mata: probstayhot(.98,.20)
mata: probstayhot(.95,.20)


mata
	frachot=.10
	staynormal=.97
	stayhot=probstayhot(staynormal,frachot)
	phit=frachot*fgh+(1-frachot)*fgn
	.5=frachot*(fgh-fgn)+fgn
	.5=frachot*d+fgn
	fgn=.5-frachot*d
	fgh=fgn+d
end

*Averge Waiting time until hot, by neg binomail
di .97/.03


*Load individual stats programs:

do individual_statistics_v6.do


cap program drop genshots_regimeswitch
program define genshots_regimeswitch
	syntax [, staynormal(real .95) frachot(real .2) nshots(integer 100) fg(real .5) dfg(real .10) ]
	*Note: If we want to keep the long run frequency of the hot hand state constant
	* 	and increase the likliness of staying in the hot state, then we have to make 
	*   the chance of trainsitioning into the hot state lower, we have to wait
	* 	longer to see a hot state nbinomial(1, 70,.01)
	* For the tranistion matrix P, let p=Pr(0|0), q=Pr(1|1)
	* The stationary pi distrubution statisfied pi0/pi1=(1-q)/(1-p)
	scalar drop _all
	scalar fg=`fg'  //Overall probability of a hit
	scalar dfg=`dfg'
	scalar pi=`frachot'
	scalar fgn=scalar(fg)-scalar(pi)*scalar(dfg) // normal state probability chosen so that they average to .5
	scalar fgh=scalar(fgn)+scalar(dfg)
	display "fgh = " fgh ", fgn = " fgn ", fg= " fg ", frachot= " pi
	assert scalar(fgn)<=1& scalar(fgn)>=0
	assert scalar(fgh)<=1& scalar(fgh)>=0
	scalar p=`staynormal'
	*from the stationary distrubtion formula pi=pi*P, we calculate the probability of staying hot
	scalar q=1-(1-scalar(p))*(1-scalar(pi))/scalar(pi) // q= 1-(1-p)*ODDSAGAINSTBEINGHOT
	
	clear
	set obs `nshots'	
	gen rand1=runiform()
	gen period=_n
	gen state=(rand1<=scalar(pi))  if _n==1  //begin according to stationary distribution
	replace state=(rand1>scalar(p) & state[_n-1]==0)+ (rand1<=scalar(q) & state[_n-1]==1)
	*sum state

	gen rand2=runiform()
	gen make=(rand2<=scalar(fgn)&state==0| rand2<=scalar(fgh) & state==1)
	gen sid=1
	gen session=1

	gen dMade3ormore2=(make[_n-1]==1&make[_n-2]==1&make[_n-3]==1) if _n>3
	gen dMissed3ormore2=(make[_n-1]==0&make[_n-2]==0&make[_n-3]==0) if _n>3
end
genshots_regimeswitch, staynormal(.96) frachot(.15) dfg(.10)

* with a two-state model you are almost certainly in state 0 if you miss 3 or more in a row
*This makes an interesting point, the conditional probability test has a trade-off in markov model
*if you restrict to Made3+ vs. Miss3+ you loss data, but on the other hand if 
*you look at Miss3+, you are almost certainly not in a hot state.
*it is not lose-lose, there is a trade off
*Clearly the hot state is less likely when you haven't made 3+
genshots_regimeswitch, staynormal(.95) frachot(.20) nshots(100) fg(.5) dfg(.40) 

tabstat state, by(dMade3ormore)
*but it is even less likely if you missed 3 or more
summ state if dMissed3ormore==1


*--------------------------------------------------
****	Positive Feedback Model
*
* Dumb 3+ Positive Feedback Model
* Everytime you hit 3+, your hit rate improves (could add with prob)
*You don't need a power simulation to see that it is underpowered, it is a pure issue of sample size
cap program drop genshots_3boost
program define genshots_3boost
	syntax [, nshots(integer 100) fgn(real .5) dfg(real .10) ]
clear
	set obs `nshots'
	gen rand1=runiform()
	scalar fgn=`fgn'
	scalar dfg=`dfg'
	scalar fghot=scalar(fgn)+scalar(dfg)
	gen period=_n
	gen make=(rand1<=scalar(fgn)) if period<4
	replace make =(rand1<=scalar(fgn)& !(make[_n-1]==1&make[_n-2]==1&make[_n-3]==1) |rand1<=scalar(fghot)&make[_n-1]==1&make[_n-2]==1&make[_n-3]==1)

	gen sid=1
	gen session=1

	gen dMade3ormore2=(make[_n-1]==1&make[_n-2]==1&make[_n-3]==1) if _n>3
	gen dMissed3ormore2=(make[_n-1]==0&make[_n-2]==0&make[_n-3]==0) if _n>3
end


*We want a fixed overall fg percentage

*Average to 60
*.58  d=.1
*.55  d=.2
*.515 d=.3
*.47  d=.4

*Average to 55
*.53  d=.1
*.51  d=.2
*.48  d=.3
*.445  d=.4

*Average to 50
*.49  d=.1
*.47  d=.2
*.45  d=.3
*.42  d=.4

*Average to 45
*.44  d=.1
*.43  d=.2
*.41  d=.3
*.39  d=.4

*Average to 40
*.395  d=.1
*.385  d=.2
*.37  d=.3
*.36  d=.4

genshots_3boost, nshots(50000) fgn(.58) dfg(.1)
summ make


	
*Begin program
cap program drop calc_diff
program define calc_diff, rclass
	syntax , command(string) method(string) [nshots(integer 300)]
	`command'
		disp "line 278"
	count if dMade3ormore ==1
	local NH3=r(N)
	count if dMissed3ormore ==1
	local NM3=r(N)
	if (`NH3'>=1 & `NM3'>=1){
	
			if "`method'"=="collapse"{
			*preserve
				keep if dMade3ormore ==1| dMissed3ormore==1
				collapse make,by(dMade3ormore dMissed3ormore )
				disp make[2]-make[1]
				return scalar diff=make[2]-make[1]
			*	restore
			
			}
			else{
	
			prtest make if dMade3ormore ==1| dMissed3ormore==1,by(dMade3ormore) // with a two-state model you are almost certainly in state 0 if you miss 3 or more in a row
			return scalar diff=r(P_2)-r(P_1)
			}
			/*
	
			*/
			
	}
	else{
			return scalar diff=.
	}
	*	
end

calc_diff, command(genshots_regimeswitch, staynormal(.95) frachot(.20) nshots(100) fg(.5) dfg(.40) )	method(collapse)

*make sure this runs fast
cap program drop testime
program define testime,rclass
	syntax, method(string)
cap: timer clear 1
timer on 1
calc_diff, command(genshots_regimeswitch, staynormal(.95) frachot(.20) nshots(100) fg(.5) dfg(.40) )	method(`method')	
timer off 1
timer list
return scalar time=r(t1)
end

simulate time=r(time), reps(1000): testime, method(a)

*test the simulation, to make sure it records well
set seed 1234
simulate diff=r(diff), reps(1000): calc_diff, command(genshots_regimeswitch, staynormal(.95) frachot(.20) nshots(100) fg(.5) dfg(.40) )	method(collapse)

*****************************************
*Regime Switching Model:
*****************************************

local fgs = "40 50 60"
local dfgs = "10 20 30 40"
local frachots= "10 15 20"

*for some reason
* `fracthot' == 10 didn't work
* "`fracthot'" == "10" din't work, and all variations, only the following worked
foreach fg of local fgs{	
	foreach dfg of local dfgs{
			foreach frachot of local frachots{
					display `frachot'
					if (`frachot' == 10) {
						local stayns = "98 99" // We can't go too low or probstayhot is to low, this give .82 and .91 respectively
					}
					else if (`frachot'== 15){
						local stayns = "97 99" // We can't go too low or probstayhot is to low, this give .83 and .94 respectively
					
					}
					else if (`frachot'== 20) {
						local stayns = "95 98" // We can't go too low or probstayhot is to low, this give .8 and .92 respectively
					}
					else {
						display "ERROR: Control flow issues"
						exit
					
					}
					foreach stayn of local stayns{
						disp "BiasAdjRegimeShift_fg`fg'-dfg`dfg'-stayn`stayn'-frachot`frachot'.dta"
						simulate  diff=r(diff) ///
						, reps(10000) saving(BiasAdjRegimeShift_fg`fg'-dfg`dfg'-stayn`stayn'-frachot`frachot'.dta, replace) : calc_diff, command(genshots_regimeswitch, staynormal(.`stayn') frachot(.`frachot') nshots(100) fg(.`fg') dfg(.`dfg') )	method(collapse)																											
					}
			}
	}
}

local fgs = "40 50 60"
local dfgs = "10 20 30 40"
local frachots= "10 15 20"

*for some reason
* `fracthot' == 10 didn't work
* "`fracthot'" == "10" din't work, and all variations, only the following worked
clear
local i=0
foreach fg of local fgs{	
	foreach dfg of local dfgs{
			foreach frachot of local frachots{
					display `frachot'
					if (`frachot' == 10) {
						local stayns = "98 99" // We can't go too low or probstayhot is to low, this give .82 and .91 respectively
					}
					else if (`frachot'== 15){
						local stayns = "97 99" // We can't go too low or probstayhot is to low, this give .83 and .94 respectively
					
					}
					else if (`frachot'== 20) {
						local stayns = "95 98" // We can't go too low or probstayhot is to low, this give .8 and .92 respectively
					}
					else {
						display "ERROR: Control flow issues"
						exit
					
					}
					foreach stayn of local stayns{
						disp "BiasAdjRegimeShift_fg`fg'-dfg`dfg'-stayn`stayn'-frachot`frachot'.dta"
						local i=`i'+1
						clear 
						use BiasAdjRegimeShift_fg`fg'-dfg`dfg'-stayn`stayn'-frachot`frachot'.dta
						collapse diff
						gen simid=`i'
						gen fg=`fg'
						gen dfg=`dfg'
						gen stayn=`stayn'
						gen frachot=`frachot'
						cap: append using BiasAdjRegimeShift
						save BiasAdjRegimeShift, replace
						*erase BiasAdjRegimeShift_fg`fg'-dfg`dfg'-stayn`stayn'-frachot`frachot'.dta
					}
					}
			}
	}

	
	**************** THESE PROGRAMS ARE PLACED IN programsSep8.do and then run via another program:
************


cap program drop simregime
program define simregime
 syntax [, fg(string) reps(integer 10)]
		local dfgs = "10 20 30 40"
		local frachots= "10 15 20"

		*for some reason
		* `fracthot' == 10 didn't work
		* "`fracthot'" == "10" din't work, and all variations, only the following worked
			foreach dfg of local dfgs{
					foreach frachot of local frachots{
							display `frachot'
							if (`frachot' == 10) {
								local stayns = "98 99" // We can't go too low or probstayhot is to low, this give .82 and .91 respectively
							}
							else if (`frachot'== 15){
								local stayns = "97 99" // We can't go too low or probstayhot is to low, this give .83 and .94 respectively
							
							}
							else if (`frachot'== 20) {
								local stayns = "95 98" // We can't go too low or probstayhot is to low, this give .8 and .92 respectively
							}
							else {
								display "ERROR: Control flow issues"
								exit
							
							}
							foreach stayn of local stayns{
								disp "BiasAdjRegimeShift_fg`fg'-dfg`dfg'-stayn`stayn'-frachot`frachot'.dta"
								simulate  diff=r(diff) ///
								, reps(`reps') saving(BiasAdjRegimeShift_fg`fg'-dfg`dfg'-stayn`stayn'-frachot`frachot'.dta, replace) : calc_diff, command(genshots_regimeswitch, staynormal(.`stayn') frachot(.`frachot') nshots(100) fg(.`fg') dfg(.`dfg') )	method(collapse)																											
							}
					}
			}
	
end


genshots_3boost, nshots(50000) fg(.58) dfg(.1)
summ make


cap program drop sim3boost
program define sim3boost

 syntax [, reps(integer 10)]
		local dfgs = "10 20 30 40"
		local fgs = "40 45 50 55 60"
		matrix matfgn=(.395, .385, .37, .36\ .44, .43, .41, .39\ .48, .47, .45, .42\ .53, .51, .48, .445\ .58, .55, .515, .47)
		
		
			local row=0
			foreach fg of local fgs{
				local row = `row'+1
				local col=0
					foreach dfg of local dfgs{
						 local col=`col'+1
							local fgn=matfgn[`row',`col']
							display "fg = `fg'" ", fgn = `fgn'" ", dfg=`dfg'"
								disp "BiasAdjRegimeShift_fg`fg'-dfg`dfg'-stayn`stayn'-frachot`frachot'.dta"
								simulate  diff=r(diff) ///
								, reps(`reps') saving(BiasAdj3boost_fg`fg'-dfg`dfg'.dta, replace) : calc_diff, command(genshots_3boost, nshots(100) fgn(`fgn') dfg(.`dfg') )	method(collapse)																											
							}
				}
	
end
sim3boost, reps(100)



		mata: probstayhot(.98,.20)
		




*************
******** GF MODEL, test the weighting:



cap program drop prop
program define prop, rclass
	syntax [, nshots(integer 20) nsids(integer 1) fg(real .5)  kstreak(integer 2) ]
	quietly{
	clear
	local obs=`nshots'
	set obs `obs'
	gen rand1=runiform()
	gen period=mod(_n-1,`nshots')+1
	gen make=(rand1<=`fg')
	gen nlast1=.
			replace nlast1=make[_n-1] if _n>1
			forvalues i=2(1)`kstreak'{
				gen nlast`i'=.
				local k=`i'-1
				replace nlast`i' = nlast`k' + make[_n-`i'] if _n>`i'
			} // end forvalues
	gen dMade`kstreak'ormore=.
		replace dMade`kstreak'ormore=(nlast`kstreak'==`kstreak') if _n>`kstreak'
			
	summ make if dMade`kstreak'ormore==1
		if r(N)==0{
			return scalar p=-99
			return scalar M=r(N)
		}
		 else{
			return scalar p=r(mean)
			return scalar M=r(N)	
		 }
	}
	
end

set seed 1234
*.415 for equally weighted T=1
*.46 for a=.5
local T=20
scalar a=2

scalar diff=1
prop
scalar tp=(min(r(M),`T'))^scalar(a)*r(p)
scalar tM=(min(r(M),`T'))^scalar(a)
scalar lastwp=scalar(tp)/scalar(tM)

forvalues i=1(1)100000 {
	prop
	 scalar tp=scalar(tp)+(min(r(M),`T'))^scalar(a)*r(p)
	 scalar tM=scalar(tM)+(min(r(M),`T'))^scalar(a)
	 scalar wp=log(scalar(tp))-log(scalar(tM))
	 disp exp(scalar(wp))
}





