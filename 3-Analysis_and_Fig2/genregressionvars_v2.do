	quietly{
capture:program drop genregressionvars
program define genregressionvars
*Regressions
	quietly{
		*Get the different sid values available for th
		sort sid session period
		levelsof sid, local(sidlevels)
			foreach i in `sidlevels' {
				levelsof session if sid==`i', local(sessionlevels`i')
		}
		*count how many sids there are
		local ndistinct=wordcount("`sidlevels'")
		
		if `ndistinct'<2 {
			noi: di  _newline  "***Warning: program genregressionvars sees only one unique sid in database***" _newline 
		}
		
		drop if make==.

		*create # of shots made in last i attempts: nlasti
		gen nlast1=.
		by sid session,sort:replace nlast1=make[_n-1] if _n>1
		forvalues i=2(1)5{
			gen nlast`i'=.
			local k=`i'-1
			by sid,sort:replace nlast`i' = nlast`k' + make[_n-`i'] if _n>`i'
		} // end forvalues

		*create indicators for making i or more shots, or exactly i previous shots
		forvalues i=1(1)5{
			gen dMade`i'ormore=.
			by sid session,sort:replace dMade`i'ormore=(nlast`i'==`i') if _n>`i'
			
			gen dMissed`i'ormore=.
			by sid session,sort:replace dMissed`i'ormore=(nlast`i'==0) if _n>`i'
			
			gen dMade`i'exactly=.
			by sid session,sort:replace  dMade`i'exactly=(nlast`i'==`i' & make[_n-`i'-1]==0) if _n>`i'+1
			
			gen dMissed`i'exactly=.
			by sid session,sort:replace  dMissed`i'exactly=(nlast`i'==0 & make[_n-`i'-1]==1) if _n>`i'+1
			} // end forvalues

} // end quietly

	
end
}
