* IV bootstraps

program ivboot, rclass
	* parse
		syntax [anything(name=0)] [if] [in] [aw fw pw iw] [,  							/*
			*/ vce(passthru) robust cluster(passthru) nulllist(numlist ascending) wald ar k j clr 	/*
			*/ wildols wildineff wildeff wildncr scorebs wcs bayesian werser wermdser 	/*
			*/ pairs resideff residineff estfun wildweight(passthru)					/*
			*/ small reps(passthru) SAving(string) trace noisily quietly Level(cilevel) *]
		marksample touse
		* parse variable list
			local n 0
			gettoken lhs 0 : 0, parse(" ,[") match(paren)
			IsStop `lhs'
			if `s(stop)' {
				error 198
			}
			while `s(stop)'==0 {
				if "`paren'"=="(" {
					local n = `n' + 1
					if `n'>1 {
						capture noi error 198
						di in red `"syntax is "(all instrumented variables = instrument variables)""'
						exit 198
					}
					gettoken p lhs : lhs, parse(" =")
					while "`p'"!="=" {
						if "`p'"=="" {
							capture noi error 198
							di in red `"syntax is "(all instrumented variables = instrument variables)""'
							di in red `"the equal sign "=" is required"'
							exit 198
						}
						local endo `endo' `p'
						gettoken p lhs : lhs, parse(" =")
					}
					* To enable Cragg HOLS estimator, allow for empty endo list
					local temp_ct  : word count `endo'
					if `temp_ct' > 0 {
						tsunab endo : `endo'
					}
					* To enable OLS estimator with (=) syntax, allow for empty exexog list
					local temp_ct  : word count `lhs'
					if `temp_ct' > 0 {
						tsunab exexog : `lhs'
					}
				}
				else {
					local inexog `inexog' `lhs'
				}
				gettoken lhs 0 : 0, parse(" ,[") match(paren)
				IsStop `lhs'
			}
			local 0 `"`lhs' `0'"'
			tsunab inexog : `inexog'
			tokenize `inexog'
			local lhs "`1'"
			local 1 " "
			local inexog `*'
			local ninexog=wordcount("`inexog'")
			local nexexog=wordcount("`exexog'")
		* vce parse
			_vce_parse `touse', argopt(cluster) opt(oim opg robust) pwallowed(cluster oim opg robust) old : [`weight'`exp'], `vce' `robust' `cluster'
			local vceopt "`r(vceopt)'"
			local N_clust=0
			if strlen("`r(cluster)'") {
				local cluster "`r(cluster)'"
				tempvar clusterindex
				qui egen `clusterindex'=group(`cluster') if `touse'
				sum `clusterindex' if `touse', meanonly
				local N_clust=r(max)
				sort `clusterindex'
				local clusterexp "cluster(`clusterindex')"
				local vceopt "vce(cluster `clusterindex')"
			}
		* parse null list (default to scalar zero)
			if !strlen("`nulllist'") {
				local nulllist=0
			}
			local nullnum=wordcount("`nulllist'")
		* user choice of tests and bootstraps
			local tests=wordcount("`wald' `ar' `k' `j' `clr'")
			if `tests'==0 {
				di as err `"must specify tests"'
				exit 198
			}
			local boots=wordcount("`wildols' `wildineff' `wildeff' `wildncr' `scorebs' `wcs' `bayesian' `werser' `wermdser' `pairs' `resideff' `residineff' `estfun'")
			if `boots'==0 {
				di as err `"must specify bootstrap"'
				exit 198
			}
			local wildboots=wordcount("`wildols' `wildineff' `wildeff' `wildncr' `scorebs' `wcs' `bayesian' `werser' `wermdser'")
			* note that Rademacher is the default wild boostrap weight (if unspecified)
			local efronboots=wordcount("`pairs' `resideff' `residineff' `estfun'")
		* parse saving
			if `"`saving'"'=="" {
				tempfile saving
				local filetmp "yes"
			}
			else {
				_prefix_saving `saving'
				local saving	`"`s(filename)'"'
				if "`double'" == "" {
					local double	`"`s(double)'"'
				}
				local every	`"`s(every)'"'
				local replace	`"`s(replace)'"'
			}
		* weight processing code and number of observations
			* fweight and aweight accepted as is
			* iweight not allowed with robust or gmm and requires a trap below when used with summarize
			* pweight is equivalent to aweight + robust
			tempvar wvar
			if "`weight'"=="fweight" | "`weight'"=="aweight" | "`weight'" == "iweight" {
				qui gen double `wvar'`exp'
			}
			if "`weight'" == "pweight" {
				qui gen double `wvar'`exp'
				local robust "robust"
			}
			if "`weight'" == "" {
				* If no weights, define neutral weight variable
				qui gen byte `wvar'=1
			}
			if "`weight'"=="fweight" | "`weight'"=="aweight" | "`weight'" == "iweight" {
				local wtexp `"[`weight'`exp']"'
			}
			else if "`weight'" == "pweight" {
				local wtexp `"[aweight`exp']"'
			}
			else {
				local wtexp ""
			}
			* Every time a weight is used, must multiply by scalar wf ("weight factor")
			* wf=1 for no weights, fw and iw, wf = scalar that normalizes sum to be N if aw or pw
			sum `wvar' if `touse' `wtexp', meanonly
			if "`weight'"=="" | "`weight'"=="fweight" | "`weight'"=="iweight" {
				* Effective number of observations is sum of weight variable.
				* If weight is "", weight var must be column of ones and N is number of rows
				local wf=1
				local N=r(sum_w)
			}
			else if "`weight'"=="aweight" | "`weight'"=="pweight" {
				local wf=r(N)/r(sum_w)
				local N=r(N)
			}
		markout `touse' `lhs' `inexog' `exexog' `endo' `cluster' `wvar', strok

	* save asymptotic test statistics
		* 2sls estimation assuming iid errors
			local nullcount=1
			foreach null of numlist `nulllist' {
				qui weakiv ivregress 2sls `lhs' `inexog' (`endo' = `exexog'), `small' null(`null')
				foreach statname in `wald' `ar' `k' `j' `clr' {
					if "`statname'"=="clr" & "`e(`statname'_stat)'"!=""		local clr_iid_`nullcount' = `e(`statname'_stat)'
					else if "`e(`statname'_chi2)'"!=""						local `statname'_iid_`nullcount' = `e(`statname'_chi2)'
					else													local `statname'_iid_`nullcount' = .b
					local p_`statname'_iid_`nullcount' = `e(`statname'_p)'
				}
				local ++nullcount
			}
		* 2sls estimation, rivtest, and store results
			qui ivregress 2sls `lhs' `inexog' (`endo' = `exexog'), `vce' `small'
			local theta=_b[`endo']
			local nullcount=1
			foreach null of numlist `nulllist' {
				qui weakiv ivregress 2sls `lhs' `inexog' (`endo' = `exexog'), `vce' `small' null(`null')
				foreach statname in `wald' `ar' `k' `j' `clr' {
					if "`statname'"=="clr" & "`e(`statname'_stat)'"!=""		local clr_rvce_`nullcount' = `e(`statname'_stat)'
					else if "`e(`statname'_chi2)'"!=""						local `statname'_rvce_`nullcount' = `e(`statname'_chi2)'
					else													local `statname'_rvce_`nullcount' = .b
					local p_`statname'_rvce_`nullcount' = `e(`statname'_p)'
				}
				local ++nullcount
			}
	* generate estimates, residuals, and predictions with and without null hypothesis imposed
		tempvar y_wb x_wb y_r y_xb yh0 yh0_r yh0_xb x_r x_xb xh0 xh0_r xh0_xb u0 e0
		* x and y with w projected out
			qui reg `lhs' `inexog' if `touse' `wtexp'
			qui predict double `y_wb' if `touse', residuals
			qui reg `endo' `inexog' if `touse' `wtexp'
			qui predict double `x_wb' if `touse', residuals
		* unrestricted second stage
			qui reg `lhs' `endo' `inexog' if `touse' `wtexp'
			qui predict double `y_r' if `touse', residuals
			qui predict double `y_xb' if `touse', xb
			
			/*
			
		* restricted second stage
			qui gen double `yh0'=`lhs'-`null'*`endo' if `touse'
			qui reg `yh0' `inexog' if `touse' `wtexp'
			qui predict double `yh0_r' if `touse', residuals
			qui predict double `yh0_xb' if `touse', xb

			*/

		* unrestricted first stage
			qui reg `endo' `exexog' `inexog' if `touse' `wtexp'
			qui predict double `x_r' if `touse', residuals
			qui predict double `x_xb' if `touse', xb
			
			/*
			
		* restricted first stage
			qui gen double `xh0'=`endo'-`yh0_r' if `touse'
			qui reg `xh0' `exexog' `inexog' if `touse' `wtexp'
			qui predict double `xh0_r' if `touse', residuals
			qui predict double `xh0_xb' if `touse', xb
			
			*/
			
			/*
			
		* vars for score bootstrap
			qui gen double `u0'=`y_wb'-`x_wb'*`null' if `touse'
			qui reg `u0' `exexog' if `touse' `wtexp', noconstant	/* constant already projected out */
			qui predict double `e0' if `touse', residuals
		* matrices from sur estimation
			tempname u0eqn xeqn
			
/* DON'T NEED TO DO THIS TWICE, SEE RIGHT ABOVE */			
			
			qui reg `u0' `exexog' if `touse' `wtexp', noconstant	/* constant already projected out */
			estimates store `u0eqn'
			* not bootstrapped, so pass from outside program
			qui reg `x_wb' `exexog' if `touse' `wtexp', noconstant	/* constant already projected out */
			estimates store `xeqn'
			* sur estimation stacking the above two models
			local names `u0eqn' `xeqn'
			tempname hcurrent V Vi b bi
			tempvar esamplei esample
			local scores
			local i 0
			foreach name of local names {
				local ++i
				nobreak {
					if "`name'" != "."		est_unhold `name' `esample'
					else					_est unhold `hcurrent'
					capture noisily break {
						GetMat `name' `bi' `Vi'
						capture drop `esamplei'
						gen byte `esamplei' = e(sample)
						// fix some irregularities in -regress-
						tempvar sc`i'_1 sc`i'_2
						quietly Fix_regress `bi' `Vi' `sc`i'_1' `sc`i'_2'
						local scoresi `sc`i'_1' `sc`i'_2'
					} // capture noisily break
					local rc = _rc
					if "`name'" != "." 		est_hold `name' `esample'
					else					_est hold `hcurrent' , restore nullok estsystem
				} // nobreak
				if (`rc') exit `rc'
				// modifies equation names into name_eq or name#
				FixEquationNames `name' `bi' `Vi'
				local neq`i' `r(neq)'
				local eqnames`i' `"`r(eqnames)'"'
				local newfullnames `"`newfullnames' `:colfullnames `bi''"'
				if `i' == 1 {
					matrix `b' = `bi'
					matrix `V' = `Vi'
				}
				else {
					// append the bi and Vi
					matrix `b' = `b' , `bi'
					local nv  = colsof(`V')
					local nvi = colsof(`Vi')
					matrix `V' = (`V', J(`nv',`nvi',0) \ J(`nvi',`nv',0), `Vi')
				}
				// score vars all models
				local scores `scores' `scoresi'
			} // loop over models
			local Stata11 = cond(c(stata_version)>=11, "version 11:", "")
			`Stata11' matrix colnames `b' = `newfullnames'
			`Stata11' matrix colnames `V' = `newfullnames'
			`Stata11' matrix rownames `V' = `newfullnames'
			_robust `scores' if `touse' `wtexp', var(`V') `clusterexp' minus(0)
		* break up vecs and mats for test components (and make small sample adjustments
			tempname btemp vtemp del pi vardel varpi vardelpi pi0
			mata: `btemp' = st_matrix("`b'")
			mata: `vtemp' = st_matrix("`V'")
			mata: `del' = `btemp'[| 1,1 \ .,`nexexog' |]
			mata: `pi' = `btemp'[| 1,`nexexog'+1 \ .,`nexexog'+`nexexog' |]
			mata: `vardel' = `vtemp'[| 1,1 \ `nexexog',`nexexog' |]
			mata: `varpi' = `vtemp'[| `nexexog'+1,`nexexog'+1 \ `nexexog'+`nexexog',`nexexog'+`nexexog' |]
			mata: `vardelpi' = `vtemp'[| `nexexog'+1,1 \ `nexexog'+`nexexog',`nexexog' |]
			mata: `pi0' = `pi'' - `vardelpi'*cholsolve(`vardel',`del'')
			
			*/
			
	* bootstraps
		if `wildboots' {
			preserve
			tempfile wildfile
			`quietly' di "Performing the following bootstraps:"
			if strlen("`wildols'")		`quietly' di "   wild OLS bootstrap"
			if strlen("`wildineff'")	`quietly' di "   inefficient wild bootstrap"
			if strlen("`wildeff'")		`quietly' di "   efficient wild bootstrap"
			if strlen("`wildncr'")		`quietly' di "   non-cluster robust wild bootstrap"
			if strlen("`scorebs'")		`quietly' di "   score bootstrap"
			if strlen("`wcs'")			`quietly' di "   wild conditional score bootstrap"
			if strlen("`bayesian'")		`quietly' di "   Bayesian conditional bootstrap"
			if strlen("`werser'")		`quietly' di "   wild efficient restricted system equation residual bootstrap"
			if strlen("`wermdser'")		`quietly' di "   wild minimum-distance efficient restricted system equation residual bootstrap"
			`quietly' simulate, `reps' saving(`wildfile', double) `trace' `noisily' nodots nolegend: wildboot if `touse' `wtexp', wvar(`wvar') wf(`wf') `vceopt' `clusterexp' nulllist(`nulllist') `wald' `ar' `k' `j' `clr' `wildols' `wildineff' `wildeff' `wildncr' `scorebs' `wcs' `bayesian' `werser' `wermdser' `wildweight' `small' x_xb(`x_xb') x_r(`x_r') y_wb(`y_wb') x_wb(`x_wb') depvar(`lhs') endo(`endo') inexog(`inexog') exexog(`exexog')
			restore
			`quietly' di
		}
		if `efronboots' {
			preserve
			tempfile efronfile
			`quietly' di "Performing the following bootstraps:"
			if strlen("`pairs'")		`quietly' di "   pairs boostrap"
			if strlen("`residineff'")	`quietly' di "   inefficient residual bootstrap"
			if strlen("`resideff'")		`quietly' di "   efficient residual bootstrap"
			if strlen("`estfun'")		`quietly' di "   estimating function bootstrap"
			if strlen("`clusterexp'") {
				tempvar bsclusterid
				local idcluster "idcluster(`bsclusterid')"
				local bootvce "vce(cluster `bsclusterid')"
				local bootclusterexp "cluster(`bsclusterid')"
			}
			else {
				local bootvce "`vceopt'"
			}
			`quietly' bootstrap, `reps' /* `clusterexp' `idcluster' */ saving(`efronfile', double) `trace' `noisily' notable nodots: efronboot if `touse' `wtexp', wvar(`wvar') wf(`wf') `vceopt' /* `bootvce' `bootclusterexp' */ nulllist(`nulllist') `wald' `ar' `k' `j' `clr' `pairs' `resideff' `residineff' `estfun' `small' x_xb(`x_xb') x_r(`x_r') y_wb(`y_wb') x_wb(`x_wb') depvar(`lhs') endo(`endo') inexog(`inexog') exexog(`exexog')
			restore
		}

	* preserve
		preserve
	* merge bootstrap files into saving
		if `wildboots' {
			use `wildfile', clear
		}
		if `efronboots' {
			if `wildboots' {
				`quietly' merge 1:1 _n using `efronfile', nogen noreport
			}
			else {
				use `efronfile', clear
			}
		}
		renpfix _b_
	* compute and print p-values
		local nullcount=1
		foreach null of numlist `nulllist' {
			foreach t in `wald' `ar' `k' `j' `clr' {
				foreach v of varlist `t'*_`nullcount' {
					
					/* BE CAREFUL HERE WHEN THE STATISTICS ARE MISSING !!! */
					
					qui gen byte r_`v' = cond(``t'_rvce_`nullcount''<=`v',1,0)
				}
			}
			local ++nullcount
		}
		`quietly' sum r_*
	* save other simulation information
		qui gen theta=`theta' in 1
		local nullcount=1
		foreach null of numlist `nulllist' {
			qui gen null_`nullcount'=`null' in 1
			foreach t in `wald' `ar' `k' `j' `clr' {
				foreach a in iid rvce {
					qui gen `t'_`a'_`nullcount' = ``t'_`a'_`nullcount'' in 1
					qui gen p_`t'_`a'_`nullcount' = `p_`t'_`a'_`nullcount'' in 1
					qui gen byte r_`t'_`a'_`nullcount' = cond(`p_`t'_`a'_`nullcount''>=1-`level'/100,1,0) in 1
				}
			}
			local ++nullcount
		}
		`quietly' save `"`saving'"', `replace'
	* restore
		restore
end



* code stolen from ivreg2

program define IsStop, sclass
				/* sic, must do tests one-at-a-time,
				 * 0, may be very large */
	version 8.2
	if `"`0'"' == "[" {
		sret local stop 1
		exit
	}
	if `"`0'"' == "," {
		sret local stop 1
		exit
	}
	if `"`0'"' == "if" {
		sret local stop 1
		exit
	}
* per official ivreg 5.1.3
	if substr(`"`0'"',1,3) == "if(" {
		sret local stop 1
		exit
	}
	if `"`0'"' == "in" {
		sret local stop 1
		exit
	}
	if `"`0'"' == "" {
		sret local stop 1
		exit
	}
	else	sret local stop 0
end

/* Programs borrowed from Stata's suest command */

program Fix_regress
/* - adds equation name "mean" to existing coefficients
   - adds an equation named "lnvar" for the log(variance)
   - returns in the two vars sc1 and sc2 the score variables
*/
	args  b V sc1 sc2
	confirm matrix `b'
	confirm matrix `V'
	tempname b0 var
	// REML estimate of variance
	scalar `var' = e(rmse)^2
	matrix `b0' = log(`var')
	matrix coln `b0' = lnvar:_cons
	local n = colsof(`b')
	matrix coleq `b' = mean
	matrix `b'  = `b', `b0'
	local names : colfullnames `b'
	matrix `V' = (`V', J(`n',1,0) \ J(1,`n',0) , 2/e(N))
	local Stata11 = cond(c(stata_version)>=11, "version 11:", "")
	`Stata11' matrix colnames `V' = `names'
	`Stata11' matrix rownames `V' = `names'
	tempvar res
	predict double `res' if e(sample), res
	gen double `sc1' = `res' / `var'		if e(sample)
	gen double `sc2' = 0.5*(`res'*`sc1' - 1)	if e(sample)
end

program GetMat
	args name b V
	local ev e(V)
	capture {
		confirm matrix e(b)
		confirm matrix `ev'
		matrix `b' = e(b)
		matrix `V' = `ev'
	}
	if _rc {
		dis as err ///
		"impossible to retrieve e(b) and e(V) in `name'"
		exit 198
	}
	if "`e(cmd)'" == "cnsreg" {
		if !missing(e(rmse)) & e(rmse) != 0 {
			matrix `V' = `V'/(e(rmse)*e(rmse))
		}
	}
end

program FixEquationNames, rclass
/* rename the equations to "name" in case of 1/0 equation, otherwise it
   prefixes "name" to equations if this yields unique equation names,
   and numbers the equations "name"_nnn otherwise.
*/
	args name b V
	if "`name'" == "." {
		local name _LAST
	}
	local qeq : coleq `b', quote
	local qeq : list clean qeq
	local eqnames : coleq `b'
	if `:length local qeq' != `:length local eqnames' {
		foreach el of local qeq {
			local new : subinstr local el " " "_", all
			local new : subinstr local new "." ",", all
			local neweq `"`neweq' `new'"'
		}
		matrix coleq `b' = `neweq'
		matrix coleq `V' = `neweq'
		matrix roweq `V' = `neweq'
		local eqnames `"`neweq'"'
	}
	local eq : list uniq eqnames
	local neq : word count `eq'
	if "`eq'" == "_" {
		local eqnames `name'
	}
	else {
		// modify equation names
		foreach e of local eq {
			local newname = substr("`name'_`e'",1,32)
			local meq `meq' `newname'
		}

		local eqmod : list uniq meq
		local neqmod : word count `eqmod'
		if `neq' == `neqmod' {
			// modified equation names are unique
			forvalues i = 1/`neq' {
				local oldname : word `i' of `eq'
				local newname : word `i' of `eqmod'
				local eqnames : subinstr local eqnames "`oldname'" "`newname'", word all
			}
		}
		else {
			// truncated modified equations not unique
			// use name_1, name_2, ...
			tokenize `eq'
			forvalues i = 1/`neq' {
				local eqnames : subinstr local eqnames "``i''" "`name'_`i'", word all
			}
		}
	}
	matrix coleq `b' = `eqnames'
	matrix roweq `V' = `eqnames'
	matrix coleq `V' = `eqnames'
	return local neq `neq'
	return local eqnames	`eq'
	return local neweqnames `eqmod'
end













/*
	// 1st level of cluster, kernel-robust OR
	// 2-level clustering, kernel-robust and time is 2nd cluster variable
			if (vcvo.kernel~="") {
				shat2=J(L,L,0)
	// First, standard cluster-robust, i.e., no lags.
				i=min(t)
				while (i<=max(t)) {  				// loop through all T clusters, adding Z'ee'Z
													// for indiv cluster in each loop
					eZ=J(1,L,0)
					svar=(t:==i)					// select obs with t=i
					if (colsum(svar)>0) {			// there are obs with t=i
						esub=select(*vcvo.e,svar)
						Zsub=select(*vcvo.Z,svar)
						wsub=select(*vcvo.wvar,svar)
						wv = esub :* wsub * vcvo.wf
						eZ = quadcross(1, wv, Zsub)		// equivalent to colsum(wv :* Zsub)
						shat2=shat2+quadcross(eZ,eZ)
					}
					i=i+vcvo.tdelta
				} // end i loop through all T clusters



	tempname delhat
	mata `delhat'=st_matrix("e(b)")
	tempvar ehat zdel
	qui predict double `ehat' if `touse', residuals
	qui predict double `zdel' if `touse', xb
	qui reg x_wb z? if `touse' `wtexp'
	tempname xeqn
	estimates store `xeqn'

	tempname zvec ehatvec wtvec zz zzinv zehat
	foreach x in zvec ehatvec wtvec zz zzinv zehat {
		mata ``x''=.
	}
	fvunab zvars : z?
	mata st_view(`zvec',.,"`zvars'","`touse'")
	mata st_view(`ehatvec',.,"`ehat'","`touse'")
	mata st_view(`wtvec',.,"`wvar'","`touse'")
	mata `zz'=quadcross(`zvec',1,`wf'*`wtvec',`zvec',1)
	mata `zzinv'=invsym(`zz')
	mata `zehat'=quadcross(`zvec',1,`wf'*`wtvec',`ehatvec',1)



	tempvar estar ustar
	qui gen double `estar'=`mult_yr'*`ww1'[1,`idcluster']*`ehat' if `touse'
	qui gen double `ustar'=`zdel'+`estar' if `touse'

	tempname zvec ustarvec estarvec wtvec zustar zestar delstar
	foreach x in zvec ustarvec estarvec wtvec zustar zestar delstar {
		mata ``x''=.
	}
	fvunab zvars : `z'
	mata st_view(`zvec',.,"`zvars'","`touse'")
	mata st_view(`ustarvec',.,"`ustar'","`touse'")
	mata st_view(`ehatvec',.,"`ehat'","`touse'")
	mata st_view(`estarvec',.,"`estar'","`touse'")
	mata st_view(`wtvec',.,"`wvar'","`touse'")
	mata `zustar'=quadcross(`zvec',1,`wf'*`wtvec',`ustarvec',1)
	mata `zestar'=quadcross(`zvec',1,`wf'*`wtvec',`estarvec',1)
	mata `delstar'=`delhat'+`zzinv'*`zestar'

	* regs and sur estimation
	tempname u0eqn
	qui reg `ustar' `z' if `touse' `wtexp', `consopt'
	estimates store `u0eqn'
	* not bootstrapped, so pass from outside program
	*qui reg `x_wb' `z' if `touse' `wtexp', `consopt'
	*estimates store `xeqn'
	* sur estimation stacking the above two models
	local names `u0eqn' `xeqn'
	tempname hcurrent V Vi b bi
	tempvar esamplei esample
	local i 0
	foreach name of local names {
		local ++i
		nobreak {
			if "`name'" != "."		est_unhold `name' `esample'
			else					_est unhold `hcurrent'
			capture noisily break {
				GetMat `name' `bi' `Vi'
				capture drop `esamplei'
				gen byte `esamplei' = e(sample)
				// fix some irregularities in -regress-
				tempvar sc`i'_1 sc`i'_2
				quietly Fix_regress `bi' `Vi' `sc`i'_1' `sc`i'_2'
				local scoresi `sc`i'_1' `sc`i'_2'
			} // capture noisily break
			local rc = _rc
			if "`name'" != "." 		est_hold `name' `esample'
			else					_est hold `hcurrent' , restore nullok estsystem
		} // nobreak
		if (`rc') exit `rc'
		// modifies equation names into name_eq or name#
		FixEquationNames `name' `bi' `Vi'
		local neq`i' `r(neq)'
		local eqnames`i' `"`r(eqnames)'"'
		local newfullnames `"`newfullnames' `:colfullnames `bi''"'
		if `i' == 1 {
			matrix `b' = `bi'
			matrix `V' = `Vi'
		}
		else {
			// append the bi and Vi
			matrix `b' = `b' , `bi'
			local nv  = colsof(`V')
			local nvi = colsof(`Vi')
			matrix `V' = (`V', J(`nv',`nvi',0) \ J(`nvi',`nv',0), `Vi')
		}
		// score vars all models
		local scores `scores' `scoresi'
	} // loop over models
	local Stata11 = cond(c(stata_version)>=11, "version 11:", "")
	`Stata11' matrix colnames `b' = `newfullnames'
	`Stata11' matrix colnames `V' = `newfullnames'
	`Stata11' matrix rownames `V' = `newfullnames'
	qui _robust `scores' if `touse' `wtexp', var(`V') `vceopt' minus(0)
* break up vecs and mats for test components (and make small sample adjustments
	tempname btemp vtemp del_z vardel pi_z var_pi_z var_pidel_z
	mata `btemp' = st_matrix("e(b)")
	mata `vtemp' = st_matrix("e(V)")

	*mata `del_z' = `btemp'[| 1,1 \ .,`nexexog' |]
	mata `del_z' = `delstar'-`delhat'

	mata `vardel' = `mult_xr' * `vtemp'[| 1,1 \ `nexexog',`nexexog' |]
	mata `pi_z' = `btemp'[| 1,`nexexog'+`ninexog'+2 \ .,`nexexog'+`ninexog'+`nexexog'+1 |]
	mata `var_pi_z' = `mult_xr' * `vtemp'[| `nexexog'+`ninexog'+2,`nexexog'+`ninexog'+2 \ `nexexog'+`ninexog'+`nexexog'+1,`nexexog'+`ninexog'+`nexexog'+1 |]
	mata `var_pidel_z' = `mult_xr' * `vtemp'[| `nexexog'+`ninexog'+2,1 \ `nexexog'+`ninexog'+`nexexog'+1,`nexexog' |]


*/















