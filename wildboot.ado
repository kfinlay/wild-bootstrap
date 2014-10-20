* wild bootstraps

program wildboot, eclass

	program wildboot
	        version 9, missing
	        local version : di "version " string(_caller()) ", missing:"

	        // <my_stuff> : <command>
	        _on_colon_parse `0'
	        local command `"`s(after)'"'
	        local 0 `"`s(before)'"'

	        syntax [anything(name=exp_list equalok)]        ///
	                [fw iw pw aw] [if] [in] [,              ///
	                        noDOTS                          ///
	                        Reps(integer -1)                ///
	                        SAving(string)                  ///
	                        DOUBle                          /// not documented
	                        noLegend                        ///
	                        Verbose                         ///
	                        SEED(string)                    ///
	                        NOIsily                         /// "prefix" options
	                        TRace                           ///
	                ]

	        if "`weight'" != "" {
	                local wgt [`weight'`exp']
	        }


	syntax [if] [in] [aw fw pw iw] [, wvar(varname) wf(real 1) 	/*
		*/ vce(passthru) nulllist(numlist ascending) wald ar /*
		*/ k j clr wildols wildineff wildeff wildncr scorebs wcs bayesian werser wermdser /*
		*/ wildweight(string) small /*
		*/ x_xb(varname) x_r(varname) y_wb(varname) x_wb(varname) zdel(varname) ehat(varname) /*
		*/ xeqn(name) delhat(name) zz(name) zzinv(name) zehat(name) /*
		*/ depvar(varname) endo(varname) inexog(varlist) exexog(varlist) *]
	marksample touse
	sum `wvar' if `touse' [`weight'`exp'], meanonly
	local N=r(sum_w)*`wf'
	`quietly' simulate, `reps' saving(`wildfile', double) `trace' `noisily' nodots nolegend: WildBootCore if `touse' `wtexp', wvar(`wvar') wf(`wf') `vceopt' `clusterexp' nulllist(`nulllist') `wald' `ar' `k' `j' `clr' `wildols' `wildineff' `wildeff' `wildncr' `scorebs' `wcs' `bayesian' `werser' `wermdser' `wildweight' `small' x_xb(`x_xb') x_r(`x_r') y_wb(`y_wb') x_wb(`x_wb') depvar(`lhs') endo(`endo') inexog(`inexog') exexog(`exexog')
	* post results
		matrix colnames `b' = `testbootlist'
		ereturn post `b' `wtexp', esample(`touse') obs(`N')
end

program WildBootCore, eclass
	* parse
		syntax [if] [in] [aw fw pw iw] [, wvar(varname) wf(real 1) 	/*
			*/ vce(passthru) nulllist(numlist ascending) wald ar /*
			*/ k j clr wildols wildineff wildeff wildncr scorebs wcs bayesian werser wermdser /*
			*/ wildweight(string) small /*
			*/ x_xb(varname) x_r(varname) y_wb(varname) x_wb(varname) zdel(varname) ehat(varname) /*
			*/ xeqn(name) delhat(name) zz(name) zzinv(name) zehat(name) /*
			*/ depvar(varname) endo(varname) inexog(varlist) exexog(varlist) *]
		marksample touse
		sum `wvar' if `touse' [`weight'`exp'], meanonly
		local N=r(sum_w)*`wf'
		* vce parse
			_vce_parse `touse', argopt(cluster) opt(oim opg robust) pwallowed(cluster oim opg robust) old : [`weight'`exp'], `vce' `robust' `cluster'
			local vceopt "`r(vceopt)'"
			local N_clust=0
			if strlen("`r(cluster)'") {
				local cluster "`r(cluster)'"
				tempvar clusterid
				qui egen `clusterid'=group(`cluster') if `touse'
				sum `clusterid' if `touse', meanonly
				local N_clust=r(max)
				sort `clusterid'
				local clusterexp "cluster(`clusterid')"
				local vceopt "vce(cluster `clusterid')"
			}
		local ninexog=wordcount("`inexog'")
		local nexexog=wordcount("`exexog'")
		local tests=wordcount("`wald' `ar' `k' `j' `clr'")
		if `tests'==0 {
			di as err `"must specify tests"'
			exit 198
		}
		local boots=wordcount("`wildols' `wildineff' `wildeff' `wildncr' `scorebs' `wcs' `bayesian' `werser' `wermdser'")
		if `boots'==0 {
			di as err `"must specify bootstrap"'
			exit 198
		}
		* parse null list (default to scalar zero)
			if !strlen("`nulllist'") {
				local nulllist=0
			}
			local nullnum=wordcount("`nulllist'")
		* results matrix
			local bindex = 0
			tempname b
			matrix `b' = J(1,`tests'*`boots'*`nullnum',.a)

	* make a finite-sample multiplier for residuals
		if "`small'"=="small" {
			local mult_yr=sqrt( (`N_clust'/(`N_clust'-1)) * ((`N'-1)/(`N'-`ninexog')) )
			local mult_xr=sqrt( (`N_clust'/(`N_clust'-1)) * ((`N'-1)/(`N'-`nexexog'-`ninexog')) )
		}
		else {
			local mult_yr=1
			local mult_xr=1
		}


/*	* make multinomial bootstrap weights
		tempname mw0 mw1
		mata:
		n = `N_clust'
		p = 1/n
		draw = J(st_nobs(),3,.)
		for(i=1; i<=rows(draw); i++) {
		         trials = uniform(1,n)
		         g1 = trials :< p1
		         g2 = trials :>= p1 :& trials :< p1 + p2
		         draw[i,.] = rowsum(g1),
		                     rowsum(g2),
		                     n - rowsum(g1) - rowsum(g2)
		}
		idx = st_addvar("long", ("G1", "G2", "G3"))
		st_store(.,idx,draw)
		end
		sum G*
*/
	* make wild cluster weights (Rademacher is the default)
		tempname ww0 ww1
		if "`wildweight'"=="mammen" {
			* Mammen weights
			local mammenp "((1+sqrt(5))/(2*sqrt(5)))"
			local mammenf "((1-sqrt(5))/2)"
			tempname mm
			mata: `mm' = runiform(1,`N_clust') :< J(1,`N_clust',`mammenp')
			mata: `ww0' = `mm':*J(1,`N_clust',`mammenp') + (I(1,`N_clust')-`mm'):*J(1,`N_clust',1-`mammenp')
		}
		else if "`wildweight'"=="liunormal" {
			* Liu normal weights
			local wildmean1 "(.5*(sqrt(17/6)+sqrt(1/6)))"
			local wildmean2 "(.5*(sqrt(17/6)-sqrt(1/6)))"
			local wildsd "(sqrt(1/2))"
			mata: `ww0' = rnormal(1,`N_clust',`wildmean1',`wildsd') :* ///
				rnormal(1,`N_clust',`wildmean2',`wildsd') - J(1,`N_clust',`wildmean1'*`wildmean2')
		}
		else if "`wildweight'"=="liugamma" {
			* Liu gamma weights - note that Stata's gamma random generator takes as input the inverse scale parameter
			mata: `ww0' = rgamma(1,`N_clust',4,.5) - J(1,`N_clust',2)
		}
		else {
			* Rademacher weights
			mata: `ww0' = 2*round(runiform(1,`N_clust')):-J(1,`N_clust',1)
		}
		* save into a stata matrix
		mata: st_matrix("`ww1'",`ww0')
	* which bootstraps?
		* things to estimate before looping over the nulls
			if "`wildineff'"=="wildineff" {
				* Wild cluster - only y restricted (inefficient)
				tempvar xstar_r1
				qui gen double `xstar_r1'=`x_xb'+`mult_xr'*`ww1'[1,`clusterid']*`x_r' if `touse'
			}
			if "`wildncr'"=="wildncr" {
				* Wild cluster - only y restricted - no cluster-robust adjustments for test stats
				tempvar xstar_iid
				qui gen double `xstar_iid'=`x_xb'+`mult_xr'*`ww1'[1,`clusterid']*`x_r' if `touse'
			}
		* now loop over the nulls
			local nullcount=1
			foreach null of numlist `nulllist' {
				* estimate null-specific elements common to multiple bootstraps here
					if strlen("`wildineff'`wildeff'`wildncr'") {
						* restricted second stage
						foreach x in yh0 yh0_r yh0_xb {
							capt drop `x'
						}
						tempvar yh0 yh0_r yh0_xb
						qui gen double `yh0'=`depvar'-`null'*`endo' if `touse'
						qui reg `yh0' `inexog' if `touse' `wtexp'
						qui predict double `yh0_r' if `touse', residuals
						qui predict double `yh0_xb' if `touse', xb
					}
					if strlen("`wildeff'") {
						* restricted first stage
						foreach x in xh0 xh0_r xh0_xb {
							capt drop `x'
						}
						tempvar xh0 xh0_r xh0_xb
						qui gen double `xh0'=`endo'-`yh0_r' if `touse'
						qui reg `xh0' `exexog' `inexog' if `touse' `wtexp'
						qui predict double `xh0_r' if `touse', residuals
						qui predict double `xh0_xb' if `touse', xb
					}
					if strlen("`scorebs'") {
						foreach x in u0 e0 {
							capt drop `x'
						}
						tempvar u0 e0
						* vars for score bootstrap
							qui gen double `u0'=`y_wb'-`x_wb'*`null' if `touse'
							qui reg `u0' `exexog' if `touse' `wtexp', noconstant	/* constant already projected out */
							qui predict double `e0' if `touse', residuals
					}
				* this is the main bootstrap code
			 		if strlen("`wildineff'") {
						* Wild cluster - only y restricted (inefficient)
						return clear
						local bootname wildineff
						foreach x in xstar_r1 ystar_r1 {
							capt drop `x'
						}
						tempvar xstar_r1
						qui gen double `xstar_r1'=`x_xb'+`mult_xr'*`ww1'[1,`clusterid']*`x_r' if `touse'
						tempvar ystar_r1
						qui gen double `ystar_r1'=`yh0_xb'+`null'*`xstar_r1'+`mult_yr'*`ww1'[1,`clusterid']*`yh0_r' if `touse'
						qui weakiv ivregress 2sls `ystar_r1' `inexog' (`xstar_r1' = `exexog') if `touse' [`weight'`exp'], `vce' `small' null(`null')
						foreach statname in `wald' `ar' `k' `j' `clr' {
							local bindex=`bindex'+1
							if "`statname'"=="clr" & "`e(`statname'_stat)'"!=""		matrix `b'[1,`bindex'] = `e(`statname'_stat)'
							else if "`e(`statname'_chi2)'"!=""						matrix `b'[1,`bindex'] = `e(`statname'_chi2)'
							else													matrix `b'[1,`bindex'] = .b
							local testbootlist "`testbootlist' `statname'_`bootname'_`nullcount'"
						}
					}
					if strlen("`wildeff'") {
						* Wild cluster - x and y restricted (efficient)
						return clear
						local bootname wildeff
						foreach x in xstar_r ystar_r {
							capt drop `x'
						}
						tempvar xstar_r ystar_r
						qui gen double `xstar_r'=`xh0_xb'+`mult_xr'*`ww1'[1,`clusterid']*`xh0_r' if `touse'
						qui gen double `ystar_r'=`yh0_xb'+`null'*`xstar_r'+`mult_yr'*`ww1'[1,`clusterid']*`yh0_r' if `touse'
						qui weakiv ivregress 2sls `ystar_r' `inexog' (`xstar_r' = `exexog') if `touse' [`weight'`exp'], `vce' `small' null(`null')
						foreach statname in `wald' `ar' `k' `j' `clr' {
							local bindex=`bindex'+1
							if "`statname'"=="clr" & "`e(`statname'_stat)'"!=""		matrix `b'[1,`bindex'] = `e(`statname'_stat)'
							else if "`e(`statname'_chi2)'"!=""						matrix `b'[1,`bindex'] = `e(`statname'_chi2)'
							else													matrix `b'[1,`bindex'] = .b
							local testbootlist "`testbootlist' `statname'_`bootname'_`nullcount'"
						}
					}
					if strlen("`wildncr'") {
						* Wild cluster - only y restricted - no cluster-robust adjustments for test stats
						return clear
						local bootname wildncr
						foreach x in ystar_iid {
							capt drop `x'
						}
						tempvar ystar_iid
						qui gen double `ystar_iid'=`yh0_xb'+`null'*`xstar_iid'+`mult_yr'*`ww1'[1,`clusterid']*`yh0_r' if `touse'
						qui weakiv ivregress 2sls `ystar_iid' `inexog' (`xstar_iid' = `exexog') if `touse' [`weight'`exp'], `small' null(`null') `kwt'
						foreach statname in `wald' `ar' `k' `j' `clr' {
							local bindex=`bindex'+1
							if "`statname'"=="clr" & "`e(`statname'_stat)'"!=""		matrix `b'[1,`bindex'] = `e(`statname'_stat)'
							else if "`e(`statname'_chi2)'"!=""						matrix `b'[1,`bindex'] = `e(`statname'_chi2)'
							else													matrix `b'[1,`bindex'] = .b
							local testbootlist "`testbootlist' `statname'_`bootname'_`nullcount'"
						}
					}
					if strlen("`scorebs'") {
						return clear
						local bootname scorebs
						tempname wvarvec clustervec e0vec Zmat scorehsum scoreV
						foreach x in wvarvec clustervec e0vec Zmat scorehsum scoreV {
							mata: ``x''=.
						}
						mata: st_view(`wvarvec',.,"`wvar'","`touse'")
						mata: st_view(`clustervec',.,"`clusterid'","`touse'")
						mata: st_view(`e0vec',.,"`e0'","`touse'")
						mata: st_view(`Zmat',.,"`exexog'","`touse'")
						mata: `scorehsum'=WildClusterCross(`Zmat',`e0vec',`clustervec',`wvarvec',`wf',`ww0',J(1,`N_clust',1))
						mata: `scoreV'=WildClusterOuterDev(`Zmat',`e0vec',(1/`N_clust')*`scorehsum',`clustervec',`wvarvec',`wf',`ww0',J(1,`N_clust',1))
						mata: ScoreTests(`scorehsum',`scoreV',`pi0',"`wald'","`ar'","`k'","`j'","`clr'")
						foreach statname in `wald' `ar' `k' `j' `clr' {
							local bindex=`bindex'+1
							if "`statname'"=="clr" & "`e(`statname'_stat)'"!=""		matrix `b'[1,`bindex'] = `e(`statname'_stat)'
							else if "`e(`statname'_chi2)'"!=""						matrix `b'[1,`bindex'] = `e(`statname'_chi2)'
							else													matrix `b'[1,`bindex'] = .b
							local testbootlist "`testbootlist' `statname'_`bootname'_`nullcount'"
						}
					}
					/* MAYBE THIS IS THE SAME AS SCOREBS ABOVE */
					if strlen("`wcs'") {
						return clear
						local bootname wcs

						foreach statname in `wald' `ar' `k' `j' `clr' {
							local bindex=`bindex'+1
							if "`statname'"=="clr" & "`e(`statname'_stat)'"!=""			matrix `b'[1,`bindex'] = `e(`statname'_stat)'
							else if "`e(`statname'_chi2)'"!=""							matrix `b'[1,`bindex'] = `e(`statname'_chi2)'
							else														matrix `b'[1,`bindex'] = .b
							local testbootlist "`testbootlist' `statname'_`bootname'_`nullcount'"
						}
					}
					if strlen("`bayesian'") {
						return clear
						local bootname bayesian
						* Dirichlet weights - note that Stata's gamma random generator takes as input the inverse scale parameter
							tempname bw0 bw1 bwvar
							* first generate vector of gammas with beta=1 and 1/alpha=1
							mata: `bw0' = rgamma(1,`N_clust',1,1/1)
							* then, divide by the sum of the gammas
							mata: `bw0' = `bw0' :/ J(1,`N_clust',sum(`bw0'))
							mata: st_matrix("`bw1'",`bw0')
							qui gen double `bwvar'=`exp'*`bw1'[1,`clusterid'] if `touse'
						* matrices from sur estimation
							qui gen double `u0'=`y_wb'-`x_wb'*`null' if `touse'
							qui reg `u0' `exexog' if `touse' `wtexp', noconstant	/* constant already projected out */
							tempname u0eqn xeqn
							qui reg `u0' `exexog' if `touse' [`weight'`bwvar'], noconstant	/* constant already projected out */
							estimates store `u0eqn'
							* not bootstrapped, so pass from outside program
							qui reg `x_wb' `exexog' if `touse' [`weight'`bwvar'], noconstant	/* constant already projected out */
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
							tempname bay_btemp bay_vtemp bay_del bay_pi bay_vardel bay_varpi bay_vardelpi bay_pi0
							mata: `bay_btemp' = st_matrix("`b'")
							mata: `bay_vtemp' = st_matrix("`V'")
							mata: `bay_del' = `bay_btemp'[| 1,1 \ .,`nexexog' |]
							mata: `bay_pi' = `bay_btemp'[| 1,`nexexog'+1 \ .,`nexexog'+`nexexog' |]
							mata: `bay_vardel' = `bay_vtemp'[| 1,1 \ `nexexog',`nexexog' |]
							mata: `bay_varpi' = `bay_vtemp'[| `nexexog'+1,`nexexog'+1 \ `nexexog'+`nexexog',`nexexog'+`nexexog' |]
							mata: `bay_vardelpi' = `bay_vtemp'[| `nexexog'+1,1 \ `nexexog'+`nexexog',`nexexog' |]
							mata: `bay_pi0' = `bay_pi'' - `bay_vardelpi'*cholsolve(`bay_vardel',`bay_del'')
						* compute tests
							tempname wvarvec clustervec e0vec Zmat scorehsum scoreV
							foreach x in wvarvec clustervec e0vec Zmat scorehsum scoreV {
								mata: ``x''=.
							}
							mata: st_view(`wvarvec',.,"`wvar'","`touse'")
							mata: st_view(`clustervec',.,"`clusterid'","`touse'")
							mata: st_view(`e0vec',.,"`e0'","`touse'")
							mata: st_view(`Zmat',.,"`exexog'","`touse'")
							mata: `scorehsum'=WildClusterCross(`Zmat',`e0vec',`clustervec',`wvarvec',`wf',J(1,`N_clust',1),`bw0')
							mata: `scoreV'=WildClusterOuterDev(`Zmat',`e0vec',(1/`N_clust')*`scorehsum',`clustervec',`wvarvec',`wf',J(1,`N_clust',1),`bw0')
							mata: ScoreTests(`scorehsum',`scoreV',`pi0',"`wald'","`ar'","`k'","`j'","`clr'")
						foreach statname in `wald' `ar' `k' `j' `clr' {
							local bindex=`bindex'+1
							if "`statname'"=="clr" & "`e(`statname'_stat)'"!=""			matrix `b'[1,`bindex'] = `e(`statname'_stat)'
							else if "`e(`statname'_chi2)'"!=""							matrix `b'[1,`bindex'] = `e(`statname'_chi2)'
							else														matrix `b'[1,`bindex'] = .b
							local testbootlist "`testbootlist' `statname'_`bootname'_`nullcount'"
						}

					}
					if strlen("`werser'") {
						return clear
						local bootname werser

						foreach statname in `wald' `ar' `k' `j' `clr' {
							local bindex=`bindex'+1
							if "`statname'"=="clr" & "`e(`statname'_stat)'"!=""			matrix `b'[1,`bindex'] = `e(`statname'_stat)'
							else if "`e(`statname'_chi2)'"!=""							matrix `b'[1,`bindex'] = `e(`statname'_chi2)'
							else														matrix `b'[1,`bindex'] = .b
							local testbootlist "`testbootlist' `statname'_`bootname'_`nullcount'"
						}
					}
					if strlen("`wermdser'") {
						return clear
						local bootname wermdser

						foreach statname in `wald' `ar' `k' `j' `clr' {
							local bindex=`bindex'+1
							if "`statname'"=="clr" & "`e(`statname'_stat)'"!=""			matrix `b'[1,`bindex'] = `e(`statname'_stat)'
							else if "`e(`statname'_chi2)'"!=""							matrix `b'[1,`bindex'] = `e(`statname'_chi2)'
							else														matrix `b'[1,`bindex'] = .b
							local testbootlist "`testbootlist' `statname'_`bootname'_`nullcount'"
						}
					}
				* end of null loop
					local ++nullcount
			}
	* post results
		matrix colnames `b' = `testbootlist'
		ereturn post `b' `wtexp', esample(`touse') obs(`N')
end

mata:
void computeivtests_robust(real matrix del_z, real matrix vardel, real matrix pi_z, real matrix var_pi_z, real matrix var_pidel_z, real matrix zz, real matrix zustar , real matrix delstar, real matrix delhat, real scalar null)
{
	// calculate matrices for test stats
		notposdef=0


		psi=zz * vardel * zz

		aux99 = cholsolve(vardel,delstar)
		if (aux1[1,1]==.) {
			notposdef = 1
			aux8 = qrsolve(vardel,delstar)
		}

		aux98 = cholsolve(vardel,delstar-delhat)
		if (aux1[1,1]==.) {
			notposdef = 1
			aux8 = qrsolve(vardel,delstar-delhat)
		}

		aux8 = cholsolve(psi,delstar)
		if (aux1[1,1]==.) {
			notposdef = 1
			aux8 = qrsolve(psi,del_z)
		}

		aux1 = cholsolve(psi,zustar)
		if (aux1[1,1]==.) {
			notposdef = 1
			aux1 = qrsolve(psi,zustar)
		}

		aux7 = cholsolve(vardel,del_z)
		if (aux1[1,1]==.) {
			notposdef = 1
			aux7 = qrsolve(vardel,del_z)
		}
		pi_beta = pi_z - var_pidel_z*aux7


		aux2 = var_pidel_z - (null)*var_pi_z
		aux3 = cholsolve(psi,aux2')
		if (aux3[1,1]==.) {
			notposdef = 1
			aux3 = qrsolve(psi,aux2')
		}

		aux4 = var_pi_z - aux2 * aux3

		rk = cholsolve(aux4, pi_beta)
		if (rk[1,1]==.) {
			notposdef = 1
			rk = qrsolve(aux3,pi_beta)
		}
		rk = pi_beta' * rk
		aux5 = cholsolve(psi,pi_beta)
		if (aux5[1,1]==.) {
			notposdef = 1
			aux5 = qrsolve(psi,pi_beta)
		}
		aux6 = cholsolve(pi_beta'*aux5,pi_beta')
		if (aux6[1,1]==.) {
			notposdef = 1
			aux6 = qrsolve(pi_beta'*aux5,pi_beta')
		}
	// calculate test stats
		wald_chi2 = delstar' * aux99
		ar_chi2 = (delstar-delhat)' * aux98
		k_chi2 = del_z' * aux5 * aux6 * aux8
		j_chi2 = ar_chi2 - k_chi2
		clr_stat = .5*(ar_chi2-rk+sqrt((ar_chi2+rk)^2 - 4*j_chi2*rk))
		if (rk[1,1]<=0)			clr_stat=.
	// return test stats in r()
		st_numscalar("r(wald_chi2)", wald_chi2[1,1])
		st_numscalar("r(ar_chi2)", ar_chi2[1,1])
		st_numscalar("r(k_chi2)", k_chi2[1,1])
		st_numscalar("r(j_chi2)", j_chi2[1,1])
		st_numscalar("r(clr_stat)", clr_stat[1,1])
		st_numscalar("r(rk)", rk[1,1])
		st_numscalar("r(notposdef)", notposdef)
}
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

mata:
void ScoreTests(real matrix efhsum, real matrix efV, real matrix pi0, string scalar wald, string scalar ar, string scalar k, string scalar j, string scalar clr)
{
	notposdef=0
	if (strlen(ar)+strlen(k)) {
		aux1 = cholsolve(efV,efhsum)
		if (aux1[1,1]==.) {
			notposdef = 1
			aux1 = qrsolve(efV,efhsum)
		}
		ar = efhsum'*aux1
		st_numscalar("r(ar_chi2)", ar[1,1])
	}
	if (strlen(k)) {
		aux2 = cholsolve(efV,pi0)
		if (aux2[1,1]==.) {
			notposdef = 1
			aux2 = qrsolve(efV,pi0)
		}
		aux3 = cholsolve(pi0'*aux2,pi0')
		if (aux3[1,1]==.) {
			notposdef = 1
			aux3 = qrsolve(pi0'*aux2,pi0')
		}
		k = efhsum' * aux2 * aux3 * aux1
		st_numscalar("r(k_chi2)", k[1,1])
	}
	st_numscalar("r(notposdef)", notposdef)
}
end

mata:
real matrix WildClusterCross(real matrix X1, real matrix X2, real vector G, real vector wvar, real scalar wf, real vector ww, real vector bw)
{
	// bw is an additional weight that is only used for the Bayesian bootstrap (Leandro's name for it)
	// otherwise, a G-dimensional vector of ones should be passed: J(1,`N_clust',1)
	shat=J(cols(X1),cols(X2),0)
	i=min(G)
	while (i<=max(G)) {  				// loop through all G clusters, adding X1'X2
										// for indiv cluster in each loop
		svar=(G:==i)					// select obs with G=i
		if (colsum(svar)>0) {			// there are obs with G=i
			X1sub=select(X1,svar)
			X2sub=select(X2,svar)
			wsub=select(wvar,svar)
			shat=shat+quadcross(X1sub,bw[1,i]*wf*wsub,ww[1,i]*X2sub)
		}
		i=i+1
	} // end i loop through all G clusters
	return(shat)
}
end

mata:
real matrix WildClusterOuterDev(real matrix X1, real matrix X2, real matrix M, real vector G, real vector wvar, real scalar wf, real vector ww, real vector bw)
{
	// bw is an additional weight that is only used for the Bayesian bootstrap (Leandro's name for it)
	// otherwise, a G-dimensional vector of ones should be passed: J(1,`N_clust',1)
	shat=J(rows(M),rows(M),0)
	i=min(G)
	while (i<=max(G)) {  				// loop through all G clusters, adding X1'X2
										// for indiv cluster in each loop
		svar=(G:==i)					// select obs with G=i
		if (colsum(svar)>0) {			// there are obs with G=i
			X1sub=select(X1,svar)
			X2sub=select(X2,svar)
			wsub=select(wvar,svar)
			h=quadcross(X1sub,bw[1,i]*wf*wsub,ww[1,i]*X2sub)-M
			shat=shat+h*h'
		}
		i=i+1
	} // end i loop through all G clusters
	return(shat)
}
end

mata:
real matrix Sandwich(real matrix Z, real matrix e, real matrix G, real matrix wvar, real scalar wf)
{
	L=cols(Z)
	shat=J(L,L,0)
	i=min(G)
	while (i<=max(G)) {  				// loop through all G clusters, adding X1'X2
										// for indiv cluster in each loop
		svar=(G:==i)					// select obs with G=i
		if (colsum(svar)>0) {			// there are obs with G=i
			esub=select(e,svar)
			Zsub=select(Z,svar)
			wsub=select(wvar,svar)
			wv = esub :* wsub * wf
			eZ = quadcross(1, wv, Zsub)		// equivalent to colsum(wv :* Zsub)
			shat=shat+quadcross(eZ,eZ)
		}
		i=i+1
	} // end i loop through all G clusters
	return(shat)
}
end













































