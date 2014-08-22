wild-bootstrap
==============

General wild bootstrap command for Stata

Goals:
* Wrapper for `simulate` and mirrors the options
* Goal is to accept all estimators
	* But start with single-equation, linear models
* Sampling reflects VCE or preset survey sampling structure
	* Menu of wild bootstrap weights
* Accepts a Stata restriction matrix to generate null hypotheses
	* But the default might accept a simpler single null hypothesis (using the `test` syntax?)
* Should largely be written in Mata to reduce run time
* Help file