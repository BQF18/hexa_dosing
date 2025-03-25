* commands for gsemmetan sj article
* 17apr2014

/* Univariate outcome meta-analysis models */

** Fixed effect meta-analysis
clear all

sjlog using gsemmetan1, replace
webuse diuretics
rename or logor
rename varor varlogor
rename std selogor
metan logor selogor, fixed nograph
sjlog close, replace


sjlog using gsemmetan2, replace
generate double weight = 1/varlogor
sem (logor <- ) [iw=weight], variance(e.logor@1) nodescribe nocnsreport nolog
sjlog close, replace


lincom _b[logor:_cons], eform


sjlog using gsemmetan3, replace
generate double invselogor = 1/selogor
generate double logortr = logor*invselogor
sem (logortr <- invselogor, noconstant), noheader nodescribe nocnsreport nolog
local Q = _b[var(e.logortr):_cons]*e(N)
local df = e(N) - 1
display "Het. test statistic = " `Q'
display "Het. test p-value = " chi2tail(`df', `Q')
display "I-squared = " (`Q' - `df')/`Q'
sjlog close, replace


** Random effects meta-analysis

sjlog using gsemmetan4, replace
mkmat varlogor, matrix(f)
matrix f = diag(f)
quietly tabulate trial, gen(tr)
gsem (logor <- M1#c.tr1@1 M2#c.tr2@1 M3#c.tr3@1 ///
	 M4#c.tr4@1 M5#c.tr5@1 M6#c.tr6@1 ///
	 M7#c.tr7@1 M8#c.tr8@1 M9#c.tr9@1) ///
	, covstructure(_LEx, fixed(f)) intmethod(laplace) nocnsreport nolog
sjlog close, replace


sjlog using gsemmetan5, replace
local setotal = sqrt(_se[logor:_cons]^2 + _b[var(e.logor):_cons])
local pilow = _b[logor:_cons] - invt(e(N) - 2, .975)*`setotal'
local piupp = _b[logor:_cons] + invt(e(N) - 2, .975)*`setotal'
display "95% Prediction interval:", `pilow', `piupp'
sjlog close, replace


metan logor selogor, random rfdist nograph
local setotal = sqrt(r(seES)^2 + r(tau2))
local pilow = r(ES) - invt(r(df) - 1, .975)*`setotal'
local piupp = r(ES) + invt(r(df) - 1, .975)*`setotal'
display "95% Prediction interval:", `pilow', `piupp'


gsem (logor <- ibn.trial#c.selogor#c.M@1), variance(M@1) nolog nocnsreport


gsem (logortr <- invselogor c.invselogor#c.M@1, noconstant), ///
	variance(e.logortr@1) latent(M) nolog nocnsreport


generate b_logor_1 = logor
generate V_logor_1_logor_1 = varlogor
mvmeta b V, vars(b_logor_1) nolog ml i2 print(bscov)
mvmeta b V, vars(b_logor_1) nolog reml i2 print(bscov)


/* Univariate outcome meta-regression models */

** Fixed effect meta-regression

sjlog using gsemmetan6, replace
use cholesterol, clear
generate weight = 1/varlogor
sem (logor <- cholreduc) [iweight=weight], variance(e.logor@1) ///
	nodescribe nolog nocnsreport
sjlog close, replace


sjlog using gsemmetan7, replace
generate double invselogor = 1/sqrt(varlogor)
generate double logortr = logor*invselogor
generate double cholreductr = cholreduc*invselogor
sem (logortr <- cholreductr invselogor, noconstant), ///
	nodescribe nolog nocnsreport
local Q = _b[var(e.logortr):_cons]*e(N)
local df = e(N) - 2
display "Het. test statistic = " `Q'
display "Het. test p-value = " chi2tail(`df', `Q')
sjlog close, replace


** Random effects meta-regression

sjlog using gsemmetan8, replace
metareg logor cholreduc, wsse(selogor)
sjlog close, replace


sjlog using gsemmetan9, replace
gsem (logor <- c.selogor#c.M@1 cholreduc), variance(M@1) intpoints(6) nolog nocnsreport
sjlog close, replace


sjlog using gsemmetan10, replace
gsem (logortr <- c.invselogor#c.M@1 invselogor cholreductr, noconstant), ///
	variance(e.logortr@1 (M, init(.0097))) latent(M) ///
	nolog nocnsreport
sjlog close, replace


generate double b_logor_1 = logor
generate double V_logor_1_logor_1 = varlogor
mvmeta b V cholreduc, vars(b_logor_1) ml nolog print(bscov)


/* Multivariate outcome meta-analysis with zero within study covariances */

** Fixed effect models

sjlog using gsemmetan11, replace
use telomerase, clear
reshape long y s, i(study) j(outcome)
generate byte y2cons = (outcome == 2)
generate weight = 1/(s^2)
sem (y <- y2cons) [iw=weight], variance(e.y@1) nocaps nodescribe nolog nocnsreport
lincom [y]_cons + [y]y2cons
sjlog close, replace


sjlog using gsemmetan12, replace
generate double invse = 1/s
generate double ytr = y*invse
generate double y2constr = y2cons*invse
quietly sem (ytr <- y2constr invse, noconstant), nocaps nodescribe nolog nocnsreport
local Q = _b[var(e.ytr):_cons]*e(N)
local df = e(N) - 2
display "Het. test statistic = " `Q'
display "Het. test p-value = " chi2tail(`df', `Q')
sjlog close, replace


use telomerase, clear
generate double S11=s1^2
generate double S22=s2^2
mvmeta y S, wscorr(0) fixed


** Random effects models

sjlog using gsemmetan13, replace
quietly use telomerase, clear
gsem (y1 <- c.s1#c.M1@1) (y2 <- c.s2#c.M2@1), ///
	covariance(M1@1 M2@1 M1*M2@0 e.y1*e.y2) latent(M1 M2) ///
	intmethod(laplace) nolog nocnsreport
sjlog close, replace


generate double y1tr = y1/s1
generate double invs1 = 1/s1
generate double y2tr = y2/s2
generate double invs2 = 1/s2
gsem (y1tr <- c.invs1#c.M1@1 invs1, noconstant) ///
	(y2tr <- c.invs2#c.M2@1 invs2, noconstant), ///
	covariance(e.y1tr@1 e.y2tr@1 e.y1tr*e.y2tr@0) ///
	latent(M1 M2) nolog nocnsreport

generate double S11=s1^2
generate double S22=s2^2
mvmeta y S, vars(y1 y2) wscorr(0) print(bscov) nolog ml


/* Multivariate outcome meta-analysis with non-zero within study covariances */

** Fixed effect models

sjlog using gsemmetan14, replace
use fscstage1, clear
forvalues i=1/31 {
	matrix V`i' = (V_Ifg_2_Ifg_2[`i'], V_Ifg_2_Ifg_5[`i'] \ ///
		V_Ifg_2_Ifg_5[`i'], V_Ifg_5_Ifg_5[`i'])
	mata V`i' = st_matrix("V`i'")
	mata invV`i' = invsym(V`i')
	mata W`i' = cholesky(invV`i')
	matrix y`i' = (b_Ifg_2[`i'] \ b_Ifg_5[`i'])
	mata y`i' = st_matrix("y`i'")
	mata ystar`i' = W`i''*y`i'
	mata x`i' = I(2)
	mata xstar`i' = W`i''*x`i'
	if `i' == 1 {
	    mata ystarstack = ystar1
	    mata xstarstack = xstar1
	}
	else {
	    mata ystarstack = (ystarstack \ ystar`i')
	    mata xstarstack = (xstarstack \ xstar`i')
	}
}
clear
getmata ystarstack (xstarstack*)=xstarstack
sem (ystarstack <- xstarstack1 xstarstack2, noconstant), ///
	variance(e.ystarstack@1) nocapslatent nolog nocnsreport nodescribe
lincom _b[ystarstack:xstarstack1], eform
lincom _b[ystarstack:xstarstack2], eform
sjlog close, replace


sjlog using gsemmetan15, replace
quietly sem (ystarstack <- xstarstack1 xstarstack2, noconstant), nocapslatent
display "var(e.ystarstack) = " _b[var(e.ystarstack):_cons]
local Q = _b[var(e.ystarstack):_cons]*e(N)
local df = e(N) - 2
display "Het. test statistic = " `Q'
display "Het. test p-value = " chi2tail(`df', `Q')
sjlog close, replace


generate study = round(_n/2)
generate outcome = (mod(_n,2)==0) + 1
reshape wide ystarstack xstarstack1 xstarstack2, i(study) j(outcome)
assert xstarstack12 == 0
sem (ystarstack1 <- xstarstack11 xstarstack21@c1) ///
	(ystarstack2 <- xstarstack22@c1), noconstant ///
	variance(e.ystarstack1@1 e.ystarstack2@1) nocaps nolog nocnsreport nodescribe


preserve
use fscstage1, clear
mvmeta b V, vars(b_Ifg_2 b_Ifg_3) fixed
restore


** Random effects models

sjlog using gsemmetan16, replace
quietly reshape long
gsem (ystarstack <- c.xstarstack1#c.M1[study]@1 ///
c.xstarstack2#c.M2[study]@1 xstarstack1 xstarstack2, nocons), ///
latent(M1 M2) nocnsreport nolog ///
cov(e.ystarstack@1 M1[study]*M2[study])
sjlog close, replace


sjlog using gsemmetan17, replace
quietly reshape wide
gsem (ystarstack1 <- c.xstarstack11#c.M1@1 c.xstarstack21#c.M2@1 ///
xstarstack11 xstarstack21@c1, nocons) ///
(ystarstack2 <- c.xstarstack22#c.M2@1 ///
xstarstack22@c1, nocons), ///
cov(e.ystarstack1@1 e.ystarstack2@1) latent(M1 M2) ///
collinear nocnsreport nolog
display "corr(M1,M2)=", _b[cov(M2,M1):_cons]/sqrt(_b[var(M1):_cons]*_b[var(M2):_cons])
sjlog close, replace


sjlog using gsemmetan18, replace
quietly use fscstage1, clear
mvmeta b V, vars(b_Ifg_2 b_Ifg_5) print(bscov) ml nolog
sjlog close, replace
