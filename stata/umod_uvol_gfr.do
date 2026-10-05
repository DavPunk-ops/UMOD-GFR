clear
use "C:\Users\dajs\OneDrive - HOPITAUX UNIVERSITAIRES DE GENEVE\recherche\UMOD copeptin SKIPOGH\stata\cross-sectional\data_all_c_25112015.dta"

// Objectives
// 1) Relationship between 24h urinary UMOD excretion and 24h urine volume
// 2) Effect of GFR on relationship 1 (continuous interaction and GFR strata)
// 3) Relationship between 24h urinary UMOD concentration and 24h urine volume
// 4) Effect of GFR on relationship 3 (continuous interaction and GFR strata)

// Covariates (same adjustment as the original analysis)
global adj "age10 sex01 diabetes ckd_epi10 kid_vol100 d_diur center"
// Same covariates without GFR, for models where GFR is already in the interaction term
global adjnogfr "age10 sex01 diabetes kid_vol100 d_diur center"






**********************
*** Main variables ***
**********************
// UMOD excretion (mg/24h)
rename umod_u24_ug UMOD24h
replace UMOD24h=UMOD24h/1000
*kdensity UMOD24h, normal
*gladder UMOD24h
gen UMOD24hsqrt=sqrt(UMOD24h)
*graph box UMOD24hsqrt

// UMOD concentration (ug/mL = mg/L)
gen UMODconc=umod_u24_ugml
*kdensity UMODconc, normal
*gladder UMODconc
gen UMODconclog=log(UMODconc)
*graph box UMODconclog
// NB: concentration = excretion / volume, so log(concentration) = log(excretion) - log(UVOL).
// A negative association between concentration and UVOL is partly expected by construction.

// UVOL (mL/24h)
rename u24_ml UVOL
*kdensity UVOL, normal
*gladder UVOL
gen UVOLlog=log(UVOL)
*graph box UVOLlog







***********************
*** Other variables ***
***********************
gen age10=age/10

gen ckd_epi10=ckd_epi/10

gen ckd=.
replace ckd=0 if ckd_epi>=60
replace ckd=1 if ckd_epi<60
tab ckd

// GFR strata: <60, 60-89, >=90 mL/min/1.73m2
gen gfrcat=.
replace gfrcat=1 if ckd_epi<60
replace gfrcat=2 if ckd_epi>=60 & ckd_epi<90
replace gfrcat=3 if ckd_epi>=90 & ckd_epi!=.
label define gfrcat 1 "eGFR <60" 2 "eGFR 60-89" 3 "eGFR >=90"
label values gfrcat gfrcat
tab gfrcat

gen kid_vol100=kid_vol/100









*************************
*** Defining outliers ***
*************************
// 1st and 99th percentiles of creatininuria, UMOD and UVOL
gen outliers=0

gen ucrt_mgkg24sqrt=sqrt(ucrt_mgkg24)
sum ucrt_mgkg24sqrt, det
replace outliers=1 if (ucrt_mgkg24sqrt<=r(p1) | ucrt_mgkg24sqrt>=r(p99)) & ucrt_mgkg24sqrt!=.

sum UMOD24hsqrt, det
replace outliers=1 if (UMOD24hsqrt<=r(p1) | UMOD24hsqrt>=r(p99)) & UMOD24hsqrt!=.

sum UVOLlog, det
replace outliers=1 if (UVOLlog<=r(p1) | UVOLlog>=r(p99)) & UVOLlog!=.









*****************
*** Flowchart ***
*****************
count if sex01!=.
count if UMOD24h!=.
count if UMODconc!=.
count if UVOL!=.
count if UMOD24h!=. & UVOL!=.











*******************************
*** Drop missing covariates ***
*******************************
drop if UMOD24h==. | UMODconc==. | UVOL==. | age==. | sex01==. | diabetes==. | ckd_epi==. | d_diur==. | kid_vol==. | center==.
count









*********************
*** Drop outliers ***
*********************
drop if outliers==1
count

// Range of UVOLlog used for predictions (1st to 99th percentile)
_pctile UVOLlog, p(1 99)
global uvlo=round(r(r1), 0.1)
global uvhi=round(r(r2), 0.1)
display "UVOLlog range for predictions: ${uvlo} to ${uvhi}"









*******************************
*** General characteristics ***
*******************************
sum UMOD24h, det
sum UMODconc, det
sum UVOL, det
sum age, det
tab sex01
tab gfrcat






*************************************
*** Table 1: by GFR strata ***
*************************************
foreach var in age bmi sbp dbp ckd_epi kid_vol UVOL UMOD24h UMODconc {
	kwallis `var', by(gfrcat)
	tabstat `var', by(gfrcat) stats(n mean sd p50 p25 p75)
}

foreach var in sex01 diabetes hta d_diur {
	tab `var' gfrcat, col chi2
}










***************************************************************
*** Objective 1: UMOD excretion (mg/24h) and UVOL (mL/24h) ***
***************************************************************
// Raw relationship
twoway (scatter UMOD24hsqrt UVOLlog, msize(small)) (lowess UMOD24hsqrt UVOLlog) (lfit UMOD24hsqrt UVOLlog)
twoway (lfitci UMOD24hsqrt UVOLlog) (scatter UMOD24hsqrt UVOLlog, msize(vsmall)), name(fig1a, replace)
*graph save "figure 1a.gph", replace

// Univariate
mixed UMOD24hsqrt UVOLlog center || treenbr:

// Fully adjusted
mixed UMOD24hsqrt UVOLlog $adj || treenbr:

margins, at(UVOLlog=(${uvlo}(0.2)${uvhi}))
marginsplot, recast(line) recastci(rarea) name(fig1b, replace)
*graph save "figure 1b.gph", replace









*******************************************************
*** Objective 2: Effect of GFR on UMOD excretion-UVOL ***
*******************************************************
// --- 2a. Continuous interaction UVOL x eGFR ---
// Univariate
mixed UMOD24hsqrt c.UVOLlog##c.ckd_epi10 center || treenbr:

// Fully adjusted
mixed UMOD24hsqrt c.UVOLlog##c.ckd_epi10 $adjnogfr || treenbr:
// p-interaction = p-value of c.UVOLlog#c.ckd_epi10

// Predicted UMOD excretion according to UVOL at eGFR 45, 75 and 105
margins, at(UVOLlog=(${uvlo}(0.2)${uvhi}) ckd_epi10=(4.5 7.5 10.5))
marginsplot, xdimension(at(UVOLlog)) noci name(fig2a, replace)
*graph save "figure 2a.gph", replace

// Slope of UMOD excretion on UVOL according to eGFR
sum ckd_epi10
margins, dydx(UVOLlog) at(ckd_epi10=(4(1)13))
marginsplot, yline(0) name(fig2b, replace)
*graph save "figure 2b.gph", replace


// --- 2b. GFR strata (<60, 60-89, >=90) ---
// Raw relationship by stratum
twoway (lfit UMOD24hsqrt UVOLlog if gfrcat==1) (lfit UMOD24hsqrt UVOLlog if gfrcat==2) (lfit UMOD24hsqrt UVOLlog if gfrcat==3), legend(order(1 "eGFR <60" 2 "eGFR 60-89" 3 "eGFR >=90")) name(fig2c, replace)

// Interaction model with categorical GFR
mixed UMOD24hsqrt c.UVOLlog##i.gfrcat center || treenbr:
testparm c.UVOLlog#i.gfrcat

mixed UMOD24hsqrt c.UVOLlog##i.gfrcat $adjnogfr || treenbr:
testparm c.UVOLlog#i.gfrcat // global p-interaction

// Slope of UMOD excretion on UVOL in each stratum
margins gfrcat, dydx(UVOLlog)

// Predicted UMOD excretion according to UVOL in each stratum
margins gfrcat, at(UVOLlog=(${uvlo}(0.2)${uvhi}))
marginsplot, xdimension(at(UVOLlog)) name(fig2d, replace)
*graph save "figure 2d.gph", replace

// Stratified models (fully adjusted)
forvalues g=1/3 {
	display "*** gfrcat = `g' ***"
	mixed UMOD24hsqrt UVOLlog $adj if gfrcat==`g' || treenbr:
}









**********************************************************************
*** Objective 3: UMOD concentration (ug/mL) and UVOL (mL/24h) ***
**********************************************************************
// Raw relationship
twoway (scatter UMODconclog UVOLlog, msize(small)) (lowess UMODconclog UVOLlog) (lfit UMODconclog UVOLlog)
twoway (lfitci UMODconclog UVOLlog) (scatter UMODconclog UVOLlog, msize(vsmall)), name(fig3a, replace)
*graph save "figure 3a.gph", replace

// Univariate
mixed UMODconclog UVOLlog center || treenbr:

// Fully adjusted
mixed UMODconclog UVOLlog $adj || treenbr:

margins, at(UVOLlog=(${uvlo}(0.2)${uvhi}))
marginsplot, recast(line) recastci(rarea) name(fig3b, replace)
*graph save "figure 3b.gph", replace









************************************************************
*** Objective 4: Effect of GFR on UMOD concentration-UVOL ***
************************************************************
// --- 4a. Continuous interaction UVOL x eGFR ---
// Univariate
mixed UMODconclog c.UVOLlog##c.ckd_epi10 center || treenbr:

// Fully adjusted
mixed UMODconclog c.UVOLlog##c.ckd_epi10 $adjnogfr || treenbr:
// p-interaction = p-value of c.UVOLlog#c.ckd_epi10

// Predicted UMOD concentration according to UVOL at eGFR 45, 75 and 105
margins, at(UVOLlog=(${uvlo}(0.2)${uvhi}) ckd_epi10=(4.5 7.5 10.5))
marginsplot, xdimension(at(UVOLlog)) noci name(fig4a, replace)
*graph save "figure 4a.gph", replace

// Slope of UMOD concentration on UVOL according to eGFR
margins, dydx(UVOLlog) at(ckd_epi10=(4(1)13))
marginsplot, yline(0) name(fig4b, replace)
*graph save "figure 4b.gph", replace


// --- 4b. GFR strata (<60, 60-89, >=90) ---
// Raw relationship by stratum
twoway (lfit UMODconclog UVOLlog if gfrcat==1) (lfit UMODconclog UVOLlog if gfrcat==2) (lfit UMODconclog UVOLlog if gfrcat==3), legend(order(1 "eGFR <60" 2 "eGFR 60-89" 3 "eGFR >=90")) name(fig4c, replace)

// Interaction model with categorical GFR
mixed UMODconclog c.UVOLlog##i.gfrcat center || treenbr:
testparm c.UVOLlog#i.gfrcat

mixed UMODconclog c.UVOLlog##i.gfrcat $adjnogfr || treenbr:
testparm c.UVOLlog#i.gfrcat // global p-interaction

// Slope of UMOD concentration on UVOL in each stratum
margins gfrcat, dydx(UVOLlog)

// Predicted UMOD concentration according to UVOL in each stratum
margins gfrcat, at(UVOLlog=(${uvlo}(0.2)${uvhi}))
marginsplot, xdimension(at(UVOLlog)) name(fig4d, replace)
*graph save "figure 4d.gph", replace

// Stratified models (fully adjusted)
forvalues g=1/3 {
	display "*** gfrcat = `g' ***"
	mixed UMODconclog UVOLlog $adj if gfrcat==`g' || treenbr:
}









*******************************************************
*** SUMMARY: answers to the 4 questions (adjusted) ***
*******************************************************
// All models fully adjusted (mixed models, random intercept for family).
// Graphs are back-transformed to original units (mg/24h, ug/mL, mL/24h):
// predictions are squared (sqrt scale) or exponentiated (log scale), i.e. medians rather than means.

capture program drop showres
program define showres
	args label expr
	quietly lincom `expr'
	display as text "`label'" as result %6.2f r(estimate) as text " (95% CI " as result %6.2f r(lb) as text " to " as result %6.2f r(ub) as text "), p = " as result %5.3f r(p)
end

// Urine volume axis labels (log scale shown as mL/24h)
global uvlab `"6.397 "600" 6.908 "1000" 7.313 "1500" 7.601 "2000" 8.006 "3000""'


// --- Question 1: UMOD excretion vs urine volume ---
quietly mixed UMOD24hsqrt UVOLlog $adj || treenbr:
display _newline as text "{hline 70}" _newline "Q1. UMOD excretion (sqrt mg/24h) vs urine volume (log mL/24h)" _newline "{hline 70}"
showres "Slope: " "UVOLlog"

quietly margins, at(UVOLlog=(${uvlo}(0.1)${uvhi})) expression(predict(xb)^2)
marginsplot, recast(line) recastci(rarea) ciopts(fcolor(%30) lwidth(none)) ///
	xlabel($uvlab) xtitle("Urine volume (mL/24h)") ytitle("UMOD excretion (mg/24h)") ///
	title("Q1. UMOD excretion vs urine volume") subtitle("Adjusted predictions") name(Q1, replace)
*graph export "Q1.png", replace


// --- Question 2: effect of GFR on UMOD excretion vs urine volume ---
display _newline as text "{hline 70}" _newline "Q2. Effect of eGFR on UMOD excretion vs urine volume" _newline "{hline 70}"
quietly mixed UMOD24hsqrt c.UVOLlog##c.ckd_epi10 $adjnogfr || treenbr:
quietly test c.UVOLlog#c.ckd_epi10
display as text "p-interaction UVOL x eGFR (continuous): " as result %5.3f r(p)

quietly mixed UMOD24hsqrt c.UVOLlog##i.gfrcat $adjnogfr || treenbr:
quietly testparm c.UVOLlog#i.gfrcat
display as text "p-interaction UVOL x eGFR strata (global): " as result %5.3f r(p)
showres "Slope eGFR <60:    " "UVOLlog"
showres "Slope eGFR 60-89:  " "UVOLlog + 2.gfrcat#c.UVOLlog"
showres "Slope eGFR >=90:   " "UVOLlog + 3.gfrcat#c.UVOLlog"

quietly margins gfrcat, at(UVOLlog=(${uvlo}(0.1)${uvhi})) expression(predict(xb)^2)
marginsplot, xdimension(at(UVOLlog)) recast(line) recastci(rarea) ciopts(fcolor(%20) lwidth(none)) ///
	xlabel($uvlab) xtitle("Urine volume (mL/24h)") ytitle("UMOD excretion (mg/24h)") ///
	title("Q2. UMOD excretion vs urine volume by eGFR") subtitle("Adjusted predictions") name(Q2, replace)
*graph export "Q2.png", replace


// --- Question 3: UMOD concentration vs urine volume ---
quietly mixed UMODconclog UVOLlog $adj || treenbr:
display _newline as text "{hline 70}" _newline "Q3. UMOD concentration (log ug/mL) vs urine volume (log mL/24h)" _newline "{hline 70}"
showres "Slope: " "UVOLlog"
display as text "(slope = -1 would mean pure dilution, i.e. UMOD excretion independent of volume)"

quietly margins, at(UVOLlog=(${uvlo}(0.1)${uvhi})) expression(exp(predict(xb)))
marginsplot, recast(line) recastci(rarea) ciopts(fcolor(%30) lwidth(none)) ///
	xlabel($uvlab) xtitle("Urine volume (mL/24h)") ytitle("UMOD concentration (ug/mL)") ///
	title("Q3. UMOD concentration vs urine volume") subtitle("Adjusted predictions") name(Q3, replace)
*graph export "Q3.png", replace


// --- Question 4: effect of GFR on UMOD concentration vs urine volume ---
display _newline as text "{hline 70}" _newline "Q4. Effect of eGFR on UMOD concentration vs urine volume" _newline "{hline 70}"
quietly mixed UMODconclog c.UVOLlog##c.ckd_epi10 $adjnogfr || treenbr:
quietly test c.UVOLlog#c.ckd_epi10
display as text "p-interaction UVOL x eGFR (continuous): " as result %5.3f r(p)

quietly mixed UMODconclog c.UVOLlog##i.gfrcat $adjnogfr || treenbr:
quietly testparm c.UVOLlog#i.gfrcat
display as text "p-interaction UVOL x eGFR strata (global): " as result %5.3f r(p)
showres "Slope eGFR <60:    " "UVOLlog"
showres "Slope eGFR 60-89:  " "UVOLlog + 2.gfrcat#c.UVOLlog"
showres "Slope eGFR >=90:   " "UVOLlog + 3.gfrcat#c.UVOLlog"

quietly margins gfrcat, at(UVOLlog=(${uvlo}(0.1)${uvhi})) expression(exp(predict(xb)))
marginsplot, xdimension(at(UVOLlog)) recast(line) recastci(rarea) ciopts(fcolor(%20) lwidth(none)) ///
	xlabel($uvlab) xtitle("Urine volume (mL/24h)") ytitle("UMOD concentration (ug/mL)") ///
	title("Q4. UMOD concentration vs urine volume by eGFR") subtitle("Adjusted predictions") name(Q4, replace)
*graph export "Q4.png", replace









*************************************************************
*** ALTERNATIVE GFR MODELLING: tertiles and heatmap (adjusted) ***
*************************************************************
// GFR tertiles (defined in the analysed sample)
xtile gfrtert=ckd_epi, nquantiles(3)
label define gfrtert 1 "eGFR T1" 2 "eGFR T2" 3 "eGFR T3"
label values gfrtert gfrtert
tabstat ckd_epi, by(gfrtert) stats(n min p50 max)

// Range of eGFR (1st to 99th percentile) for the heatmaps
_pctile ckd_epi10, p(1 99)
global gfrlo=round(r(r1), 0.5)
global gfrhi=round(r(r2), 0.5)


// --- Q2 with GFR tertiles ---
display _newline as text "{hline 70}" _newline "Q2 (tertiles). Effect of eGFR on UMOD excretion vs urine volume" _newline "{hline 70}"
quietly mixed UMOD24hsqrt c.UVOLlog##i.gfrtert $adjnogfr || treenbr:
quietly testparm c.UVOLlog#i.gfrtert
display as text "p-interaction UVOL x eGFR tertiles (global): " as result %5.3f r(p)
showres "Slope eGFR T1: " "UVOLlog"
showres "Slope eGFR T2: " "UVOLlog + 2.gfrtert#c.UVOLlog"
showres "Slope eGFR T3: " "UVOLlog + 3.gfrtert#c.UVOLlog"

quietly margins gfrtert, at(UVOLlog=(${uvlo}(0.1)${uvhi})) expression(predict(xb)^2)
marginsplot, xdimension(at(UVOLlog)) recast(line) recastci(rarea) ciopts(fcolor(%20) lwidth(none)) ///
	xlabel($uvlab) xtitle("Urine volume (mL/24h)") ytitle("UMOD excretion (mg/24h)") ///
	title("Q2. UMOD excretion vs urine volume by eGFR tertile") subtitle("Adjusted predictions") name(Q2tert, replace)
*graph export "Q2_tertiles.png", replace


// --- Q4 with GFR tertiles ---
display _newline as text "{hline 70}" _newline "Q4 (tertiles). Effect of eGFR on UMOD concentration vs urine volume" _newline "{hline 70}"
quietly mixed UMODconclog c.UVOLlog##i.gfrtert $adjnogfr || treenbr:
quietly testparm c.UVOLlog#i.gfrtert
display as text "p-interaction UVOL x eGFR tertiles (global): " as result %5.3f r(p)
showres "Slope eGFR T1: " "UVOLlog"
showres "Slope eGFR T2: " "UVOLlog + 2.gfrtert#c.UVOLlog"
showres "Slope eGFR T3: " "UVOLlog + 3.gfrtert#c.UVOLlog"

quietly margins gfrtert, at(UVOLlog=(${uvlo}(0.1)${uvhi})) expression(exp(predict(xb)))
marginsplot, xdimension(at(UVOLlog)) recast(line) recastci(rarea) ciopts(fcolor(%20) lwidth(none)) ///
	xlabel($uvlab) xtitle("Urine volume (mL/24h)") ytitle("UMOD concentration (ug/mL)") ///
	title("Q4. UMOD concentration vs urine volume by eGFR tertile") subtitle("Adjusted predictions") name(Q4tert, replace)
*graph export "Q4_tertiles.png", replace


// --- Heatmaps: continuous UVOL x continuous eGFR interaction ---
// Predictions from the fully adjusted continuous interaction model,
// on a grid of UVOL and eGFR, with other covariates fixed at their sample mean.
capture program drop gfrheatmap
program define gfrheatmap
	args outcome backtransform ztitle gname gtitle
	quietly mixed `outcome' c.UVOLlog##c.ckd_epi10 $adjnogfr || treenbr:
	foreach v of global adjnogfr {
		quietly sum `v'
		local m_`v'=r(mean)
	}
	preserve
	clear
	local nuv=round((${uvhi}-${uvlo})/0.05)+1
	local ngfr=round((${gfrhi}-${gfrlo})/0.25)+1
	quietly set obs `=`nuv'*`ngfr''
	gen UVOLlog=${uvlo}+mod(_n-1, `nuv')*0.05
	gen ckd_epi10=${gfrlo}+floor((_n-1)/`nuv')*0.25
	foreach v of global adjnogfr {
		gen `v'=`m_`v''
	}
	quietly predict xb, xb
	gen pred=`backtransform'
	gen eGFR=ckd_epi10*10
	twoway (contour pred eGFR UVOLlog, levels(15)), ///
		xlabel($uvlab) xtitle("Urine volume (mL/24h)") ytitle("eGFR (mL/min/1.73m{superscript:2})") ///
		clegend(title("`ztitle'", size(small))) title("`gtitle'") subtitle("Adjusted predictions") name(`gname', replace)
	restore
end

gfrheatmap UMOD24hsqrt "xb^2" "UMOD excretion (mg/24h)" Q2heat "Q2. UMOD excretion by urine volume and eGFR"
*graph export "Q2_heatmap.png", replace

gfrheatmap UMODconclog "exp(xb)" "UMOD concentration (ug/mL)" Q4heat "Q4. UMOD concentration by urine volume and eGFR"
*graph export "Q4_heatmap.png", replace
