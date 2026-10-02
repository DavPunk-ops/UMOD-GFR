clear
use "C:\Users\dajs\OneDrive - HOPITAUX UNIVERSITAIRES DE GENEVE\recherche\UMOD copeptin SKIPOGH\stata\cross-sectional\data_all_c_25112015.dta"






**********************
*** Main variables ***
**********************
// POSM (mOsm/kg) (calculated)
gen POSM=2*(na+k)+ure+glu
*kdensity POSM, normal
*graph box POSM
sum osmc osmo2 osmo1 osmo osmolality POSM // Comparison between plasma osmolarity variables

// UOSM (mOsm/kg) (calculated)
rename uosm_u24 UOSMcalculated
*kdensity UOSMcalculated, normal
*gladder UOSMcalculated
gen UOSMcalculatedsqrt=sqrt(UOSMcalculated)
*graph box UOSMcalculatedsqrt

// UOSM (mOsm/kg) (measured)
rename uosm_u24_od UOSMmeasured
*kdensity UOSMmeasured, normal
*gladder UOSMmeasured
gen UOSMmeasuredsqrt=sqrt(UOSMmeasured)
*graph box UOSMmeasuredsqrt

// Comparison between calculated and measured UOSM
sum UOSMcalculated UOSMmeasured
*twoway (scatter UOSMcalculated UOSMmeasured) (lowess UOSMcalculated UOSMmeasured) (lfit UOSMcalculated UOSMmeasured)

sum UOSMcalculatedsqrt UOSMmeasuredsqrt
*twoway (scatter UOSMcalculatedsqrt UOSMmeasuredsqrt) (lowess UOSMcalculatedsqrt UOSMmeasuredsqrt) (lfit UOSMcalculatedsqrt UOSMmeasuredsqrt)

// Copeptin (pmol/L)
*kdensity copeptin, normal
*gladder copeptin
gen copeptinlog=log(copeptin)
*graph box copeptinlog

// UMOD (mg/24h)
rename umod_u24_ug UMOD24h
replace UMOD24h=UMOD24h/1000
*kdensity UMOD24h, normal
*gladder UMOD24h
gen UMOD24hsqrt=sqrt(UMOD24h)
*graph box UMOD24hsqrt

// UVOL (mL/24h)
rename u24_ml UVOL
*kdensity UVOL, normal
*gladder UVOL
gen UVOLlog=log(UVOL)
*graph box UVOLlog

// UOSM (mOsm/24h) (calculated)
gen UOSM24calculated=(UOSMcalculated*UVOL)/1000
*kdensity UOSM24calculated, normal
*gladder UOSM24calculated
gen UOSM24calculatedsqrt=sqrt(UOSM24calculated)
*graph box UOSM24calculatedsqrt

// UOSM (mOsm/24h) (measured)
gen UOSM24measured=(UOSMmeasured*UVOL)/1000
*kdensity UOSM24measured
*gladder UOSM24measured
gen UOSM24measuredsqrt=sqrt(UOSM24measured)
*graph box UOSM24measuredsqrt







***********************
*** Other variables ***
***********************
*kdensity age, normal
*gladder age
gen agesqrt=sqrt(age)
*graph box agesqrt

gen age10=age/10

*kdensity ckd_epi, normal
*graph box ckd_epi

gen ckd_epi10=ckd_epi/10

gen ckd=.
replace ckd=0 if ckd_epi>=60
replace ckd=1 if ckd_epi<60
tab ckd

*kdensity kid_length, normal
*gladder kid_length
gen kid_lengthlog=log(kid_length)
*graph box kid_lengthlog

*kdensity kid_vol, normal
*gladder kid_vol
gen kid_vollog=log(kid_vol)
*graph box kid_vollog

gen kid_vol100=kid_vol/100

*kdensity bmi, normal
*graph box bmi
*gladder bmi
gen bmisqrt=sqrt(bmi)
*graph box bmisqrt

*kdensity sbp, normal
*graph box sbp

*kdensity dbp, normal
*graph box dbp









*************************
*** Defining outliers ***
*************************
// 99th percentiles of creatininuria, UMOD, copeptin, UVOL, UOSM, POSM
gen outliers=0

*gladder ucrt_mgkg24
gen ucrt_mgkg24sqrt=sqrt(ucrt_mgkg24)
sum ucrt_mgkg24sqrt, det
replace outliers=1 if (ucrt_mgkg24sqrt<=r(p1) | ucrt_mgkg24sqrt>=r(p99)) & ucrt_mgkg24sqrt!=.

sum UMOD24hsqrt, det
replace outliers=1 if (UMOD24hsqrt<=r(p1) | UMOD24hsqrt>=r(p99)) & UMOD24hsqrt!=.

sum copeptinlog, det
replace outliers=1 if (copeptinlog<=r(p1) | copeptinlog>=r(p99)) & copeptinlog!=.

sum UVOLlog, det
replace outliers=1 if (UVOLlog<=r(p1) | UVOLlog>=r(p99)) & UVOLlog!=.

sum UOSMmeasuredsqrt, det
replace outliers=1 if (UOSMmeasuredsqrt<=r(p1) | UOSMmeasuredsqrt>=r(p99)) & UOSMmeasuredsqrt!=.









*****************
*** Flowchart ***
*****************
count if sex01!=.
count if UMOD24h!=.
count if copeptin!=.
count if UMOD24h!=. & copeptin!=.











*******************************
*** Drop missing covariates ***
*******************************
drop if UMOD24h==. | copeptin==. | UOSMmeasured==. | UVOL==. | age==. | sex01==. | diabetes==. | ckd_epi==. | d_diur==. | kid_vol==. | center==.
count if UMOD24h!=. & copeptin!=.








*********************
*** Drop outliers ***
*********************
drop if outliers==1
count if UMOD24h!=. & copeptin!=.











*************************
*** Defining tertiles ***
*************************
// UMOD
xtile UMODtertiles=UMOD24hsqrt, nquantiles(3)
tab UMODtertiles
by UMODtertiles, sort: egen UMODmedian=median(UMOD24hsqrt)
tab UMODmedian UMODtertiles // 4.90; 6.35 and 7.81

// UVOL
xtile UVOLtertiles=UVOLlog, nquantiles(3)
tab UVOLtertiles
by UVOLtertiles, sort: egen UVOLmedian=median(UVOLlog)
tab UVOLmedian UVOLtertiles // 6.95; 7.37 and 7.78








*******************************
*** General characteristics ***
*******************************
sum UMOD24h, det
sum age, det
tab sex01
tab ckd






***************
*** Table 1 ***
***************
oneway UMOD24h UMODtertiles // Bartlett rejeté
kwallis UMOD24h, by(UMODtertiles)
sum UMOD24h if UMODtertiles==1, det
sum UMOD24h if UMODtertiles==2, det
sum UMOD24h if UMODtertiles==3, det

oneway age UMODtertiles // Bartlett rejeté
kwallis age, by(UMODtertiles)
sum age if UMODtertiles==1, det
sum age if UMODtertiles==2, det
sum age if UMODtertiles==3, det

tab sex01 UMODtertiles, col chi2

oneway bmi UMODtertiles // Bartlett OK
sum bmi if UMODtertiles==1, det
sum bmi if UMODtertiles==2, det
sum bmi if UMODtertiles==3, det

tab diabetes UMODtertiles, col chi2

tab hta UMODtertiles, col chi2

tab d_diur UMODtertiles, col chi2

oneway sbp UMODtertiles // Bartlett rejeté
kwallis sbp, by(UMODtertiles)
sum sbp if UMODtertiles==1, det
sum sbp if UMODtertiles==2, det
sum sbp if UMODtertiles==3, det

oneway dbp UMODtertiles // Bartlett OK
sum dbp if UMODtertiles==1, det
sum dbp if UMODtertiles==2, det
sum dbp if UMODtertiles==3, det

oneway ckd_epi UMODtertiles // Bartlett rejeté
kwallis ckd_epi, by(UMODtertiles)
sum ckd_epi if UMODtertiles==1, det
sum ckd_epi if UMODtertiles==2, det
sum ckd_epi if UMODtertiles==3, det

tab ckd UMODtertiles, col chi2

oneway POSM UMODtertiles // Bartlett OK
sum POSM if UMODtertiles==1, det
sum POSM if UMODtertiles==2, det
sum POSM if UMODtertiles==3, det

oneway UOSMmeasured UMODtertiles // Bartlett OK
sum UOSMmeasured if UMODtertiles==1, det
sum UOSMmeasured if UMODtertiles==2, det
sum UOSMmeasured if UMODtertiles==3, det

oneway UVOL UMODtertiles // Bartlett OK
sum UVOL if UMODtertiles==1, det
sum UVOL if UMODtertiles==2, det
sum UVOL if UMODtertiles==3, det

oneway copeptin UMODtertiles // Bartlett rejeté
kwallis copeptin, by(UMODtertiles)
sum copeptin if UMODtertiles==1, det
sum copeptin if UMODtertiles==2, det
sum copeptin if UMODtertiles==3, det

oneway kid_length UMODtertiles // Bartlett rejeté
kwallis kid_length, by(UMODtertiles)
sum kid_length if UMODtertiles==1, det
sum kid_length if UMODtertiles==2, det
sum kid_length if UMODtertiles==3, det

oneway kid_vol UMODtertiles // Bartlett OK
sum kid_vol if UMODtertiles==1, det
sum kid_vol if UMODtertiles==2, det
sum kid_vol if UMODtertiles==3, det










************************
*** Basic assumption ***
************************
// Copeptin and POSM ==> OK
twoway (scatter copeptinlog POSM) (lowess copeptinlog POSM) (lfit copeptinlog POSM) (qfit copeptinlog POSM)
twoway (lfitci copeptinlog POSM)
*graph save "supplementary figure S2a.gph", replace
regress copeptinlog POSM

// Copeptin and UOSM ==> OK
twoway (scatter UOSMmeasuredsqrt copeptinlog) (lowess UOSMmeasuredsqrt copeptinlog) (lfit UOSMmeasuredsqrt copeptinlog)
twoway (lfitci UOSMmeasuredsqrt copeptinlog)
*graph save "supplementary figure S2b.gph", replace
regress UOSMmeasuredsqrt copeptinlog

// Copeptin and UVOL ==> OK
twoway (scatter UVOLlog copeptinlog) (lowess UVOLlog copeptinlog) (lfit UVOLlog copeptinlog)
twoway (lfitci UVOLlog copeptinlog)
*graph save "supplementary figure S2c.gph", replace
regress UVOLlog copeptinlog

// Combine graphs
*graph combine "supplementary figure S2a.gph" "supplementary figure S2b.gph" "supplementary figure S2c.gph", row(1) ysize(2) xsize(6) imargin (5 5 5)
*graph save "supplementary figure S2abc.gph", replace









****************************************************
*** Hypothesis 1: Determinants of UMOD excretion ***
****************************************************
// Univariate
mixed UMOD24hsqrt UVOLlog center || treenbr:

mixed UMOD24hsqrt center || treenbr:
predict r, residuals
egen meanUMOD=mean(UMOD24hsqrt)
gen adjustedUMOD=meanUMOD+r
twoway (lfitci adjustedUMOD UVOLlog)
*graph save "figure 1a.gph", replace
drop r meanUMOD adjustedUMOD


// Univariate
mixed UMOD24hsqrt copeptinlog center || treenbr:

mixed UMOD24hsqrt center || treenbr:
predict r, residuals
egen meanUMOD=mean(UMOD24hsqrt)
gen adjustedUMOD=meanUMOD+r
twoway (lfitci adjustedUMOD copeptinlog)
*graph save "figure 1b.gph", replace
drop r meanUMOD adjustedUMOD


// Univariate
mixed UMOD24hsqrt age10 center || treenbr:


// Univariate
mixed UMOD24hsqrt sex01 center || treenbr:


// Univariate
mixed UMOD24hsqrt diabetes center || treenbr:


// Univariate
mixed UMOD24hsqrt ckd_epi10 center || treenbr:


// Univariate
mixed UMOD24hsqrt kid_vol100 center || treenbr:


// Univariate
mixed UMOD24hsqrt d_diur center || treenbr:



// Partially adjusted
mixed UMOD24hsqrt c.UVOLlog##c.copeptinlog center || treenbr:

sum copeptinlog
margins, at(copeptinlog=(0(0.4)2.9) UVOLlog=(6.95 7.37 7.78))
marginsplot, xlabel (0(0.4)2.9)
*graph save "figure 1c.gph", replace


// Combine graphs
*graph combine "figure 1a.gph" "figure 1b.gph" "figure 1c.gph", row(1) ysize(2) xsize(6) imargin (5 5 5)
*graph save "figure 1abc.gph", replace


// Fully adjusted
mixed UMOD24hsqrt c.UVOLlog##c.copeptinlog age10 sex01 diabetes ckd_epi10 kid_vol100 d_diur center || treenbr:

sum copeptinlog
margins, at(copeptinlog=(0(0.4)2.9) UVOLlog=(6.95 7.37 7.78))
marginsplot, xlabel (0(0.4)2.9)
*graph save "figure 2a.gph", replace

sum UVOLlog
margins, at(copeptinlog=(0(0.4)2.9) UVOLlog=(6.2(0.4)8.2)) saving(predictions, replace)
use predictions.dta, clear
describe
twoway (contour _margin _at2 _at1, xlabel(6.2(0.4)8.2) ylabel(0(0.4)2.9) zlabel(#10) levels(500))
*graph save "figure 2b.gph", replace


// Combine graphs
*graph combine "figure 2a.gph" "figure 2b.gph", row(1) ysize(2) xsize(5) imargin (5 5)
*graph save "figure 2ab.gph", replace


// Colinearity of key variables
regress UMOD24hsqrt copeptinlog UVOLlog
estat vif





**********************************************************
*** Hypothesis 2: Effect of UMOD on water reabsorption ***
**********************************************************
// Univariate
mixed UOSMmeasuredsqrt UMOD24hsqrt center || treenbr:

mixed UOSMmeasuredsqrt center || treenbr:
predict r, residuals
egen meanUOSM=mean(UOSMmeasuredsqrt)
gen adjustedUOSM=meanUOSM+r
twoway (lfitci adjustedUOSM UMOD24hsqrt)
drop r meanUOSM adjustedUOSM
*graph save "figure 3a.gph", replace


// Univariate
mixed UOSMmeasuredsqrt copeptinlog center || treenbr:


// Univariate
mixed UOSMmeasuredsqrt UVOLlog center || treenbr:


// Univariate
mixed UOSMmeasuredsqrt age10 center || treenbr:


// Univariate
mixed UOSMmeasuredsqrt sex01 center || treenbr:


// Univariate
mixed UOSMmeasuredsqrt diabetes center || treenbr:


// Univariate
mixed UOSMmeasuredsqrt ckd_epi10 center || treenbr:


// Univariate
mixed UOSMmeasuredsqrt kid_vol100 center || treenbr:


// Univariate
mixed UOSMmeasuredsqrt d_diur center || treenbr:




// Partially adjusted
mixed UOSMmeasuredsqrt c.UMOD24hsqrt##c.copeptinlog center || treenbr:

sum copeptinlog
margins, at(copeptinlog=(0(0.4)2.9) UMOD24hsqrt=(4.90 6.35 7.81))
marginsplot, xlabel (0(0.4)2.9)
*graph save "figure 3b.gph", replace


// Combine graphs
*graph combine "figure 3a.gph" "figure 3b.gph", row(1) ysize(2) xsize(6) imargin (5 5)
*graph save "figure 3ab.gph", replace


// Fully adjsuted
mixed UOSMmeasuredsqrt c.UMOD24hsqrt##c.copeptinlog UVOLlog age10 sex01 diabetes ckd_epi10 kid_vol100 d_diur center || treenbr:

sum copeptinlog
margins, at(copeptinlog=(0(0.4)2.9) UMOD24hsqrt=(4.90 6.35 7.81))
marginsplot, xlabel (0(0.4)2.9)
*graph save "figure 4a.gph", replace

sum UMOD24hsqrt
margins, at(copeptinlog=(0(0.4)2.9) UMOD24hsqrt=(2.3(1)10.3)) saving(predictions, replace)
use predictions.dta, clear
describe
twoway (contour _margin _at1 _at2, xlabel(0(0.4)2.9) ylabel(2.3(1)10.3) zlabel(#10) levels(500))
*graph save "figure 4b.gph", replace


// Combine graphs
*graph combine "figure 4a.gph" "figure 4b.gph", row(1) ysize(2) xsize(5) imargin (5 5)
*graph save "figure 4ab.gph", replace


* Colinearity of key variables *
regress UOSMmeasuredsqrt UMOD24hsqrt copeptinlog UVOLlog
estat vif












****************************
*** Sensitivity analyses ***
****************************
// With outliers: Run code without excluding outliers

* Hypothesis 1 *
mixed UMOD24hsqrt c.UVOLlog##c.copeptinlog age10 sex01 diabetes ckd_epi10 kid_vol100 d_diur center || treenbr:

* Hypothesis 2 *
mixed UOSMmeasuredsqrt c.UMOD24hsqrt##c.copeptinlog UVOLlog age10 sex01 diabetes ckd_epi10 kid_vol100 d_diur center || treenbr:






// Adjusting for creatininuria
* Hypothesis 1 *
mixed UMOD24hsqrt c.UVOLlog##c.copeptinlog age10 sex01 diabetes ckd_epi10 kid_vol100 d_diur ucrt_mgkg24sqrt center || treenbr:

* Hypothesis 2 *
mixed UOSMmeasuredsqrt c.UMOD24hsqrt##c.copeptinlog UVOLlog age10 sex01 diabetes ckd_epi10 kid_vol100 d_diur ucrt_mgkg24sqrt center || treenbr:






// Without CKD patients
drop if ckd==1

* Hypothesis 1 *
mixed UMOD24hsqrt c.UVOLlog##c.copeptinlog age10 sex01 diabetes ckd_epi10 kid_vol100 d_diur center || treenbr:

* Hypothesis 2 *
mixed UOSMmeasuredsqrt c.UMOD24hsqrt##c.copeptinlog UVOLlog age10 sex01 diabetes ckd_epi10 kid_vol100 d_diur center || treenbr:





// With calculated UOSM instead of measured UOSM
* ==> Separate dofile



// Adjusting for UNa24h as a proxy for salt consumption
* Hypothesis 1 *
mixed UMOD24hsqrt c.UVOLlog##c.copeptinlog age10 sex01 diabetes ckd_epi10 kid_vol100 d_diur una_u24 center || treenbr:

* Hypothesis 2 *
mixed UOSMmeasuredsqrt c.UMOD24hsqrt##c.copeptinlog UVOLlog age10 sex01 diabetes ckd_epi10 kid_vol100 d_diur una_u24 center || treenbr:
