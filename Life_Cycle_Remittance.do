** Life Cycle Remittance of Rural-Urban Migrants
** Author: Weifeng Dai
** Please do not circulate this file without author's permission.

clear all
set scheme s1mono
set more off

****************************** Directory Setting *******************************


global file "/Users/weifengdai/Documents/"
global code "$file/Undergrad/2020Summer/MigrantsRemittancePaper" 
global stata "$code/STATA"  // file for STATA file
global matlab "$code/Matlab"	// file for calibration and counterfactural
global chip "$file/Database/CHIP" // file for CHIP database
global draft "$code/Draft" // file for latex
global chip2007Migrant  "$chip/2007/CHIP2007-流动人口数据"  // migrant
global chip2007Urban	"$chip/2007/CHIP2007-城镇数据"	// urban
global chip2007Rural 	"$chip/2007/CHIP2007-农村数据"	// rural
global chip2008Rural 	"$chip/2008/CHIP2008-农村数据"	// rural(2008)


********************************************************************************
********************************************************************************
*** 								Rural Data 								 ***
********************************************************************************
********************************************************************************

cd $chip2007Rural
use "RHS_w1_abc.dta",replace
keep if a02==1 // household head
keep a05_1 c18 c21
rename a05_1 birth_year
rename c18 income_1
rename c21 income_2
// change to annual income
gen log_income_1 = log(12*income_1)  // change to annual income
gen log_income_2 = max(log(12*income_2),log_income_1)
// delete the outliers 
qui summ log_income_1, detail
drop if log_income_1 < r(p25) - 1.5*(r(p75)-r(p25)) | log_income_1 > r(p75) + 1.5*(r(p75)-r(p25))
qui summ log_income_2, detail
drop if log_income_2 < r(p25) - 1.5*(r(p75)-r(p25)) | log_income_2 > r(p75) + 1.5*(r(p75)-r(p25))
gen age = 2007 - birth_year
keep if age>=25 &age<=55

************************* Mean of Logged Rural Income **************************
summ log_income_1 // moment for calibration

************************** Age Distribution (Rural) ****************************
keep if age>=25 &age<=55
gen count =1
collapse (sum) no_at_age = count, by(age)
cd $code/Matlab
export delim using "age_distribution_rural.csv",replace


** Since there is no land info in 2007 survey, 
** I pick 2008 survey data as approximation.
cd $chip2008Rural
use "RHS_w2_hijk.dta",clear
keep h01 h01_1
rename h01 land_size
rename h01_1 irrigated_land_size
gen log_land_size = log(land_size+1) // deal with 0
gen log_irrigated_land_size = log(irrigated_land_size+1)

********** Initial Land Endowment Distribution for Rural Residient *************
count if log_land_size ==0
count if log_irrigated_land_size ==0
summ log_land_size if log_land_size>0
summ log_irrigated_land_size if log_irrigated_land_size>0
histogram log_land_size
histogram log_irrigated_land_size,percent



********************************************************************************
********************************************************************************
*** 								Urban Data 								 ***
********************************************************************************
********************************************************************************

cd $chip2007Urban
use "UHS_w1_abc.dta",replace
keep if a02==1 // household head
keep a05_1 c17 c20
rename a05_1 birth_year
rename c17 income_1
rename c20 income_2
gen log_income_1 = log(12*income_1)
gen log_income_2 = max(log(12*income_2),log_income_1)
qui summ log_income_1, detail
drop if log_income_1 < r(p25) - 1.5*(r(p75)-r(p25)) | log_income_1 > r(p75) + 1.5*(r(p75)-r(p25))
qui summ log_income_2, detail
drop if log_income_2 < r(p25) - 1.5*(r(p75)-r(p25)) | log_income_2 > r(p75) + 1.5*(r(p75)-r(p25))
gen age = 2007 - birth_year
keep if age>=25 &age<=55
************************* Mean of Logged Urban Income **************************
summ log_income_1 // moment for counterfactural



********************************************************************************
********************************************************************************
*** 								Migrant Data 							 ***
********************************************************************************
********************************************************************************

cd $chip2007Migrant
use "MHS_w1_fghijls.dta",clear

keep hhcode g100 g108 g122 g136 g139 j112_1 j113
rename j112_1 dummy_land
rename j113	land_size
rename g100 income
rename g108 asset_return
rename g122 consum
rename g136 remit
rename g139 saving
foreach var of varlist income consum remit saving{
gen log_`var' = log(1+`var')
}

preserve

tempfile 2007owner
cd $chip2007Migrant
** owner data
use "MHS_w1_abc.dta",clear
keep hhcode a02 a05 a051 a15 a16 ///
	 b102 /// 
	 c165_1 c165_2 c188
keep if a051 == 2 // age >= 16
keep if a02 == 1 // home owner info
drop a02 a051
rename a05 age
rename a15 hukou
rename a16 hukou_time
rename b102 education
rename c165_1 migrant_year
rename c165_2 migrant_month
rename c188 rural_monthly_income
keep if age>=25 & age<=55
gen rural_income = 12*rural_monthly_income
gen log_rural_income = log(rural_income)
qui summ log_rural_income, detail
drop if log_rural_income < r(p25) - 1.5*(r(p75)-r(p25)) | log_rural_income > r(p75) + 1.5*(r(p75)-r(p25))
save `2007owner'

restore
merge 1:1 hhcode using `2007owner'
drop _merge // Checked that "Not Matched" due to the truncation of age
keep if age>=25 & age<=55
******************* Variance of Income at Age 25 for Migrants *****************
summ log_income if age==25

******************* Life Cycle of Logged Income for Migrants *******************
preserve
keep age log_income log_saving
sort age
collapse (mean) mean_log_income = log_income ///
		 (mean) mean_log_saving = log_saving, ///
		 by(age)
tw (scatter mean_log_income age)(qfit mean_log_income age), ///
ytitle("Mean of Logged Income") legend(off)
cd $draft
graph export "wage profile.png",replace
restore

************** Initial Land Endowment Distribution for Migrants  ***************
gen log_land_size = log(1+land_size)
count if log_land_size ==0
summ log_land_size if log_land_size>0

************ Life Cycle of Number of Parent and Logged Remitttance *************
cd $chip2007Migrant
use "MHS_w1_e3.dta",clear
rename e3 parent_code
rename e304 dummy_living
rename e320 remit
cd $code/STATA

* adjust the dummy for the convenience of summation up number of existing parents
replace dummy_living = 0 if dummy_living == 2 
replace remit = 0 if dummy_living ==0

* keep household head (owner) info
keep if parent_code == 1 | parent_code == 2
keep hhcode parent_code dummy_living remit 
reshape wide dummy_living remit , i(hhcode) j(parent_code)
gen no_parent = dummy_living1+dummy_living2 // number of parent
gen remit = remit1+remit2
merge 1:1 hhcode using `2007owner'
drop _merge
keep if age>=25 & age<=55
gen log_remit = log(1+remit)
qui summ log_remit, detail
drop if log_remit < r(p25) - 1.5*(r(p75)-r(p25)) | log_remit > r(p75) + 1.5*(r(p75)-r(p25))

************** Mean of Logged Remitttance (By Number of Parents) ***************
summ log_remit if no_parent ==2
summ log_remit if no_parent ==1
histogram(log_remit) if no_parent == 2, ytitle("Percent")
histogram(log_remit) if no_parent == 1, ytitle("Percent")

*********************  Life Cycle of Remittance (Fact 2)  **********************
preserve
tempfile life_cycle_remit
keep log_remit age
collapse (mean) mean_log_remit = log_remit, by(age)
cd $draft
twoway (scatter mean_log_remit age)(qfit mean_log_remit age), ///
xtitle("Age") ytitle("Mean of Logged Remittance") legend(off)
graph export "remittance profile.png",replace
restore

*************************  Age Distribution of Migrant  ************************
gen count =1
collapse (sum) no_at_age = count, by(age)
cd $code/Matlab
export delim using "age_distribution_migrant.csv",replace

