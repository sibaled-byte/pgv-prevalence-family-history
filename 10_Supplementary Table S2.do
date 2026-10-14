* Supplementary Table S2: Variant-level summary and population allele frequencies

clear all
set more off

use "familyhistory_pv_hpc.dta", clear

keep if sampleqc == "QC clear"

count
display "QC-clear HERPACC participants = " r(N)

gen strL V_APC = APC_H
replace V_APC = APC_P if V_APC == "" & APC_P != ""
replace V_APC = V_APC + ";" + APC_P if V_APC != "" & APC_P != "" & V_APC != APC_P

gen strL V_ATM = ATM_H
replace V_ATM = ATM_P if V_ATM == "" & ATM_P != ""
replace V_ATM = V_ATM + ";" + ATM_P if V_ATM != "" & ATM_P != "" & V_ATM != ATM_P

gen strL V_BARD1 = BARD1_H
replace V_BARD1 = BARD1_P if V_BARD1 == "" & BARD1_P != ""
replace V_BARD1 = V_BARD1 + ";" + BARD1_P if V_BARD1 != "" & BARD1_P != "" & V_BARD1 != BARD1_P

gen strL V_BMPR1A = BMPR1A_H

gen strL V_BRCA1 = BRCA1_H
replace V_BRCA1 = BRCA1_P if V_BRCA1 == "" & BRCA1_P != ""
replace V_BRCA1 = V_BRCA1 + ";" + BRCA1_P if V_BRCA1 != "" & BRCA1_P != "" & V_BRCA1 != BRCA1_P

gen strL V_BRCA2 = BRCA2_H
replace V_BRCA2 = BRCA2_P if V_BRCA2 == "" & BRCA2_P != ""
replace V_BRCA2 = V_BRCA2 + ";" + BRCA2_P if V_BRCA2 != "" & BRCA2_P != "" & V_BRCA2 != BRCA2_P

gen strL V_BRIP1 = BRIP1_H
replace V_BRIP1 = BRIP1_P if V_BRIP1 == "" & BRIP1_P != ""
replace V_BRIP1 = V_BRIP1 + ";" + BRIP1_P if V_BRIP1 != "" & BRIP1_P != "" & V_BRIP1 != BRIP1_P

gen strL V_CDH1 = CDH1_H
replace V_CDH1 = CDH1_P if V_CDH1 == "" & CDH1_P != ""
replace V_CDH1 = V_CDH1 + ";" + CDH1_P if V_CDH1 != "" & CDH1_P != "" & V_CDH1 != CDH1_P

gen strL V_CDK4 = CDK4_H
replace V_CDK4 = CDK4_P if V_CDK4 == "" & CDK4_P != ""
replace V_CDK4 = V_CDK4 + ";" + CDK4_P if V_CDK4 != "" & CDK4_P != "" & V_CDK4 != CDK4_P

gen strL V_CDKN2A = CDKN2A_H
replace V_CDKN2A = CDKN2A_P if V_CDKN2A == "" & CDKN2A_P != ""
replace V_CDKN2A = V_CDKN2A + ";" + CDKN2A_P if V_CDKN2A != "" & CDKN2A_P != "" & V_CDKN2A != CDKN2A_P

gen strL V_CHEK2 = CHEK2_H
replace V_CHEK2 = CHEK2_P if V_CHEK2 == "" & CHEK2_P != ""
replace V_CHEK2 = V_CHEK2 + ";" + CHEK2_P if V_CHEK2 != "" & CHEK2_P != "" & V_CHEK2 != CHEK2_P

gen strL V_EPCAM = EPCAM_H
replace V_EPCAM = EPCAM_P if V_EPCAM == "" & EPCAM_P != ""
replace V_EPCAM = V_EPCAM + ";" + EPCAM_P if V_EPCAM != "" & EPCAM_P != "" & V_EPCAM != EPCAM_P

gen strL V_HOXB13 = HOXB13_H

gen strL V_MLH1 = MLH1_H
replace V_MLH1 = MLH1_P if V_MLH1 == "" & MLH1_P != ""
replace V_MLH1 = V_MLH1 + ";" + MLH1_P if V_MLH1 != "" & MLH1_P != "" & V_MLH1 != MLH1_P

gen strL V_MSH2 = MSH2_H
replace V_MSH2 = MSH2_P if V_MSH2 == "" & MSH2_P != ""
replace V_MSH2 = V_MSH2 + ";" + MSH2_P if V_MSH2 != "" & MSH2_P != "" & V_MSH2 != MSH2_P

gen strL V_MSH6 = MSH6_H
replace V_MSH6 = MSH6_P if V_MSH6 == "" & MSH6_P != ""
replace V_MSH6 = V_MSH6 + ";" + MSH6_P if V_MSH6 != "" & MSH6_P != "" & V_MSH6 != MSH6_P

gen strL V_MUTYH = MUTYH_P

gen strL V_NBN = NBN_H
replace V_NBN = NBN_P if V_NBN == "" & NBN_P != ""
replace V_NBN = V_NBN + ";" + NBN_P if V_NBN != "" & NBN_P != "" & V_NBN != NBN_P

gen strL V_NF1 = NF1_H
replace V_NF1 = NF1_P if V_NF1 == "" & NF1_P != ""
replace V_NF1 = V_NF1 + ";" + NF1_P if V_NF1 != "" & NF1_P != "" & V_NF1 != NF1_P

gen strL V_PALB2 = PALB2_H
replace V_PALB2 = PALB2_P if V_PALB2 == "" & PALB2_P != ""
replace V_PALB2 = V_PALB2 + ";" + PALB2_P if V_PALB2 != "" & PALB2_P != "" & V_PALB2 != PALB2_P

gen strL V_PMS2 = PMS2_H
replace V_PMS2 = PMS2_P if V_PMS2 == "" & PMS2_P != ""
replace V_PMS2 = V_PMS2 + ";" + PMS2_P if V_PMS2 != "" & PMS2_P != "" & V_PMS2 != PMS2_P

gen strL V_PTEN = PTEN_H
replace V_PTEN = PTEN_P if V_PTEN == "" & PTEN_P != ""
replace V_PTEN = V_PTEN + ";" + PTEN_P if V_PTEN != "" & PTEN_P != "" & V_PTEN != PTEN_P

gen strL V_RAD51C = RAD51C_H
replace V_RAD51C = RAD51C_P if V_RAD51C == "" & RAD51C_P != ""
replace V_RAD51C = V_RAD51C + ";" + RAD51C_P if V_RAD51C != "" & RAD51C_P != "" & V_RAD51C != RAD51C_P

gen strL V_RAD51D = RAD51D_H
replace V_RAD51D = RAD51D_P if V_RAD51D == "" & RAD51D_P != ""
replace V_RAD51D = V_RAD51D + ";" + RAD51D_P if V_RAD51D != "" & RAD51D_P != "" & V_RAD51D != RAD51D_P

gen strL V_SMAD4 = ""

gen strL V_STK11 = STK11_H
replace V_STK11 = STK11_P if V_STK11 == "" & STK11_P != ""
replace V_STK11 = V_STK11 + ";" + STK11_P if V_STK11 != "" & STK11_P != "" & V_STK11 != STK11_P

gen strL V_TP53 = TP53_H
replace V_TP53 = TP53_P if V_TP53 == "" & TP53_P != ""
replace V_TP53 = V_TP53 + ";" + TP53_P if V_TP53 != "" & TP53_P != "" & V_TP53 != TP53_P

keep acc V_*

reshape long V_, i(acc) j(gene) string
rename V_ variant

replace variant = trim(variant)
drop if variant == ""

split variant, parse(";") generate(piece)
drop variant

reshape long piece, i(acc gene) j(variant_number)
rename piece variant

replace variant = trim(variant)
drop if variant == ""

duplicates drop acc gene variant, force

egen carrier_tag = tag(acc)
count if carrier_tag
display "QC-clear HERPACC participants with at least one recorded variant = " r(N)
drop carrier_tag

count
display "Participant-variant records after deduplication = " r(N)

bysort gene variant: gen herpacc_carrier_n = _N
bysort gene variant: keep if _n == 1

count
display "Unique variants among QC-clear HERPACC participants = " r(N)

assert strpos(variant, ";") == 0
assert herpacc_carrier_n >= 1

keep gene variant herpacc_carrier_n
gsort gene -herpacc_carrier_n variant

save "HERPACC_variant_counts.dta", replace

export excel ///
    gene variant herpacc_carrier_n ///
    using "HERPACC_variant_counts.xlsx", ///
    firstrow(variables) replace

list gene variant herpacc_carrier_n in 1/30, noobs


	
	
	
	
*this study 	
	
clear all
set more off

use "familyhistory_pv_hpc.dta", clear

keep if sampleqc == "QC clear"

merge 1:1 acc using "ACC_PV_FH_controls_20260803.dta", ///
    keep(match) nogen

count
display "Final analytic cohort = " r(N)
assert r(N) == 11811


gen strL V_APC = APC_H
replace V_APC = APC_P if V_APC == "" & APC_P != ""
replace V_APC = V_APC + ";" + APC_P ///
    if V_APC != "" & APC_P != "" & V_APC != APC_P

gen strL V_ATM = ATM_H
replace V_ATM = ATM_P if V_ATM == "" & ATM_P != ""
replace V_ATM = V_ATM + ";" + ATM_P ///
    if V_ATM != "" & ATM_P != "" & V_ATM != ATM_P

gen strL V_BARD1 = BARD1_H
replace V_BARD1 = BARD1_P if V_BARD1 == "" & BARD1_P != ""
replace V_BARD1 = V_BARD1 + ";" + BARD1_P ///
    if V_BARD1 != "" & BARD1_P != "" & V_BARD1 != BARD1_P

gen strL V_BMPR1A = BMPR1A_H

gen strL V_BRCA1 = BRCA1_H
replace V_BRCA1 = BRCA1_P if V_BRCA1 == "" & BRCA1_P != ""
replace V_BRCA1 = V_BRCA1 + ";" + BRCA1_P ///
    if V_BRCA1 != "" & BRCA1_P != "" & V_BRCA1 != BRCA1_P

gen strL V_BRCA2 = BRCA2_H
replace V_BRCA2 = BRCA2_P if V_BRCA2 == "" & BRCA2_P != ""
replace V_BRCA2 = V_BRCA2 + ";" + BRCA2_P ///
    if V_BRCA2 != "" & BRCA2_P != "" & V_BRCA2 != BRCA2_P

gen strL V_BRIP1 = BRIP1_H
replace V_BRIP1 = BRIP1_P if V_BRIP1 == "" & BRIP1_P != ""
replace V_BRIP1 = V_BRIP1 + ";" + BRIP1_P ///
    if V_BRIP1 != "" & BRIP1_P != "" & V_BRIP1 != BRIP1_P

gen strL V_CDH1 = CDH1_H
replace V_CDH1 = CDH1_P if V_CDH1 == "" & CDH1_P != ""
replace V_CDH1 = V_CDH1 + ";" + CDH1_P ///
    if V_CDH1 != "" & CDH1_P != "" & V_CDH1 != CDH1_P

gen strL V_CDK4 = CDK4_H
replace V_CDK4 = CDK4_P if V_CDK4 == "" & CDK4_P != ""
replace V_CDK4 = V_CDK4 + ";" + CDK4_P ///
    if V_CDK4 != "" & CDK4_P != "" & V_CDK4 != CDK4_P

gen strL V_CDKN2A = CDKN2A_H
replace V_CDKN2A = CDKN2A_P if V_CDKN2A == "" & CDKN2A_P != ""
replace V_CDKN2A = V_CDKN2A + ";" + CDKN2A_P ///
    if V_CDKN2A != "" & CDKN2A_P != "" & V_CDKN2A != CDKN2A_P

gen strL V_CHEK2 = CHEK2_H
replace V_CHEK2 = CHEK2_P if V_CHEK2 == "" & CHEK2_P != ""
replace V_CHEK2 = V_CHEK2 + ";" + CHEK2_P ///
    if V_CHEK2 != "" & CHEK2_P != "" & V_CHEK2 != CHEK2_P

gen strL V_EPCAM = EPCAM_H
replace V_EPCAM = EPCAM_P if V_EPCAM == "" & EPCAM_P != ""
replace V_EPCAM = V_EPCAM + ";" + EPCAM_P ///
    if V_EPCAM != "" & EPCAM_P != "" & V_EPCAM != EPCAM_P

gen strL V_HOXB13 = HOXB13_H

gen strL V_MLH1 = MLH1_H
replace V_MLH1 = MLH1_P if V_MLH1 == "" & MLH1_P != ""
replace V_MLH1 = V_MLH1 + ";" + MLH1_P ///
    if V_MLH1 != "" & MLH1_P != "" & V_MLH1 != MLH1_P

gen strL V_MSH2 = MSH2_H
replace V_MSH2 = MSH2_P if V_MSH2 == "" & MSH2_P != ""
replace V_MSH2 = V_MSH2 + ";" + MSH2_P ///
    if V_MSH2 != "" & MSH2_P != "" & V_MSH2 != MSH2_P

gen strL V_MSH6 = MSH6_H
replace V_MSH6 = MSH6_P if V_MSH6 == "" & MSH6_P != ""
replace V_MSH6 = V_MSH6 + ";" + MSH6_P ///
    if V_MSH6 != "" & MSH6_P != "" & V_MSH6 != MSH6_P

gen strL V_MUTYH = MUTYH_P

gen strL V_NBN = NBN_H
replace V_NBN = NBN_P if V_NBN == "" & NBN_P != ""
replace V_NBN = V_NBN + ";" + NBN_P ///
    if V_NBN != "" & NBN_P != "" & V_NBN != NBN_P

gen strL V_NF1 = NF1_H
replace V_NF1 = NF1_P if V_NF1 == "" & NF1_P != ""
replace V_NF1 = V_NF1 + ";" + NF1_P ///
    if V_NF1 != "" & NF1_P != "" & V_NF1 != NF1_P

gen strL V_PALB2 = PALB2_H
replace V_PALB2 = PALB2_P if V_PALB2 == "" & PALB2_P != ""
replace V_PALB2 = V_PALB2 + ";" + PALB2_P ///
    if V_PALB2 != "" & PALB2_P != "" & V_PALB2 != PALB2_P

gen strL V_PMS2 = PMS2_H
replace V_PMS2 = PMS2_P if V_PMS2 == "" & PMS2_P != ""
replace V_PMS2 = V_PMS2 + ";" + PMS2_P ///
    if V_PMS2 != "" & PMS2_P != "" & V_PMS2 != PMS2_P

gen strL V_PTEN = PTEN_H
replace V_PTEN = PTEN_P if V_PTEN == "" & PTEN_P != ""
replace V_PTEN = V_PTEN + ";" + PTEN_P ///
    if V_PTEN != "" & PTEN_P != "" & V_PTEN != PTEN_P

gen strL V_RAD51C = RAD51C_H
replace V_RAD51C = RAD51C_P if V_RAD51C == "" & RAD51C_P != ""
replace V_RAD51C = V_RAD51C + ";" + RAD51C_P ///
    if V_RAD51C != "" & RAD51C_P != "" & V_RAD51C != RAD51C_P

gen strL V_RAD51D = RAD51D_H
replace V_RAD51D = RAD51D_P if V_RAD51D == "" & RAD51D_P != ""
replace V_RAD51D = V_RAD51D + ";" + RAD51D_P ///
    if V_RAD51D != "" & RAD51D_P != "" & V_RAD51D != RAD51D_P

gen strL V_SMAD4 = ""

gen strL V_STK11 = STK11_H
replace V_STK11 = STK11_P if V_STK11 == "" & STK11_P != ""
replace V_STK11 = V_STK11 + ";" + STK11_P ///
    if V_STK11 != "" & STK11_P != "" & V_STK11 != STK11_P

gen strL V_TP53 = TP53_H
replace V_TP53 = TP53_P if V_TP53 == "" & TP53_P != ""
replace V_TP53 = V_TP53 + ";" + TP53_P ///
    if V_TP53 != "" & TP53_P != "" & V_TP53 != TP53_P


keep acc V_*

reshape long V_, i(acc) j(gene) string
rename V_ variant

replace variant = trim(variant)
drop if variant == ""

split variant, parse(";") generate(piece)
drop variant

reshape long piece, i(acc gene) j(variant_number)
rename piece variant

replace variant = trim(variant)
drop if variant == ""

duplicates drop acc gene variant, force


egen carrier_tag = tag(acc)

count if carrier_tag
display "Analytic-cohort participants with at least one variant = " r(N)
assert r(N) == 220

drop carrier_tag


count
display "Participant-variant records after deduplication = " r(N)


bysort gene variant: gen study_carrier_n = _N

bysort gene variant: keep if _n == 1


count
display "Unique variants in analytic cohort = " r(N)


gen study_carrier_frequency_percent = ///
    100 * study_carrier_n / 11811

format study_carrier_frequency_percent %9.4f

assert strpos(variant, ";") == 0
assert study_carrier_n >= 1
assert study_carrier_n <= 11811


keep gene variant study_carrier_n study_carrier_frequency_percent

gsort gene -study_carrier_n variant


save "Study_variant_counts.dta", replace

export excel ///
    gene ///
    variant ///
    study_carrier_n ///
    study_carrier_frequency_percent ///
    using "Study_variant_counts.xlsx", ///
    firstrow(variables) replace


list gene variant study_carrier_n ///
    study_carrier_frequency_percent ///
    in 1/30, noobs
	
	
	
	
	*combine both outputs 
use "Study_variant_counts.dta", clear

merge 1:1 gene variant using "HERPACC_variant_counts.dta", nogen


replace study_carrier_n = 0 if missing(study_carrier_n)
replace study_carrier_frequency_percent = 0 ///
    if missing(study_carrier_frequency_percent)
	
	

	
save "Variant_summary_combined.dta", replace

order gene variant herpacc_carrier_n ///
    study_carrier_n study_carrier_frequency_percent

gsort gene -study_carrier_n -herpacc_carrier_n variant

export excel ///
    gene variant herpacc_carrier_n ///
    study_carrier_n study_carrier_frequency_percent ///
    using "Variant_summary_combined.xlsx", ///
    firstrow(variables) replace	
	
	
	
	
use "Variant_summary_combined.dta", clear	
	* Keep only variants observed in the final analytic cohort
keep if study_carrier_n > 0

* Create a frequency column that will display correctly in Excel
gen study_carrier_frequency_display = ///
    string(study_carrier_frequency_percent, "%9.4f")

* Rename columns for a clearer supplementary table
rename gene Gene
rename variant Genomic_variant
rename herpacc_carrier_n HERPACC_carriers
rename study_carrier_n Study_carriers

label variable Gene "Gene"
label variable Genomic_variant "Genomic variant"
label variable HERPACC_carriers "QC-clear HERPACC carriers, No."
label variable Study_carriers "Analytic cohort carriers, No."
label variable study_carrier_frequency_display ///
    "Carrier frequency in analytic cohort (%)"

* Add blank columns for the annotations that will be added next
gen strL HGVSc = ""
gen strL HGVSp = ""
gen strL Consequence = ""
gen strL ClinVar_classification = ""
gen str20 Japanese_AF = ""
gen str20 East_Asian_AF = ""
gen str20 Non_Finnish_European_AF = ""

* Arrange and sort the table
order Gene Genomic_variant HGVSc HGVSp Consequence ///
    ClinVar_classification HERPACC_carriers Study_carriers ///
    study_carrier_frequency_display Japanese_AF ///
    East_Asian_AF Non_Finnish_European_AF

gsort Gene -Study_carriers -HERPACC_carriers Genomic_variant

* Confirm expected number of variants
count
assert r(N) == 139

save "Supplementary_variant_summary.dta", replace

* Export the clean table
export excel ///
    Gene ///
    Genomic_variant ///
    HGVSc ///
    HGVSp ///
    Consequence ///
    ClinVar_classification ///
    HERPACC_carriers ///
    Study_carriers ///
    study_carrier_frequency_display ///
    Japanese_AF ///
    East_Asian_AF ///
    Non_Finnish_European_AF ///
    using "Supplementary_variant_summary.xlsx", ///
    firstrow(varlabels) replace
	

	
	
 
* CREATE VEP INPUT VARIABLES FROM GENOMIC VARIANT
 

use "Supplementary_variant_summary.dta", clear

* Remove variables that may remain from a previous run
capture drop v1 v2 v3 v4
capture drop chromosome position reference alternate
capture drop pos start end allele variant_type

 
* SPLIT GENOMIC VARIANT: chromosome_position_REF_ALT
 

split Genomic_variant, parse("_") gen(v)

rename v1 chromosome
rename v2 position
rename v3 reference
rename v4 alternate

replace chromosome = trim(chromosome)
replace position   = trim(position)
replace reference  = trim(reference)
replace alternate  = trim(alternate)

* Check that all four components were successfully extracted
assert !missing(chromosome)
assert !missing(position)
assert !missing(reference)
assert !missing(alternate)

 
* CREATE NUMERIC GENOMIC POSITION
 

gen double pos = real(position)

assert !missing(pos)

* Check that numeric conversion preserved the original coordinate
assert string(pos, "%12.0f") == trim(position)

 
* CLASSIFY VARIANT TYPE
 

gen str20 variant_type = ""

* Single nucleotide variant
replace variant_type = "SNV" ///
    if length(reference)==1 & length(alternate)==1

* Anchored insertion:
* ALT is longer than REF and begins with the REF sequence
replace variant_type = "insertion" ///
    if length(alternate)>length(reference) ///
    & substr(alternate,1,length(reference))==reference

* Anchored deletion:
* REF is longer than ALT and begins with the ALT sequence
replace variant_type = "deletion" ///
    if length(reference)>length(alternate) ///
    & substr(reference,1,length(alternate))==alternate

* Variants not represented by the simple rules above
replace variant_type = "complex" ///
    if variant_type==""

tab variant_type

 
* CREATE VEP COORDINATES AND ALLELES
 

gen double start = .
gen double end   = .
gen strL allele  = ""

 
* SNVs
 

replace start = pos if variant_type=="SNV"
replace end   = pos if variant_type=="SNV"

replace allele = reference + "/" + alternate ///
    if variant_type=="SNV"

 
* DELETIONS
 

* Remove the shared left anchor.
* Example: REF=AT, ALT=A becomes T/-
replace start = pos + length(alternate) ///
    if variant_type=="deletion"

replace end = pos + length(reference) - 1 ///
    if variant_type=="deletion"

replace allele = ///
    substr(reference,length(alternate)+1,.) + "/-" ///
    if variant_type=="deletion"

 
* INSERTIONS
 

* VEP represents an insertion between bases using start = end + 1.
* Remove the shared left anchor from ALT.
replace start = pos + length(reference) ///
    if variant_type=="insertion"

replace end = start - 1 ///
    if variant_type=="insertion"

replace allele = "-/" + ///
    substr(alternate,length(reference)+1,.) ///
    if variant_type=="insertion"

 
* IDENTIFY COMPLEX VARIANTS
 

count if variant_type=="complex"

list Gene Genomic_variant chromosome pos reference alternate ///
    if variant_type=="complex", noobs

 
* FINAL QUALITY CHECKS
 

* Every non-complex variant should now have valid VEP fields
assert !missing(start) if variant_type!="complex"
assert !missing(end)   if variant_type!="complex"
assert allele!=""      if variant_type!="complex"

* Check for unexpected variant classifications
tab variant_type, missing

 
* ARRANGE VARIABLES
 

order Gene Genomic_variant ///
      chromosome pos start end ///
      reference alternate allele variant_type

 
* SAVE STATA VERSION
 

save "Supplementary_variant_summary_VEP.dta", replace

 
* EXPORT EXCEL VERSION FOR CHECKING
 

export excel ///
    Gene Genomic_variant ///
    chromosome pos start end ///
    reference alternate allele variant_type ///
    HERPACC_carriers Study_carriers ///
    study_carrier_frequency_display ///
    using "Supplementary_variant_summary_VEP.xlsx", ///
    firstrow(variables) replace

	
	
	
	
 
* BUILD VARIANT ANNOTATION TABLE - SLOT 1
 

import excel "Variant_annotations.xlsx", ///
    firstrow clear

keep chr_pos_ref_alt1 Gene1 Annotation1 ///
    Clinvar_202106191 HGVSc1 HGVSp1

drop if missing(chr_pos_ref_alt1)

rename chr_pos_ref_alt1 Genomic_variant
rename Gene1 Gene_DrK
rename Annotation1 Consequence_DrK
rename Clinvar_202106191 ClinVar_DrK
rename HGVSc1 HGVSc_DrK
rename HGVSp1 HGVSp_DrK

save "Variant_annotations.dta", replace


 
* BUILD VARIANT ANNOTATION TABLE - SLOT 2
 

import excel "Variant_annotations.xlsx", ///
    firstrow clear

keep chr_pos_ref_alt2 Gene2 Annotation2 ///
    Clinvar_202106192 HGVSc2 HGVSp2

drop if missing(chr_pos_ref_alt2)

rename chr_pos_ref_alt2 Genomic_variant
rename Gene2 Gene_DrK
rename Annotation2 Consequence_DrK
rename Clinvar_202106192 ClinVar_DrK
rename HGVSc2 HGVSc_DrK
rename HGVSp2 HGVSp_DrK

save "Variant_annotations_2.dta", replace


 
* COMBINE BOTH VARIANT ANNOTATION SLOTS
 

use "Variant_annotations.dta", clear

append using "Variant_annotations_2.dta"

* Remove exact duplicate annotation rows
duplicates drop Genomic_variant Consequence_DrK ///
    ClinVar_DrK HGVSc_DrK HGVSp_DrK, force

* Check whether any genomic variant still appears more than once
duplicates tag Genomic_variant, gen(dup)

tab dup

list Genomic_variant Gene_DrK Consequence_DrK ///
    ClinVar_DrK HGVSc_DrK HGVSp_DrK ///
    if dup > 0, noobs

drop dup

* If inspection shows repeated rows contain the same annotation,
* keep one row per genomic variant
duplicates drop Genomic_variant, force

save "Variant_annotations_all.dta", replace



 
* MERGE VARIANT ANNOTATIONS INTO THE 139-VARIANT MASTER TABLE
 

use "Supplementary_variant_summary.dta", clear

* Confirm master table has exactly 139 variants
assert _N == 139
isid Genomic_variant

merge m:1 Genomic_variant using "Variant_annotations_all.dta", ///
    keep(master match)

* Check merge
tab _merge


* Fill annotations only where current table is blank
replace Consequence = Consequence_DrK ///
    if missing(Consequence) & !missing(Consequence_DrK)

replace ClinVar_classification = ClinVar_DrK ///
    if missing(ClinVar_classification) & !missing(ClinVar_DrK)

replace HGVSc = HGVSc_DrK ///
    if missing(HGVSc) & !missing(HGVSc_DrK)

replace HGVSp = HGVSp_DrK ///
    if missing(HGVSp) & !missing(HGVSp_DrK)

* Remove temporary variables
drop Consequence_DrK ClinVar_DrK HGVSc_DrK HGVSp_DrK
drop _merge

* Final safety checks
assert _N == 139
isid Genomic_variant


 
* SAVE CLEAN 139-VARIANT ANNOTATED TABLE
 

save "Supplementary_variant_summary_annotated.dta", replace

export excel using ///
    "Supplementary_variant_summary_annotated.xlsx", ///
    firstrow(variables) replace


 
* IDENTIFY COMPLETELY UNANNOTATED VARIANTS
 

count if missing(Consequence) & ///
         missing(HGVSc) & ///
         missing(HGVSp) & ///
         missing(ClinVar_classification)

* Expected = 49


 
* CREATE CHROMOSOME / POSITION / REF / ALT
 

capture drop chromosome position reference alternate
capture drop v1 v2 v3 v4
capture drop pos start end allele variant_type

split Genomic_variant, parse("_") gen(v)

rename v1 chromosome
rename v2 position
rename v3 reference
rename v4 alternate

replace chromosome = trim(chromosome)
replace position   = trim(position)
replace reference  = trim(reference)
replace alternate  = trim(alternate)

assert !missing(chromosome)
assert !missing(position)
assert !missing(reference)
assert !missing(alternate)

gen double pos = real(position)
assert !missing(pos)


 
* CLASSIFY VARIANT TYPE
 

gen str20 variant_type = ""

* SNV
replace variant_type = "SNV" ///
    if length(reference)==1 & length(alternate)==1

* Anchored insertion
replace variant_type = "insertion" ///
    if length(alternate)>length(reference) & ///
       substr(alternate,1,length(reference))==reference

* Anchored deletion
replace variant_type = "deletion" ///
    if length(reference)>length(alternate) & ///
       substr(reference,1,length(alternate))==alternate

* Anything not represented by simple rules
replace variant_type = "complex" ///
    if variant_type==""

tab variant_type


 
* CREATE VEP COORDINATES
 

gen double start = .
gen double end   = .
gen strL allele = ""


* SNVs
replace start = pos ///
    if variant_type=="SNV"

replace end = pos ///
    if variant_type=="SNV"

replace allele = reference + "/" + alternate ///
    if variant_type=="SNV"


* DELETIONS
* Remove common left anchor
replace start = pos + length(alternate) ///
    if variant_type=="deletion"

replace end = pos + length(reference) - 1 ///
    if variant_type=="deletion"

replace allele = ///
    substr(reference,length(alternate)+1,.) + "/-" ///
    if variant_type=="deletion"


* INSERTIONS
replace start = pos + length(reference) ///
    if variant_type=="insertion"

replace end = start - 1 ///
    if variant_type=="insertion"

replace allele = "-/" + ///
    substr(alternate,length(reference)+1,.) ///
    if variant_type=="insertion"


 
* CHECK THE 49 UNANNOTATED VARIANTS
 

count if missing(Consequence) & ///
         missing(HGVSc) & ///
         missing(HGVSp) & ///
         missing(ClinVar_classification)

* Expected = 49


 
* CHECK HOW MANY CAN GO DIRECTLY TO VEP
 

count if missing(Consequence) & ///
         missing(HGVSc) & ///
         missing(HGVSp) & ///
         missing(ClinVar_classification) & ///
         variant_type!="complex"

* Expected = 48


 
* IDENTIFY THE COMPLEX VARIANT
 

list Gene Genomic_variant reference alternate ///
    if missing(Consequence) & ///
       missing(HGVSc) & ///
       missing(HGVSp) & ///
       missing(ClinVar_classification) & ///
       variant_type=="complex", noobs


 
* INSPECT THE 48 STANDARD VEP VARIANTS
 

list Gene Genomic_variant chromosome start end allele variant_type ///
    if missing(Consequence) & ///
       missing(HGVSc) & ///
       missing(HGVSp) & ///
       missing(ClinVar_classification) & ///
       variant_type!="complex", noobs


 
* EXPORT THE 48 STANDARD VARIANTS FOR VEP
 

gen str1 strand = "+"

export delimited chromosome start end allele strand Genomic_variant ///
    if missing(Consequence) & ///
       missing(HGVSc) & ///
       missing(HGVSp) & ///
       missing(ClinVar_classification) & ///
       variant_type!="complex" ///
    using "VEP_48_unmatched_FIXED.txt", ///
    delimiter(tab) novarnames replace
	
	
	* Export variant subset for quality-control review
export excel ///
    Gene Genomic_variant chromosome start end ///
    reference alternate allele variant_type ///
    if missing(Consequence) & ///
       missing(HGVSc) & ///
       missing(HGVSp) & ///
       missing(ClinVar_classification) & ///
       variant_type!="complex" ///
    using "VEP_unmatched.xlsx", ///
    firstrow(variables) replace
