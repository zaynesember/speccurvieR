#' CalCOFI Bottle Data
#'
#' A subset of data from the California Cooperative Oceanic Fisheries
#' Investigations. Each observation describes a sample of ocean water
#' collected.
#'
#' @format ## `bottles`
#' A data frame with 500 rows and 62 columns:
#' \describe{
#'   \item{Cst_Cnt}{Cast count}
#'   \item{Btl_Cnt}{Bottle Count}
#'   \item{Sta_ID}{Line and Station}
#'   \item{Depth_ID}{Depth ID}
#'   \item{Depthm}{Bottle depth in meters}
#'   \item{T_degC}{Water temperature in degrees Celsius}
#'   \item{Salnty}{Salinity (Practical Salinity Scale 1978)}
#'   \item{O2ml_L}{Milliliters of oxygen per liter of seawater}
#'   \item{STheta}{Potential Density (Sigma Theta), Kg/M³}
#'   \item{O2Sat}{Oxygen percent saturation}
#'   \item{Oxy_µmol/Kg}{Oxygen micromoles per kilogram seawater}
#'   \item{BtlNum}{Niskin bottle sample was collected from}
#'   \item{RecInd}{Record Indicator}
#'   \item{T_prec}{Temperature Precision}
#'   \item{T_qual}{Quality Code}
#'   \item{S_prec}{Salinity Precision}
#'   \item{S_qual}{Quality Code}
#'   \item{P_qual}{Quality Code}
#'   \item{O_qual}{Quality Code}
#'   \item{SThtaq}{Quality Code}
#'   \item{O2Satq}{Quality Code}
#'   \item{ChlorA}{Micrograms Chlorophyll-a per liter seawater}
#'   \item{Chlqua}{Quality Code}
#'   \item{Phaeop}{Micrograms Phaeopigment per liter seawater}
#'   \item{Phaqua}{Quality Code}
#'   \item{PO4uM}{Micromoles Phosphate per liter of seawater}
#'   \item{PO4q}{Quality Code}
#'   \item{SiO3uM}{Micromoles Silicate per liter of seawater}
#'   \item{SiO3qu}{Quality Code}
#'   \item{NO2uM}{Micromoles Nitrite per liter of seawater}
#'   \item{NO2q}{Quality Code}
#'   \item{NO3uM}{Micromoles Nitrate per liter of seawater}
#'   \item{NO3q}{Quality Code}
#'   \item{NH3uM}{Micromoles Ammonia per liter of seawater}
#'   \item{NH3q}{Quality Code}
#'   \item{C14As1}{14C Assimilation of Replicate 1}
#'   \item{C14A1p}{Precision of 14C Assimilation of Replicate 1}
#'   \item{C14A1q}{Quality Code}
#'   \item{C14As2}{14C Assimilation of Replicate 2}
#'   \item{C14A2p}{Precision of 14C Assimilation of Replicate 2}
#'   \item{C14A2q}{Quality Code}
#'   \item{DarkAs}{14C Assimilation of Dark/Control Bottle}
#'   \item{DarkAp}{Precision of 14C Assimilationof Dark/Control Bottle}
#'   \item{Darkaq}{Quality Code}
#'   \item{MeanAs}{Mean 14C Assimilation of Replicates 1 and 2}
#'   \item{MeanAp}{Precision of Mean 14C Assimilation of Replicates 1 and 2}
#'   \item{MeanAq}{Quality Code}
#'   \item{IncTim}{Elapsed incubation time of primary productivity experiment}
#'   \item{LightP}{Light intensities of the incubation tubes}
#'   \item{R_Depth}{Reported Depth (from pressure) in meters}
#'   \item{R_Temp}{Reported (Potential) Temperature in degrees Celsius}
#'   \item{R_Sal}{Reported Salinity (from Specific Volume Anomoly, M³/Kg)}
#'   \item{R_DYNHT}{Reported Dynamic Height in units of dynamic meters}
#'   \item{R_Nuts}{Reported Ammonium concentration}
#'   \item{R_Oxy_µmol/Kg}{Reported Oxygen micromoles/kilogram}
#'   \item{DIC1}{Dissolved Inorganic Carbon micromoles per kilogram solution}
#'   \item{DIC2}{Dissolved Inorganic Carbon on a replicate sample}
#'   \item{TA1}{Total Alkalinity micromoles per kilogram solution}
#'   \item{TA2}{Total Alkalinity on a replicate sample}
#'   \item{pH1}{pH (the degree of acidity/alkalinity of a solution)}
#'   \item{pH2}{pH on a replicate sample}
#'   \item{DIC Quality Comment}{Quality Comment}
#' }
#' @source <https://calcofi.org/data/oceanographic-data/bottle-database/>
"bottles"
