load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

function get_title_position(plot)
local vpy, vpx, vph
begin
  getvalues plot
      "vpYF"       : vpy
      "vpXF"       : vpx
      "vpHeightF"  : vph
      "vpWidthF"   : vpw
  end getvalues
  pos = (/vpx,vpy, vpw, vph/)
  return(pos)
end

begin
;============================================================================
; sst 
ncdir0  = "//Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SST/Had/"
fil0    = systemfunc("ls " + ncdir0 + "anom_seasonly_sst_had.nc")
f0      = addfile(fil0, "r")
sst_djf = f0->sst_djf
sst_son = f0->sst_son
sst_jja = f0->sst_jja
sst_mam = f0->sst_mam
sst_annual = f0->sst_ann

EP_el = (/1982,1997,2015/)
el_year = dimsizes(EP_el)

nlat = dimsizes(sst_djf&lat)
nlon = dimsizes(sst_djf&lon)

el_ta = new((/el_year, nlat, nlon/), float)
el_ta(0,:,:) = sst_djf((EP_el(0)-1948),:,:)
el_ta(1,:,:) = sst_djf((EP_el(1)-1948),:,:)
el_ta(2,:,:) = sst_djf((EP_el(2)-1948),:,:)

ep_ta = dim_avg_n_Wrap(el_ta, 0)
printVarSummary(ep_ta)
ep_ta@lat2d = conform(ep_ta, ep_ta&lat,0)
ep_ta@lon2d = conform(ep_ta, ep_ta&lon,1)
ep_ta = where(ep_ta@lat2d .ge. (-20) .and. ep_ta@lat2d .le. (15) .and. ep_ta@lon2d .ge. (150) .and. ep_ta@lon2d .le. 300, ep_ta, ep_ta@_FillValue)
ep_ta = where(ep_ta .ge. 0.0, ep_ta, ep_ta@_FillValue)
; read rsi
ncdir00  = "//Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SST/Had/"
fil00    = systemfunc("ls " + ncdir00 + "rsi_seasonally_sst_had.nc")
f00      = addfile(fil00, "r")
sst_rsi_djf = f00->rsi_djf
sst_rsi_son = f00->rsi_son
sst_rsi_jja = f00->rsi_jja
sst_rsi_mam = f00->rsi_mam
sst_rsi_annual = f00->rsi_mam
;------------select sst regions-----------------------
; region 1  Central-North Pac
sst1_djf   = sst_djf(:, {45}, {200})
sst_rsi1_djf   = sst_rsi_djf(:, {45}, {200})
; region 2  KOE Pac
sst2_djf   = sst_djf(:, {25}, {150})
sst_rsi2_djf   = sst_rsi_djf(:, {25}, {150})
; region 3 south Pac
sst3_djf   = sst_djf(:, {-55}, {197.5})
sst_rsi3_djf   = sst_rsi_djf(:, {-55}, {197.5})
; region 4 south Pac
 sst4_djf   = sst_son(:, {-35}, {262.5})
sst_rsi4_djf   = sst_rsi_son(:, {-35}, {262.5})

sst4_son   = sst_son(:, {-42.5}, {202.5})
sst_rsi4_son   = sst_rsi_son(:, {-42.5}, {202.5})
; southeastern Indian ocean 
sst5_son   = sst_son(:, {-30}, {107.5})
sst_rsi5_son   = sst_rsi_son(:, {-30}, {107.5})
; north atalantic 
sst1_mam   = sst_mam(:, {42.5}, {317.5})
sst_rsi1_mam   = sst_rsi_mam(:, {42.5}, {317.5})
; Gulf of Mexico
sst2_annual   = sst_djf(:, {25}, {270})
sst_rsi2_annual  = sst_rsi_djf(:, {25}, {270})
;=============================================================
ncdir1  = "//Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SAT/NCEP/"
fil1    = systemfunc("ls " + ncdir1 + "anom_seasonly_sat_ncep.nc")
f1      = addfile(fil1, "r")
sat_djf = f1->sat_djf
sat_son = f1->sat_son
sat_jja = f1->sat_jja
sat_mam = f1->sat_mam
sat_annual = f1->sat_ann

ncdir11  = "//Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SAT/NCEP/"
fil11    = systemfunc("ls " + ncdir11 + "rsi_seasonally_sat_ncep.nc")
f11      = addfile(fil11, "r")
sat_rsi_djf = f11->rsi_djf
sat_rsi_son = f11->rsi_son
sat_rsi_jja = f11->rsi_jja
sat_rsi_mam = f11->rsi_mam
sat_rsi_annual = f11->rsi_ann
; region 1  South atlantic
sat1_djf   = sat_djf(:, {-22.5}, {337.5})
sat_rsi1_djf   = sat_rsi_djf(:,  {-22.5}, {337.5})

; region2 SWP 
sat2_djf   = sat_djf(:, {-12.5}, {162.5})
sat_rsi2_djf   = sat_rsi_djf(:,{-12.5}, {162.5})
;print(sat_rsi2_djf)
; eastern europe 
sat1_mam   = sat_mam(:, {60}, {40})
sat_rsi1_mam   = sat_rsi_mam(:,{60}, {40})

; eastern africa
sat2_mam   = sat_mam(:, {5}, {30})
sat_rsi2_mam   = sat_rsi_mam(:,{5}, {30})

; south asia 
sat3_mam   = sat_mam(:, {42.5}, {80})
sat_rsi3_mam   = sat_rsi_mam(:,{42.5}, {80})
print(sat_rsi3_mam)

; north america 
sat1_ann   = sat_annual(:, {72.5}, {290})
sat_rsi1_ann   = sat_rsi_annual(:,{72.5}, {290})
print(sat_rsi1_ann)
;=================================================================
;read soil 
ncdir2  = "//Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/CPC/"
fil2    = systemfunc("ls " + ncdir2 + "anom_seasonly_soilw_cpc.nc")
f2      = addfile(fil2, "r")
soil_djf = f2->soilw_djf
soil_son = f2->soilw_son
soil_jja = f2->soilw_jja
soil_mam = f2->soilw_mam
soil_annual = f2->soilw_ann

ncdir22  = "//Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/CPC/"
fil21    = systemfunc("ls " + ncdir22 + "cpc_soilw_rsi_djf.nc")
fil22    = systemfunc("ls " + ncdir22 + "cpc_soilw_rsi_son.nc")
fil23    = systemfunc("ls " + ncdir22 + "cpc_soilw_rsi_jja.nc")
fil24    = systemfunc("ls " + ncdir22 + "cpc_soilw_rsi_mam.nc")
fil25    = systemfunc("ls " + ncdir22 + "cpc_soilw_rsi_annual.nc")

f21      = addfile(fil21, "r")
f22      = addfile(fil22, "r")
f23      = addfile(fil23, "r")
f24      = addfile(fil24, "r")
f25      = addfile(fil25, "r")

soil_rsi_djf = f21->rsi
soil_rsi_son = f22->rsi
soil_rsi_jja = f23->rsi
soil_rsi_mam = f24->rsi
soil_rsi_annual = f25->rsi
; region1 central austrilia
soil1_annual = soil_annual(:,{-27.5},{127.5})
soil1_rsi_annual = soil_rsi_annual(:,{-27.5},{127.5})
print(soil1_rsi_annual)

; region2 eastern europe
soil2_annual = soil_annual(:,{45},{70})
soil2_rsi_annual = soil_rsi_annual(:,{45},{70})
print(soil2_rsi_annual)

; region3 south asia 
soil3_annual = soil_jja(:,{30},{92.5})
soil3_rsi_annual = soil_rsi_jja(:,{30},{92.5})
print(soil3_rsi_annual)

; region4 amazon 
soil4_annual = soil_annual(:,{-5},{305})
soil4_rsi_annual = soil_rsi_annual(:,{-5},{305})
print(soil4_rsi_annual)

soil1_jja = soil_jja(:,{17.5},{17.5})
soil1_rsi_jja = soil_rsi_jja(:,{17.5},{17.5})
;;=================================================================
;read SST data and rsi
ncdir1  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SST/time_series/"
fil1    = systemfunc("ls " + ncdir1 + "HadISST_DJF_sst_rsi_time_series_single_points.nc")
f1      = addfile(fil1, "r")
sst1_DJF  =f1->sst_ECP_DJF 
sst2_DJF  =f1->sst_NTA_DJF 
sst3_DJF  =f1->sst_CSA_DJF
sst4_DJF  =f1->sst_CNP_DJF 
sst5_DJF  =f1->sst_WNP_DJF 
sst6_DJF  =f1->sst_EAT_DJF 

rsi1_DJF =f1->rsi_ECP_DJF                        
rsi2_DJF =f1->rsi_NTA_DJF 
rsi3_DJF =f1->rsi_CSA_DJF  
rsi4_DJF =f1->rsi_CNP_DJF 
rsi5_DJF =f1->rsi_WNP_DJF 
rsi6_DJF =f1->rsi_EAT_DJF
;===========================================================
fil2    = systemfunc("ls " + ncdir1 + "HadISST_SON_sst_rsi_time_series_single_points.nc")
f2      = addfile(fil2, "r")
sst1_SON = f2->sst_CNP_SON    
sst2_SON = f2->sst_SP_SON     
sst3_SON = f2->sst_EAT_SON    
sst4_SON = f2->sst_SEI_SON  
sst5_SON = f2->sst_CSA_SON  
sst6_SON = f2->sst_SA_SON  

rsi1_SON = f2->rsi_CNP_SON    
rsi2_SON = f2->rsi_SP_SON     
rsi3_SON = f2->rsi_EAT_SON    
rsi4_SON = f2->rsi_SEI_SON  
rsi5_SON = f2->rsi_CSA_SON  
rsi6_SON = f2->rsi_SA_SON  
;=========================================================
fil3    = systemfunc("ls " + ncdir1 + "HadISST_JJA_sst_rsi_time_series_single_points.nc")
f3      = addfile(fil3, "r")

sst1_JJA = f3->sst_WNP_JJA      
sst2_JJA = f3->sst_CNP_JJA      
sst3_JJA = f3->sst_SEI_JJA      
sst4_JJA = f3->sst_EAT_JJA      
sst5_JJA = f3->sst_NEP_JJA      
sst6_JJA = f3->sst_SA_JJA       
sst7_JJA = f3->sst_SP_JJA       
sst8_JJA = f3->sst_SEP_JJA  

rsi1_JJA = f3->rsi_WNP_JJA      
rsi2_JJA = f3->rsi_CNP_JJA      
rsi3_JJA = f3->rsi_SEI_JJA      
rsi4_JJA = f3->rsi_EAT_JJA      
rsi5_JJA = f3->rsi_NEP_JJA      
rsi6_JJA = f3->rsi_SA_JJA       
rsi7_JJA = f3->rsi_SP_JJA       
rsi8_JJA = f3->rsi_SEP_JJA  
;============================================================
fil4    = systemfunc("ls " + ncdir1 + "HadISST_MAM_sst_rsi_time_series_single_points.nc")
f4      = addfile(fil4, "r")

sst1_MAM = f4->sst_CNP_MAM  
sst2_MAM = f4->sst_SAT_MAM  
sst3_MAM = f4->sst_EAT_MAM  
sst4_MAM = f4->sst_SWP_MAM  
sst5_MAM = f4->sst_SEI_MAM 


rsi1_MAM = f4->rsi_CNP_MAM  
rsi2_MAM = f4->rsi_SAT_MAM  
rsi3_MAM = f4->rsi_EAT_MAM  
rsi4_MAM = f4->rsi_SWP_MAM  
rsi5_MAM = f4->rsi_SEI_MAM 
;=================================================================
ncdir2  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SAT/time_series/"
fil01    = systemfunc("ls " + ncdir2 + "SAT_DJF_ts_rsi_time_series_single_points.nc")
f01      = addfile(fil01, "r")

ts1_DJF  = f01->ts_EAF_DJF 
ts2_DJF  = f01->ts_SAS_DJF 
ts3_DJF  = f01->ts_MC_DJF  
ts4_DJF  = f01->ts_CNP_DJF 
ts5_DJF  = f01->ts_SSP_DJF 
ts6_DJF  = f01->ts_CAM_DJF 
ts7_DJF  = f01->ts_AMA_DJF
ts8_DJF  = f01->ts_SAM_DJF
ts9_DJF  = f01->ts_SAT_DJF
ts10_DJF = f01->ts_SWP_DJF 

rsi01_DJF  = f01->rsi_EAF_DJF 
rsi02_DJF  = f01->rsi_SAS_DJF 
rsi03_DJF  = f01->rsi_MC_DJF  
rsi04_DJF  = f01->rsi_CNP_DJF 
rsi05_DJF  = f01->rsi_SSP_DJF 
rsi06_DJF  = f01->rsi_CAM_DJF 
rsi07_DJF  = f01->rsi_AMA_DJF
rsi08_DJF  = f01->rsi_SAM_DJF
rsi09_DJF  = f01->rsi_SAT_DJF
rsi010_DJF = f01->rsi_SWP_DJF 
;================================================================
fil02    = systemfunc("ls " + ncdir2 + "SAT_SON_ts_rsi_time_series_single_points.nc")
f02      = addfile(fil02, "r")

ts1_SON = f02->ts_EAF_SON 
ts2_SON = f02->ts_SAS_SON 
ts3_SON = f02->ts_MC_SON
ts4_SON = f02->ts_CNP_SON 
ts5_SON = f02->ts_SEI_SON 
ts6_SON = f02->ts_CAM_SON 
ts7_SON = f02->ts_AMA_SON 
ts8_SON = f02->ts_SAM_SON 
ts9_SON = f02->ts_SAT_SON 
ts10_SON = f02->ts_SWP_SON 

rsi01_SON  = f02->rsi_EAF_SON 
rsi02_SON  = f02->rsi_SAS_SON 
rsi03_SON  = f02->rsi_MC_SON
rsi04_SON  = f02->rsi_CNP_SON 
rsi05_SON  = f02->rsi_SEI_SON 
rsi06_SON  = f02->rsi_CAM_SON 
rsi07_SON  = f02->rsi_AMA_SON 
rsi08_SON  = f02->rsi_SAM_SON 
rsi09_SON  = f02->rsi_SAT_SON 
rsi010_SON = f02->rsi_SWP_SON 
;===================================================================
fil03    = systemfunc("ls " + ncdir2 + "SAT_JJA_ts_rsi_time_series_single_points.nc")
f03      = addfile(fil03, "r")

ts1_JJA  = f03->ts_EAF_JJA     
ts2_JJA  = f03->ts_SAS_JJA     
ts3_JJA  = f03->ts_MC_JJA    
ts4_JJA  = f03->ts_CNP_JJA     
ts5_JJA  = f03->ts_SEI_JJA     
ts6_JJA  = f03->ts_CAM_JJA     
ts7_JJA  = f03->ts_AMA_JJA   
ts8_JJA  = f03->ts_SAM_JJA   
ts9_JJA  = f03->ts_SAT_JJA   
ts10_JJA = f03->ts_SWP_JJA  

rsi01_JJA  = f03->rsi_EAF_JJA     
rsi02_JJA  = f03->rsi_SAS_JJA     
rsi03_JJA  = f03->rsi_MC_JJA    
rsi04_JJA  = f03->rsi_CNP_JJA     
rsi05_JJA  = f03->rsi_SEI_JJA     
rsi06_JJA  = f03->rsi_CAM_JJA     
rsi07_JJA  = f03->rsi_AMA_JJA   
rsi08_JJA  = f03->rsi_SAM_JJA   
rsi09_JJA  = f03->rsi_SAT_JJA   
rsi010_JJA = f03->rsi_SWP_JJA  
;================================================================
fil04    = systemfunc("ls " + ncdir2 + "SAT_MAM_ts_rsi_time_series_single_points.nc")
f04      = addfile(fil04, "r")

ts1_MAM  = f04->ts_EAF_MAM  
ts2_MAM  = f04->ts_NAF_MAM  
ts3_MAM  = f04->ts_MC_MAM  
ts4_MAM  = f04->ts_CNP_MAM  
ts5_MAM  = f04->ts_SEI_MAM  
ts6_MAM  = f04->ts_NTA_MAM  
ts7_MAM  = f04->ts_AMA_MAM 
ts8_MAM  = f04->ts_SAM_MAM 
ts9_MAM  = f04->ts_SAT_MAM 
ts10_MAM = f04->ts_SWP_MAM  

rsi01_MAM  = f04->rsi_EAF_MAM  
rsi02_MAM  = f04->rsi_NAF_MAM  
rsi03_MAM  = f04->rsi_MC_MAM  
rsi04_MAM  = f04->rsi_CNP_MAM  
rsi05_MAM  = f04->rsi_SEI_MAM  
rsi06_MAM  = f04->rsi_NTA_MAM  
rsi07_MAM  = f04->rsi_AMA_MAM 
rsi08_MAM  = f04->rsi_SAM_MAM 
rsi09_MAM  = f04->rsi_SAT_MAM 
rsi010_MAM = f04->rsi_SWP_MAM  
;===========================================================
ncdir3  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/time_series//"
fil001    = systemfunc("ls " + ncdir3 + "soil_DJF_ts_rsi_time_series_single_points.nc")
f001      = addfile(fil001, "r")

 soil1_DJF  = f001->ts_SAF_DJF  
 soil2_DJF  = f001->ts_NAF_DJF  
 soil3_DJF  = f001->ts_SAS_DJF  
 soil4_DJF  = f001->ts_CAS_DJF  
 soil5_DJF  = f001->ts_NAS_DJF  
 soil6_DJF  = f001->ts_CAU_DJF  
 soil7_DJF  = f001->ts_NAM_DJF
 soil8_DJF  = f001->ts_WAZ_DJF
 soil9_DJF  = f001->ts_WGL_DJF
soil10_DJF = f001->ts_ICP_DJF    
soil11_DJF = f001->ts_EAM_DJF

 rsi001_DJF  = f001->rsi_SAF_DJF  
 rsi002_DJF  = f001->rsi_NAF_DJF  
 rsi003_DJF  = f001->rsi_SAS_DJF  
 rsi004_DJF  = f001->rsi_CAS_DJF  
 rsi005_DJF  = f001->rsi_NAS_DJF  
 rsi006_DJF  = f001->rsi_CAU_DJF  
 rsi007_DJF  = f001->rsi_NAM_DJF
 rsi008_DJF  = f001->rsi_WAZ_DJF
 rsi009_DJF  = f001->rsi_WGL_DJF
 rsi0010_DJF  = f001->rsi_ICP_DJF    
 rsi0011_DJF  = f001->rsi_EAM_DJF
;===========================================================
fil002    = systemfunc("ls " + ncdir3 + "soil_SON_ts_rsi_time_series_single_points.nc")
f002      = addfile(fil002, "r")
 
 soil1_SON  = f002->ts_SAF_SON  
 soil2_SON  = f002->ts_NAF_SON  
 soil3_SON  = f002->ts_SAS_SON  
 soil4_SON  = f002->ts_CAS_SON  
 soil5_SON  = f002->ts_NAS_SON  
 soil6_SON  = f002->ts_CAU_SON  
 soil7_SON  = f002->ts_NAM_SON 
 soil8_SON  = f002->ts_WAZ_SON 
 soil9_SON  = f002->ts_WGL_SON 
soil10_SON  = f002->ts_ICP_SON  
soil11_SON  = f002->ts_CAM_SON      
soil12_SON  = f002->ts_EAM_SON

rsi001_SON     = f002->rsi_SAF_SON  
rsi002_SON     = f002->rsi_NAF_SON  
rsi003_SON     = f002->rsi_SAS_SON  
rsi004_SON     = f002->rsi_CAS_SON  
rsi005_SON     = f002->rsi_NAS_SON  
rsi006_SON     = f002->rsi_CAU_SON  
rsi007_SON     = f002->rsi_NAM_SON 
rsi008_SON     = f002->rsi_WAZ_SON 
rsi009_SON     = f002->rsi_WGL_SON 
rsi0010_SON    = f002->rsi_ICP_SON  
rsi0011_SON    = f002->rsi_CAM_SON      
rsi0012_SON    = f002->rsi_EAM_SON
;============================================================
fil003    = systemfunc("ls " + ncdir3 + "soil_JJA_ts_rsi_time_series_single_points.nc")
f003      = addfile(fil003, "r")

 soil1_JJA  = f003->ts_SAF_JJA  
 soil2_JJA  = f003->ts_NAF_JJA  
 soil3_JJA  = f003->ts_SAS_JJA  
 soil4_JJA  = f003->ts_CAS_JJA  
 soil5_JJA  = f003->ts_NAS_JJA  
 soil6_JJA  = f003->ts_CAU_JJA  
 soil7_JJA  = f003->ts_NAM_JJA
 soil8_JJA  = f003->ts_CAZ_JJA
 soil9_JJA  = f003->ts_WGL_JJA
soil10_JJA  = f003->ts_ICP_JJA 
soil11_JJA  = f003->ts_CAM_JJA      
soil12_JJA  = f003->ts_EAM_JJA  

 rsi001_JJA   = f003->rsi_SAF_JJA  
 rsi002_JJA   = f003->rsi_NAF_JJA  
 rsi003_JJA   = f003->rsi_SAS_JJA  
 rsi004_JJA   = f003->rsi_CAS_JJA  
 rsi005_JJA   = f003->rsi_NAS_JJA  
 rsi006_JJA   = f003->rsi_CAU_JJA  
 rsi007_JJA   = f003->rsi_NAM_JJA
 rsi008_JJA   = f003->rsi_CAZ_JJA
 rsi009_JJA   = f003->rsi_WGL_JJA
 rsi0010_JJA  = f003->rsi_ICP_JJA 
 rsi0011_JJA  = f003->rsi_CAM_JJA      
 rsi0012_JJA  = f003->rsi_EAM_JJA 
;============================================
fil004    = systemfunc("ls " + ncdir3 + "soil_MAM_ts_rsi_time_series_single_points.nc")
f004      = addfile(fil004, "r")

   soil1_MAM  = f004->ts_SAF_MAM     
   soil2_MAM  = f004->ts_NAF_MAM     
   soil3_MAM  = f004->ts_SAS_MAM     
   soil4_MAM  = f004->ts_CAS_MAM     
   soil5_MAM  = f004->ts_NAS_MAM     
   soil6_MAM  = f004->ts_CAU_MAM     
   soil7_MAM  = f004->ts_NAM_MAM   
   soil8_MAM  = f004->ts_WAZ_MAM   
   soil9_MAM  = f004->ts_WGL_MAM   
  soil10_MAM  = f004->ts_ICP_MAM        
  soil11_MAM  = f004->ts_EAM_MAM    

  rsi001_MAM  = f004->rsi_SAF_MAM     
  rsi002_MAM  = f004->rsi_NAF_MAM     
  rsi003_MAM  = f004->rsi_SAS_MAM     
  rsi004_MAM  = f004->rsi_CAS_MAM     
  rsi005_MAM  = f004->rsi_NAS_MAM     
  rsi006_MAM  = f004->rsi_CAU_MAM     
  rsi007_MAM  = f004->rsi_NAM_MAM   
  rsi008_MAM  = f004->rsi_WAZ_MAM   
  rsi009_MAM  = f004->rsi_WGL_MAM   
  rsi0010_MAM  = f004->rsi_ICP_MAM        
  rsi0011_MAM  = f004->rsi_EAM_MAM    
;===========================================================
  ;=========================================================
  ; read SST 
  ; CNP
  data1 = new((/2,75/), float)
  data1(0,:) = sst_rsi1_djf
  data1(1,:) = sst1_djf
  ; KOE 
  data2 = new((/2,75/), float)
  data2(0,:) = sst_rsi2_djf
  data2(1,:) = sst2_djf
  ;South Pac 
  data3 = new((/2,75/), float)
  data3(0,:) = rsi4_SON
  data3(1,:) = sst4_SON
  ;north atlantic
  data4 = new((/2,75/), float)
  data4(0,:) = sst_rsi4_djf
  data4(1,:) = sst4_djf

  ;Gulf of Mexico
  data04 = new((/2,75/), float)
  data04(0,:) = sst_rsi2_annual
  data04(1,:) = sst2_annual

  ;=============================================
  ; read SAT 
  ; DJF south atlantic
  data5 = new((/2,75/), float)
  data5(0,:) = rsi01_DJF
  data5(1,:) = ts1_DJF
  ; southwestern Pac
  data6 = new((/2,75/), float)
  data6(0,:) = rsi03_SON
  data6(1,:) = ts3_SON
  print(data6(0,:))
  ; East europe
  data7 = new((/2,75/), float)
 ;data7(0,:) = sat_rsi1_ann
 ;data7(1,:) = sat1_ann
  data7(0,:) = sat_rsi1_mam
  data7(1,:) = sat1_mam
  ; SWP
  data8 = new((/2,75/), float)
  data8(0,:) = rsi010_JJA
  data8(1,:) = ts10_JJA

 ; read soil 
 ; annual central austrilia
  data9 = new((/2,75/), float)
  data9(0,:) = soil1_rsi_annual
  data9(1,:) = soil1_annual
  ; central Asia
  data10 = new((/2,75/), float)
  data10(0,:) = soil2_rsi_annual
  data10(1,:) = soil2_annual
  ; south asia
  data11 = new((/2,75/), float)
 ;data11(0,:) = soil3_rsi_annual
 ;data11(1,:) = soil3_annual
 data11(0,:) = rsi009_SON
 data11(1,:) = soil9_SON
  ; amazon
  data12 = new((/2,75/), float)
  data12(0,:) = soil4_rsi_annual
  data12(1,:) = soil4_annual

 data1(0,71:74) = sst_rsi1_djf@_FillValue
 data2(0,71:74) = sst_rsi1_djf@_FillValue
 data3(0,71:74) = sst_rsi1_djf@_FillValue
 data4(0,71:74) = sst_rsi1_djf@_FillValue
 data4(0,71:74) = sst_rsi1_djf@_FillValue
 data04(0,71:74) = sst_rsi1_djf@_FillValue
;========================================
 data5(0,71:74) = sst_rsi1_djf@_FillValue
 data6(0,71:74) = sst_rsi1_djf@_FillValue
 data7(0,71:74) = sst_rsi1_djf@_FillValue
 data8(0,71:74) = sst_rsi1_djf@_FillValue
;========================================
 data9(0,71:74) = sst_rsi1_djf@_FillValue
 data10(0,71:74) = sst_rsi1_djf@_FillValue
 data11(0,71:74) = sst_rsi1_djf@_FillValue
 data12(0,71:74) = sst_rsi1_djf@_FillValue

  data1(0,:)   = where(data1(0,:).eq.sst_rsi1_djf@_FillValue, data1@_FillValue,data1(0,:))
  data2(0,:)   = where(data2(0,:).eq.sst_rsi1_djf@_FillValue, data2@_FillValue,data2(0,:))
  data3(0,:)   = where(data3(0,:).eq.sst_rsi1_djf@_FillValue, data3@_FillValue,data3(0,:))
  data4(0,:)   = where(data4(0,:).eq.sst_rsi1_djf@_FillValue, data4@_FillValue,data4(0,:))
  data04(0,:)   = where(data04(0,:).eq.sst_rsi1_djf@_FillValue, data04@_FillValue,data04(0,:))
;=========================================================================================
  data5(0,:)   = where(data5(0,:).eq.sst_rsi1_djf@_FillValue, data5@_FillValue,data5(0,:))
  data6(0,:)   = where(data6(0,:).eq.sst_rsi1_djf@_FillValue, data6@_FillValue,data6(0,:))
  data7(0,:)   = where(data7(0,:).eq.sst_rsi1_djf@_FillValue, data7@_FillValue,data7(0,:))
  data8(0,:)   = where(data8(0,:).eq.sst_rsi1_djf@_FillValue, data8@_FillValue,data8(0,:))
  ;=======================================================================================
  data9(0,:)   = where(data9(0,:).eq.sst_rsi1_djf@_FillValue, data9@_FillValue,data9(0,:))
  data10(0,:) = where(data10(0,:).eq.sst_rsi1_djf@_FillValue, data10@_FillValue,data10(0,:))
  data11(0,:) = where(data11(0,:).eq.sst_rsi1_djf@_FillValue, data11@_FillValue,data11(0,:))
  data12(0,:) = where(data12(0,:).eq.sst_rsi1_djf@_FillValue, data12@_FillValue,data12(0,:))

  data1(0,:)   = abs(data1(0,:))
  data2(0,:)   = abs(data2(0,:))
  data3(0,:)   = abs(data3(0,:))
  data4(0,:)   = abs(data4(0,:))
  data04(0,:)   = abs(data04(0,:))
;======================================
  data5(0,:)   = abs(data5(0,:))
  data6(0,:)   = abs(data6(0,:))
  data7(0,:)   = abs(data7(0,:))
  data8(0,:)   = abs(data8(0,:))
;======================================
  data9(0,:)   = abs(data9(0,:))
  data10(0,:)  = abs(data10(0,:))
  data11(0,:)  = abs(data11(0,:))
  data12(0,:)  = abs(data12(0,:))

  cond1 = ind(data1(0,:).ne.0)
  cond2 = ind(data2(0,:).ne.0)
  cond3 = ind(data3(0,:).ne.0)
  cond4 = ind(data4(0,:).ne.0)
  cond41 = ind(data04(0,:).ne.0)
;=================================
  cond5 = ind(data5(0,:).ne.0)
  cond6 = ind(data6(0,:).ne.0)
  cond7 = ind(data7(0,:).ne.0)
  cond8 = ind(data8(0,:).ne.0)
  cond9 = ind(data9(0,:).ne.0)
  cond10 = ind(data10(0,:).ne.0)
  cond11 = ind(data11(0,:).ne.0)
  cond12 = ind(data12(0,:).ne.0)
  years = (/1981-1948,1982-1948,1983-1948,1984-1948,1997-1948,1996-1948,1997-1948,1998-1948,1999-1948,2014-1948,2015-1948,2016-1948,2017-1948/)
  print(years)

  cond01   = cond1(ind(.not.ismissing(get1Dindex(years,cond1))))
  cond02   = cond2(ind(.not.ismissing(get1Dindex(years,cond2))))
  cond03   = cond3(ind(.not.ismissing(get1Dindex(years,cond3))))
  cond04   = cond4(ind(.not.ismissing(get1Dindex(years,cond4))))
  cond041   = cond41(ind(.not.ismissing(get1Dindex(years,cond41))))

  cond05   = cond5(ind(.not.ismissing(get1Dindex(years,cond5))))
  cond06   = cond6(ind(.not.ismissing(get1Dindex(years,cond6))))
  cond07   = cond7(ind(.not.ismissing(get1Dindex(years,cond7))))
  cond08   = cond8(ind(.not.ismissing(get1Dindex(years,cond8))))
  cond09   = cond9(ind(.not.ismissing(get1Dindex(years,cond9))))
  cond010  = cond10(ind(.not.ismissing(get1Dindex(years,cond10))))
  cond011  = cond11(ind(.not.ismissing(get1Dindex(years,cond11)))) 
  cond012  = cond12(ind(.not.ismissing(get1Dindex(years,cond12)))) 
  print(cond01)
  n1  = dimsizes(cond01)
  n2  = dimsizes(cond02)
  n3  = dimsizes(cond03)
  n4  = dimsizes(cond04)
  n04  = dimsizes(cond041)
  n5  = dimsizes(cond05)
  n6  = dimsizes(cond06)
  n7  = dimsizes(cond07)
  n8  = dimsizes(cond08)
  n9  = dimsizes(cond09)
  n10 = dimsizes(cond010)
  n11 = dimsizes(cond011)
  n12 = dimsizes(cond012)

wks  = gsn_open_wks("pdf","/Users/xueaoyun/Documents/work/regime_shift/New/Obs/figures/map_pattern_regime_shift_evidence")
gsn_define_colormap(wks,"MPL_RdBu")
;----------------------------------------------------------------------
; Set some resources common to all map projections.
;----------------------------------------------------------------------
    res                      = True
   ; res@gsnMaximize          = True          ; maximize plot-frame
    res@vpWidthF             = 0.40
    res@vpHeightF            = 0.2
    res@vpXF                 = 0.3
    res@vpYF                 = 0.78
    res@gsnDraw              = False 
    res@gsnFrame             = False
    res@mpProjection           = "Robinson"
    res@mpFillOn               = True
   ; res@mpProjection = projections(1)
    res@mpMinLatF           = -90.0         ; map area
    res@mpMaxLatF           =  90.0         ; latitudes
    res@mpMinLonF           =   0.0         ; 
    res@mpMaxLonF           = 360.0         ; longitudes
    res@mpCenterLonF        = 180.0
   ; res@mpDataBaseVersion    = "MediumRes"    ; better map outlines
   ; res@mpGridAndLimbOn      = True           ; turn on limb and grid lines
  ;  res@mpGridLineColor      = -1             ; ...but don't draw grid lines...
    res@mpGridAndLimbOn     = True  
    res@mpPerimOn            = False          ; turn off map perimeter
    res@mpFillDrawOrder      = "PreDraw"     ; draw map fill last
    res@mpLandFillColor      = "gray95"
    res@mpOceanFillColor     = "White"
            ; turn on lat/lon lines
  ; res@mpGridLatSpacingF =  90               ; change latitude  line spacing
  ; res@mpGridLonSpacingF = 180.              ; change longitude line spacing
    res@mpGridLineColor   = "transparent"     ; trick ncl into drawing
    res@mpGridMaskMode    = "MaskLand"
    res@mpGeophysicalLineColor      = "black"       ; color of cont. outlines
    res@mpGeophysicalLineThicknessF = 1.5          ; thickness of outlines
    res@gsnRightString       = ""             ; turn off special titles
    res@gsnLeftString        = ""
    res@tiMainFontHeightF    = 0.05           ; change size of main title
    res@tmXTOn   = False
    res@tmYROn   = False
    res@gsnLeftString = ""
    res@gsnRightString = ""

    res@cnFillOn             = True           ; turn on contour fill
    res@cnLinesOn            = False          ; turn off contour lines
    res@cnLineLabelsOn       = False          ; turn off contour labels
    res@lbLabelBarOn         = False          ; turn off labelbar
    res@lbOrientation         = "Horizontal"
    res@lbBoxEndCapStyle      = "TriangleBothEnds"
    res@pmLabelBarWidthF      = 0.20               ; default is shorter
    res@pmLabelBarHeightF     = 0.05
    res@pmLabelBarOrthogonalPosF   = -0.01
    res@pmLabelBarParallelPosF   = 0.0 

    res@cnLevelSelectionMode = "ExplicitLevels"	
    res@cnLevels =(/-2.5,-2.0,-1.5,-1.0,0,1.0,1.5,2.0,2.5/)
    res@cnFillColors =(/120,110,100,90,0,0,50,40,30,20/)
    res@tmXTOn   = False
    res@tmYROn   = False
    res@gsnLeftString = ""
    res@gsnRightString = ""
    res@cnCellFillMissingValEdgeColor  = "paleturquoise"
  ; total numbers for four seasons 
  map  = gsn_csm_contour_map(wks,ep_ta,res)
  ; SAT 
  ;lon_new1 = (/190,180,130,35,100,270,290,85,305,345/)
  ;lat_new1 = (/35,-35,0,5,-37.5,25,-30,30,0,-40/)
  ;text1 = (/"CNP","SWP","MC","EAF","SEI","CAM","SAM","SAS","AMA","SAT"/)
 ; lon_new1 = (/180,130, 35,295,270/)
 ; lat_new1 = (/-35,  0, 5, -30,22/)
 ; text1 = (/"SWP","MC","EAF","SAM","CAM"/)
  lon_new1 = (/180,140, 35,40/)
  lat_new1 = (/-35,  0, 5,60/)
  text1 = (/"SWP","MC","EAF","NEU"/)
  ; soil 
 ; lon_new2 = (/20,90,100,135,245/)
 ; lat_new2 = (/-10,30,55,-25,50/)
 ; text2 = (/"SAF","SAS","NAS","CAU","WNA"/)
  lon_new2 = (/310,70,135,305/)
  lat_new2 = (/70,40,-25,-5/)
  text2 = (/"WGL","CAS","CAU","AMA"/)
  ; SST 
 ; lon_new3 = (/190,155,100,345,90/)
 ; lat_new3 = (/ 35,30,-37.5,-40,-10/)
 ; text3 = (/"CNP","KOE","SEI","SAT","CIO"/)
   lon_new3 = (/190,150,260,270/)
   lat_new3 = (/ 40,25,-35,25/)
 ; text3 = (/"CNP","SEI","WNP","SEP","GM"/)
  text3 = (/"CNP","WNP","SEP","GM"/)

  pmid = new(12,graphic)
    gsres = True
    gsres@gsMarkerIndex = 16          ; Use filled dots for markers.
    do i = 0, 3
        gsres@gsMarkerColor      = "lightskyblue"
        gsres@gsMarkerSizeF      = 0.02
        gsres@gsFillOpacityF     = 0.0
        pmid(i) = gsn_add_polymarker(wks,map,lon_new1(i),lat_new1(i),gsres)
    end do

    do i = 4,7
    gsres@gsMarkerColor      = "tan"
    gsres@gsMarkerSizeF      = 0.02
  ;  gsres@gsFillOpacityF     = 1.0
    pmid(i) = gsn_add_polymarker(wks,map,lon_new2(i-4),lat_new2(i-4),gsres)
    end do
    
    do i = 8, 11
    gsres@gsMarkerColor      = "coral"
    gsres@gsMarkerSizeF      = 0.02
 ;  gsres@gsFillOpacityF     = 0.2
    pmid(i) = gsn_add_polymarker(wks,map,lon_new3(i-8),lat_new3(i-8),gsres)
    end do


    txres                      = True                ; text mods desired
    txres@txFontHeightF        = 0.007               ; make smaller
    txres@txFontColor          = "black"
    txres@txFontThicknessF     = 6.0

    dum = new(12,graphic)
    do i=0,3
      dum(i) = gsn_add_text(wks,map,text1(i),lon_new1(i),lat_new1(i),txres)
    end do
    do i=4,7
      dum(i) = gsn_add_text(wks,map,text2(i-4),lon_new2(i-4),lat_new2(i-4),txres)
    end do
    do i=8,11
      dum(i) = gsn_add_text(wks,map,text3(i-8),lon_new3(i-8),lat_new3(i-8),txres)
    end do

    labels = (/"SAT","SST","Soil Moisture"/)
    gsres1               = True
    gsres1@gsMarkerIndex = 16          ; Use filled dots for markers.
    gsres1@gsMarkerColor      = "lightskyblue"
    gsres1@gsMarkerSizeF      = 0.015

    txres1                      = True             
    txres1@txFontHeightF        = 0.007               ; make smaller
    txres1@txFontColor          = "black"
    txres1@txFontThicknessF     = 10.0

    lon_new = (/100.,180.,260./)
    lat_new = (/-60.,-60.,-60./)
    colors = (/"lightskyblue","coral","tan"/)
    txt = new(3,graphic)
    dumm = new(3, graphic)
  do i = 0, 1
      gsres1@gsMarkerColor      = colors(i)
    txt(i) =  gsn_add_polymarker(wks,map,lon_new(i),lat_new(i),gsres1)
    dumm(i) = gsn_add_text(wks,map,labels(i),lon_new(i)+25,lat_new(i),txres1)
  end do
  gsres1@gsMarkerColor      = colors(2)
   txt(2) =  gsn_add_polymarker(wks,map,lon_new(2),lat_new(2),gsres1)
  dumm(2) = gsn_add_text(wks,map,labels(2),lon_new(2)+45,lat_new(2),txres1)
draw(map)
;=============================================================================================

triangle = NhlNewMarker(wks, "u",  34, 0.0, 0.0, 1.2, 1.0, 0.)

x0=ispan(0,74,1)
; draw the bar and time series
aresL               = True 
aresL@gsnFrame      = False
aresL@gsnDraw       = False
aresL@trXMinF       = 0
aresL@trXMaxF       = 75
aresL@vpWidthF      = 0.16
aresL@vpHeightF     = 0.09
aresL@tmXBMode   = "Explicit"          
aresL@tmXBValues      = ispan(2,72,10)
aresL@tmXBLabels      = (/"50","60","70","80","90","00","10","20"/)
aresL@tmXBMinorValues       = ispan(2,72,10)
aresL@tmXBLabelAngleF = 0
aresL@tiXAxisString = ""
aresL@tiXAxisFontHeightF       = 0.012
aresL@tmXBLabelFontHeightF     = 0.006
aresL@tmYLLabelFontHeightF     = 0.006
;aresL@tmXBLabelFont            = 4
;aresL@tmYLLabelFont            = 4
aresL@tmXBLabelFontThicknessF  = 1.0
aresL@tmYLLabelFontThicknessF  = 1.0 
aresL@tmBorderThicknessF       = 1.5
aresL@tmXBMajorLengthF         = 0.005
aresL@tmYLMajorLengthF         = 0.005
aresL@tmXBMinorLengthF         = 0.006
aresL@tmYLMinorLengthF         = 0.006 
aresL@tmXBMajorThicknessF      = 1.0
aresL@tmYLMajorThicknessF      = 1.0
aresL@tmXBMinorThicknessF      = 0.5
aresL@tmYLMinorThicknessF      = 0.5 
aresL@tmXBLabelDeltaF          = -0.2
aresL@tmYLLabelDeltaF          = -0.2
;aresL@xyDashPattern    = 0          ; Make curves all solid
;aresL@xyLineThicknessF = 2          ; make 2nd lines thicker
;aresL@xyLineColor      =(/"red"/)   
aresL@gsnYRefLine            = 0.              ; reference line
aresL@gsnYRefLineDashPattern = (/0/)
aresL@gsnXRefLine            =cond01             ; reference line
aresL@gsnXRefLineColors      = (/"black"/)
aresL@gsnXRefLineDashPatterns = (/0/)
aresL@gsnXRefLineThicknesses  = (/1.5/)
aresL@trYMinF       = -2
aresL@trYMaxF       = 2
aresL@tiYAxisString=""
aresL@tiYAxisFontHeightF     = 0.01
aresL@tiMainString          =" "
aresL@tiMainStringFontHeightF =0.0075

;aresL@tiMainStringFontColor  = "black"
aresL@tmYLMode        = "Explicit"
aresL@tmYLValues      = fspan(-2,2,5)
aresL@tmYLLabels      = (/"-2","-1","0","1","2"/)
aresL@tmYLLabelFontColor     = "blue"
aresL@tiYAxisFontColor       = "blue"
;aresL@gsnXYBarChart          =  True 
;aresL@gsnAboveYRefLineBarColors = "orchid1"
;aresL@gsnBelowYRefLineBarColors = "skyblue"
aresL@gsnAboveYRefLineColor = "orchid1"
aresL@gsnBelowYRefLineColor = "skyblue"
aresL@tmXTOn                    = False 
aresL@tmYROn                    = False 
aresL@vpXF      = 0.09
aresL@vpYF      = 0.91
aresL@tmXBLabels      = (/"50","60","70","80","90","00","10","20"/)
aresL@tiMainString  = "~F22~(A) SST-DJF-CNP (Central-North Pac)"
aresL@tiMainFontAspectF = 1.3
aresL@tiMainFontColor = "red"
;plot0   = gsn_csm_xy2(wks,x0,data1(1,:),data1(0,:),aresL, aresR)
plot1   = gsn_csm_xy(wks,x0,data1(1,:),aresL)
  dum1 = new(dimsizes(cond01),graphic)
  pmres = True
  do i=0,dimsizes(cond01)-1

    xmarkers1 = cond01(i)
    ymarkers1 = -1.8

    pmres@gsMarkerIndex = triangle 
    pmres@gsMarkerSizeF    = 0.015
    pmres@gsMarkerColor    = "blue"
    dum1(i) = gsn_add_polymarker(wks, plot1, xmarkers1, ymarkers1, pmres)
  end do
;=======================================================================
delete(aresL@tmYLValues)
delete(aresL@tmYLLabels)
delete(aresL@gsnXRefLine)
delete(aresL@gsnXRefLineColors)
delete(aresL@gsnXRefLineDashPatterns)
aresL@gsnXRefLine            = (/cond02,97-48/)             ; reference line
aresL@vpXF      = 0.31
aresL@vpYF      = 0.91
aresL@trYMinF       = -1
aresL@trYMaxF       = 1
aresL@tmYLValues      = (/-1,-0.5,0,0.5,1/)
aresL@tmYLLabels      = (/"-1","-0.5","0","0.5","1"/)
aresL@tmYROn                    = False 
aresL@tiMainString  = "~F22~(B) SST-DJF-WNP (Western North Pac)"
;plot1   = gsn_csm_xy2(wks,x0,data2(1,:),data2(0,:),aresL, aresR)
plot2   = gsn_csm_xy(wks,x0,data2(1,:),aresL)
dum2 = new(dimsizes(cond02)+1,graphic)
pmres = True
cond02_new = (/cond02,97-48/)
do i=0,dimsizes(cond02_new)-1

  xmarkers2 = cond02_new(i)
  ymarkers2 = -0.9

  pmres@gsMarkerIndex = triangle
  pmres@gsMarkerSizeF    = 0.015
  pmres@gsMarkerColor    = "blue"
  dum2(i) = gsn_add_polymarker(wks, plot2, xmarkers2, ymarkers2, pmres)
end do
;=====================================================================
delete(aresL@gsnXRefLine)
delete(aresL@gsnXRefLineColors)
delete(aresL@gsnXRefLineDashPatterns)
aresL@gsnXRefLine            = cond041       
aresL@vpXF      = 0.75
aresL@vpYF      = 0.91
delete(aresL@tmYLValues)
delete(aresL@tmYLLabels)
aresL@trYMinF       = -2
aresL@trYMaxF       = 2
aresL@tmYLValues      = (/-2,-1,0,1,2/)
aresL@tmYLLabels      = (/"-2","-1","0","1","2"/)
aresL@tiMainString  = "~F22~(D) SST-DJF-GM (Gulf of Mexico)"
;plot2   = gsn_csm_xy2(wks,x0,data3(1,:),data3(0,:),aresL, aresR)
plot3   = gsn_csm_xy(wks,x0,data04(1,:),aresL)
dum3 = new(dimsizes(cond041),graphic)
pmres = True
do i=0,dimsizes(cond041)-1

  xmarkers3 = cond041(i)
  ymarkers3 = -1.8

  pmres@gsMarkerIndex = triangle
  pmres@gsMarkerSizeF    = 0.015
  pmres@gsMarkerColor    = "blue"
  dum3(i) = gsn_add_polymarker(wks, plot3, xmarkers3, ymarkers3, pmres)
end do
;========================================================================
delete(aresL@gsnXRefLine)
delete(aresL@gsnXRefLineColors)
delete(aresL@gsnXRefLineDashPatterns)
print(cond04)
aresL@gsnXRefLine            = (/cond04/)  
delete(aresL@tmYLValues)
delete(aresL@tmYLLabels)
aresL@vpXF      = 0.53
aresL@vpYF      = 0.91
aresL@trYMinF       = -1
aresL@trYMaxF       = 1
aresL@tmYLValues      = (/-1,-0.5,0,0.5,1/)
aresL@tmYLLabels      = (/"-1","-0.5","0","0.5","1"/)
aresL@tmXBLabels      = (/"50","60","70","80","90","00","10","20"/)
aresL@tiMainString  = "~F22~(C) SST-SON-SEP (Southeastern Pac)"
;plot3   = gsn_csm_xy2(wks,x0,data4(1,:),data4(0,:),aresL, aresR)
plot4   = gsn_csm_xy(wks,x0,data4(1,:),aresL)
dum4 = new(dimsizes(cond04),graphic)
pmres = True
do i=0,dimsizes(cond04)-1
  xmarkers4 = cond04(i)
  ymarkers4 = -0.9

  pmres@gsMarkerIndex = triangle
  pmres@gsMarkerSizeF    = 0.015
  pmres@gsMarkerColor    = "blue"
  dum4(i) = gsn_add_polymarker(wks, plot4, xmarkers4, ymarkers4, pmres)
end do
;===============================================================
; plot SAT
delete(aresL@gsnXRefLine)
delete(aresL@gsnXRefLineColors)
delete(aresL@gsnXRefLineDashPatterns)
aresL@gsnXRefLine            = cond05  
delete(aresL@tmYLValues)
delete(aresL@tmYLLabels)
aresL@vpXF      = 0.75
aresL@vpYF      = 0.78
aresL@tiMainStringFontColor="blue"
aresL@tiMainString  = "~F22~(G) SAT-DJF-EAF (Eastern Africa)"
aresL@tmXBLabels      = (/"50","60","70","80","90","00","10","20"/)
aresL@trYMinF       = -2
aresL@trYMaxF       = 2
aresL@tmYLValues      = fspan(-2,2,5)
aresL@tmYLLabels      = (/"-2","-1","0","1","2"/)
aresL@tiMainFontColor = "blue"
;plot4   = gsn_csm_xy2(wks,x0,data5(1,:),data5(0,:),aresL, aresR)
plot5   = gsn_csm_xy(wks,x0,data5(1,:),aresL)
dum5    = new(dimsizes(cond05),graphic)
pmres = True
do i=0,dimsizes(cond05)-1
  xmarkers5 = cond05(i)
  ymarkers5 = -1.8

  pmres@gsMarkerIndex = triangle
  pmres@gsMarkerSizeF    = 0.015
  pmres@gsMarkerColor    = "blue"
  dum5(i) = gsn_add_polymarker(wks, plot5, xmarkers5, ymarkers5, pmres)
end do
;===============================================================

delete(aresL@gsnXRefLine)
delete(aresL@gsnXRefLineColors)
delete(aresL@gsnXRefLineDashPatterns)
aresL@gsnXRefLine            = cond06  
delete(aresL@tmYLValues)
delete(aresL@tmYLLabels)

aresL@vpXF      = 0.09
aresL@vpYF      = 0.65
aresL@trYMinF       = -1
aresL@trYMaxF       = 1
aresL@tmYLValues      = (/-1,-0.5,0,0.5,1/)
aresL@tmYLLabels      = (/"-1","-0.5","0","0.5","1"/)
aresL@tiMainString  = "~F22~(F) SAT-SON-MC (Maritime Continent)"
;plot6   = gsn_csm_xy2(wks,x0,data6(1,:),data6(0,:),aresL, aresR)
plot6   = gsn_csm_xy(wks,x0,data6(1,:),aresL)
dum6    = new(dimsizes(cond06),graphic)
pmres = True
do i=0,dimsizes(cond06)-1
  xmarkers6 = cond06(i)
  ymarkers6 = -0.9

  pmres@gsMarkerIndex = triangle 
  pmres@gsMarkerSizeF    = 0.015
  pmres@gsMarkerColor    = "blue"
  dum6(i) = gsn_add_polymarker(wks, plot6, xmarkers6, ymarkers6, pmres)
end do
;==============================================================
delete(aresL@tmYLValues)
delete(aresL@tmYLLabels)
delete(aresL@gsnXRefLine)
delete(aresL@gsnXRefLineColors)
delete(aresL@gsnXRefLineDashPatterns)
aresL@gsnXRefLine            = cond07  
aresL@vpXF      = 0.75
aresL@vpYF      = 0.65
aresL@trYMinF       = -5
aresL@trYMaxF       = 5
aresL@tmYLValues      = (/-5,-4,-3,-2,-1,0,1,2,3,4,5/)
aresL@tmYLLabels      = (/"-5","-4","-3","-2","-1","0","1","2","3","4","5"/)
aresL@tiMainString  = "~F22~(H) SAT-DJF-NEU (Northern Europe)"
;plot7   = gsn_csm_xy2(wks,x0,data7(1,:),data7(0,:),aresL, aresR)
plot7   = gsn_csm_xy(wks,x0,data7(1,:),aresL)
dum7    = new(dimsizes(cond07),graphic)

pmres = True
do i=0,dimsizes(cond07)-1
  xmarkers7 = cond07(i)
  ymarkers7 = -4.5

  pmres@gsMarkerIndex = triangle
  pmres@gsMarkerSizeF    = 0.015
  pmres@gsMarkerColor    = "blue"
  dum7(i) = gsn_add_polymarker(wks, plot7, xmarkers7, ymarkers7, pmres)
end do
;=====================================================================
delete(aresL@tmYLValues)
delete(aresL@tmYLLabels)
delete(aresL@gsnXRefLine)
delete(aresL@gsnXRefLineColors)
delete(aresL@gsnXRefLineDashPatterns)
aresL@gsnXRefLine            = cond08
aresL@vpXF      = 0.09
aresL@vpYF      = 0.78
aresL@trYMinF       = -2
aresL@trYMaxF       = 2
aresL@tmYLValues      = (/-1,-0.5,0,0.5,1/)*2
aresL@tmYLLabels      = (/"-2","-1","0","1","2"/)
aresL@tiMainString  = "~F22~(E) SAT-JJA-SWP (Southwestern Pac)"
;plot8   = gsn_csm_xy2(wks,x0,data8(1,:),data8(0,:),aresL, aresR)
plot8   = gsn_csm_xy(wks,x0,data8(1,:),aresL)
dum8    = new(dimsizes(cond08),graphic)
pmres = True
do i=0,dimsizes(cond08)-1
  xmarkers8 = cond08(i)
  ymarkers8 = -0.9*2

  pmres@gsMarkerIndex = triangle
  pmres@gsMarkerSizeF    = 0.015
  pmres@gsMarkerColor    = "blue"
  dum8(i) = gsn_add_polymarker(wks, plot8, xmarkers8, ymarkers8, pmres)
end do
;===================================================================
;plot soil moisture
delete(aresL@gsnXRefLine)
delete(aresL@gsnXRefLineColors)
delete(aresL@gsnXRefLineDashPatterns)
aresL@gsnXRefLine            = cond09
aresL@vpXF      = 0.09
aresL@vpYF      = 0.52
aresL@trYMinF       = -120
aresL@trYMaxF       = 120
aresL@tmYLValues      = fspan(-120,120,5)
aresL@tiMainStringFontColor="sienna"
aresL@tmXBLabels      = (/"50","60","70","80","90","00","10","20"/)
aresL@tmYLLabels      = (/"-120","-60","0","60","120"/)
aresL@tiMainString  = "~F22~(I) Soil-JJA-CAU (Central Australia)"
aresL@tiMainFontColor = "sienna"
;plot9   = gsn_csm_xy2(wks,x0,data9(1,:),data9(0,:),aresL, aresR)
plot9   = gsn_csm_xy(wks,x0,data9(1,:),aresL)
dum9    = new(dimsizes(cond09),graphic)
pmres = True
do i=0,dimsizes(cond09)-1
  xmarkers9 = cond09(i)
  ymarkers9 = -0.9*120

  pmres@gsMarkerIndex = triangle
  pmres@gsMarkerSizeF    = 0.015
  pmres@gsMarkerColor    = "blue"
  dum9(i) = gsn_add_polymarker(wks, plot9, xmarkers9, ymarkers9, pmres)
end do
;=================================================================
delete(aresL@gsnXRefLine)
delete(aresL@gsnXRefLineColors)
delete(aresL@gsnXRefLineDashPatterns)
aresL@gsnXRefLine            = cond010  
aresL@vpXF      = 0.31
aresL@vpYF      = 0.52
aresL@tmXBLabels      = (/"50","60","70","80","90","00","10","20"/)
aresL@tiMainString  = "~F22~(J) Soil-annual-CAS (Central Asia)"
;plot10   = gsn_csm_xy2(wks,x0,data10(1,:),data10(0,:),aresL, aresR)
plot10   = gsn_csm_xy(wks,x0,data10(1,:),aresL)
dum10    = new(dimsizes(cond010),graphic)
pmres = True
do i=0,dimsizes(cond010)-1
  xmarkers10 = cond010(i)
  ymarkers10 = -0.9*120

  pmres@gsMarkerIndex = triangle 
  pmres@gsMarkerSizeF    = 0.015
  pmres@gsMarkerColor    = "blue"
  dum10(i) = gsn_add_polymarker(wks, plot10, xmarkers10, ymarkers10, pmres)
end do
;========================================================================
delete(aresL@gsnXRefLine)
delete(aresL@gsnXRefLineColors)
delete(aresL@gsnXRefLineDashPatterns)
print(cond011)
aresL@gsnXRefLine            = (/35,69,97-48/)
aresL@vpXF      = 0.53
aresL@vpYF      = 0.52
aresL@tiMainString  = "~F22~(K) Soil-SON-WGL ((West Greenland))"
;plot10  = gsn_csm_xy2(wks,x0,data11(1,:),data11(0,:),aresL, aresR)
plot11   = gsn_csm_xy(wks,x0,data11(1,:),aresL)
cond011_new = (/35,69,97-48/)
dum11    = new(dimsizes(cond011_new),graphic)
pmres = True
do i=0,dimsizes(cond011_new)-1
  xmarkers11 = cond011_new(i)
  ymarkers11 = -0.9*120

  pmres@gsMarkerIndex = triangle
  pmres@gsMarkerSizeF    = 0.015
  pmres@gsMarkerColor    = "blue"
  dum11(i) = gsn_add_polymarker(wks, plot11, xmarkers11, ymarkers11, pmres)
end do
;===============================================================================
delete(aresL@gsnXRefLine)
delete(aresL@gsnXRefLineColors)
delete(aresL@gsnXRefLineDashPatterns)
delete(aresL@tmYLLabels)
delete(aresL@tmYLValues)
aresL@gsnXRefLine            = cond012  
aresL@vpXF      = 0.75
aresL@vpYF      = 0.52
aresL@trYMinF       = -150
aresL@trYMaxF       = 150
aresL@tmYLValues      = fspan(-150,150,7)
aresL@tmYLLabels      = (/"-150","-100","-50","0","50","100","150"/)
aresL@tiMainString  = "~F22~(L) Soil-annual-amazon (Amazon)"
;plot12  = gsn_csm_xy2(wks,x0,data12(1,:),data12(0,:),aresL, aresR)
plot12   = gsn_csm_xy(wks,x0,data12(1,:),aresL)
dum12    = new(dimsizes(cond012),graphic)
pmres = True
do i=0,dimsizes(cond012)-1
  xmarkers12 = cond012(i)
  ymarkers12 = -0.9*150

  pmres@gsMarkerIndex = triangle
  pmres@gsMarkerSizeF    = 0.015
  pmres@gsMarkerColor    = "blue"
  dum12(i) = gsn_add_polymarker(wks, plot12, xmarkers12, ymarkers12, pmres)
end do
;=================================================================================
draw(plot1)
draw(plot2)
draw(plot3)
draw(plot4)
draw(plot5)
draw(plot6)
draw(plot7)
draw(plot8)
draw(plot9)
draw(plot10)
draw(plot11)
draw(plot12)
;===================================================================
end



  
