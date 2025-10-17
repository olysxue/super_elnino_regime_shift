load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

function get_title_position1(plot)
local vpy, vpx, vph
begin
  getvalues plot
      "vpYF"       : vpy
      "vpXF"       : vpx
      "vpHeightF"  : vph
  end getvalues
  pos = (/vpx+0.15/2,vpy+0.03/)
  return(pos)
end



begin
;============================================================================================

  sat_y1  = new((/2,5/), float)
  sat_y1(0,0) = 1.416089
  sat_y1(0,1) = 1.447912
  sat_y1(0,2) = 1.467572
  sat_y1(0,3) = 1.487168
  sat_y1(0,4) = 1.510134
  
  sat_y1(1,0) = 1.661704 
  sat_y1(1,1) = 1.728422
  sat_y1(1,2) = 1.77369
  sat_y1(1,3) = 1.819598
  sat_y1(1,4) = 1.875973 
;============================================================================================

  sst_y1  = new((/2,5/), float)
  sst_y1(0,0) = 1.35543 
  sst_y1(0,1) = 1.374752
  sst_y1(0,2) = 1.389183
  sst_y1(0,3) = 1.403979
  sst_y1(0,4) = 1.426149 

  sst_y1(1,0) = 1.563915 
  sst_y1(1,1) = 1.609235
  sst_y1(1,2) = 1.642808
  sst_y1(1,3) = 1.675877
  sst_y1(1,4) = 1.723884 
;============================================================================================
; Soil
  soil_y1  = new((/2,5/), float)
  soil_y1(0,0) =  0.273682
  soil_y1(0,1) =  0.2758267
  soil_y1(0,2) =  0.2776575
  soil_y1(0,3) =  0.279511
  soil_y1(0,4) =  0.2823261

  soil_y1(1,0) =  0.2901407
  soil_y1(1,1) =  0.2980953
  soil_y1(1,2) =  0.3030504
  soil_y1(1,3) =  0.3080258
  soil_y1(1,4) =  0.3161311
;==============================================================================================

sat_y2  = new((/2,5/), float)
sat_y2(0,0) =  1.524771
sat_y2(0,1) =  1.546949
sat_y2(0,2) =  1.566438
sat_y2(0,3) =  1.58616
sat_y2(0,4) =  1.614545

sat_y2(1,0) =  1.80505
sat_y2(1,1) =  2.090765
sat_y2(1,2) =  2.271671
sat_y2(1,3) =  2.464148
sat_y2(1,4) =  2.747574

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SST

sst_y2  = new((/2,5/), float)
sst_y2(0,0) =  1.60914
sst_y2(0,1) =  1.637146
sst_y2(0,2) =  1.655188
sst_y2(0,3) =  1.673348
sst_y2(0,4) =  1.698775

sst_y2(1,0) =  1.836562
sst_y2(1,1) =  2.052951
sst_y2(1,2) =  2.212148
sst_y2(1,3) =  2.37871
sst_y2(1,4) =  2.603869

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
soil_y2  = new((/2,5/), float)
soil_y2(0,0) = 0.3637861
soil_y2(0,1) = 0.3671262
soil_y2(0,2) = 0.3696403
soil_y2(0,3) = 0.3721274
soil_y2(0,4) = 0.375254

soil_y2(1,0) = 0.3199122+0.01
soil_y2(1,1) = 0.3397308+0.01
soil_y2(1,2) = 0.3584305+0.01
soil_y2(1,3) = 0.3774258+0.01
soil_y2(1,4) = 0.3997416+0.01
;===========================================================================
x = (/1,2/)
wks                        = gsn_open_wks ("pdf", get_script_prefix_name)       ; open workstation
; plot boxplot
  xyres                    = True                 ; plot mods desired         
  xyres@trYMinF            = 1.0                    ; bring bars down to zero
  xyres@trYMaxF            = 3.0                  ; bring bars down to zero
  xyres@trXMinF            = 0.5                 ; adds space on either end
  xyres@trXMaxF            = 2.5                 ; of the 1st and last bars 
  xyres@vpWidthF           = 0.28
  xyres@vpHeightF          = 0.15
  xyres@vpXF               = 0.25 
  xyres@vpYF               = 0.90
  
  xyres@tmXBOn                = True                 ; turn off tickmarks at bot
  xyres@tmYLOn                = True 
  xyres@tmYROn                = False                 ; turn off tickmarks at bot
  xyres@tmXTOn                = False                 ; turn off tickmarks at bot
  xyres@tmXTBorderOn          = True
  xyres@tmYRBorderOn          = True
 ; xyres@tmXBMode              = "Explicit"         ; explicit labels
 ; xyres@tmXBValues            = (/1,2/)
  xyres@tmXBLabels            = (/"~F22~non-Super El Ni"+"n~H-13V2F35~D~FV-2H3~"+"~F22~o", "~F22~Super El Ni"+"n~H-13V2F35~D~FV-2H3~"+"~F22~o"/)
  xyres@tmXBLabelAngleF       = 0
 ; xyres@tmXBLabelJust         = "CenterRight"
  xyres@tmYLMode              = "Explicit"         ; explicit labels
  xyres@tmYLValues            = (/1.0, 1.5, 2.0, 2.5, 3.0/)
  xyres@tmYLLabels            = (/"1.0","1.5","2.0","2.5","3.0"/)
  xyres@tmYLMinorOn           = False
  xyres@tmYLMinorValues       = fspan(1.0,3.0,21)

  xyres@tmXBLabelFontHeightF = 0.007
  xyres@tmYLLabelFontHeightF = 0.007
  xyres@tmXBLabelFontThicknessF  = 1.5
  xyres@tmYLLabelFontThicknessF  = 1.5 
  xyres@tmBorderThicknessF       = 1.5
  xyres@tmXBMajorLengthF         = 0.006
  xyres@tmYLMajorLengthF         = 0.006
  xyres@tmXBMajorOutwardLengthF  = 0.006
  xyres@tmYLMajorOutwardLengthF  = 0.006
  xyres@tmXBMinorLengthF         = 0.004
  xyres@tmYLMinorLengthF         = 0.004
  xyres@tmXBMinorOutwardLengthF  = 0.004
  xyres@tmYLMinorOutwardLengthF  = 0.004 
  xyres@tmXBMajorThicknessF      = 1.5
  xyres@tmYLMajorThicknessF      = 1.5
  xyres@tmXBMinorThicknessF      = 1.5
  xyres@tmYLMinorThicknessF      = 1.5 
 ;xyres@tmXBLabelDeltaF          = -0.3
 ;xyres@tmYLLabelDeltaF          = -0.3 
  xyres@tmXTOn                   = False
  xyres@tmYROn                   = False

  xyres@tiMainString             = " "
  xyres@tiMainFontHeightF        = 0.008
  xyres@tiMainFont               = 22
  xyres@tiMainOffsetYF           = -0.005
; lres
  lres                           = True
  lres@gsLineThicknessF          = 1
  lres@gsLineDashPattern         = 0
  ; opti 
   opti                          = True
   opti@boxWidth                 = 0.2
   opti@boxColors                = (/"blue", "blue"/)
   plot1 = boxplot(wks, x, sat_y1, opti, xyres, lres)
   opti@boxColors                = (/"red", "red"/)
   plot01 = boxplot(wks, x, sat_y2, opti, xyres, lres)
   overlay(plot1, plot01)
   xcoords1 = 2-opti@boxWidth/2-0.05
   ycoords1 = sat_y1(1,2)
  
   xcoords2 = 2-opti@boxWidth/2-0.05
   ycoords2 = sat_y2(1,2) 

 cmap = (/                         \  
     (/ 1.0, 1.0, 1.0 /),    \  ; color index 0 - white
     (/ 0.0, 0.0, 0.0 /),    \  ; color index 1 - black
     (/ 1.0, 0.6, 0.6 /),    \  ; color index 2 - red 
     (/ 0.5, 0.7, 1.0/),    \  ; color index 3 - blue
     (/ 0.0, 1.0, 0.0 /),    \  ; color index 4 - green
     (/ 0.9, 0.9, 0.9 /),    \  ; color index 5 - gray
     (/ 1.0, 1.0, 1.0 /)     \
  /)

  gsn_define_colormap(wks, cmap)  
  
  size = 0.031         ;  applies to all subsequent arrows.
  wmsetp("ARS - arrow size",size)
  wmsetp("AWC - color index for interior of arrow",3)
  wmsetp("ARD - arrow direction",13.)
  wmsetp("AOC - color index for arrow outline",3)
  wmsetp("ARL - scale factor for length of arrow tail",3.35)
  wmlabs(wks,xcoords1,ycoords1,"Arrow")

      
  size = 0.031         ;  applies to all subsequent arrows.
  wmsetp("ARS - arrow size",size)
  wmsetp("AWC - color index for interior of arrow",2)
  wmsetp("ARD - arrow direction",28.)
  wmsetp("AOC - color index for arrow outline",2)
  wmsetp("ARL - scale factor for length of arrow tail",3.65)
  wmlabs(wks,xcoords2,ycoords2,"Arrow")

  txres                    = True
  txres@txFontHeightF      = 0.0065
  txres@txAngleF           = 13.
  txres@txFontColor        = "blue"

dum1  = gsn_add_text(wks, plot1, "~F22~Present Day",   1.5, sat_y1(0,0)+0.08, txres)
  txres@txFontColor        = "red"
  txres@txAngleF      = 28.
dum2  = gsn_add_text(wks, plot1, "~F22~Warm Climate",  1.5-0.1, sat_y2(0,4)+0.40, txres)

txres1 = True
txres1@txFontColor        = "black"
txres1@txFontHeightF      =  0.010
txres1@txAngleF      = 0.
gsn_text_ndc(wks, "~F22~(a)  SAT", 0.28, 0.91, txres1)
;=======================================================================
draw(plot1)
   ; draw SST 
   xyres@trYMinF            = 1.0                    ; bring bars down to zero
   xyres@trYMaxF            = 3.0                  ; bring bars down to zero
   xyres@trXMinF            = 0.5                 ; adds space on either end
   xyres@trXMaxF            = 2.5                 ; of the 1st and last bars 
   xyres@vpWidthF           = 0.28
   xyres@vpHeightF          = 0.15
   xyres@vpXF               = 0.25 
   xyres@vpYF               = 0.70
   xyres@tiMainString             = " "
   opti                          = True
   opti@boxWidth                 = 0.2
   opti@boxColors                = (/"blue", "blue"/)
   plot2 = boxplot(wks, x, sst_y1, opti, xyres, lres)
   opti@boxColors                = (/"red", "red"/)
   plot02 = boxplot(wks, x, sst_y2, opti, xyres, lres)
   overlay(plot2, plot02)

   xcoords3 = 2-opti@boxWidth/2-0.05
   ycoords3 = sst_y1(1,2)
  
   xcoords4 = 2-opti@boxWidth/2-0.05
   ycoords4 = sst_y2(1,2) 

  cmap = (/                         \  
      (/ 1.0, 1.0, 1.0 /),    \  ; color index 0 - white
      (/ 0.0, 0.0, 0.0 /),    \  ; color index 1 - black
      (/ 1.0, 0.6, 0.6 /),    \  ; color index 2 - red 
      (/ 0.5, 0.7, 1.0/),    \  ; color index 3 - blue
      (/ 0.0, 1.0, 0.0 /),    \  ; color index 4 - green
      (/ 0.9, 0.9, 0.9 /),    \  ; color index 5 - gray
      (/ 1.0, 1.0, 1.0 /)     \
   /)
   gsn_define_colormap(wks,cmap)
  
  size = 0.031         ;  applies to all subsequent arrows.
  wmsetp("ARS - arrow size",size)
  wmsetp("AWC - color index for interior of arrow",3)
  wmsetp("ARD - arrow direction",10.5)
  wmsetp("AOC - color index for arrow outline",3)
  wmsetp("ARL - scale factor for length of arrow tail",3.3)
  wmlabs(wks,xcoords3,ycoords3,"Arrow")

      
  size = 0.031         ;  applies to all subsequent arrows.
  wmsetp("ARS - arrow size",size)
  wmsetp("AWC - color index for interior of arrow",2)
  wmsetp("ARD - arrow direction",22.)
  wmsetp("AOC - color index for arrow outline",2)
  wmsetp("ARL - scale factor for length of arrow tail",3.5)
  wmlabs(wks,xcoords4,ycoords4,"Arrow")


  txres                    = True
  txres@txFontHeightF      = 0.0065
  txres@txAngleF           = 10.5
  txres@txFontColor        = "blue"

dum01  = gsn_add_text(wks, plot2, "~F22~Present Day",   1.5, sst_y1(0,0)+0.04, txres)
  txres@txFontColor        = "red"
  txres@txAngleF      = 22.
dum02  = gsn_add_text(wks, plot2, "~F22~Warm Climate",  1.4, sst_y2(0,4)+0.35, txres)

txres1 = True
txres1@txFontColor        = "black"
txres1@txFontHeightF      =  0.010
txres1@txAngleF      = 0.
gsn_text_ndc(wks, "~F22~(b)  SST", 0.28, 0.71, txres1)
   delete(xyres@tmYLValues)
   delete(xyres@tmYLLabels)
draw(plot2)
;=====================================================================
    ; draw Soil
    xyres@trYMinF            = 0.2                    ; bring bars down to zero
    xyres@trYMaxF            = 0.5                  ; bring bars down to zero
    xyres@tmYLMode           = "Explicit"         ; explicit labels
    xyres@tmYLValues         = (/0.2, 0.3, 0.4, 0.5/)
    xyres@tmYLLabels         = (/"0.2","0.3","0.4","0.5"/)
    xyres@trXMinF            = 0.5                 ; adds space on either end
    xyres@trXMaxF            = 2.5                 ; of the 1st and last bars 
    xyres@vpWidthF           = 0.28
    xyres@vpHeightF          = 0.15
    xyres@vpXF               = 0.25 
    xyres@vpYF               = 0.50
    xyres@tiMainString             = " "

    opti                          = True
    opti@boxWidth                 = 0.2
    opti@boxColors                = (/"blue", "blue"/)
    plot3 = boxplot(wks, x, soil_y1, opti, xyres, lres)
    opti@boxColors                = (/"red", "red"/)
    plot03 = boxplot(wks, x, soil_y2, opti, xyres, lres)
    overlay(plot3, plot03)
     
    xcoords5 = 2-opti@boxWidth/2-0.05
    ycoords5 = soil_y1(1,2)
   
    xcoords6 = 2-opti@boxWidth/2-0.05
    ycoords6 = soil_y2(1,2)
   
   cmap = (/                         \  
       (/ 1.0, 1.0, 1.0 /),    \  ; color index 0 - white
       (/ 0.0, 0.0, 0.0 /),    \  ; color index 1 - black
       (/ 1.0, 0.6, 0.6 /),    \  ; color index 2 - red 
       (/ 0.5, 0.7, 1.0/),    \  ; color index 3 - blue
       (/ 0.0, 1.0, 0.0 /),    \  ; color index 4 - green
       (/ 0.9, 0.9, 0.9 /),    \  ; color index 5 - gray
       (/ 1.0, 1.0, 1.0 /)     \
    /)
    gsn_define_colormap(wks,cmap)
   
   size = 0.031         ;  applies to all subsequent arrows.
   wmsetp("ARS - arrow size",size)
   wmsetp("AWC - color index for interior of arrow",3)
   wmsetp("ARD - arrow direction",7.)
   wmsetp("AOC - color index for arrow outline",3)
   wmsetp("ARL - scale factor for length of arrow tail",3.25)
   wmlabs(wks,xcoords5,ycoords5,"Arrow")

   size = 0.031         ;  applies to all subsequent arrows.
   wmsetp("ARS - arrow size",size)
   wmsetp("AWC - color index for interior of arrow",2)
   wmsetp("ARD - arrow direction",-0.5)
   wmsetp("AOC - color index for arrow outline",2)
   wmsetp("ARL - scale factor for length of arrow tail",3.25)
   wmlabs(wks,xcoords6,ycoords6,"Arrow")

   txres                    = True
   txres@txFontHeightF      = 0.0065
   txres@txAngleF           = 7.
   txres@txFontColor        = "blue"
 
 dum001  = gsn_add_text(wks, plot3, "~F22~Present Day",   1.5, soil_y1(0,0)-0.002, txres)
   txres@txFontColor        = "red"
   txres@txAngleF      = -2.9
 dum002  = gsn_add_text(wks, plot3, "~F22~Warm Climate",  1.45, soil_y2(0,4)+0.012, txres)

 txres1 = True
txres1@txFontColor        = "black"
txres1@txFontHeightF      =  0.010
txres1@txAngleF      = 0.
gsn_text_ndc(wks, "~F22~(c) Soil Moisture", 0.307, 0.51, txres1)

draw(plot3)

;system("/Users/xueaoyun/ImageMagick-7.0.10/bin/magick convert -density 720 /Users/xueaoyun/Documents/work/regime_shift/New/model/scripts/Figure4_new.pdf /Users/xueaoyun/Documents/work/regime_shift/New/model/scripts/Figure4_new.jpg")
;system("/Users/xueaoyun/ImageMagick-7.0.10/bin/magick convert -trim /Users/xueaoyun/Documents/work/regime_shift/New/model/scripts/Figure4_new.jpg /Users/xueaoyun/Documents/work/regime_shift/New/model/scripts/Figure4_new.jpg")

end 
