load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

begin

  ncdir1   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/hist/sat/"
  ncdir2   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/hist/sst/"
  ncdir3   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/hist/soil/soilw/"  
  ncdir4   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/nino/"
  L        = 10

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  fili     = systemfunc("cd "+ncdir1+" ; ls rsi_seasonally_sat_cesm2_hist_ensemble*.nc")
  nfili    = dimsizes(fili)
  print("nfili ="+nfili)

  do i         = 0,nfili-1
 ;   print("i   = "+i)
 ;   print("file= "+fili(i))
    f          = addfile(ncdir1+fili(i),"r")
    rsi_son    = f->rsi_son
    rsi_djf    = f->rsi_djf
    rsi_mam    = f->rsi_mam    
    rsi_jja    = f->rsi_jja
    nlat       = dimsizes(f->lat)
    nlon       = dimsizes(f->lon)
    ntime      = dimsizes(f->time)
    rsi_son    = where(rsi_son.ne.0, 1.0, rsi_son)
    rsi_djf    = where(rsi_djf.ne.0, 1.0, rsi_djf)
    rsi_mam    = where(rsi_mam.ne.0, 1.0, rsi_mam)
    rsi_jja    = where(rsi_jja.ne.0, 1.0, rsi_jja)
    rsi_sum    = rsi_son
    rsi_sum    = ((/rsi_son/)+(/rsi_djf/)+(/rsi_mam/)+(/rsi_jja/))    

    if(i.eq.0) then 
      rsi_tot      = new((/nfili,ntime,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot!0    = "ensemble"
      rsi_tot!1    = "time"
      rsi_tot!2    = "lat"
      rsi_tot!3    = "lon"
      rsi_tot&time = rsi_sum&time      
      rsi_tot&lat  = rsi_sum&lat    
      rsi_tot&lon  = rsi_sum&lon
      rsi_tot_sen      = new((/nfili,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot_sen!0    = "ensemble"
      rsi_tot_sen!1    = "lat"
      rsi_tot_sen!2    = "lon"     
      rsi_tot_sen&lat  = rsi_sum&lat    
      rsi_tot_sen&lon  = rsi_sum&lon
      rsi_tot_all      = rsi_tot_sen

      rsi_tot_ren      = new((/nfili,nlat,nlon/),"float",rsi_sum@_FillValue) ; ren regular el nino
      rsi_tot_ren!0    = "ensemble"
      rsi_tot_ren!1    = "lat"
      rsi_tot_ren!2    = "lon"     
      rsi_tot_ren&lat  = rsi_sum&lat    
      rsi_tot_ren&lon  = rsi_sum&lon

      rsi_tot_rla      = new((/nfili,nlat,nlon/),"float",rsi_sum@_FillValue); rla regular la Nina
      rsi_tot_rla!0    = "ensemble"
      rsi_tot_rla!1    = "lat"
      rsi_tot_rla!2    = "lon"     
      rsi_tot_rla&lat  = rsi_sum&lat    
      rsi_tot_rla&lon  = rsi_sum&lon

    end if
    rsi_tot(i,:,:,:) = (/rsi_sum(:,:,:)/)
    masks            = ntime-(L/2)
    maske            = ntime-1
    rsi_tot(i,masks:maske,:,:) = 0

    f1          = addfile(ncdir4+"nino_cesm2_hist.nc","r")
    nino        = f1->nino(i,:)
    nino(0:L-1)         = 0  
    nino(masks-1:maske) = 0
    sen1        = ind(nino.ge.2.0)
    sen2        = sen1+1
    sen         = array_append_record(sen1, sen2, 0)
    ren         = ind(nino.ge.0.5 .and. nino.le.1.5)
    rla         = ind(nino.le.-0.5 .and. nino.ge.-1.0)

    rsi_tot_sen(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,sen,:,:),0)
    rsi_tot_ren(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,ren,:,:),0)
    rsi_tot_rla(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,rla,:,:),0)
    rsi_tot_all(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,:,:,:),0)
    delete([/sen1,sen2,sen, ren, rla/])
  end do
  printVarSummary(rsi_tot_sen)
  delete([/nino/])

  f2       = addfile(ncdir4+"nino_cesm2_hist.nc","r")
  nino     = f2->nino(0:nfili-1,L:masks-2)
  nino1D   = ndtooned(nino)
  sen      = ind(nino1D.ge.2.0)
  ren      = ind(nino1D.ge.0.5 .and. nino1D.le.1.5)
  rla      = ind(nino1D.le.-0.5  .and. nino1D.ge.-1.0)
  nsen     = dimsizes(sen)*2
  nren     = dimsizes(ren)
  nrla     = dimsizes(rla)
  nall     = nfili*(ntime-(L+(L/2)))

  rsp_all   = dim_sum_n_Wrap(rsi_tot_all,0)
  rsp_all   = rsp_all/(4.0*nall)*100.0
  rsp_sen   = dim_sum_n_Wrap(rsi_tot_sen,0)
  rsp_sen   = rsp_sen/(4.0*nsen)*100.0
  rsp_ren   = dim_sum_n_Wrap(rsi_tot_ren,0)
  rsp_ren   = rsp_ren/(4.0*nren)*100.0
  rsp_rla   = dim_sum_n_Wrap(rsi_tot_rla,0)
  rsp_rla   = rsp_rla/(4.0*nrla)*100.0

  rsp_sen_all = (rsp_sen-rsp_all)
  rsp_ren_all = (rsp_ren-rsp_all)
  rsp_rla_all = (rsp_rla-rsp_all)
  copy_VarCoords(rsp_sen, rsp_sen_all)
  copy_VarCoords(rsp_ren, rsp_ren_all)
  copy_VarCoords(rsp_rla, rsp_rla_all)
printMinMax(rsp_sen_all, 0)
  wgty1                   = cos(get_d2r(1.0)*f->lat)
  gm_ratio_all            = wgt_areaave_Wrap(rsp_all(:,:), wgty1, 1.0, 0)
  gm_ratio_sen            = wgt_areaave_Wrap(rsp_sen(:,:), wgty1, 1.0, 0)
  gm_ratio_sen_all        = wgt_areaave_Wrap(rsp_sen_all(:,:), wgty1, 1.0, 0)
  print("gm_ratio_all     = "+gm_ratio_all)
  print("gm_ratio_sen     = "+gm_ratio_sen)
  print("gm_ratio_sen_all = "+gm_ratio_sen_all)
  print("*************************************")

  delete([/rsi_son, rsi_djf, rsi_mam, rsi_jja, rsi_sum, rsi_tot/])
  delete([/rsi_tot_ren, rsi_tot_rla, rsi_tot_sen, rsi_tot_all/])
;************************************************
; plot 

  wks                        = gsn_open_wks ("pdf", get_script_prefix_name)       ; open workstation
  gsn_define_colormap(wks,"GMT_polar")
  gsn_define_colormap(wks,"MPL_bwr")       
  res                        = True 
  res@gsnDraw                = False              
  res@gsnFrame               = False 
  res@gsnAddCyclic           = True
   
  res@mpShapeMode            = "FreeAspect" 
  res@vpHeightF              = 0.10
  res@vpWidthF               = 0.20
  res@vpXF                   = 0.04                      ; Move slightly to left
  res@vpYF                   = 0.95

  res@cnFillOn               = True                 ; turn on color fill
  res@cnLevelSelectionMode   = "ExplicitLevels"
 ;res@cnLevels               = (/1,4,7,10,13,16,19,22,25,28/)
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.25
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
 ; res@cnLevels               = fspan(0.1,5.0,50)
 ; res@cnFillColors           = ispan(66,116,1)
  res@cnLevels                 = fspan(1.1,4.0,30)
  res@cnFillColors           = ispan(66,126,2)
 ;res@cnLevels               = (/1,3,5,7,9,11,13,15/)  
 ;res@cnFillColors           = (/0,13,14,15,16,17,18,19,20/)
 ;res@cnLevels               = (/-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6/)
 ;res@cnFillColors           = (/3,5,6,8,9,10,11,12,13,14,15,17,18,20/)
  res@cnLinesOn              = False             ; turn off contour lines
  res@cnLineColor            = "white"              ; turn off contour lines
  res@cnLineThicknessF       = 3.5              ; turn off contour lines  
  res@cnFillDrawOrder        = "PreDraw"          ; Make sure map fill happens
  res@cnLineDrawOrder        = "PostDraw"          ; Make sure map fill happens
  res@cnSmoothingOn          = False
  res@cnSmoothingTensionF    = -1
  res@cnInfoLabelOn          = False 
  res@cnInfoLabelOn          = False
  res@cnLabelMasking         = True
  res@cnLineLabelsOn            = True
  res@cnLineLabelPlacementMode  = "constant"
  res@cnLineLabelFontThicknessF = 3.5
  res@cnLineLabelFormat         = "@^sg" 
  res@cnLineLabelFontHeightF    = 0.01
  res@cnLineLabelFontColor      = "black"  
  res@cnLineLabelFontAspectF    = 1.4
  res@cnLineLabelPerimColor     = "white"
  res@cnLineLabelBackgroundColor= "white"
 ;res@cnLineLabelInterval       = 1  
  res@cnSmoothingOn             = False
 ;res@cnSmoothingDistanceF      = -0.01  

  res@mpFillOn               = False
  res@mpLandFillColor        = "white"
  res@mpCenterLonF           = 180
  res@mpMaxLatF              = 90             ; zoom in over australia
  res@mpMinLatF              = -90
  res@mpMaxLonF              = 360
  res@mpMinLonF              = 0
  res@mpOutlineOn            = True
  res@mpOutlineDrawOrder     = "Draw"
  res@mpProjection           = "WinkelTripel"
  res@mpProjection           = "Robinson"  
  res@mpGeophysicalLineColor      = "black"
  res@mpGeophysicalLineThicknessF = 0.5 
  res@mpGridAndLimbOn        = True              ; turn on lat/lon lines
  res@mpGridLineColor        = "Transparent"     ; line color
  res@mpGridLineThicknessF   = 1.5
  res@mpGridLineDashPattern  = 2
  res@mpGridAndLimbDrawOrder = "PostDraw"     ; draw lines first, so filled land is on top
  res@mpLimbLineThicknessF   = 1.0
  res@mpLimbLineColor        = "black"
  res@mpPerimOn              = False       ; turn off map perimeter
 ;res@mpPerimLineThicknessF  = 1.0
 ;res@mpPerimDrawOrder       = "PostDraw"

  res@tmXBMajorLengthF         = 0.008
  res@tmYLMajorLengthF         = 0.008
  res@tmXBMinorLengthF         = 0.005
  res@tmYLMinorLengthF         = 0.005 
  res@tmXBMajorThicknessF      = 3.5
  res@tmYLMajorThicknessF      = 3.5
  res@tmXBMinorThicknessF      = 3.0
  res@tmYLMinorThicknessF      = 3.0 
  res@tmXBLabelDeltaF          = -0.3
  res@tmYLLabelDeltaF          = -0.3
  res@tmXBLabelFontHeightF     = 0.015
  res@tmYLLabelFontHeightF     = 0.015
  res@tmXBLabelFontThicknessF  = 3.5
  res@tmYLLabelFontThicknessF  = 3.5
  res@tmBorderThicknessF       = 4.0

  res@tiYAxisString            = ""
  res@gsnLeftString            = ""
  res@gsnRightString           = "" 
  res@tiYAxisFontHeightF       = 0.015
  res@tiXAxisFontHeightF       = 0.015
  res@tiYAxisFontThicknessF    = 2
  res@tiXAxisFontThicknessF    = 2
  res@tiXAxisString            = ""
  res@tiYAxisString            = ""
  res@tiMainString             = "(a) Mean SAT Prob"
  res@tiMainFont               = 22 
  res@tiMainFontHeightF        = 0.008
 ;res@tiMainOffsetYF           = -0.005

  res@lbLabelBarOn          = True
 ;res@lbOrientation         = "Vertical"
  res@lbBoxLinesOn          = False 
  res@lbBoxEndCapStyle      = "TriangleBothEnds" 
  res@lbLabelFontHeightF    = 0.007  
  res@lbLabelFontThicknessF = 2.0
  res@pmLabelBarWidthF      = 0.20               ; default is shorter
  res@pmLabelBarHeightF     = 0.04 
  res@pmLabelBarOrthogonalPosF = 0.03
  res@pmLabelBarParallelPosF   = 0.5
  plot0                        = gsn_csm_contour_map(wks, rsp_all, res)

  res@vpXF                     = 0.26
  res@lbLabelBarOn          =  False
  delete([/res@cnLevels,res@cnFillColors/])
 ;res@cnLevels                 = (/-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50/)/5.0*1
 ;res@cnFillColors             = (/2,3,4,5,6,7,8,9,10,11,0,12,13,14,15,16,17,18,19,20,21/)
 ;res@cnLevels                 = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.5
 ;res@cnFillColors             = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels                 = fspan(0.1,3.0,30)
  res@cnFillColors           = ispan(66,126,2)
  res@tiMainString             = "(b) Regular El Ni"+"n~H-13V2F35~D~FV-2H3~"+"o SAT "+" ~F8~D~F22~Prob"
  plot1                        = gsn_csm_contour_map(wks, rsp_ren_all, res)

  res@vpXF                     = 0.48
    res@lbLabelBarOn          =  True 
   res@pmLabelBarWidthF      = 0.30               ; default is shorter
  res@pmLabelBarHeightF     = 0.04 
  res@pmLabelBarOrthogonalPosF = 0.03
  res@pmLabelBarParallelPosF   =  0.5
  res@cnLevels                 = fspan(0.1,3.0,30)
  res@cnFillColors           = ispan(66,126,2)
  res@tiMainString             = "(c) La Ni"+"n~H-13V2F35~D~FV-2H3~"+"a SAT "+" ~F8~D~F22~Prob"
  plot2                        = gsn_csm_contour_map(wks, rsp_rla_all, res) 

  res@vpXF                     = 0.70
  res@lbLabelBarOn          = False
  res@cnLevels                 = fspan(0.1,3.0,30)
  res@cnFillColors           = ispan(66,126,2)
  res@tiMainString             = "(d) Super El Ni"+"n~H-13V2F35~D~FV-2H3~"+"o SAT "+" ~F8~D~F22~Prob"
  plot3                        = gsn_csm_contour_map(wks, rsp_sen_all, res) 

  gres                   = True
  gres@gsLineColor       = "black"
  gres@gsLineDashPattern = 2
  gres@gsLineThicknessF  = 0.5
  gres@tfPolyDrawOrder   = "Draw"  ; this can be used for polylines, polymarkers, or polygons

  lat4                 = new((/360/),"float")
  lat4(:)              = 0
  lon4                 = rsp_all&lon
  dum0                 = gsn_add_polyline(wks,plot0,lon4,lat4,gres) 
  dum1                 = gsn_add_polyline(wks,plot1,lon4,lat4,gres)
  dum2                 = gsn_add_polyline(wks,plot2,lon4,lat4,gres)
  dum3                 = gsn_add_polyline(wks,plot3,lon4,lat4,gres)

  draw(plot0)
  draw(plot1)
  draw(plot2)
  draw(plot3)


  txres               = True         ; Text resources
  txres@txFontHeightF = 0.012
  txres@txFontColor  = "black"
  txres@txFontThicknessF = 4.0

 ; text1 = (/"72/73", "82/83", "97/98", "15/16"/)
 text1 = (/ "~F22~CESM2-LE Historical", "~F22~CESM2-LE SSP370"/)

    gsn_text_ndc(wks, text1(0), 0.47, 0.98, txres)
    ;gsn_text_ndc(wks, text1(1), 0.75, 0.98, txres)
    delete([/fili,nino,ren, rla, sen,nino1D/])
  delete([/rsp_all, rsp_sen, rsp_sen_all, rsp_ren_all, rsp_rla_all/])
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SST_HAD
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  fili     = systemfunc("cd "+ncdir2+" ; ls rsi_seasonally_sst_cesm2_hist_ensemble*.nc")
  nfili    = dimsizes(fili)
  print("nfili ="+nfili)

  do i         = 0,nfili-1
  ; print("i   = "+i)
  ; print("file= "+fili(i))
    f          = addfile(ncdir2+fili(i),"r")
    rsi_son    = f->rsi_son
    rsi_djf    = f->rsi_djf
    rsi_mam    = f->rsi_mam    
    rsi_jja    = f->rsi_jja
    nlat       = dimsizes(f->lat)
    nlon       = dimsizes(f->lon)
    ntime      = dimsizes(f->time)
    rsi_son    = where(rsi_son.ne.0, 1.0, rsi_son)
    rsi_djf    = where(rsi_djf.ne.0, 1.0, rsi_djf)
    rsi_mam    = where(rsi_mam.ne.0, 1.0, rsi_mam)
    rsi_jja    = where(rsi_jja.ne.0, 1.0, rsi_jja)
    rsi_sum    = rsi_son 
    rsi_sum    = ((/rsi_son/)+(/rsi_djf/)+(/rsi_mam/)+(/rsi_jja/))    

    if(i.eq.0) then 
      rsi_tot      = new((/nfili,ntime,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot!0    = "ensemble"
      rsi_tot!1    = "time"
      rsi_tot!2    = "lat"
      rsi_tot!3    = "lon"
      rsi_tot&time = rsi_sum&time      
      rsi_tot&lat  = rsi_sum&lat    
      rsi_tot&lon  = rsi_sum&lon
      rsi_tot_sen      = new((/nfili,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot_sen!0    = "ensemble"
      rsi_tot_sen!1    = "lat"
      rsi_tot_sen!2    = "lon"     
      rsi_tot_sen&lat  = rsi_sum&lat    
      rsi_tot_sen&lon  = rsi_sum&lon
      rsi_tot_all      = rsi_tot_sen

      rsi_tot_ren      = new((/nfili,nlat,nlon/),"float",rsi_sum@_FillValue) ; ren regular el nino
      rsi_tot_ren!0    = "ensemble"
      rsi_tot_ren!1    = "lat"
      rsi_tot_ren!2    = "lon"     
      rsi_tot_ren&lat  = rsi_sum&lat    
      rsi_tot_ren&lon  = rsi_sum&lon

      rsi_tot_rla      = new((/nfili,nlat,nlon/),"float",rsi_sum@_FillValue); rla regular la Nina
      rsi_tot_rla!0    = "ensemble"
      rsi_tot_rla!1    = "lat"
      rsi_tot_rla!2    = "lon"     
      rsi_tot_rla&lat  = rsi_sum&lat    
      rsi_tot_rla&lon  = rsi_sum&lon

    end if
    rsi_tot(i,:,:,:) = (/rsi_sum(:,:,:)/)
    masks            = ntime-(L/2)
    maske            = ntime-1
    rsi_tot(i,masks:maske,:,:) = 0

    f1          = addfile(ncdir4+"nino_cesm2_hist.nc","r")
    nino        = f1->nino(i,:)
    nino(0:L-1)         = 0  
    nino(masks-1:maske) = 0
    sen1        = ind(nino.ge.2.0)
    sen2        = sen1+1
    sen         = array_append_record(sen1, sen2, 0)
    nsen     = dimsizes(sen)*2
    ren      = ind(nino.ge.0.5 .and. nino.le.1.5)
    rla      = ind(nino.le.-0.5 .and. nino.ge.-1.0)
     nren     = dimsizes(ren)
     nrla     = dimsizes(rla)
    rsi_tot_sen(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,sen,:,:),0)
    rsi_tot_ren(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,ren,:,:),0)
    rsi_tot_rla(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,rla,:,:),0)
    rsi_tot_all(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,:,:,:),0)
    delete([/sen1,sen2,sen, ren, rla/])
  end do
  printVarSummary(rsi_tot_sen)
  delete([/nino/])

  f2       = addfile(ncdir4+"nino_cesm2_hist.nc","r")
  nino     = f2->nino(0:nfili-1,L:masks-2)
  nino1D   = ndtooned(nino)
  sen      = ind(nino1D.ge.2.0)
  ren      = ind(nino1D.ge.0.5 .and. nino1D.le.1.5)
  rla       = ind(nino1D.le.-0.5 .and. nino1D.ge.-1.0)
  nsen     = dimsizes(sen)*2
  nren     = dimsizes(ren)
  nrla     = dimsizes(rla)
  nall     = nfili*(ntime-(L+(L/2)))

  rsp_all   = dim_sum_n_Wrap(rsi_tot_all,0)
  rsp_all   = rsp_all/(4.0*nall)*100.0
  rsp_sen   = dim_sum_n_Wrap(rsi_tot_sen,0)
  rsp_sen   = rsp_sen/(4.0*nsen)*100.0
  rsp_ren   = dim_sum_n_Wrap(rsi_tot_ren,0)
  rsp_ren   = rsp_ren/(4.0*nren)*100.0
  rsp_rla   = dim_sum_n_Wrap(rsi_tot_rla,0)
  rsp_rla   = rsp_rla/(4.0*nrla)*100.0

  rsp_sen_all = (rsp_sen-rsp_all)
  rsp_ren_all = (rsp_ren-rsp_all)
  rsp_rla_all = (rsp_rla-rsp_all)
  copy_VarCoords(rsp_sen, rsp_sen_all)
  copy_VarCoords(rsp_ren, rsp_ren_all)
  copy_VarCoords(rsp_rla, rsp_rla_all)
printMinMax(rsp_sen_all, 0)
  wgty1                   = cos(get_d2r(1.0)*f->lat)
  gm_ratio_all            = wgt_areaave_Wrap(rsp_all(:,:), wgty1, 1.0, 0)
  gm_ratio_sen            = wgt_areaave_Wrap(rsp_sen(:,:), wgty1, 1.0, 0)
  gm_ratio_sen_all        = wgt_areaave_Wrap(rsp_sen_all(:,:), wgty1, 1.0, 0)
  print("gm_ratio_all     = "+gm_ratio_all)
  print("gm_ratio_sen     = "+gm_ratio_sen)
  print("gm_ratio_sen_all = "+gm_ratio_sen_all)
  print("*************************************")
  delete([/rsi_tot_ren, rsi_tot_rla, rsi_tot_sen, rsi_tot_all/])

;************************************************
; plot 

  delete([/res@cnLevels,res@cnFillColors/])
  res@vpXF                   = 0.04                      ; Move slightly to left
  res@vpYF                   = 0.78
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.25
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels                 = fspan(1.1,4.0,30)
  res@cnFillColors           = ispan(66,126,2)
  res@tiMainString           = "(e) Mean SST Prob"
 ;res@mpLandFillColor        = 8
  res@lbLabelBarOn          = True
     res@pmLabelBarWidthF      = 0.20               ; default is shorter
  plot4                      = gsn_csm_contour_map(wks, rsp_all, res)

  res@vpXF                   = 0.26
    res@lbLabelBarOn          = False 
  delete([/res@cnLevels,res@cnFillColors/])
 ;res@cnLevels               = (/-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50/)/5.0*1
 ;res@cnFillColors           = (/2,3,4,5,6,7,8,9,10,11,0,12,13,14,15,16,17,18,19,20,21/)
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.5
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels               = fspan(0.1,3.0,30)
  res@cnFillColors           = ispan(66,126,2)
  res@tiMainString           = "(f) Regular El Ni"+"n~H-13V2F35~D~FV-2H3~"+"o SST "+" ~F8~D~F22~Prob"
  plot5                      = gsn_csm_contour_map(wks, rsp_ren_all, res)
  res@vpXF                     = 0.48
  res@lbLabelBarOn          =  True 
   res@pmLabelBarWidthF      = 0.30               ; default is shorter
  res@pmLabelBarHeightF     = 0.04 
  res@pmLabelBarOrthogonalPosF = 0.03
  res@pmLabelBarParallelPosF   =  0.5    
  delete([/res@cnLevels,res@cnFillColors/])
 ;res@cnLevels               = (/-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50/)/5.0*1
 ;res@cnFillColors           = (/2,3,4,5,6,7,8,9,10,11,0,12,13,14,15,16,17,18,19,20,21/)
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.5
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels               = fspan(0.1,3.0,30)
  res@cnFillColors           = ispan(66,126,2)
  res@tiMainString           = "(g) La Ni"+"n~H-13V2F35~D~FV-2H3~"+"a SST "+" ~F8~D~F22~Prob"
  plot6                      = gsn_csm_contour_map(wks, rsp_rla_all, res)

  res@vpXF                   = 0.70
  res@lbLabelBarOn          = False
  delete([/res@cnLevels,res@cnFillColors/])
 ;res@cnLevels               = (/-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50/)/5.0*1
 ;res@cnFillColors           = (/2,3,4,5,6,7,8,9,10,11,0,12,13,14,15,16,17,18,19,20,21/)
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.5
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels               = fspan(0.1,3.0,30)
  res@cnFillColors           = ispan(66,126,2)
  res@tiMainString           = "(h) Super El Ni"+"n~H-13V2F35~D~FV-2H3~"+"o SST "+" ~F8~D~F22~Prob"
  plot7                      = gsn_csm_contour_map(wks, rsp_sen_all, res)

  lat4                 = new((/360/),"float")
  lat4(:)              = 0
  lon4                 = rsp_all&lon
  dum4                 = gsn_add_polyline(wks,plot4,lon4,lat4,gres) 
  dum5                 = gsn_add_polyline(wks,plot5,lon4,lat4,gres)
  dum6                 = gsn_add_polyline(wks,plot6,lon4,lat4,gres)
  dum7                 = gsn_add_polyline(wks,plot7,lon4,lat4,gres)

  draw(plot4)
  draw(plot5)
  draw(plot6)
  draw(plot7)
    delete([/fili,nino,sen,ren, rla,rsi_tot,nino1D/])
    delete([/rsp_all, rsp_ren, rsp_rla, rsp_sen, rsp_sen_all, rsp_ren_all, rsp_rla_all/])
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; soil
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  fili     = systemfunc("cd "+ncdir3+" ; ls rsi_seasonally_soilw_cesm2_hist_ensemble*.nc")
  nfili    = dimsizes(fili)
  print("nfili ="+nfili)

  do i         = 0,nfili-1
  ;  print("i   = "+i)
  ;  print("file= "+fili(i))
    f          = addfile(ncdir3+fili(i),"r")
    rsi_son    = f->rsi_son
    rsi_djf    = f->rsi_djf
    rsi_mam    = f->rsi_mam    
    rsi_jja    = f->rsi_jja
    nlat       = dimsizes(f->lat)
    nlon       = dimsizes(f->lon)
    ntime      = dimsizes(f->time)
    rsi_son    = where(rsi_son.ne.0, 1.0, rsi_son)
    rsi_djf    = where(rsi_djf.ne.0, 1.0, rsi_djf)
    rsi_mam    = where(rsi_mam.ne.0, 1.0, rsi_mam)
    rsi_jja    = where(rsi_jja.ne.0, 1.0, rsi_jja)
    rsi_sum    = rsi_son 
    rsi_sum    = ((/rsi_son/)+(/rsi_djf/)+(/rsi_mam/)+(/rsi_jja/))    

    if(i.eq.0) then 
      rsi_tot      = new((/nfili,ntime,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot!0    = "ensemble"
      rsi_tot!1    = "time"
      rsi_tot!2    = "lat"
      rsi_tot!3    = "lon"
      rsi_tot&time = rsi_sum&time      
      rsi_tot&lat  = rsi_sum&lat    
      rsi_tot&lon  = rsi_sum&lon
      rsi_tot_sen      = new((/nfili,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot_sen!0    = "ensemble"
      rsi_tot_sen!1    = "lat"
      rsi_tot_sen!2    = "lon"     
      rsi_tot_sen&lat  = rsi_sum&lat    
      rsi_tot_sen&lon  = rsi_sum&lon
      rsi_tot_ren      = new((/nfili,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot_ren!0    = "ensemble"
      rsi_tot_ren!1    = "lat"
      rsi_tot_ren!2    = "lon"     
      rsi_tot_ren&lat  = rsi_sum&lat    
      rsi_tot_ren&lon  = rsi_sum&lon
      rsi_tot_rla      = new((/nfili,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot_rla!0    = "ensemble"
      rsi_tot_rla!1    = "lat"
      rsi_tot_rla!2    = "lon"     
      rsi_tot_rla&lat  = rsi_sum&lat    
      rsi_tot_rla&lon  = rsi_sum&lon

      rsi_tot_all      = rsi_tot_sen
    end if
    rsi_tot(i,:,:,:) = (/rsi_sum(:,:,:)/)
    masks            = ntime-(L/2)
    maske            = ntime-1
    rsi_tot(i,masks:maske,:,:) = 0

   ;f1          = addfile(ncdir4+"nino_cesm2_hist.nc","r")
   ;nino        = f1->nino(i,:)


    f1          = addfile(ncdir4+"nino_djf_cesm2_for_soilw_tsoilw_hist.nc","r")
    nino        = f1->nino(i,:)
    nino(0:L-1)         = 0  
    nino(masks-1:maske) = 0
    sen1        = ind(nino.ge.2.0)
    sen2        = sen1+1
    sen         = array_append_record(sen1, sen2, 0)
    ren      = ind(nino.ge.0.5 .and. nino.le.1.5)
    rla      = ind(nino.le.-0.5 .and. nino.ge.-1.0)
    rsi_tot_sen(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,sen,:,:),0)
    rsi_tot_all(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,:,:,:),0)
    rsi_tot_ren(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,ren,:,:),0)
    rsi_tot_rla(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,rla,:,:),0)
    delete([/sen1,sen2,sen, ren, rla/])
  end do
  printVarSummary(rsi_tot_sen)
  delete([/nino/])

  f2       = addfile(ncdir4+"nino_djf_cesm2_for_soilw_tsoilw_hist.nc","r")
  nino     = f2->nino(0:nfili-1,L:masks-2)
  nino1D   = ndtooned(nino)
  sen      = ind(nino1D.ge.2.0)
  ren      = ind(nino1D.ge.0.5 .and. nino1D.le.1.5)
  rla       = ind(nino1D.le.-0.5 .and. nino1D.ge.-1.0)
  nsen     = dimsizes(sen)*2
  nren     = dimsizes(ren)
  nrla     = dimsizes(rla)
  nall     = nfili*(ntime-(L+(L/2)))


  rsp_all   = dim_sum_n_Wrap(rsi_tot_all,0)
  rsp_all   = rsp_all/(4.0*nall)*100.0
  rsp_sen   = dim_sum_n_Wrap(rsi_tot_sen,0)
  rsp_sen   = rsp_sen/(4.0*nsen)*100.0
  rsp_ren   = dim_sum_n_Wrap(rsi_tot_ren,0)
  rsp_ren   = rsp_ren/(4.0*nren)*100.0
  rsp_rla   = dim_sum_n_Wrap(rsi_tot_rla,0)
  rsp_rla   = rsp_rla/(4.0*nrla)*100.0

  rsp_sen_all = (rsp_sen-rsp_all)
  rsp_ren_all = (rsp_ren-rsp_all)
  rsp_rla_all = (rsp_rla-rsp_all)
  copy_VarCoords(rsp_sen, rsp_sen_all)
  copy_VarCoords(rsp_ren, rsp_ren_all)
  copy_VarCoords(rsp_rla, rsp_rla_all)
printMinMax(rsp_sen_all, 0)
  wgty1                   = cos(get_d2r(1.0)*f->lat)
  gm_ratio_all            = wgt_areaave_Wrap(rsp_all(:,:), wgty1, 1.0, 0)
  gm_ratio_sen            = wgt_areaave_Wrap(rsp_sen(:,:), wgty1, 1.0, 0)
  gm_ratio_sen_all        = wgt_areaave_Wrap(rsp_sen_all(:,:), wgty1, 1.0, 0)
  print("gm_ratio_all     = "+gm_ratio_all)
  print("gm_ratio_sen     = "+gm_ratio_sen)
  print("gm_ratio_sen_all = "+gm_ratio_sen_all)

  print("*************************************")
  delete([/rsi_tot_ren, rsi_tot_rla, rsi_tot_sen, rsi_tot_all/])
;************************************************
; plot 

  delete([/res@cnLevels,res@cnFillColors/])
  res@vpXF                   = 0.04                     ; Move slightly to left
  res@vpYF                   = 0.61
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.25
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels                 = fspan(1.1,4.0,30)
  res@cnFillColors           = ispan(66,126,2)
  res@tiMainString           = "(i) Mean SOIL Prob"
 ;res@mpFillOn               = False
  res@lbLabelBarOn           = True
     res@pmLabelBarWidthF      = 0.20               ; default is shorter
  plot8                      = gsn_csm_contour_map(wks, rsp_all, res)

  res@vpXF                   = 0.26
  delete([/res@cnLevels,res@cnFillColors/])
    res@lbLabelBarOn          = False 
 ;res@cnLevels               = (/-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50/)/5.0*1
 ;res@cnFillColors           = (/2,3,4,5,6,7,8,9,10,11,0,12,13,14,15,16,17,18,19,20,21/)
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.5
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels               = fspan(0.1,2.0,20)
  res@cnFillColors           = ispan(66,126,3)
  res@tiMainString           = "(j) Regular El Ni"+"n~H-13V2F35~D~FV-2H3~"+"o SOIL "+" ~F8~D~F22~Prob"
  plot9                      = gsn_csm_contour_map(wks, rsp_ren_all, res)

 res@vpXF                   = 0.48
  delete([/res@cnLevels,res@cnFillColors/])
   res@lbLabelBarOn          =  True 
   res@pmLabelBarWidthF      = 0.30               ; default is shorter
  res@pmLabelBarHeightF     = 0.04 
  res@pmLabelBarOrthogonalPosF = 0.03
  res@pmLabelBarParallelPosF   =  0.5
 ;res@cnLevels               = (/-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50/)/5.0*1
 ;res@cnFillColors           = (/2,3,4,5,6,7,8,9,10,11,0,12,13,14,15,16,17,18,19,20,21/)
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.5
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels               = fspan(0.1,2.0,20)
  res@cnFillColors           = ispan(66,126,3)
  res@tiMainString           = "(k) La Ni"+"n~H-13V2F35~D~FV-2H3~"+"a SOIL "+" ~F8~D~F22~Prob"
  plot10                      = gsn_csm_contour_map(wks, rsp_rla_all, res)


   res@vpXF                   = 0.70
   res@lbLabelBarOn          = False
  delete([/res@cnLevels,res@cnFillColors/])
 ;res@cnLevels               = (/-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50/)/5.0*1
 ;res@cnFillColors           = (/2,3,4,5,6,7,8,9,10,11,0,12,13,14,15,16,17,18,19,20,21/)
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.5
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels               = fspan(0.1,2.0,20)
  res@cnFillColors           = ispan(66,126,3)
  res@tiMainString           = "(l) Super El Ni"+"n~H-13V2F35~D~FV-2H3~"+"o SOIL "+" ~F8~D~F22~Prob"
  plot11                      = gsn_csm_contour_map(wks, rsp_sen_all, res)


  lat4                 = new((/360/),"float")
  lat4(:)              = 0
  lon4                 = rsp_all&lon
  dum8                 = gsn_add_polyline(wks,plot8,lon4,lat4,gres) 
  dum9                 = gsn_add_polyline(wks,plot9,lon4,lat4,gres)
  dum10                = gsn_add_polyline(wks,plot10,lon4,lat4,gres)
  dum11                = gsn_add_polyline(wks,plot11,lon4,lat4,gres)

  draw(plot8)
  draw(plot9)
  draw(plot10)
  draw(plot11)

  delete([/rsi_son, rsi_djf, rsi_mam, rsi_jja, rsi_sum, rsi_tot/])
  delete([/rsp_all, rsp_sen, rsp_sen_all, rsp_ren_all, rsp_rla_all/])
  delete([/fili,nino,sen,ren, rla,nino1D/])
;============================================================================================
end