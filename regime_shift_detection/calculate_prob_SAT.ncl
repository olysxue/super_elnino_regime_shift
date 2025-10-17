load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

begin

  ncdir1   = "/glade/derecho/scratch/aoyunxue/regime_shift/obs/data/soil/ERA5/"
  opdir1   = "/glade/derecho/scratch/aoyunxue/regime_shift/obs/data/soil/ERA5/"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SAT_era5
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  f1        = addfile(ncdir1+"rsp_total_soil_ERA5_25degree.nc", "r")
  lat      = (f1->lat)*-1.0
  lon      = f1->lon 
  time     = f1->time
  rsp_tot1  = f1->rsp_tot
  nlat     = dimsizes(lat)
  nlon     = dimsizes(lon)
  ntime    = dimsizes(time)
  wgty     = tofloat(cos(get_d2r(1.0)*lat))
  do   i   = 0,dimsizes(lat)-1
    rsp_tot1(:,i,:) = (/rsp_tot1(:,i,:)/)*wgty(i)
  end do
  printVarSummary(rsp_tot1)

  L        = 10
  masks    = ntime-(L/2)
  maske    = ntime-1
  rsp_tot1(masks:maske,:,:) = 0
  rsi_tot1  = dim_sum_n_Wrap(rsp_tot1,0)  ;;spatial distribution
  printVarSummary(rsp_tot1)

  year     = ispan(1948,2022,1)
  sen1     = ind((year.eq.1982).or.(year.eq.1997).or.(year.eq.2015))
  sen2     = ind((year.eq.1982).or.(year.eq.1997).or.(year.eq.2015))+1
  sen      = array_append_record(sen1, sen2, 0)
  nor      = ind((year.ne.1982.and.year.ne.1983).and.(year.ne.1997.and.year.ne.1998).and.(year.ne.2015.and.year.ne.2016))
  nsen     = dimsizes(sen)
  nnor     = dimsizes(nor)-(L+(L/2))
  printVarSummary(sen)

  rsp_all   = dim_sum_n_Wrap(rsp_tot1(:,:,:),0)
  rsp_all   = rsp_all/(4.0*(ntime-(L+(L/2))))*100.0
  rsp_sen   = dim_sum_n_Wrap(rsp_tot1(sen,:,:),0)
  rsp_sen   = rsp_sen/(4.0*nsen)*100.0
  rsp_nor   = dim_sum_n_Wrap(rsp_tot1(nor,:,:),0)
  rsp_nor   = rsp_nor/(4.0*nnor)*100.0
  ratio_sen_nor = rsp_sen 
  ratio_sen_nor = (rsp_sen-rsp_nor)
  ratio_sen_all = (rsp_sen-rsp_all)
  
  fname       = "prob_soil_ERA5.nc"
  system("/bin/rm -f "+opdir1+fname)
  a1          = addfile(opdir1+fname,"c")    ; write netCDF file   
  a1->prob_all = rsp_all
  a1->prob_sen = rsp_sen
  a1->prob_nor = rsp_nor
  a1->prob_increase = ratio_sen_nor
  a1->prob_increase2 = ratio_sen_all

end 