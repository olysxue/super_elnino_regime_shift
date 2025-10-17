 %   stars_Precipitation.m
 %
 %   A demonstration of the STARS change point (regime shift) algorithm.
 %   
 %
 %   written by olysxue
 %   <olysxue@gmail.com>
 %   2022-pohang
 %
  % -------------------------------------------------------------------------
  %   Initialization
  % -------------------------------------------------------------------------
  rng(0)                          % seed the RNG for reproducibility
  
  L = 10;                         % cutoff        |
  p = 0.05;                       % significance  |   -> stars parameters
  h = 1.0;                        % huber weight  |
  addpath('/glade/derecho/scratch/aoyunxue/regime_shift/function/');
  tic;
  mydir = ['/glade/derecho/scratch/aoyunxue/regime_shift/obs/data/soil_moisture/ERA5/ERA5_land/'];
  myfile = 'ERA5_land_anom_soilmoisture1.nc';
  outfile = 'rsi_ERA5_land_anom_soilmoisture1.nc';
  
     lon     = ncread(strcat(mydir,myfile),'lon');
     lat     = ncread(strcat(mydir,myfile),'lat');
     time    = ncread(strcat(mydir,myfile),'valid_time');
     nlon    = length(lon);
     nlat    = length(lat);
     ntime   = length(time);
     sm     = zeros([nlon nlat ntime]);
     sm     = ncread(strcat(mydir,myfile),'sm_son');
    % -------------------------------------------------------------------------
    % Calculate change points (RSI series)
    % -------------------------------------------------------------------------
  rsi = zeros(nlon, nlat, 74, 'single');
   for m = 1:nlon;
       for n = 1:nlat;
             rsi(m,n,:) = stars(squeeze(sm(m,n,:)), 'L', L, 'p', p, 'h', h);
       end
   end
  varname = 'rsi';
  fillvalue = -999;
  rsi(isnan(rsi)) = fillvalue;
  outdir = ['/glade/derecho/scratch/aoyunxue/regime_shift/obs/data/soil_moisture/ERA5/ERA5_land/'];
  fname0  = strcat(outdir,outfile);
  fname   = char(fname0);
                                                                                                                                                                                                                                                                                       
     nccreate(fname,'lat', 'Dimensions',{'lat',length(lat)},'Datatype','double','Format','64bit');
     ncwrite(fname, 'lat', lat);
     nccreate(fname,'lon', 'Dimensions',{'lon',length(lon)},'Datatype','double','Format','64bit');
     ncwrite(fname, 'lon', lon);
     nccreate(fname,'time', 'Dimensions',{'time',length(time)},'Datatype','double','Format','64bit');
     ncwrite(fname, 'time', time);
  
     nccreate(fname, varname,'Dimensions',{'lon' nlon 'lat' nlat 'time' ntime},'Datatype','single','Format','64bit');
     ncwrite(fname, varname, rsi);
     ncwriteatt(fname, varname, '_FillValue', fillvalue);
     ncwriteatt(fname,'lat','units','degrees_north');
     ncwriteatt(fname,'lat','long_name','latitude');
     ncwriteatt(fname,'lon','units','degrees_east');
     ncwriteatt(fname,'lon','long_name','longitude');
  
     ncwriteatt(fname,'time','units','years since 1950-01-01 00:00:00');
     ncwriteatt(fname,'time','long_name','time');
   toc;
