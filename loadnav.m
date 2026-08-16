%======================================================================
%                    L O A D N A V . M 
%                    doc: Thu Jun 17 18:01:50 2004
%                    dlm: Tue Jan 28 13:24:01 2020
%                    (c) 2004 ladcp@
%                    uE-Info: 43 64 NIL 0 0 72 2 2 8 NIL ofnI
%======================================================================

% MODIFICATIONS BY ANT:
%   Jun 26, 2004: - totally re-wrote the file-reading, which reduces
%		    run time by over factor 90(!!!) --- 5.7 instead of
%		    527 seconds in the test case --- for deep casts
%   Jul  1, 2004: - added support for time bases
%		  - made input definition flexible
%		  - removed ipos argument
%   Dec 18, 2004: - BUG: header_lines > 0 did not work
%   Jun  9, 2005: - improved time-base comments
%   Aug  9, 2006: - added support for elapsed-time
%   Jan  5, 2007: - removed LADDER-1 specific code for IX_4 distribution
%		  - added Dan Torres' code to interpolate irregular
%		    GPS time series
%		  - added nav file layout into p structure
%   Jan 17, 2007: - added support for new [geomag.m]
%   Jan 26, 2007: - BUG: file layout default was in p (rather than f) structure
%   Jun 30, 2008: - adapted to calculate magdev using external program geomag60
%   Jul  2, 2008: - BUG: missing "" in geomag command to allow spcs in 
%			 cof-file path
%   Jul 17, 2008: - moved some code to [loadctd.m] to fix bug associated with
%		    adjusting of start/end positions in case of significant
%		    ADCP vs GPS/CTD clock offset (unclear whether this happened
%		    only when elapsed time was used)
%   Jul 27, 2008: - nanmean() -> meannan()
%   Oct 15, 2008: - replaced mean by median to get lat/lon (bad outliers in L1 data set)
%   Dec  1, 2009: - BUG: geomag date check was wrong (Dec 1 2009 resulted in a date >= 2010)
%   Jan 22, 2010: - adapted to Eric Firing's much simplified magdec utility
%   Jan  3, 2011: - changed IGRF11 validity to end of 2015 (from 2010)
%   Feb 18, 2016: - BUG: geomag year range check bombed in 2016
%		  - added p.interp_missing_GPS using code provided by Jay Hooper
%   Sep 14, 2018: - BUG: 2008 code move broke working with single GPS file for entire cruise
%   Sep  4, 2019: - adapted to GK new magdev.m
%		  - added p.magdec_source
%   Jan 28, 2020: - 2008 code move is required because at the loadnav stage p.time_start and
%		    end reflect the LADCP turn-on and -off times

function [d,p]=loadnav(f,d,p)
% function [d,p]=loadnav(f,d,p)

% This routine works for generic ASCII files, containing GPS time series,
% with fields time, lat and lon.

%====================
% TWEAKABLES
%====================

% PRE-AVERAGING OF GPS TIME (IN DAYS)
p = setdefv(p,'navtime_av',2/60/24);

% INTERPOLATE IRREGULAR NAV TIME SERIES (Dan Torres)
p=setdefv(p,'interp_nav_times',0);

% INTERPOLATE MISSING GPS VALUES (Jay Hooper)
p=setdefv(p,'interp_missing_GPS',1);

% MAGNETIC DECLINATION
%    There are three different sources for magnetic declination
%	- specify p.drot manually
%	- p.magdec_source = 1	% use external magdec program, if available, otherwise ...
%	- p.magdec_source = 2	% use [magdev.m] from Gerd Krahmann
p=setdefv(p,'magdec_source',2);

% FILE LAYOUT
f = setdefv(f,'nav_header_lines',0);
f = setdefv(f,'nav_fields_per_line',3);
f = setdefv(f,'nav_time_field',1);
f = setdefv(f,'nav_lat_field',2);
f = setdefv(f,'nav_lon_field',3);

% TIME BASE
% 	0 for elapsed time in seconds
% 	1 for year-day (1.0 = Jan 1, 00:00)
% 	2 for Visbeck's Gregorian (see gregoria.m)
f = setdefv(f,'nav_time_base',0);

%======================================================================

% MODIFICATIONS BY ANT:
%   Jun 26, 2004: - totally re-wrote the file-reading, which reduces
%		    run time by over factor 90(!!!) --- 5.7 instead of
%		    527 seconds in the test case --- for deep casts
%   Jul  1, 2004: - added support for time bases
%		  - made input definition flexible
%		  - removed ipos argument
%   Dec 18, 2004: - BUG: header_lines > 0 did not work
%   Jun  9, 2005: - improved time-base comments
%   Aug  9, 2006: - added support for elapsed-time
%   Jan  5, 2007: - removed LADDER-1 specific code for IX_4 distribution
%		  - added Dan Torres' code to interpolate irregular
%		    GPS time series
%		  - added nav file layout into p structure
%   Jan 17, 2007: - added support for new [geomag.m]
%   Jan 26, 2007: - BUG: file layout default was in p (rather than f) structure
%   Jan  7, 2009: - tightened use of exist()

disp(['LOADNAV: load NAV time series ',f.nav])
if ~exist(f.nav,'file')
 disp([' can not find ',f.nav])
 return
end

% construct input format
cur_field = 1; input_format = '';
for i=1:f.nav_fields_per_line
  switch i,
    case f.nav_time_field,
      i_time = cur_field; cur_field = cur_field + 1;
      input_format = [input_format ' %g'];
    case f.nav_lat_field
      i_lat = cur_field; cur_field = cur_field + 1;
      input_format = [input_format ' %g'];
    case f.nav_lon_field
      i_lon = cur_field; cur_field = cur_field + 1;
      input_format = [input_format ' %g'];
    otherwise
      input_format = [input_format ' %*g'];
  end
end
if cur_field ~= 4
  error('File format definition error');
end

% open input & skip header
header_lines = f.nav_header_lines;
fp=fopen(f.nav);
while header_lines > 0
  fgets(fp);
  header_lines = header_lines - 1;
end

% read time series
[A,nread] = fscanf(fp,input_format,[3,inf]);

% close file
fclose(fp);

% NAV time
d.navtime_jul=A(i_time,:)';
switch f.nav_time_base
  case 0 % elapsed time in seconds
    d.navtime_jul = d.navtime_jul/24/3600 + julian(p.time_start);
  case 1 % year-day
    d.navtime_jul = d.navtime_jul + julian([p.time_start(1) 1 0 0 0 0]);
end

disp([' number of NAV scans: ',int2str(length(d.navtime_jul)),...
       '  delta t : ',num2str(median(diff(d.navtime_jul))*24*3600),' seconds'])

%----------------------------------------
% interpolate to regular time series
%	code provided by Dan Torres
%----------------------------------------

if p.interp_nav_times
  min_t = min(d.navtime_jul);
  max_t = max(d.navtime_jul);
  delta_t = median(diff(d.navtime_jul));
  data = interp1q(d.navtime_jul,A([i_lat i_lon],:)',[min_t:delta_t:max_t]');
  d.navtime_jul = [min_t:delta_t:max_t]';
  disp(sprintf(' interpolated to %d NAV scans; delta_t = %.2f seconds',...
		length(d.navtime_jul),median(diff(d.navtime_jul))*24*3600));
else
  data=A([i_lat i_lon],:)';
end

p.navdata = 1;
d.slat = data(:,1);
d.slon = data(:,2);

%----------------------------------------
% interpolate missign GPS values
%	code provided by Jay Hooper
%----------------------------------------

if p.interp_missing_GPS
    bad_lon = find(d.slon == -9.990e-29);
    if ~isempty(bad_lon),
	good_lon = find(d.slon ~= -9.990e-29);
	d.slon = interp1(d.navtime_jul(good_lon),d.slon(good_lon),d.navtime_jul);
    end
    bad_lat = find(d.slat == -9.990e-29);
    if ~isempty(bad_lat),
	good_lat = find(d.slat ~= -9.990e-29);
	d.slat = interp1(d.navtime_jul(good_lat),d.slat(good_lat),d.navtime_jul);
    end
end

%----------------------------------------------------------------------
% The following code was moved into [loadctd.m] in July 2008 to fix a
% bug that most likely only occurs with elapsed times. In 2018, it was
% discovered that the bug fix introduced another bug that only affects
% data sets with a single GPS time series from which the per-station
% information must be extracted. Therefore, the code was moved back
% here but only called when not working with elapsed time.
% In Jan 2020 I realized that the code cannot be here at all, because 
% p.time_start and end are still the LADCP turn-on and -off times.
% Therefore, I disabled the code here and re-enabled it for all
% time bases in [loadctd.m]. I do not remember what bug this now causes.
%----------------------------------------------------------------------

if 0
  if min(d.navtime_jul)>max(d.time_jul) | max(d.navtime_jul)<min(d.time_jul)
    disp('NAV timeseries does not overlap WRONG STATION????')
    disp(['  NAV-data : ',datestrj(d.navtime_jul(1)),'	to ',datestrj(d.navtime_jul(end))])
    disp([' ADCP-data : ',datestrj(d.time_jul(1)),'  to ',datestrj(d.time_jul(end))])
    disp(' will ignore nav data')
    p.navdata = 0;
    d.navtime_jul = d.time_jul; % make sure same length
    d.slat = d.navtime_jul*NaN;
    d.slon = d.slat;
  else
    if d.navtime_jul(1) > d.time_jul(1)
      disp('NAV timeseries starts after ADCP timeseries: used first NAV value to patch ')
      disp(['  NAV-data : ',datestrj(d.navtime_jul(1)),'  to ',datestrj(d.navtime_jul(end))])
      disp([' ADCP-data : ',datestrj(d.time_jul(1)),'  to ',datestrj(d.time_jul(end))])
      d.navtime_jul(1) = d.time_jul(1);
    end
    
    if d.navtime_jul(end) < d.time_jul(end)
      disp('NAV timeseries ends before ADCP timeseries: used last NAV value to patch ')
      disp(['  NAV-data : ',datestrj(d.navtime_jul(1)),'  to ',datestrj(d.navtime_jul(end))])
      disp([' ADCP-data : ',datestrj(d.time_jul(1)),'  to ',datestrj(d.time_jul(end))])
      d.navtime_jul(end) = d.time_jul(end);
    end
  
    % find valid
    ii=find(diff(d.navtime_jul)>0);
    ii=[ii;length(d.navtime_jul)];
  
    % average over p.navtime_av days
    dt2m=[-p.navtime_av:(1/3600/24):p.navtime_av]';
    slon=medianan(interp1q(d.navtime_jul(ii),d.slon(ii),julian(p.time_start)+dt2m));
    slat=medianan(interp1q(d.navtime_jul(ii),d.slat(ii),julian(p.time_start)+dt2m));
keyboard
    p.nav_start=[fix(slat), (slat-fix(slat))*60, fix(slon), (slon-fix(slon))*60];
    elon=medianan(interp1q(d.navtime_jul(ii),d.slon(ii),julian(p.time_end)+dt2m));
    elat=medianan(interp1q(d.navtime_jul(ii),d.slat(ii),julian(p.time_end)+dt2m));
    p.nav_end=[fix(elat), (elat-fix(elat))*60, fix(elon), (elon-fix(elon))*60];
  
    % interpolate on RDI data
    % this also shortens vectors to length of d.time_jul, which may be the
    % only thing that is really needed
    d.slon=interp1q(d.navtime_jul(ii),d.slon(ii),d.time_jul')';
    d.slat=interp1q(d.navtime_jul(ii),d.slat(ii),d.time_jul')';

    p.poss = [NaN NaN NaN NaN];
    p.pose = [NaN NaN NaN NaN];
  end

  p=setdefv(p,'poss',[NaN NaN NaN NaN]);
  [slat,slon] = pos2str(p.poss(1)+p.poss(2)/60,p.poss(3)+p.poss(4)/60);
  disp([' update start pos  from:',slat,'  ',slon])
  p.poss=p.nav_start;
  [slat,slon] = pos2str(p.poss(1)+p.poss(2)/60,p.poss(3)+p.poss(4)/60);
  disp(['			to:',slat,'  ',slon])
  
  p=setdefv(p,'pose',[NaN NaN NaN NaN]);
  [slat,slon] = pos2str(p.pose(1)+p.pose(2)/60,p.pose(3)+p.pose(4)/60);
  disp([' update end pos    from:',slat,'  ',slon])
  p.pose=p.nav_end;
  [slat,slon] = pos2str(p.pose(1)+p.pose(2)/60,p.pose(3)+p.pose(4)/60);
  disp(['			to:',slat,'  ',slon])

end % if p.navdata

% =================================================================
% - at this point nav data is in d.navtime_jul, d.slon, d.slat
%   and p.navdata is set to 1
% - time shifting & extraction of begin/end position is handled
%   in [loadctd.m]
% =================================================================

if ~isfinite(p.drot)						% set magdecl
 if p.magdec_source == 1					% external magdec program
   [s,o] = system('magdec');
   if s == 1
     p.drot = geomag(f,meannan(d.navtime_jul),medianan(d.slat),medianan(d.slon));
   else    
     warn = sprintf('"magdec" not found; using magdev Matlab code');
     disp(['WARNING: ' warn]);
     p.warn(size(p.warn,1)+1,1:length(warn))=warn;
   end
 end

 if ~isfinite(p.drot)
   p.drot = magdev(medianan(d.slat),medianan(d.slon),0,p.time_start(1));
 end

 [d.ru,d.rv]=uvrot(d.ru,d.rv,p.drot);
 [d.bvel(:,1),d.bvel(:,2)]=uvrot(d.bvel(:,1),d.bvel(:,2),p.drot);
 disp(sprintf(' corrected for magnetic declination of %.1f deg',p.drot));
end

% ================================================================

function depth=p2z(p,lat)
% !!!!!! USES Z=0 AT P=0  (I.E. NOT 1ATM AT SEA SURFACE)
%	pressure to depth conversion using
%	saunders&fofonoff's method (deep sea res.,
%	1976,23,109-111)
%	formula refitted for alpha(p,t,s) = eos80
%	units:
%		depth         z        meter
%		pressure      p        dbars  (original in bars, but below
%                                              division by 10 is included)
%		latitude      lat      deg
%	checkvalue:
%		depth =       9712.654  m
%	for
%		p     =         1000.     bars
%		lat   =           30.     deg
%	real lat,p
        if nargin < 2, lat=54; end
        p=p/10.;
	x=sin(lat/57.29578);
	x=x*x;
	gr=9.780318*(1.0+(5.2788e-3+2.36e-5*x)*x)+1.092e-5*p;
	depth=(((-1.82e-11*p+2.279e-7).*p-2.2512e-3).*p+97.2659).*p;
	depth=depth./gr;
%

function a=datestrj(b)
% 
% print julian date string
%
a=datestr(b-julian([0 1 0 0 0 0]));

%
%==============================================================
function [hours]=hms2h(h,m,s);
%HMS2H converts hours, minutes, and seconds to hours
%
%  Usage:  [hours]=hms2h(h,m,s);   or [hours]=hms2h(hhmmss);
%
if nargin== 1,
   hms=h;
   h=floor(hms/10000);
   ms=hms-h*10000;
   m=floor(ms/100);
   s=ms-m*100;
   hours=h+m/60+s/3600;
else
   hours=h+(m+s/60)/60;
end

%===============================================================
function  dev=geomag(f,date,lat,lon);
% function  dev=geomag(f,date,lat,lon);
% 
% call SOEST magdec to compute magnetic deviation

% INSTALLATION INSTRUCTIONS:
%	- src available as "geomag" at http://currents.soest.hawaii.edu/hg/hgwebdir.cgi
%	- on UNIX systems (I tested Linux, MacOSX & FreeBSD), 
%		1) compile by typing "make" or "gmake" in source directory
%		2) install by typing "make install" or "gmake install" as root
%		3) test by executing matlab command "system('magdec')"
%			if this test produces an error and a return value of 127, 
%			the path is not set correctly

dstr = gregoria(date);					% convert date (approx)
year = dstr(1); month = dstr(2); day = dstr(3);
 if (year < 1980 || year > 2031)
	error(sprintf('year = %d out of range',year));
end
							% execute external program
CMD = sprintf('magdec %g %g %d %d %d',lon,lat,year,month,day);
disp(sprintf('executing %s',CMD));
[status,work] = system(CMD);
if status ~= 0
	error(['cannot execute <' CMD '>']);
end

vals = sscanf(work,'%g');				% parse output
if length(vals) ~= 4
	error(['unexpected output from <' CMD '>']);
end

dev = vals(1);						% return result
