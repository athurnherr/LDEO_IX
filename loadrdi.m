function [d,p,de]=loadrdi(f,p)
% function [d,p,de]=loadrdi(f,p)
% LADCP-2 software 
%
% Martin Visbeck and Gerd Krahmann , LDEO April-2000
% read RDI raw data
% thanks to Christian Mertens who contributed most of the functions
%
% Added pg_save variable which saves percent good data, 6/2003
% Added added up-looking beam coordiante data, 6/2004
%

%======================================================================
%                    L O A D R D I . M 
%                    doc: Fri Jun 18 18:21:56 2004
%                    dlm: Tue Jun 25 10:46:50 2024
%                    (c) 2004 ladcp@
%                    uE-Info: 55 74 NIL 0 0 72 2 2 8 NIL ofnI
%======================================================================

% CHANGES BY ANT
%  Jun 18, 2004: - added p.mask_dn_bins, p.mask_up_bins (later moved)
%  Jun 21, 2004: - clarified large-velocity warning message
%  Jun 22, 2004: - changed large-velocity warning to check only
%		   central hour of cast, as large velocities were found
%	           to occur commonly in uplooker near beginning and end
%		   of cast, but not in middle. The max number of allowed
%		   large velocities before warning is issued was reduced
%		   from 25 to zero, on the other hand.
%		 - removed logging messages that are not useful in general
%		   (e.g. number of bytes per ensemble)
%		 - removed `removed' messages if no ensembles were removed
%  Jun 24, 2004: - BUG: new large-velocity warning code failed for casts
%		        shorter than one hour
%  Jul  4, 2004: - BUG: new large-velocity warning code failed because
%		        in some cases l.tim contains NaNs
%  Jul 21, 2004: - moved bin masking to [edit_data.m]
%  Nov 17, 2007: - bad error message in b2earth()
%		 - removed a lot of commented-out code from b2earth()
%  Nov 18, 2007: - added code for p.allow_3beam_solutions, p.ignore_beam
%  Jan  7, 2009: - tightened use of exist()
%  Oct 28, 2009: - modified p.ignore_beam for dual-head systems
%  Jun 27, 2011: - l.blen removed because bin lenght can be different for UL & DL
%		 - apparently unused z-variable commented out
%  Jun 30, 2011: - buggy bin-remapping disabled
%  Aug 18, 2011: - added comment to coord-transformation code (gimbal pitch)
%  Jun 24, 2013: - blen re-added but separately for DL/UL
%		 - added separate nbin, blnk, dist for DL/UL to p struct
%  Jan 23, 2015: - made updown() bomb when UL file is not found
%  Apr 15, 2015: - modified ambiguity-velocity warning as suggested by Diana Cardoso
%  May 27, 2015: - clarified time-related warnings
%  Feb 23, 2016: - clarified header id error message
%  Jan 27, 2020: - moved magdev call further back, where start time is known
%		   (with GK's new magdev code, this howto is miraculously correct again - I think)
%  Jun 25, 2024: - BUG: RDI BT data were used when processing only UL data

% p=setdefv(p,'pg_save',[1 2 3 4]);
% Default =3 for loadctd_whoi.
p=setdefv(p,'drot',nan);
% how many db should the last bin be below bin 1
p=setdefv(p,'ts_signal_min',-5);

p=setdefv(p,'ignore_beam',[nan nan]);

if nargin<2, p.name='unknown'; end

if existf(f,'ladcpdo')==0
 error([' need file name f.ladcpdo  !!!!! '])
end

if ~exist(f.ladcpdo,'file') | length(f.ladcpdo)==sum(f.ladcpdo==' ')
 error([' can not find ADCP data file : ',f.ladcpdo])
end

f=setdefv(f,'ladcpup',' ');

disp('LOADRDI:')

%
% set some defaults, in case they have not been set yet
%
p=setdefv(p,'pglim',0); 
p=setdefv(p,'elim',0.5); 
p=setdefv(p,'vlim',2.5); 
p=setdefv(p,'wlim',0.2);
p=setdefv(p,'timoff',0);
p=setdefv(p,'drot',NaN);
p=setdefv(p,'orig',0);
p=setdefv(p,'weighbin1',1);
p=setdefv(p,'tiltmax',[22 4]);
p=setdefv(p,'name','unknown');
p=setdefv(p,'maxbinrange',0);
p=setdefv(p,'ts_save',0);
p=setdefv(p,'cm_save',0);
p=setdefv(p,'pg_save',3);
p=setdefv(p,'xmv_min',0);

%
% load RDI data
%
[l,message,le]=updown(f.ladcpdo,f.ladcpup,p.pglim,p.elim,...
                      p.maxbinrange,p.ts_save,p.cm_save,p.pg_save,p);
p.warn=l.warn;


% store original data
%
if p.orig
 d.l = l;
end

% check for data set
%
de=0;
if le==1
  disp(' !!!  NO DATA !!! ')
  disp(message)
  de=1;
  return
end


%
% remember which bin come from which instrument (up-down) configuration
% get instrument configuration
%
d.izd		= 1:length(l.zd);
d.zd		= l.zd;
d.bbadcp	= l.bbadcp;

p.serial_cpu_d	= l.serial_cpu_d;
p.nping_total	= l.npng_d*l.nens_d;
p.instid(1)	= prod(p.serial_cpu_d+1)+sum(p.serial_cpu_d);
p.blen_d 	= l.blen_d;
p.nbin_d 	= l.nbin_d;
p.blnk_d 	= l.blnk_d;
p.dist_d 	= l.dist_d;

[dummy,d.down]=rditype(f.ladcpdo);
if d.down.Up
  warn=(' up looking instrument detected in do-file');
  p.warnp(size(p.warnp,1)+1,1:length(warn))=warn;
  disp(warn)
  d.zd=-d.zd;
end
 
if existf(l,'zu')
   d.izu		= fliplr(1:length(l.zu));
   d.izd		= d.izd+length(d.izu);
   d.zu			= l.zu;

   p.serial_cpu_u	= l.serial_cpu_u;
   p.instid(2)		= prod(p.serial_cpu_u+1)+sum(p.serial_cpu_u);
   p.nping_total(2)	= l.npng_u*l.nens_u;
   p.blen_u		= l.blen_u;
   p.nbin_u		= l.nbin_u;
   p.blnk_u		= l.blnk_u;
   p.dist_u     	= l.dist_u;
   
   [dummy,d.up]=rditype(f.ladcpup);
   if d.up.Up==0
    warn=(' down looking instrument detected in up-file');
    p.warnp(size(p.warnp,1)+1,1:length(warn))=warn;
    disp(warn)
    d.zu=-d.zu;
   end
else
   d.izu=[];
   d.zu=[];
end

%
% apply w velocity threshold
%
d.wrange=5;
izr=d.izd(1:d.wrange);
w(d.izd,:)=meshgrid(medianan(l.w(izr,:),1),d.izd,2);
if existf(l,'zu') 
  izr=[izr,d.izu(1:d.wrange)]; 
  w(d.izu,:)=meshgrid(medianan(l.w(d.izu(1:d.wrange),:),1),d.izu,2);
end
p=setdefv(p,'wizr',izr);

% to normal velocity data
j = find(abs(l.w-w) > p.wlim);
l.u(j) = NaN;
l.v(j) = NaN;
l.w(j) = NaN;
if p.orig
 d.l.problem(j) = d.l.problem(j)+1;
end

% Estimate single ping velocity error from std(W)

nmax=min([length(d.izd),6]);
sw=stdnan(l.w(d.izd(2:nmax),:));
ii=find(sw>0); sw=medianan(sw(ii));
d.down.Single_Ping_Err=sw/tan(d.down.Beam_angle*pi/180)/...
                        sqrt(d.down.Pings_per_Ensemble);
p.beamangle=d.down.Beam_angle;
if existf(l,'zu') 
 nmax=min([length(d.izu),6]);
 sw=stdnan(l.w(d.izu(2:nmax),:));
 ii=find(sw>0); sw=medianan(sw(ii));
 d.up.Single_Ping_Err=sw/tan(d.up.Beam_angle*pi/180)/...
                        sqrt(d.up.Pings_per_Ensemble);
end

% to bottom track velocity data
j = find(abs(l.wb-w(d.izd(1),:)) > p.wlim);
l.ub(j) = NaN;
l.vb(j) = NaN;
l.wb(j) = NaN;
if p.orig
 d.l.problemb(j) = d.l.problemb(j)+1;
end

% Horizontal Velocity limit
%
% to normal velocity data

vel=sqrt(l.u.^2+l.v.^2);
j = find(vel > p.vlim);
l.u(j) = NaN;
l.v(j) = NaN;
l.w(j) = NaN;
if length(j) > 0
  disp(sprintf(' removed %d values because of horizontal speed > %g m/s',length(j),p.vlim));
end

% only warn if large velocities occur during middle hour of cast; this
% excludes near-surface effects when large velocities can be common.
% However, reduce number of allowed large velocities before warning is
% issued from 25 to 10.

nens = size(l.u,2);
enstime = 86400 * (max(l.tim(1,:))-min(l.tim(1,:))) / nens;
jstart = floor(nens/2-1800/enstime);
jend = ceil(nens/2+1800/enstime);
if (jstart < 1), jstart = 1; end
if (jend > length(vel)), jend = length(vel); end
jj = find(vel(jstart:jend) > p.vlim);

skipnens = 1200 / enstime;
[j1,j2] = find(vel > p.vlim);
jj = find(j2>skipnens & j2<size(l.u,2)-skipnens);

if length(jj)>10
 warn = sprintf('** found %d (%.1f%% of total) velocity measurements > %g m/s',...
		length(jj),length(jj)/length(vel)*100,p.vlim);
 disp(warn);
 if length(jj)>100
  p.warn(size(p.warn,1)+1,1:length(warn))=warn;
 end
 disp('** WARNING  check ambiguity velocity setting in CMD-file   ** ')
end
if p.orig
 d.l.problem(j) = d.l.problem(j)+1;
end

% to bottom track velocity data
vel=sqrt(l.ub.^2+l.vb.^2);
j = find(vel > p.vlim);
l.ub(j) = NaN;
l.vb(j) = NaN;
l.wb(j) = NaN;
if p.orig
 d.l.problemb(j) = d.l.problemb(j)+1;
end

%
% apply a time offset, if given
%
if p.timoff~=0
 disp([' WARNING adjusted ADCP time by ',num2str(p.timoff),' days']),
 l.tim=l.tim+p.timoff;
end

% cut time range and apply a time offset, if given

if existf(p,'time_start')==0
  disp(' using whole profile since no start time was given')
  it=find(isfinite(l.tim(1,:)));
  p.time_start=gregoria(l.tim(1,it(1)));
  p.time_end=gregoria(l.tim(1,it(end)));
  d.time_jul=l.tim(1,it);
else
  % fix time for NB-ADCP data
  if l.bbadcp==0
   dum=gregoria(l.tim(1));
   l.tim=l.tim-julian(dum)+julian([p.time_start(1) dum(2:end)]);
   disp(' adjust year for NB-ADCP using given start time ')
  end
  d.time_jul=l.tim(1,:);
  it=find(julian(p.time_start)<=d.time_jul & julian(p.time_end)>=d.time_jul);
  it2=find(julian(p.time_start)>d.time_jul | julian(p.time_end)<d.time_jul);
  if p.orig
   d.l.problem(:,it2) = d.l.problem(:,it2)+10;
   d.l.problemb(:,it2) = d.l.problemb(:,it2)+10;
  end
  d.time_jul = d.time_jul(it);
  disp([' extracting ',int2str(length(it)),' ensembles as profile'])
end


%
% check whether the given time ranges have lead to a profile
%
if length(it)==0
  disp(' ')
  disp(' given time range resulted in empty extracted array!')
  disp(' ')
  disp(' check times')
  disp([' given start time :'])
  disp(p.time_start)
  disp([' internal LADCP start time :'])
  disp(gregoria(l.tim(1,1)))
  disp([' given end time :'])
  disp(p.time_end)
  disp([' internal LADCP end time :'])
  disp(gregoria(maxnan(l.tim(1,:))))
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
  disp('Will try to use all data instead')
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
  pause(2)
  it=find(isfinite(l.tim(1,:)));
  d.time_jul=l.tim(1,it);
%  p.time_start=gregoria(l.tim(it(1)));
%  p.time_end=gregoria(l.tim(it(end)));
end

if existf(p,'poss')
 if isfinite(sum(p.poss))
  drot=magdev(p.poss(1)+p.poss(2)/60, p.poss(3)+p.poss(4)/60,0,p.time_start(1));
  if existf(p,'drot')
   if isfinite(p.drot)
    disp([' found drot:',num2str(p.drot),' should be ',num2str(drot)])
   else
    p.drot=drot;
   end
  else
   p.drot=drot;
  end
 end
end
 
%
% rotate for magnetic deviation
%

d.soundc=0;
if isfinite(p.drot)
 [d.ru,d.rv]=uvrot(l.u(:,it),l.v(:,it),p.drot);
 [ub,vb]=uvrot(l.ub(it),l.vb(it),p.drot);
 disp(' apply magnetic deviation, rotate bottom track and water velocities')
else
 d.ru=l.u(:,it);
 d.rv=l.v(:,it);
 ub=l.ub(it);
 vb=l.vb(it);
end

%
% save all bottom track data together and remove bad data
%
d.bvel=[ub',vb',l.wb(it)',l.eb(it)'];
ii=find(d.bvel<-30); 
d.bvel(ii)=NaN; 
d.hbot=l.hb(it);
d.hbot4=l.hb4(:,it);
p.btrk_used=l.btrk_used;
if existf(l,'hs')==1
 d.hsurf=l.hs(it);
end



%
% extract the profile from all recorded data
%
d.firstlastindx = [it(1),it(end)];
d.rw=l.w(:,it);
d.re=l.e(:,it);
d.ts=l.ts(:,it);
if p.ts_save(1)~=0
 d.ts_all_d=l.ts_all_d(it,:,:);
 if size(l.pit,1)==2
  d.ts_all_u=l.ts_all_u(it,:,:);
 end
end
if p.cm_save(1)~=0
 d.cm_all_d=l.cm_all_d(it,:,:);
 if size(l.pit,1)==2
  d.cm_all_u=l.cm_all_u(it,:,:);
 end
end
if p.pg_save(1)~=0
 d.pg_all_d=l.pg_all_d(it,:,:);
 if size(l.pit,1)==2
  d.pg_all_u=l.pg_all_u(it,:,:);
 end
end
d.hdg=l.hdg(:,it);
d.xmc=l.xmc(:,it);
d.xmv=l.xmv(:,it);
d.tint=l.tint(:,it);
d.sv=l.sv(:,it);
d.temp=l.t(:,it);
d.weight=l.cm(:,it);
d.weight=d.weight./medianan(maxnan(d.weight));
d.pit = l.pit(:,it);
d.rol = l.rol(:,it);
%d.tilt=sqrt(l.pit(1,it).^2 + l.rol(1,it).^2);
% more accurate calculation
d.tilt=real(asin(sqrt(sin(l.pit(1,it)/180*pi).^2 +...
 sin(l.rol(1,it)/180*pi).^2)))/pi*180;

% compute tilt difference
rold=mean(abs(diff([0,d.rol(1,:);d.rol(1,:),0]'))');
pitd=mean(abs(diff([0,d.pit(1,:);d.pit(1,:),0]'))');
d.tiltd=sqrt(rold.^2+pitd.^2);

% reduce weight for strong tilt difference
ii=find(d.tilt>p.tiltmax(1));
if length(ii) > 0
  disp([' removed ',num2str(length(ii)),...
        ' profiles due to tilt > ',num2str(p.tiltmax(1)) ' degrees'])
end
d.weight(:,ii)=NaN;
if length(ii)>length(d.tilt)*0.1;
  warn=([' ',int2str(length(ii)*100/length(d.tilt)),...
         '%  tilt > ',int2str(p.tiltmax(1)),' ']);
  disp(warn)
  p.warn(size(p.warn,1)+1,1:length(warn))=warn;
  p.warnp(size(p.warnp,1)+1,1:length(warn))=warn;

end

ii=find(d.tiltd>p.tiltmax(2));
if length(ii) > 0
  disp([' removed ',num2str(length(ii)),...
        ' profiles due to tilt derivative > ',num2str(p.tiltmax(2)) ' degrees'])
end        
d.weight(:,ii)=NaN;

% reduce weight for strong echos possibly from crosstalk or bottom
d.tsw=d.weight*0+1;
for i=1:size(d.tsw,1)
 tsmed=median(d.ts(i,:));
 ts=d.ts(i,:)-tsmed;
 ii=find(ts>0);
 if length(ii)>0
  d.weight(i,ii)=d.weight(i,ii).*(1-(ts(ii)/max(ts)).^1.5);
  d.tsw(i,ii)=(1-(ts(ii)/max(ts)).^1.5);
 end
end

% save transmit current volt and internal temperature
for j=1:size(d.xmv,1)
 p.xmc(j)=medianan(d.xmc(j,:),size(d.xmc,2)/4);
 p.xmv(j)=medianan(d.xmv(j,:),size(d.xmv,2)/4);
 p.tint(j)=medianan(d.tint(j,:),size(d.tint,2)/4);
end

% warn if battery low
if p.xmv(1)<p.xmv_min
  warn=([' median Xmit-volt ',num2str(p.xmv(1)),' < ',...
           num2str(p.xmv_min),' BATTERY weak ? ']);
  disp(warn)
  p.warn(size(p.warn,1)+1,1:length(warn))=warn;
end

% check for dead instrument (no pings just listen)
%
% remove suspect ensembles
drw=medianan(abs(diff(d.rw(d.izd,:))));
dru=medianan(abs(diff(d.rv(d.izd,:))));
drv=medianan(abs(diff(d.ru(d.izd,:))));
nbad=find(abs(drw)<0.005 & abs(dru)<0.005 & abs(drv)<0.005);
if length(nbad) > 0.2*length(it)
  warn=([' down looker ',int2str(length(nbad)),' ensembles  ',...
           ' have no flow gradient. ']);
  disp(warn)
  p.warn(size(p.warn,1)+1,1:length(warn))=warn;
  warn=(['DOWN LOOKER NOT PINGING ?']);
  disp(warn)
  p.warn(size(p.warn,1)+1,1:length(warn))=warn;
  disp(' WARNING WARNING ')
end
if length(nbad)>0
 d.weight(d.izd,nbad)=nan;
 disp([' removed ',int2str(length(nbad)),...
        ' suspect non pinging? (low velocity gradient) ensembles from down-looker'])
end

if length(d.izu)>1
 % remove suspect ensembles
 drw=medianan(abs(diff(d.rw(d.izu,:))));
 dru=medianan(abs(diff(d.rv(d.izu,:))));
 drv=medianan(abs(diff(d.ru(d.izu,:))));
 nbad=find(abs(drw)<0.005 & abs(dru)<0.005 & abs(drv)<0.005);
 if length(nbad) > 0.2*length(it)
   warn=(['   up looker ',int2str(length(nbad)),' ensembles  ',...
            ' have no flow gradient.']);
   disp(warn)
   p.warn(size(p.warn,1)+1,1:length(warn))=warn;
   warn=(['UP LOOKER NOT PINGING ?']);
   disp(warn)
   p.warn(size(p.warn,1)+1,1:length(warn))=warn;
   disp(' WARNING WARNING ')
 end
 if length(nbad)>0
  d.weight(d.izu,nbad)=nan;
  disp([' removed ',int2str(length(nbad)),...
        ' suspect non pinging? (low velocity gradient) ensembles from up-looker'])
 end
end

%
% reduce certainty in bin 1
%
idb1=d.izd(1);
if length(idb1)==1
  d.weight(idb1,:)=d.weight(idb1,:)*p.weighbin1;
end
if length(d.izu)>1
  iub1=find(d.izu==1);
  if length(iub1)==1
    d.weight(iub1,:)=d.weight(iub1,:)*p.weighbin1;
  end
end
if p.weighbin1~=1
 disp([' multiply weight of bin 1 by ',num2str(p.weighbin1)])
end


%
% prepare array
%
d.izm=d.weight+NaN;

%
% save mean correlation and echo amp profile
%
if length(size(l.tsd_m))==2
  d.tsd_m=reshape(l.tsd_m,length(l.tsd_m)/4,4);
else
  d.tsd_m=squeeze(l.tsd_m);
end
if length(size(l.cmd_m))==2
  d.cmd_m=reshape(l.cmd_m,length(l.cmd_m)/4,4);
else
  d.cmd_m=squeeze(l.cmd_m);
end

t=d.tsd_m(1,:);

if abs(min(t)-median(t))>15
  disp('!!!! WARNING one beam might be broken !!!!!!!!!')
  [m,it]=min(t);
  warn=(['  broken down looking beam ',int2str(it)]);
  disp(warn)
  p.warn(size(p.warn,1)+1,1:length(warn))=warn;
elseif abs(min(t)-median(t))>5
  disp('WARNING one beam might be weak !')
  [m,it]=min(t);
  warn=(['  weak down looking beam ',int2str(it)]);
  disp(warn)
  p.warn(size(p.warn,1)+1,1:length(warn))=warn;
end

if existf(l,'tsu_m')
  if length(size(l.tsu_m))==2
    d.tsu_m=reshape(l.tsu_m,length(l.tsu_m)/4,4);
  else
    d.tsu_m=squeeze(l.tsu_m);
  end
  if length(size(l.cmu_m))==2
    d.cmu_m=reshape(l.cmu_m,length(l.cmu_m)/4,4);
  else
    d.cmu_m=squeeze(l.cmu_m);
  end

  t=d.tsu_m(1,:);
  if abs(min(t)-median(t))>10
    disp('!!!! WARNING one beam might be broken !!!!!!!!!')
    [m,it]=min(t);
    warn=(['  broken up looking beam ',int2str(it)]);
    disp(warn)
    p.warn(size(p.warn,1)+1,1:length(warn))=warn;
  elseif abs(min(t)-median(t))>5
    disp('WARNING one beam might be weak ! ')
    [m,it]=min(t);
    warn=(['  weak up looking beam ',int2str(it)]);
    disp(warn)
    p.warn(size(p.warn,1)+1,1:length(warn))=warn;
  end

end

[p.nbins,p.nt]=size(d.ru);

% check for outlier within the whole data set
[d,p]=outlier(d,p);

%-----------------------------------------------------------
function [l,message,le] = updown(fdown,fup,pglim,elim,...
                                 bmax,tssave,cmsave,pgsave,p)
%UPDOWN Load and merge upward and downward looking ADCP raw data.
%  L = UPDOWN('filedown','fileup') reads ADCP raw data from the specified
%  files. L is a structure array with the following fields:
%
%    blen: bin length		%%% ANT: REMOVED 2011/06/28 BECAUSE IT CAN BE DIFFERENT FOR UL/DL
%    nbin: number of bins
%    blnk: blank after transmit
%    dist: distance of bin 1 from transducer
%     tim: time axis
%     pit: pitch
%     rol: roll
%     hdg: heading
%       s: salinity
%       t: temperature
%      sv: sound velocity
%       u: east velocity
%       v: north velocity
%       w: vertical velocity
%       e: error velocity
%      ts: target strength
%      cm: correlation
%      hb: bottom track distance
%      ub: bottom track east velocity
%      vb: bottom track north velocity
%      wb: bottom track vertical velocity
%      eb: bottom track error velocity
%
%  [L,MESSAGE] = UPDOWN('filedown','fileup') returns a system dependent error
%  message if the opening of 'filedown' is not successful. In this case -1 is
%  returned for L.

%  Christian Mertens, IfM Kiel

% default editing parameters
if nargin < 5
  bmax = 0;
end
if nargin < 4
  elim = 0.5;
end
if nargin < 3
  pglim = 0.3;
end

% fixed leader
f.nbin = 1;   % number of depth cells
f.npng = 2;   % pings per ensemble
f.blen = 3;   % depth cell length
f.blnk = 4;   % blank after transmit
f.dist = 5;   % distance to the middle of the first depth cell
f.plen = 6;   % transmit pulse length
f.serial = '7:14'; % serial number of CPU board

% variable leader
v.tim = 1;    % true time (Julian days)
v.pit = 2;    % pitch
v.rol = 3;    % roll
v.hdg = 4;    % heading
v.t   = 5;    % temperature
v.s   = 6;    % salinity
v.sv  = 7;    % sound velocity
v.xmc = 8;    % transmit current
v.xmv = 9;    % transmit volt
v.tint = 10;  % internal temperature

% load downward looking ADCP
[fid,message] = fopen(fdown,'r','l');
le=0;
if fid == -1
  le = 1;
  message = sprintf('%s: %s',fdown,message);
  disp(' LOADRDI problem with down looking RDI file ')
  disp(message)
  error('terminate LADCP processing')
end
disp([' loading down-data ',fdown])

% check if BB data
if isbb(fid)
 [fd,vd,veld,cmd,ead,pgd,btd] = rdread(fid);
 l.bbadcp=1;
else
 fclose(fid);
 fid = fopen(fdown,'r','b');
 [fd,vd,veld,swd,ead,pgd] = nbread(fid);
 cmd=ead*0+100;
 btd = NaN*ones(size(vd,1),1,16);
 
 ok = double(prod(veld,3)~=sum(veld,3));
 ii=find(ok==0);
 ok(ii)=NaN;
 disp([' removed ',int2str(length(ii)),...
 ' values because of 0 in nbdata'])
 for k = 1:4
  veld(:,:,k) = veld(:,:,k).*ok;
 end
 l.bbadcp=0;
end

% check for beam coordinates
[dummy,dd]=rditype(fdown);
if dd.Coordinates==0
 disp(' DETECTED BEAM coordinates: rotating to EARTH coordinates')
 veld=b2earth(veld,vd,dd,p,p.ignore_beam(1));
end
%

fclose(fid);
l1=size(veld,1);
l2=size(veld,2);
disp([' read ',int2str(l1),' ensembles with ',int2str(l2),' bins each']) 

% remove extra (?) bottom track dimension
btd = squeeze(btd);
if ndims(btd)>2
  warning(' removal of extra bottom track dimension failed !!!')
end


% median echoamplitude and correlation
%ead = targs(mean(ead,3)',z(:))';
ead_m=medianan(ead);
cmd_m=medianan(cmd);
if tssave(1)~=0
 ead_all=ead(:,:,tssave);
end
if cmsave(1)~=0
 cmd_all=cmd(:,:,cmsave);
end
if pgsave(1)~=0
 pgd_all=pgd(:,:,pgsave);
end
ead=median(ead,3);
cmd=median(cmd,3);

% prepare vector containing info why a value has been discarded
% b is the same for the bottom track data
%
% 1st digit : w deviates more than p.wlim from median w (checked later)
% 2nd digit : out of time range    (checked later)
% 3rd digit : below percent good threshold
% 4th digit : above error velocity threshold
% 5th digit : 3 beam solution
% 6th digit : no vel
%
%GK
dproblem = repmat(0,[size(veld,1),size(veld,2)]);
problemb = repmat(0,[size(btd,1),size(btd,2)]);


% grep 3-beam solution and no velocities at all
i = find( isnan(veld(:,:,4)) & ~isnan(veld(:,:,3)) );
dproblem(i) = dproblem(i) + 10000;
i = find( isnan(veld(:,:,1)) );
dproblem(i) = dproblem(i) + 100000;

l.warn=('LADCP WARNINGS');
l.btrk_used = 0;
if sum(isfinite(btd(:)))>0
  [dummy,db]=rditype(fdown);
  if db.Up
    warn=(' Warning: ignoring RDI bottom track data from upward-looking instrument');
    disp(warn)
    l.warn(size(l.warn,1)+1,1:length(warn))=warn;
  else    
    l.btrk_used = 1;
  end
end

% transform to earth coordinates
if l.btrk_used == 1
 [dummy,db]=rditype(fdown);
 if db.Coordinates==0
  for i=1:4
   velb(:,1,i)=btd(:,4+i);
   velb(:,2,i)=NaN;
  end
  disp(' DETECTED BEAM bottom track coordinates!')
  db.use_binremap=0;
  velb=b2earth(velb,vd,db,p,p.ignore_beam(1));
  for i=1:4
   btd(:,4+i)=velb(:,1,i);
  end
 end
end

% apply percent-good threshold
pgd = pgd(:,:,4);
i = pgd < pglim;
pgd(i) = NaN;
pgd(~i) = 1;
if length(find(i)) > 0
	disp(sprintf(' removed %d downlooker values because of percent good < %g',...
		length(find(i)),pglim));
end
for k = 1:4
  veld(:,:,k) = veld(:,:,k).*pgd;
end
dproblem(i) = dproblem(i) + 100;
 
% load upward looking ADCP
up = nargin>1;

if strcmp(fup,'') || strcmp(fup,' ')
  up = 0;
end

if up
  fid = fopen(fup,'r','l');
  if fid == -1
    error(sprintf('%s: no such file or directory',fup));
  end

  if length(bmax)<2, bmax(2)=bmax(1); end
  disp([' loading up-data ',fup])
  [fu,vu,velu,cmu,eau,pgu,btu] = rdread(fid);

  % check for beam coordinates
  [dummy,du]=rditype(fup);
  if du.Coordinates==0
   disp(' DETECTED BEAM coordinates: rotating to EARTH coordinates')
   velu=b2earth(velu,vu,du,p,p.ignore_beam(2));
  end
  %

  l1=size(velu,1);
  l2=size(velu,2);
  btu = squeeze(btu);
  disp([' read ',int2str(l1),' ensemble and ',int2str(l2),' bins ']) 
  % median echoamplitude and correlation
  %eau = targs(mean(eau,3)',z(:))';
  eau_m=medianan(eau);
  cmu_m=medianan(cmu);
  if tssave(1)~=0
   eau_all=eau(:,:,tssave);
  end
  if cmsave(1)~=0
   cmu_all=cmu(:,:,cmsave);
  end
  if pgsave(1)~=0
   pgu_all=pgu(:,:,pgsave);
  end
  eau=median(eau,3);
  cmu=median(cmu,3);
  fclose(fid);

  % prepare vector containing info why a value has been discarded
  % b is the same for the bottom track data
  uproblem = repmat(0,[size(velu,1),size(velu,2)]);

  % grep 3-beam solution and no velocities at all
  i = find( velu(:,:,4) == 0 & velu(:,:,3)~=0 );
  uproblem(i) = uproblem(i) + 10000;
  i = find( velu(:,:,1) == 0 );
  uproblem(i) = uproblem(i) + 100000;

  % apply percent-good threshold
  pgu = pgu(:,:,4);
  i = pgu < pglim;
  if length(find(i)) > 0
  	  disp(sprintf(' removed %d uplooker values because of percent good < %g',...
	 	  length(find(i)),pglim));
  end
  pgu(i) = NaN;
  pgu(~i) = 1;
  for k = 1:4
    velu(:,:,k) = velu(:,:,k).*pgu;
  end
  uproblem(i) = uproblem(i) + 100; 

end

% distance vectors
%%% z = fd(f.dist) + fd(f.blen)*([1:fd(f.nbin)] - 1);	% unused?! ANT 2011/06/28
if bmax(1)>0, fd(f.nbin)=bmax(1); end
idb=1:fd(f.nbin);

if up
   if bmax(2)>0, fu(f.nbin)=bmax(2); end
   iub=1:fu(f.nbin);

  % check if ping rate is the same for both instruments
  timd = vd(:,1,v.tim)';
  timu = vu(:,1,v.tim)';
  pingconst=abs(maxnan(diff(timd))+maxnan(-diff(timd))) > (0.05/24/3600);
  pingdiff=abs(medianan(diff(timd))-medianan(diff(timu))) > (0.05/24/3600);
  ii=min(length(timd),length(timu));
  timediff=abs((timd(1)-timd(ii)-(timu(1)-timu(ii)))) > (5/24/3600);
  if pingdiff | timediff | pingconst
   disp(' WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING ')
   if pingconst
    warn=(' Warning: non-constant ping rate in downlooker data (staggered pinging?)');
    disp(warn)
    disp(['  min down ping rate :',num2str(-24*3600*maxnan(-diff(timd))),...
         '  max down ping rate :',num2str(24*3600*maxnan(diff(timd)))])
   end
   if pingdiff
    warn=(' Warning: mean ping rates differ in downlooker/uplooker data ');
    disp(warn)
    l.warn(size(l.warn,1)+1,1:length(warn))=warn;
    disp(['  mean down ping rate :',num2str(24*3600*meannan(diff(timd))),...
         '  mean up ping rate :',num2str(24*3600*meannan(diff(timu)))])
   end
   if timediff
    warn=(' Warning: cast duration differs in downlooker/uplooker data ');
    disp(warn)
    l.warn(size(l.warn,1)+1,1:length(warn))=warn;
    disp(['  down dt for common ping number:',num2str((timd(ii)-timd(1))*24),...
          '  up dt :',num2str((timu(ii)-timu(1))*24),' hours '])
   end
   iu=1:length(timd);
   ii=find(iu>length(timu));
   iu(ii)=length(timu);
   disp(' find best time match of up-looking ADCP to down looking ADCP')
   for i=find(isfinite(timd))
    [m,iu(i)]=min(abs(timu-timd(i)));
   end
   ilast=min(length(iu),length(timu));
   disp([' up instrument is different by ',num2str(iu(ilast)-ilast),...
          ' ensembles'])
   disp(' WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING ')
   id=1:length(timd);

  else

   id=1:length(timd);
   iu=1:length(timu);

  end
 
   % find best lag to match up vertical velocity 
   wu=squeeze(velu(iu,:,3));
   wd=squeeze(veld(id,:,3));
   wb2u=medianan(wu');
   wb2d=medianan(wd');

   maxlag=20;
   [lag,iiu,id,co]=bestlag(wb2u,wb2d,maxlag);
   disp([' try to shift timeseries by lag: ',num2str(lag),...
         ' correlation: ',num2str(co)])
   if abs(lag)==maxlag | co<0.9, 
    disp(' best lag not obvious!  use time to match up-down looking ADCP')
    id=1:length(timd);
    iu=1:length(timd);
    ii=find(iu>length(timu));
    iu(ii)=length(timu);
    for i=find(isfinite(timd))
      [m,iu(i)]=min(abs(timu-timd(i)));
    end
    lag=mean(iu-id);
    disp([' mean lag is ',num2str(lag),' ensembles']);
   else
    disp([' shift ADCP timeseries by ',num2str(lag),' ensembles']);
    iu=iu(iiu);
   end

  disp([' number of joint ensembles is : ',num2str(length(iu))]);

  % merge upward and downward
  l.zu=[0:(fu(f.nbin)-1)]*fu(f.blen)+fu(f.dist);
  l.zd=[0:(fd(f.nbin)-1)]*fd(f.blen)+fd(f.dist);
  eval(['l.serial_cpu_u=fu(',f.serial,');']);
  eval(['l.serial_cpu_d=fd(',f.serial,');']);
  l.npng_u = fu(f.npng);
  l.npng_d = fd(f.npng);
  l.nens_u = size(vu,1);
  l.nens_d = size(vd,1);
  l.blen_u = fu(f.blen);
  l.blen_d = fd(f.blen);
  l.nbin_u = fu(f.nbin);
  l.nbin_d = fd(f.nbin);
  l.blnk_u = fu(f.blnk);
  l.blnk_d = fd(f.blnk);
  l.dist_u = fu(f.dist);
  l.dist_d = fd(f.dist);
  l.tim = [vd(id,1,v.tim),vu(iu,1,v.tim)]';
  l.pit = [vd(id,1,v.pit),vu(iu,1,v.pit)]';
  l.rol = [vd(id,1,v.rol),vu(iu,1,v.rol)]';
  l.hdg = [vd(id,1,v.hdg),vu(iu,1,v.hdg)]';
  l.s = [vd(id,1,v.s),vu(iu,1,v.s)]';
  l.t = [vd(id,1,v.t),vu(iu,1,v.t)]';
  l.sv = [vd(id,1,v.sv),vu(iu,1,v.sv)]';
  l.xmc = [vd(id,1,v.xmc),vu(iu,1,v.xmc)]';
  l.xmv = [vd(id,1,v.xmv),vu(iu,1,v.xmv)]';
  l.tint = [vd(id,1,v.tint),vu(iu,1,v.tint)]';
  l.u = [fliplr(velu(iu,iub,1)) veld(id,idb,1)]';
  l.v = [fliplr(velu(iu,iub,2)) veld(id,idb,2)]';
  l.w = [fliplr(velu(iu,iub,3)) veld(id,idb,3)]';
  l.e = [fliplr(velu(iu,iub,4)) veld(id,idb,4)]';
  l.ts = [fliplr(eau(iu,iub)) ead(id,idb)]';
  l.cm = [fliplr(cmu(iu,iub)) cmd(id,idb)]';
% No reason to keep this since pgu and pgd don't mean much anymore  
%l.pg = [fliplr(pgu(iu,iub)) pgd(id,idb)]';
  if tssave(1)~=0
   l.ts_all_u = eau_all(iu,iub,:);
   l.ts_all_d = ead_all(id,idb,:);
  end
  if cmsave(1)~=0
   l.cm_all_u = cmu_all(iu,iub,:);
   l.cm_all_d = cmd_all(id,idb,:);
  end
  if pgsave(1)~=0
   l.pg_all_u = pgu_all(iu,iub,:);
   l.pg_all_d = pgd_all(id,idb,:);
  end
% distance to surface
  hs = median(btu(iu,1:4),2)';
  if sum(isfinite(hs))>1
   l.hs = hs;
  else
% try to use targestength to find surface
   if sum(isfinite(eau))>1
    disp(' use target strength of up looking to find surface ')
    eaum=medianan(eau);
    eaua=eau-meshgrid(eaum,eau(:,1));
    [eam,hsb]=max(eaua(iu,:)');
%    l.hs=l.zu(hsb)+(l.zu(2)-l.zu(1))/2;
    l.hs=l.zu(hsb);
    ii=find(eam<20);
    l.hs(ii)=NaN;
    ii=find(hsb==1 | hsb==size(eau,2));
    l.hs(ii)=NaN;
   end
  end
  l.hb = median(btd(id,1:4),2)';
  l.hb4 = btd(id,1:4)';
  l.ub = btd(id,5)';
  l.vb = btd(id,6)';
  l.wb = btd(id,7)';
  l.eb = btd(id,8)';
  l.tsd_m=ead_m;
  l.cmd_m=cmd_m;
  l.tsu_m=eau_m;
  l.cmu_m=cmu_m;
  l.problem = [fliplr(uproblem(iu,:)) dproblem(id,:)]';
  l.problemb = problemb(id,:)';

else % single instrument

  l.zd=[0:(fd(f.nbin)-1)]*fd(f.blen)+fd(f.dist);
  eval(['l.serial_cpu_d=fd(',f.serial,');']);
  l.npng_d = fd(f.npng);
  l.nens_d = length(vd(v.tim));
  l.blen_d = fd(f.blen);
  l.nbin_d = fd(f.nbin);
  l.blnk_d = fd(f.blnk);
  l.dist_d = fd(f.dist);
  l.tim = vd(:,1,v.tim)';
  l.pit = vd(:,1,v.pit)';
  l.rol = vd(:,1,v.rol)';
  l.hdg = vd(:,1,v.hdg)';
  l.s = vd(:,1,v.s)';
  l.t = vd(:,1,v.t)';
  l.sv = vd(:,1,v.sv)';
  l.xmc = vd(:,1,v.xmc)';
  l.xmv = vd(:,1,v.xmv)';
  l.tint = vd(:,1,v.tint)';
  l.hdg = vd(:,1,v.hdg)';
  l.u = veld(:,idb,1)';
  l.v = veld(:,idb,2)';
  l.w = veld(:,idb,3)';
  l.e = veld(:,idb,4)';
  l.ts = ead(:,idb)';
  if tssave(1)~=0
   l.ts_all_d = ead_all(:,idb,:);
  end
  l.cm = cmd(:,idb)';
  if cmsave(1)~=0
   l.cm_all_d = cmd_all(:,idb,:);
  end
  l.pg = pgd(:,idb)';
  if pgsave(1)~=0
   l.pg_all_d = pgd_all(:,idb,:);
  end
  % fix to reduce funny bottom track dimension
  id=1:length(l.tim);
  l.hb = median(btd(id,1:4),2)';
  l.hb4 =btd(id,1:4)';
  l.ub = btd(id,5)';
  l.vb = btd(id,6)';
  l.wb = btd(id,7)';
  l.eb = btd(id,8)';
  l.tsd_m=ead_m;
  l.cmd_m=cmd_m;
  l.problem = dproblem';
  l.problemb = problemb';

end

if l.btrk_used == 1
 good = find(isfinite(l.wb));
 disp([' found ',int2str(length(good)),' finite RDI bottom velocities'])
end

%GK
% discard dummy error velocities
bad = find(l.eb==-32.768);
if ~isempty(bad)
  disp([' found ',int2str(length(bad)),' NaN bottom error velocities and discarded them'])
  l.eb(bad) = nan;
end

% check for 3-beam solution
jok = cumprod(size(find(~isnan(l.w))));
j = cumprod(size(find(isnan(l.e) & ~isnan(l.w))));
if j/jok > 0.2
 disp(['!!!!!!!!!!!!! WARNING  WARNING  WARNING !!!!!!!!!!!!!!'])
 warn=([' detected  ',int2str(j*100/jok),' %  3 BEAM solutions ']);
 disp(warn)
 l.warn(size(l.warn,1)+1,1:length(warn))=warn;
 disp(['!!!!!!!!!!!!! WARNING  WARNING  WARNING !!!!!!!!!!!!!!'])
end


% apply error velocity threshold
j = find(abs(l.e) > elim);
if length(j) > 0
	disp(sprintf(' removed %d values because of error velocity > %g m/s',...
		length(j),elim));
end
l.u(j) = NaN;
l.v(j) = NaN;
l.w(j) = NaN;
l.problem(j) = l.problem(j) + 1000;
j = find(abs(l.eb) > elim);
if length(j) > 0
	disp(sprintf(' removed %d bottom-track values because of error velocity > %g m/s',...
		length(j),elim));
end
l.ub(j) = NaN;
l.vb(j) = NaN;
l.wb(j) = NaN;
l.problemb(j) = l.problemb(j) + 1000;


% ---------------------------------------------------------
function varargout = rdread(fid)
%RDREAD Read RDI BB data.
%  RDREAD(FID)

%  Christian Mertens, IfM Kiel

% rewind to beginning of file and read header to get the number of bytes
% in each ensemble, the number of data types, and the address offsets
status = fseek(fid,0,'bof');
[nbytes,dtype,offset] = rdhead(fid);
% disp([' raw data has ',int2str(nbytes),'+2 bytes per ensemble'])
ntypes = length(dtype);

% get the number of ensembles from file size; each ensemble has nbytes
% plus two bytes for the checksum
status = fseek(fid,0,'eof');
m = floor(ftell(fid)/(nbytes + 2));
status = fseek(fid,0,'bof');

% number of bins is the offset difference (minus 2 bytes for the ID code)
% between velocity data and correlation magnitude devided by 4 beams
% times 2 bytes
n = (offset(4) - offset(3) - 2)/(2*4);

% data parameters
scale = [NaN,NaN,0.001,1,0.45,1,0.001];
precision = {'','','int16','uint8','uint8','uint8',''};
varid = [0,128,256,512,768,1024,1536];
bad = scale.*[NaN,NaN,-32768,0,NaN,NaN,-32768];

% initialize output variables
for k = 1:length(varid)
  if k == 1
    varargout{k} = NaN*ones(1,1,14);
  elseif k == 2
    varargout{k} = NaN*ones(m,1,10);
  elseif (k >= 3 & k <= 6)
    varargout{k} = NaN*ones(m,n,4);
  elseif k == 7
    varargout{k} = NaN*ones(m,1,16);
  end
end

% read fixed leader data
status = fseek(fid,offset(1)+2,'bof');
varargout{1} = rdflead(fid);

icheck=0;

for i = 1:m
  % read ensemble to verify the checksum
  status = fseek(fid,(i-1)*(nbytes+2),'bof');
  buffer = fread(fid,nbytes,'uint8');
  checksum = fread(fid,1,'uint16');

  % read ensemble if checksum is ok
  if checksum == rem(sum(buffer),2^16);
    for kk = 2:length(dtype)
      k = dtype(kk);
      % set file pointer to beginning of data
      status = fseek(fid,(i-1)*(nbytes+2)+offset(kk)+2,'bof');
      switch varid(k)
        case varid(2)
          % variable leader data
          varargout{k}(i,1,:) = rdvlead(fid);
        case varid(7)
          % bottom track data
          varargout{k}(i,1,:) = rdbtrack(fid);
        otherwise
          % velocity, correlation, echo intensity, or percent-good data
          a = fread(fid,4*n,precision{k});
          varargout{k}(i,:,:) = scale(k)*reshape(a,4,n)';
      end
    end
  else
   icheck=icheck+1;
  end
end

if icheck > m*0.01 
 disp([' WARNING  found ',int2str(icheck),' ensembles with bad checksum '])
end

% check for bad values
i = find(varargout{3} == bad(3));
varargout{3}(i) = NaN;
i = find(varargout{4} == bad(4));
varargout{4}(i) = NaN;
% bottom track
i = find(varargout{7}(:,1,5) == bad(7));
varargout{7}(i,1,1:8) = NaN;
i = find(varargout{7}(:,1,6) == bad(7));
varargout{7}(i,1,1:8) = NaN;
i = find(varargout{7}(:,1,7) == bad(7));
varargout{7}(i,1,1:8) = NaN;
i = find(varargout{7}(:,1,8) == bad(7));
varargout{7}(i,1,1:8) = NaN;


%-------------------------------------------------------------------------------

function [nbytes,dtype,offset] = rdhead(fid)
%RDHEAD Read the header data from a raw ADCP data file.
%  [NBYTES,DTYPE,OFFSET] = RDHEAD(FID)

hid = 127;  % header identification byte
sid = 127;  % data source identification byte

% get file position pointer
fpos = ftell(fid);

% check header and data source identification bytes
[id,n] = fread(fid,2,'uint8');
if  (n < 2 | feof(fid))
  error('Unexpected end of file.')
end
if (id(1) ~= hid | id(2) ~= sid)
  error(sprintf('Header identification byte not found (%02x %02x).',id(1),id(2)))
end

% read the number of bytes
nbytes = fread(fid,1,'uint16');

% skip spare byte
fseek(fid,1,'cof');

% read the number of data types
ndt = fread(fid,1,'uint8');
if ndt >= 8; ndt=7; end;	%%% DT bug fix 2009-01-07

% read address offsets for data types
offset = fread(fid,ndt,'uint16');

% read variable identifiers
varid = [0 128 256 512 768 1024 1536];
for i = 1:ndt
  fseek(fid,fpos+offset(i),'bof');
  id = fread(fid,1,'uint16');
  dtype(i) = find(id == varid);
end

% rewind to the beginning of the ensemble
fseek(fid,fpos,'bof');


%-------------------------------------------------------------------------------

function fl = rdflead(fid);
%RDFLEAD Read the fixed leader data from a raw ADCP data file.
%  FL = RDFLEAD(FID)

fseek(fid,7,'cof');

% number of depth cells
fl(1) = fread(fid,1,'uint8');

% pings per ensemble, depth cell length in cm, blank after transmit
fl(2:4) = fread(fid,3,'uint16');
fl(3) = 0.01*fl(3);
fl(4) = 0.01*fl(4);

fseek(fid,16,'cof');

% Bin 1 distance, xmit pulse length
fl(5:6) = 0.01*fread(fid,2,'ushort');

fseek(fid,6,'cof');
% Serial Number of CPU board
fl(7:14) = fread(fid,8,'uint8');



%-------------------------------------------------------------------------------

function vl = rdvlead(fid)
%RDVLEAD Read the variable leader data from a raw ADCP data file.
%  VL = RDVLEAD(FID)

fseek(fid,2,'cof');

% time of ensemble
c = fread(fid,7,'uint8');
c(1)=y2k(c(1));
vl(1) = julian(c(1),c(2),c(3),c(4)+c(5)/60+c(6)/3600+c(7)/360000);
fseek(fid,3,'cof');

% speed of sound (EC)
vl(7) = fread(fid,1,'uint16');
fseek(fid,2,'cof');

% heading (EH)
vl(4) = 0.01*fread(fid,1,'uint16');

% pitch (EP) and roll (ER)
vl(2:3) = 0.01*fread(fid,2,'int16');

% salinity (ES)
vl(6) = 0.001*fread(fid,1,'uint16');

% temperature (ET)
vl(5) = 0.01*fread(fid,1,'int16');
fseek(fid,6,'cof');

% ADC channels
% Transmit Current
vl(8) = fread(fid,1,'uint8');
% Transmit Volt
vl(9) = fread(fid,1,'uint8');
% Internal Temperature
vl(10) = fread(fid,1,'uint8');



%-------------------------------------------------------------------------------

function bt = rdbtrack(fid)
%RDBTRACK Read the bottom track data from a raw ADCP data file.
%  BT = RDBTRACK(FID)

fseek(fid,14,'cof');

% range
bt(1:4) = 0.01*fread(fid,4,'uint16');

% velocity
bt(5:8) = 0.001*fread(fid,4,'int16');

% correlation magnitude and percent good
bt(9:16) = fread(fid,8,'uint8');

%-------------------------------------------------------------------------------
function i = isbb(fid)
%ISBB True if broad-band ADCP.

% check header and data source identification bytes
hid = 127;
sid = 127;
id = fread(fid,2,'uint8');
if length(id)<2
 err('ISBB: ****** can not read file id *****')
else
 i = id(1) == hid & id(2) == sid;
 % if i, disp('ISBB: BB-data '), end
end

% rewind file
fseek(fid,0,'bof');

%-------------------------------------------------------------------------------


function varargout = nbread(fid)
%NBREAD

% rewind to beginning of file and read header to get the number of bytes
% in each ensemble, the number of data types, and the address offsets
status = fseek(fid,0,'bof');
[nbytes,dtype,offset] = nbhead(fid);
ntypes = length(dtype);

% get the number of ensembles from file size; each ensemble has nbytes
% plus two bytes for the checksum
status = fseek(fid,0,'eof');
ftell(fid);
m = floor(ftell(fid)/(nbytes + 2));
status = fseek(fid,0,'bof');

% number of bins
n = (offset(4) - offset(3))/6;

% data parameters
varid = [1:7];
scale = [NaN,NaN,0.0025,1,1,1,1];
precision = {'','','bit12','uint8','uint8','uint8','bit4'};
bad = scale.*[NaN,NaN,NaN,NaN,NaN,NaN,NaN];

% initialize output variables
for k = 1:length(varid)
  if k == 1
    varargout{k} = NaN*ones(1,1,14);
  elseif k == 2
    varargout{k} = NaN*ones(m,1,13);
  else
    varargout{k} = NaN*ones(m,n,4);
  end
end

% read fixed leader
status = fseek(fid,offset(1),'bof');
varargout{1} = nbflead(fid);

for i = 1:m

  % read ensemble to verify the checksum
  status = fseek(fid,(i-1)*(nbytes+2),'bof');
  buffer = fread(fid,nbytes,'uint8');
  checksum = fread(fid,1,'uint16');

  % read ensemble if checksum is ok
  if checksum == rem(sum(buffer),65536)
    for k = dtype(2:end)
      % set file pointer to beginning of data
      status = fseek(fid,(i-1)*(nbytes+2)+offset(k),'bof');
      switch varid(k)
        case varid(2)
          % variable leader data
          varargout{k}(i,1,:) = nbvlead(fid);
        otherwise
          % velocity, spectral width, amplitude, percent-good, or status data
          a = fread(fid,4*n,precision{k});
          varargout{k}(i,:,:) = scale(k)*reshape(a,4,n)';
      end
    end
  end

end


% scale pitch, roll, and heading
varargout{2}(:,1,2:4) = varargout{2}(:,1,2:4)*360/65536;

% scale temperature
varargout{2}(:,1,5) = 45 - varargout{2}(:,1,5)*50/4096;


%-------------------------------------------------------------------------------

function [nbytes,dtype,offset] = nbhead(fid)
%NBHEAD Read header data from raw narrow-band ADCP data file.
%  [NBYTES,DTYPE,OFFSET] = nbhead(FID)

h = fread(fid,7,'uint16')';

% number of bytes
nbytes = h(1);

% address offsets and data types
varid = [1:7];
offset = [14 h(2:end)];
dtype = varid(offset ~= 0);
offset = [14 cumsum(offset(1:end-1))];


%-------------------------------------------------------------------------------

function fl = nbflead(fid)
%NBFLEAD Read fixed leader data from raw narrow-band ADCP data file.

fl = zeros(1,6);
f = cos(20*pi/180)/cos(30*pi/180);
fseek(fid,8,'cof');

% pings per ensemble
fl(2) = fread(fid,1,'uint16');

% bins per ping
fl(1) = fread(fid,1,'uint8');

% bin length
fl(3) = fread(fid,1,'uint8');
fl(3) = f*2^fl(3);

% transmit interval
fl(6) = fread(fid,1,'uint8');

% delay
fl(4) = f*fread(fid,1,'uint8');

% bin 1 distance
fl(5) = fl(4) + fl(3)/2;

% attenuation
fl(7) = 0.039;

% source level;
fl(8) = 100;

% serial number;
fl(9:14)=[3 4 5 6 7 8];

%-------------------------------------------------------------------------------

function vl = nbvlead(fid)
%NBVLEAD Read variable leader data from raw narrow-band ADCP data file.

vl = zeros(1,13);

% time of ensemble (mm/dd hh:mm:ss)
a = fread(fid,5,'uint8');
a = dec2hex(a);
c(1)=1900;
c(2:6) = str2num(a);
vl(1) = julian(c(1),c(2),c(3),c(4)+c(5)/60+c(6)/3600);

% pitch, roll, heading, and temperature
fseek(fid,16,'cof');
vl(2:5) = fread(fid,4,'uint16');
vl(2:3) = vl(2:3) - floor(vl(2:3)/(183*180))*360*182;

vl(6) = 35;
vl(7) = 1536;
vl(8:10)=[nan nan nan];
%-------------------------------------------------------------------------------

function bt = nbbtrack(fid)
%RDBTRACK Read the bottom track data from a raw ADCP data file.
%  BT = NBBTRACK(FID)

fseek(fid,14,'cof');

% range
%bt(1:4) = 0.01*fread(fid,4,'uint16');

% velocity
%bt(5:8) = 0.001*fread(fid,4,'int16');

% correlation magnitude and percent good
%bt(9:16) = fread(fid,8,'uint8');

% not implemented
bt(1:16)=NaN;

%-------------------------------------------------------------------------------


function d=y2k(d)
% fix date string
if d<80, d=2000+d; end
if d<100, d=1900+d; end


%-------------------------------------------------------------------------------

function [vele]=b2earth(velb,v,a,p,ignore_beam)
% 
% convert beam ADCP data to earth velocities
%
% input velb:  beam coordinates
%          v:  attitude vector
%          a:  ADCP information
%	   p:  global p structure
% 	ignore_beam: nan or beam number to ignore
% 
% output vele: earth coordinates
%
% hard wired for LADCP systems
% M. Visbeck  Jan 2004

if a.Coordinates~=0
 disp('Data are not in beam coordinates!')
 vele=velb;
 return
end

a=setdefv(a,'use_tilt',1);
a=setdefv(a,'use_heading',1);
a=setdefv(a,'use_binremap',0);			%%% CODE IS BUGGY! DO NOT USE
a=setdefv(a,'beams_up',a.Up);
a=setdefv(a,'beamangle',a.Beam_angle);

p=setdefv(p,'allow_3beam_solutions',1);

a.sensor_config=1;
a.convex=1;

N_3beam = 0;
N_4beam = 0;

% Written by Marinna Martini for the 
% U.S. Geological Survey
% Branch of Atlantic Marine Geology
% Thanks to Al Pluddeman at WHOI for helping to identify the 
% tougher bugs in developing this algorithm


% precompute some constants
d2r=pi/180; % conversion from degrees to radians
C30=cos(a.beamangle*d2r);
S30=sin(a.beamangle*d2r);

if a.beams_up == 1, % for upward looking
  ZSG = [+1, -1, +1, -1];
else % for downward looking
  ZSG = [+1, -1, -1, +1];
end

% size of problem
nb=size(velb,2);
ne=size(velb,1);

vele=velb*nan;

%big loop over profiles
for ii=1:ne

roll=v(ii,1,3);
pitch=v(ii,1,2);
head=v(ii,1,4);
beam=squeeze(velb(ii,:,:));

% Step 1 - determine rotation angles from sensor readings
% fixed sensor case
% make sure everything is expressed in radians for MATLAB
RR=roll.*d2r;
KA=sqrt(1.0 - (sin(pitch.*d2r).*sin(roll.*d2r)).^2);
PP=asin(sin(pitch.*d2r).*cos(roll.*d2r)./KA);
%% NB: The preceding two lines could be replaced with
%% 		PP=atan(tan(pitch.*d2r) * cos(roll.*d2r));
%%     which is the expression given by RDI in the coord-
%%     trans manual. I have tried this with a single
%%     file from DIMES UK2 and the max velocity differences
%%     are 1e-13 m/s, i.e. they look consistent with
%%     roundoff errors.

HH=head.*d2r;

% Step 2 - calculate trig functions and scaling factors
if a.use_tilt
 CP=cos(PP); CR=cos(RR); 
 SP=sin(PP); SR=sin(RR); 
else
 CP=1; CR=1;
 SP=0; SR=0;
end

if a.use_heading
 CH=cos(HH); SH=sin(HH);
else
 CH=1; SH=0; 
end

% fixed sensor case
M(1) = -SR.*CP;
M(2) = SP;
M(3) = CP.*CR;

% compute scale factor for each beam to transform depths
% in a tilted frame to depths in a fixed frame
SC(1) = (M(3).*C30 + ZSG(1).*M(1).*S30);
SC(2) = (M(3).*C30 + ZSG(2).*M(1).*S30);
SC(3) = (M(3).*C30 + ZSG(3).*M(2).*S30);
SC(4) = (M(3).*C30 + ZSG(4).*M(2).*S30);

SSCOR = 1;
% my version of Al's scaling constant, using RDI's
% convention for theta as beam angle from the vertical
VXS = SSCOR/(2.0*S30);
VYS = VXS;
VZS = SSCOR/(4.0*C30);
VES = VZS;

[NBINS, n]=size(beam);
earth=zeros(size(beam));
clear n;
J=zeros(1,4);

for IB=1:NBINS,
 % Step 3:  correct depth cell index for pitch and roll
 for i=1:4, 
  if a.use_binremap
%   J(i)=fix(IB.*SC(i)+0.5);
   J(i)=IB; %%%
  else
   J(i)=IB;
  end
 end

 % Step 4:  ADCP coordinate velocity components
 if all(J > 0) & all(J <= NBINS),

  if ~isnan(ignore_beam)
    beam(:,ignore_beam) = nan;
  end

  this_3beam = 0;
  if (p.allow_3beam_solutions) && ...
      (isnan(beam(J(1),1)) + isnan(beam(J(2),2)) + ...
       isnan(beam(J(3),3)) + isnan(beam(J(4),4)) == 1)
    N_3beam = N_3beam + 1;
    this_3beam = 1;
    if isnan(beam(J(1),1))
      beam(J(1),1) = -beam(J(2),2) + beam(J(3),3) + beam(J(4),4);
    elseif isnan(beam(J(2),2))
      beam(J(2),2) = -beam(J(1),1) + beam(J(3),3) + beam(J(4),4);
    elseif isnan(beam(J(3),3))
      beam(J(3),3) = beam(J(1),1) + beam(J(2),2) - beam(J(4),4);
    else
      beam(J(4),4) = beam(J(1),1) + beam(J(2),2) - beam(J(3),3);
    end
  elseif isnan(beam(J(1),1)) + isnan(beam(J(2),2)) + ...	
         isnan(beam(J(3),3)) + isnan(beam(J(4),4)) == 0
    N_4beam = N_4beam + 1;
  end
  
  if isnan(beam(J(1),1)) || isnan(beam(J(2),2)) || ...
     isnan(beam(J(3),3)) || isnan(beam(J(4),4)),
    earth(IB,:)=ones(size(beam(IB,:))).*NaN;
  else
    if a.beams_up ,
     % for upward looking convex
     VX = VXS.*(-beam(J(1),1)+beam(J(2),2));
     VY = VYS.*(-beam(J(3),3)+beam(J(4),4));
     VZ = VZS.*(-beam(J(1),1)-beam(J(2),2)-beam(J(3),3)-beam(J(4),4));
     VE = VES.*(+beam(J(1),1)+beam(J(2),2)-beam(J(3),3)-beam(J(4),4));
    else
     % for downward looking convex
     VX = VXS.*(+beam(J(1),1)-beam(J(2),2));
     VY = VYS.*(-beam(J(3),3)+beam(J(4),4));
     VZ = VZS.*(+beam(J(1),1)+beam(J(2),2)+beam(J(3),3)+beam(J(4),4));
     VE = VES.*(+beam(J(1),1)+beam(J(2),2)-beam(J(3),3)-beam(J(4),4));
    end
    if this_3beam && abs(VE) > 1e-9
     error('3-beam code assertion failed');
    end

   % Step 5: convert to earth coodinates
   VXE =  VX.*(CH*CR + SH*SR*SP) + VY.*SH.*CP + VZ.*(CH*SR - SH*CR*SP);
   VYE = -VX.*(SH*CR - CH*SR*SP) + VY.*CH.*CP - VZ.*(SH*SR + CH*SP*CR);
   VZE = -VX.*(SR*CP)            + VY.*SP     + VZ.*(CP*CR);
   earth(IB,:) = [VXE, VYE, VZE, VE];
  end % end of if any(isnan(beam(IB,:))),
 else
  earth(IB,:)=ones(size(beam(IB,:))).*NaN;
 end % end of if all(J > 0) && all(J < NBINS),
end % end of IB = 1:NBINS

% save results
vele(ii,:,:)=earth;

end % Big Loop

if p.allow_3beam_solutions
  disp(sprintf(' %d 3-beam solutions calculated (%d%% of total)',...
	N_3beam,round(100*N_3beam/(N_3beam+N_4beam))));
end
