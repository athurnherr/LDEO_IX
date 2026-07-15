%======================================================================
%                    G E T D P T H I . M 
%                    doc: Wed Jan  7 16:25:26 2009
%                    dlm: Fri Nov 19 12:24:59 2021
%                    (c) 2009 A.M. Thurnherr
%                    uE-Info: 13 0 NIL 0 0 72 0 2 4 NIL ofnI
%======================================================================

% CHANGES BY ANT:
%   Jan  7, 2009: - tightened use of exist()
%	Aug 30, 2019: - BUG: missing pressure values cause problem in output
%	Nov 19, 2021: - added error message

function [d,p]=getdpthi(d,p)
% function [d,p]=getdpthi(d,p)
% LADCP-2 processing software v 7.0
%
%  -- make depth from raw data 
%  use inverse approach
% Matin Visbeck 
%  December 2002, LDEO

disp('GETDPTHI: Depth from vertical velocity inverse method')

% set default start, deepst  and end depth
p=setdefv(p,'zpar',[10 NaN 10]);
p=setdefv(p,'cut',15);
p=setdefv(p,'dzbelow',[2 -1]*medianan(abs(diff(d.zd))));
p=setdefv(p,'guessbottom',NaN);
p=setdefv(p,'zbottom',NaN);
p=setdefv(p,'ctddepth',0);
p=setdefv(p,'navtime_av',2/60/24);
p.zpar=abs(p.zpar);

% remove empty profiles at beginning and end
ii=ones(size(d.rw(1,:)));
jj=find(~isnan(meannan(d.rw)));
ii(1:jj(1))=0;
ii(jj(end):end)=0;
d=cutstruct(d,ii);

% get time difference for w-values
dt=diff(d.time_jul)*24*3600;
dt=mean([dt([1,1:end]);dt([1:end,end])]);

% initialize a matrix where below bottom and above surface data
% get flagged as bad
d.izmflag = d.rw*0;

figure(4)
orient tall

% two sweeps throught the data
n2=2;

% only one sweep if ctd pressure data exist
if p.ctddepth>0, n2=1; end

for n=1:n2
  
  disp([' starting run ',int2str(n),' to get LADCP depth'])
  if n==1 
%  if n==1 | (~exist('ddoz') & ~exist('ddoz'))
%  first sweep used intergal of median w 
   ibad=find(p.wizr<1 | p.wizr>size(d.rw,1));
   p.wizr(ibad)=[];
   if length(p.wizr)>1
    dw=medianan(d.rw(p.wizr,:),1);
   elseif length(p.wizr)==1
    dw=d.rw(p.wizr,:);
   else
    dw=medianan(d.rw,1);
   end
%  set non finite w to zero
   ii=find(~isfinite(dw));
   dw(ii)=0;
%  integrate result
   zz=cumsum(dw.*dt);
%  make sure that start and end depth are as wanted
   zz=zz-linspace(-maxnan([0 p.zpar(1)]),-maxnan([0 p.zpar(3)])+zz(end),length(zz));
%  set maxdepth to be as requested
   if isfinite(p.zpar(2))
    zz=zz/max(zz)*p.zpar(2);
   end

  else

%  second sweep used inverse method to get z including bottom and
%  surface track data
%  "erase" subplots
  subplot(321), hold off
  subplot(322), hold off
  subplot(312), hold off

% set base matrix
   dw=d.rw(p.wizr,:)+d.izmflag(p.wizr,:);
   [d1,A1,ibot]=dinset(dw,dt);

   % set boundary conditions for inversion
   if isfinite(p.zpar(1))
     d1=[d1;p.zpar(1)];
   else 
     d1=[d1;10];
   end
   A1(length(d1),1)=1;

   if isfinite(p.zpar(2))
     d1=[d1;p.zpar(2)*10];
     A1(length(d1),ibot)=10;
   end

   if isfinite(p.zpar(3))
     d1=[d1;p.zpar(3)];
   else
     d1=[d1;10];
   end
   A1(length(d1),end)=1;

   % add surface/bottom reflections if present

   if exist('ddoz','var')
    [ld,lz]=size(A1);
    d1=[d1;ddoz*0.1];
    A2=sparse(1:length(ddoi),ddoi,0.1);
    A2(1,lz)=0;
    A1=[A1;A2];
   end
 
   if exist('dupz','var')
    [ld,lz]=size(A1);
    d1=[d1;dupz*0.1];
    A2=sparse(1:length(dupi),dupi,0.1);
    A2(1,lz)=0;
    A1=[A1;A2];
   end

   if exist('dbotdz','var') 
    bfac=0.1;
    [ld,lz]=size(A1);
    d1=[d1;dbotdz*bfac];
    ix=[1:length(dboti)]';
    A2=sparse([ix;ix],[dboti;ix*0+ibot],[ix*0-bfac;ix*0+bfac]);
    A2(1,lz)=0;
    A1=[A1;A2];
   end

% require z to be smooth (needed otherwise ill constrained)
   [A1,d1]=dismoo(A1,d1,0.01);

% solve for best depth time series
   zz=lesqchol(d1,A1)';

  end

  % set depth below which to delete all data 
  dzbelow=p.dzbelow(n);


  if n==1
  % get sound speed time series
   if ~existf(d,'ss')
    if existf(d,'ctd_ss')
      disp(' take soundspeed from CTD time series')
      d.ss=d.ctd_ss;
    elseif existf(d,'ctdprof_ss')
     disp(' take soundspeed from CTD profile ')
     zctd=d.ctdprof_z;
     zctd(1)=-1e5;
     zctd(end)=1e5;
     d.ss=interp1q(zctd,d.ctdprof_ss,zz')';
    else
      if existf(d,'ctd_temp')
        disp(' make soundspeed based on CTD pressure and temp')
        pp=press(abs(d.z));
        d.ss=sounds(pp,d.ctd_temp,34.5);
      else
        disp(' make soundspeed based on pressure and ADCP temp')
        pp=press(abs(zz));
        d.ss=sounds(pp,d.temp(1,:),34.5);
      end
    end
   end
   % sound speed correction
   if d.soundc==0
    disp(' correct velocities for sound speed ')
    sc=meshgrid(d.ss./d.sv(1,:),d.izd);
    d.ru(d.izd,:)=d.ru(d.izd,:).*sc;
    d.rv(d.izd,:)=d.rv(d.izd,:).*sc;
    d.rw(d.izd,:)=d.rw(d.izd,:).*sc;
    if existf(d,'hbot')
     d.hbot=d.hbot.*sc(1,:);
     d.bvel(:,1:3)=d.bvel(:,1:3).*sc(1:3,:)';
     if existf(d,'bvel_rdi')
      d.bvel_rdi(:,1:3)=d.bvel_rdi(:,1:3).*sc(1:3,:)';
     end
     if existf(d,'bvel_own')
      d.bvel_own(:,1:3)=d.bvel_own(:,1:3).*sc(1:3,:)';
     end
    end
    if length(d.zd)~=length(d.ru(:,1))
     sc=meshgrid(d.ss./d.sv(2,:),d.izu);
     d.ru(d.izu,:)=d.ru(d.izu,:).*sc;
     d.rv(d.izu,:)=d.rv(d.izu,:).*sc;
     d.rw(d.izu,:)=d.rw(d.izu,:).*sc;
    end
    if existf(d,'hsurf')
     d.hsurf=d.hsurf.*sc(1,:);
    end
    d.soundc=1;
   else
    disp(' will not correct for sound speed twice')
   end

  end


  if p.ctddepth==0
   % save results only if CTD-depth was not available
   d.z=-zz;
   p.ladcpdepth=2;
   disp(' use LADCP depth from constrained integrated W ')
  else
   d.z_ladcp=-zz;
   dz=d.z_ladcp-d.z;
   ii=find(isfinite(dz));
   p.ladcpr_CTD_depth_std=[mean(dz(ii)), std(dz(ii))];
   disp(' use CTD time series depth, will not do depth inversion ')
   disp([' LADCP minus CTD depth mean: ',num2str(p.ladcpr_CTD_depth_std(1)),...
          '  std: ',num2str(p.ladcpr_CTD_depth_std(2))]);
   if ~isfinite(p.ladcpr_CTD_depth_std(1))
	error('non-numeric result, try reprocessing with p.getdepth = 1');
   end
   p.ladcpdepth=0;
  end

  [p.maxdepth,ibottom]=max(-d.z);

  % plot near surface LADCP data 
  [dum,ibot]=max(-d.z);
  % for very shallow stations turn of surface detection
  if dum<100, p.surfdist=0; disp(' shallow station no surface detection '),end
  % first down cats
  ii=1:ibot;
  iok=ii(find(d.z(ii)>-200 & d.z(ii)<-30 ));
  if length(iok)>2
   iok2=1:iok(end);
   subplot(321)
   plot(iok2,d.z(iok2),'.r',iok2,iok2*0,'-k'),
   hold on
   title('start depth (.red) (.blue) Surface distance')
   xlabel('time in ensembles')
   ylabel('depth in meter')
   axis tight
   ax=axis;
   ax(4)=max(ax(4),20);
   axis(ax);
  else
   disp('no surface reflections down-cast')
  end

  % then up cast
  ii=ibot:length(d.z);
  iok=ii(find(d.z(ii)>-200 & d.z(ii)<-30 ));
  if length(iok)>2
   iok2=iok(1):length(d.z);
   subplot(322)
   plot(iok2,d.z(iok2),'.r',iok2,iok2*0,'-k'),
   hold on
   title('end depth')
   xlabel('time in ensembles')
   ylabel('depth in meter')
   axis tight
   ax=axis;
   ax(4)=max(ax(4),20);
   axis(ax);
  else
   disp('no surface reflections up-cast')
  end

  % extract surface distance to find start depth
  if existf(d,'hsurf') & p.surfdist
   if sum(isfinite(d.hsurf))>10
    [zmax,ibot]=max(-d.z);
   % start depth
    ii=find(isfinite(d.hsurf(1:ibot)));
   % limit to surface detections within 30m of answer
    iok=ii(find(d.z(ii)>-200 & d.z(ii)<-40 & abs(d.hsurf(ii)+d.z(ii))<30));
    if ~isempty(iok)
     ddoz=d.hsurf(iok)';
     ddoi=iok';
    % temporary plot of surface detection
      iok2=1:iok(end);
      subplot(321)
      plot(iok,d.hsurf(iok)+d.z(iok),'.',iok2,d.z(iok2),'.r'),
      ax=axis;
      ax(4)=max(ax(4),20);
      axis(ax);
      axis tight
    else
     disp(' no surface reflections down-cast')
    end

   % surface distance to find end depth
    ii=find(isfinite(d.hsurf(ibot:end)))+ibot-1;
   % limit to surface detections within 30m of answer
    iok=ii(find(d.z(ii)>-200 & d.z(ii)<-40 & abs(d.hsurf(ii)+d.z(ii))<30));
    if ~isempty(iok)
     dupz=d.hsurf(iok)';
     dupi=iok';
    % temporary plot of surface detection
      iok2=iok(1):length(d.z);
      subplot(322)
      plot(iok,d.hsurf(iok)+d.z(iok),'.',iok2,d.z(iok2),'.r'),
      ax=axis;
      ax(4)=max(ax(4),20);
      axis(ax);
      axis tight
    else
     disp(' no surface reflections up-cast')
    end
   end 
  end 

  disp([' maximum depth from int W is :',int2str(max(zz))])
  disp([' should be                   :',int2str(p.zpar(2))])

  % bottom depth
  subplot(312)
  iok=find((max(-d.z)+d.z)<200);
  iok2=iok(1):iok(end);
  plot(iok2,d.z(iok2),'.r'),
  if ~isfinite(p.zbottom)
   if sum(isfinite(d.hbot))>10
    % look for bottom only close to deepest CTD depth
    iok=find((max(-d.z)+d.z)<200 & d.hbot>0 );
    if length(iok)>2
      % fit polynom to bottom depth time series
      if n==n2
      % use deepest point to set bottom depth
       if isfinite(p.guessbottom)
	  zbottom=p.guessbottom;
       else
          zbottom= median(d.hbot(iok)-d.z(iok));
       end
       zbottomerr= zbottom-(d.hbot(iok)-d.z(iok)) ;
       [dum,is]=sort(abs(zbottomerr));
       is=is(1:fix(length(is)/2));
       c=polyfit(iok(is),d.hbot(iok(is))-d.z(iok(is)),1);
       zbottomerr= polyval(c,iok)-(d.hbot(iok)-d.z(iok)) ;
       [dum,is]=sort(abs(zbottomerr));
       is=is(1:fix(length(is)/2));
       iok=iok(find(abs(zbottomerr)<2*std(zbottomerr(is)) | abs(zbottomerr)<30 ));
       c=polyfit(iok,d.hbot(iok)-d.z(iok),2);
       zbottomerr= polyval(c,iok)-(d.hbot(iok)-d.z(iok)) ;
       zbottom=polyval(c,ibottom);
       p.zbottom=zbottom;
      else
       if isfinite(p.guessbottom)
	  zbottom=p.guessbottom;
       else
          zbottom= median(-d.hbot(iok)-d.z(iok));
       end
       zbottomerr= zbottom-(d.hbot(iok)-d.z(iok)) ;
       [dum,is]=sort(abs(zbottomerr));
       is=is(1:fix(length(is)/2));
       zbottom=medianan(d.hbot(iok(is))-d.z(iok(is)));
       zbottomerr= zbottom-(d.hbot(iok)-d.z(iok)) ;
       [dum,is]=sort(abs(zbottomerr));
       is=is(1:fix(length(is)/2));
       iok=iok(find(abs(zbottomerr)<2*std(zbottomerr(is)) | abs(zbottomerr)<50 ));
       c=polyfit(iok,d.hbot(iok)-d.z(iok),1);
       % save bottom distances for inversion
       dbotdz=(d.hbot(iok)-polyval(c,iok)-d.z(ibottom))';
       dboti=iok';
      end
      p.zbottomerror = medianan(abs(zbottomerr)); 
    % temporary plot of bottom detection
      iok2=iok(1):iok(end);
      plot(iok,-d.hbot(iok)+d.z(iok),'.',iok2,d.z(iok2),'.r'),
      hold on
      plot(iok2,-zz(iok2),'-k')
      hold on, plot(iok2,iok2*0-p.zbottom,'--k') 
      hold on, plot(iok2,-polyval(c,iok2),'-b')
      title('bottom (--k) LADCP depth (-k) bottom distance (.b)')
      xlabel('time in ensembles')
      ylabel('depth in meter')
      if n==n2
 % remove suspicious bottom track data
       axis tight
       ax=axis;
       text(ax(1)+abs(diff(ax(1:2)))*.15,ax(3)+abs(diff(ax(3:4)))*.8,...
       ['bottom at: ',int2str(p.zbottom),' [m]   ADCP was ',...
         int2str(p.zbottom-max(-d.z)),' m above bottom'])
       ibad=1:length(d.hbot);
       % good data are
       ibad(iok)=[];
       d.hbot(ibad)=NaN;
       d.bvel(ibad,:)=NaN;
      else
      end
    else
      zbottom=NaN;
      p.zbottomerror = nan;
    end
   else
    zbottom=NaN;
    p.zbottomerror=NaN;
   end
  else
    disp(['  bottom preset at ',int2str(p.zbottom)])
    zbottom=p.zbottom;
    p.zbottomerror=0;
  end
   % check if bottom is shallower that maxctd-depth an
    if ((zbottom-p.maxdepth<-(p.maxdepth*0.01+10) & isfinite(zbottom)) |...
       	p.zbottomerror > 20 )
      disp('  no bottom found')
      disp(['   given maximum profile depth : ',int2str(p.maxdepth)])
      disp(['   extracted bottom depth      : ',int2str(zbottom)])
      disp(['        bottom depth error     : ',int2str(p.zbottomerror)])
      p.zbottom=NaN;
    elseif n==n2
      p.zbottom=zbottom;
      disp(['  bottom found at ',int2str(p.zbottom),' +/- ',...
                                 int2str(p.zbottomerror),' m'])
      if (p.zbottom<p.maxdepth)
        disp('  extracted bottom within 20m above given maximum profile depth')
      end
    end
  pause(0.1)

  % assign a depth to each bin
  [izm1,izm]=meshgrid([fliplr(d.zu),-d.zd],d.z);
  izm1=izm1';
  if d.soundc==1
  % make sound speed correction for depth vector
   sc=meshgrid(d.ss./d.sv(1,:),d.izd);
   izm1(d.izd,:)=izm1(d.izd,:).*sc;
   if length(d.zu)>0
    sc=meshgrid(d.ss./d.sv(2,:),d.izu);
    izm1(d.izu,:)=izm1(d.izu,:).*sc;
   end
   disp(' correct bin length for sound speed')
  end
  % add the two parts
  d.izm=izm'+izm1;

  % flag all data below bottom as bad
  if ~isnan(p.zbottom)
    ii = find(d.izm<-p.zbottom-dzbelow);
    d.izmflag(ii)=NaN; 
  end

  % flag all data close to the surface as bad
  if length(d.zu)>0
   ii = find(d.izm>-(d.zu(2)-d.zu(1))/2);
   d.izmflag(ii)=NaN; 
  end
  
end 


% set velocities deeper than bottom to NaN
bad = find( isnan(d.izmflag) & isfinite(d.ru) );
if ~isempty(bad) & isfinite(p.zbottom)
  disp([' removing ',int2str(length(bad)),...
	' values below bottom'])

 %d.ru = d.ru+d.izmflag;
 %d.rv = d.rv+d.izmflag;
 %d.rw = d.rw+d.izmflag;
 d.weight = d.weight + d.izmflag;
end

% compute pressure from depth
d.p=press(abs(d.z));

if existf(d,'wctd')==1
 d.wm=d.wctd;
else
 d.wm=-gradient(d.z)./gradient(d.time_jul*24*3600);
end

if d.z(1)<-50, 
 warn=[' first LADCP depth is ',int2str(d.z(1))];
 disp(warn)
 p.warnp(size(p.warnp,1)+1,1:length(warn))=warn;
end

if d.z(end)<-50, 
 warn=[' last LADCP depth is ',int2str(d.z(end))];
 disp(warn)
 p.warnp(size(p.warnp,1)+1,1:length(warn))=warn;
end

% cut raw data to only include profile
i1=find(d.z(1:ibottom)>-p.cut);
i2=find(d.z(ibottom:end)>-p.cut)+ibottom-1;
if length(i1)==0
 i1=1;
end
if length(i2)==0
 i2=length(d.z);
end
ii=d.z*0;
ic=i1(end):i2(1);
ii(ic)=1;
if (sum(ii)~=length(ii)) & p.cut>0 
 disp(' remove data at begining and end of cast')
 disp(' adjust start and end time ')
 d=cutstruct(d,ii);
 p.zpar([1 3])=p.cut;
 subplot(321), ax=axis; plot([1 1]*ic(1),ax(3:4),'--k')
 subplot(322), ax=axis; plot([1 1]*ic(end),ax(3:4),'--k')
 p.time_start=gregoria(d.time_jul(1));
 p.time_end=gregoria(d.time_jul(end));
 if existf(d,'slon')
  % average over first p.navtime_av days
  ii=find(d.time_jul<(d.time_jul(1)+p.navtime_av));
  slon=median(d.slon(ii));
  slat=median(d.slat(ii));
  if isfinite(slon+slat)
   p.poss=[fix(slat), (slat-fix(slat))*60, fix(slon), (slon-fix(slon))*60];
  end
  ii=find(d.time_jul>(d.time_jul(end)-p.navtime_av));
  elon=median(d.slon(ii));
  elat=median(d.slat(ii));
  if isfinite(elon+elat)
   p.pose=[fix(elat), (elat-fix(elat))*60, fix(elon), (elon-fix(elon))*60];
  end
  disp(' find new start end position')
 elseif existf(p,'poss')
  % don't have time variable navigation data, adjust end position
  % to account for shorter cast time
  slat=p.poss(1)+p.poss(2)/60; 
  slon=p.poss(3)+p.poss(4)/60; 
  elat=p.pose(1)+p.pose(2)/60; 
  elon=p.pose(3)+p.pose(4)/60; 
  elon=slon+(elon-slon)*sum(ii)/length(ii);
  elat=slat+(elat-slat)*sum(ii)/length(ii);
  disp(' adjust end position for shorter cast time')
  if isfinite(elon+elat)
   p.pose=[fix(elat), (elat-fix(elat))*60, fix(elon), (elon-fix(elon))*60];
  end
 end
	  
 % save start and end depth
 p.zpar(1)=-d.z(1);
 p.zpar(3)=-d.z(end);
end

streamer([p.name,'   Figure 4']);
pause(0.01)

%-----------------------------------------------------------------
function [d,A,izbot]=dinset(dw,dt)
% function [d,A]=dinset(dw,dt)
% set up sparse Matrix for depth inversion
%
[nb,nt]=size(dw);

% find bottom roughly and devide cast in down and up trace
if nb>1
 wm=medianan(dw);
else
 wm=dw;
end
ii=find(~isfinite(wm));
wm(ii)=0;
zz=cumsum(wm.*dt);

[zbot,izbot]=maxnan(zz);
disp([' bottom:',int2str(zbot),' @ ',int2str(izbot)])

ido=1:izbot;
iup=(izbot+1):length(dt);

% 
dtm=repmat(dt,[nb 1]);
izm=repmat(1:nt,[nb 1]);

d=reshape(dw.*dtm,nb*nt,1);
izv=reshape(izm,nb*nt,1);

ibad=find((izv-1)<1 | (izv+1)>nt);
d(ibad)=[];
izv(ibad)=[];

iweak=find(~isfinite(d) | (izv-1)<1 | (izv+1)>nt);

it=[1:length(d)]';
i1=it*0+0.5;
d(iweak)=0;
i1(iweak)=0.01;

A=sparse([it;it],[izv+1;izv-1],[i1;-i1]);
A(1,nt)=0;


%==============================================================
function a=cutstruct(a,ii)
% reduce array size in structure
lz=length(ii);
iok=find(ii==1);
if existf(a,'cutindx')
 a.cutindx=a.cutindx(1)-1+[iok(1) iok(end)];
else
 a.cutindx=[iok(1) iok(end)];
end
if isstruct(a)
  fnames = fieldnames(a);
  for n=1:size(fnames,1)
   dummy = getfield(a,fnames{n});
   [ly,lx]=size(dummy);
   if ly==lz
    a=setfield(a,fnames{n},dummy(iok,:));
   elseif lx==lz
    a=setfield(a,fnames{n},dummy(:,iok));
   end
  end
end


%-------------------------------------------------------------------
function [A,d]=dismoo(A,d,fs0,cur);
% function [A,d]=dismoo(A,d,fs0,cur);
%
% smooth results by minimizing curvature
% also smooth if elements are not constrained
%
if nargin<3, fs0=1; end
if nargin<4, cur=[-1 2 -1]; end

[ld,ls]=size(A);
fs=full(sum(abs(A)));
fsm=max(median(fs),0.1);
% find ill constrained data
ibad=find(fs<fsm*0.1);

% increase weight for poorly constrained data
fs=max(fs,fsm*0.1);
fs=fsm^2./fs * fs0(1);


if length(ibad)>0
% set ill constrainded data to a minimum weight
 fs(ibad)=max(fs(ibad),0.5);
 if fs0==0
  disp([' found ',int2str(length(ibad)),' ill constrained elements will smooth '])
 else
  disp([' found ',int2str(length(ibad)),' ill constrained elements'])
 end
end

if sum(fs>0)>0

 cur=cur-mean(cur);

 lc=length(cur);
 lc2=fix(lc/2);
 fs2=fs((lc2+1):(end-lc2));
 inc=[1:length(cur)]-lc2;

 ii=find(fs2>0);
 % find how many smooth constraints to apply

 if length(ii)>0 
  [i1,i2]=meshgrid(inc,ii+lc2-1);
  [curm,fsm]=meshgrid(cur,fs2(ii));

  As=sparse(i2,i1+i2,curm.*fsm);
  [lt,lm]=size(A);
  if size(As,2)<lm 
   As(1,lm)=0;
  end
  A=[A;As];

 end

% smooth start and end of vector
 for j=1:lc2
  j0=j-1;
  [lt,lm]=size(A);
  if fs(1+j0)>0
   A(lt+1,[1:2]+j0)=[2 -2]*fs(1+j0);
  end
  if fs(end-j0)>0
   A(lt+2,end-[1,0]-j0)=[-2 2]*fs(end-j0);
  end
 end

 [lt,lm]=size(A);
 d(lt)=0;
else
 disp(' no smoothness constraint applied ')
end

%-------------------------------------------------------------------
function [m,dm,c]=lesqchol(d,g)
% function [m,dm,c]=lesqcholw(d,g)

% fit least squares method to linear problem 
% Use Cholesky transform
% 
%input parameters:
%  d:= data vector ;  g:= model matrix 
% output parameters:
% m=model factors; dm= model data, c=correlation 

n=length(d);
[i,j]=size(g);
if i~=n; disp(' wrong arguments'),return,end
[r,b] = chol( g.' * g);
if b~=0, m=g(1,:)'+NaN; dm=d+NaN; c=NaN; return, end
y = forwardsub(r.' , g.' * d);
m  = backsub(r,y);
if nargout<2, return, end
dm = g * m;
if nargout<3, return, end
co = cov([d,dm]);
c  = co(1,2) / sqrt( co(1,1)*co(2,2) );
 
%-------------------------------------------------------------------
function  X = backsub(A,B)

% X = BACKSUB(A,B)  Solves AX=B where A is upper triangular.
% A is an nxn upper-triangular matrix (input)
% B is an nxp matrix (input)
% X is an nxp matrix (output)

[n,p] = size(B);
X = zeros(n,p);

X(n,:) = B(n,:)/A(n,n);

for i = n-1:-1:1,
  X(i,:) = (B(i,:) - A(i,i+1:n)*X(i+1:n,:))/A(i,i);
end
%-------------------------------------------------------------------
function  X = forwardsub(A,B)

% X = FORWARDSUB(A,B))  Solves AX=B where A is lower triangular.
% A is an nxn lower-triangular matrix, input.
% B is an nxp matrix (input)
% X is an nxp matrix (output)

[n,p] = size(B);
X = zeros(n,p);

X(1,:) = B(1,:)/A(1,1);

for i = 2:n,
  X(i,:) = (B(i,:) - A(i,1:i-1)*X(1:i-1,:))/A(i,i);
end

