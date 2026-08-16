function p=plotraw(d,p)
% function p=plotraw(d,p)
% plot some results
%

%======================================================================
%                    P L O T R A W . M 
%                    doc: Fri Jan  5 15:38:43 2007
%                    dlm: Wed Apr 15 11:30:50 2015
%                    (c) 2007 M. Visbeck with contribs from A. Thurnherr
%                    uE-Info: 19 76 NIL 0 0 72 0 2 4 NIL ofnI
%======================================================================

% MODIFICATIONS BY ANT:
%	Jan  5, 2007: - fixed checkbeam() as suggested by B. Huber
%	Jun  3, 2013: - BUG: top panel of Fig. 2 was wrong for dual-headed
%						 LADCPs with different UL/DL bin sizes
%	Apr 15, 2015: - BUG: ax() is always used but only sometimes defined()
%						 reported by Achim Randelhoff (yahoo email 04/23/14)

pmax=200;
orient tall

if length(d.time_jul)>pmax
 ii=fix(linspace(2,length(d.time_jul)-2,pmax));
else
 ii=1:length(d.time_jul);
end

clf
subplot(411)
z=-d.izm(:,1)+d.z(1);
zz=[];
rw=[];

% find three beam solutions
% a profile is considered 3-beam is half of the good data are 3-beam
iz=fliplr(d.izu);
n3bu=0;
n3bd=0;
if length(iz)>1
 zz=[zz;-z(iz)];
 rw=[rw;d.rw(iz,ii)+d.weight(iz,ii)*0];
 % check for 3-beam solutions
 iw=sum(~isnan(d.rw(iz,ii)));
 ie=sum(~isnan(d.re(iz,ii)));
 i3bu=find(iw>(2*ie));
 n3bu=length(i3bu)/length(ii)*100;
 rw=[rw;rw(1,:)*NaN];
 zz=[zz;0];
end
iz=d.izd;
if length(iz)>1
 zz=[zz;-z(iz)];
 rw=[rw;d.rw(iz,ii)+d.weight(iz,ii)*0];
 % check for 3-beam solutions
 iw=sum(~isnan(d.rw(iz,ii)));
 ie=sum(~isnan(d.re(iz,ii)));
 i3bd=find(iw>(2*ie));
 n3bd=length(i3bd)/length(ii)*100;
end

% plot vertical velocities
col=jet(128);
col=([[1 1 1]; col]);
colormap(col)

pcolorn(ii,zz,rw) % this is neccessary even for dual-headed data sets

if (length(d.izu) > 0)
	pcolorn(ii,zz(1:length(d.izu)),rw(1:length(d.izu),:))
	hold on, shading flat
	ax = axis;
	ax(1) = ii(1);
	ax(2) = ii(end);
	ax(4) = -zz(end);
	axis(ax);
	pcolorn(ii,zz(length(d.izu)+1:end),rw(length(d.izu)+1:end,:))
else
	ax = axis;	% alternative solution to bug reported by Achim Randelhoff yahoo email 04/23/2014
end

hold on
plot([ii(1),ii(end)],[0 0],'-k')
colorbar('horiz')

% mark 3-beam solutions
if n3bu>10
 l3b=NaN*ii;
 l3b(i3bu)=ax(3);
 plot(ii,l3b,'-r','linewidth',8)
 text(mean(ax(1:2)),ax(4),['found ',int2str(n3bu),'% profiles 3 beam solutions'],...
 'VerticalAlignment','top','HorizontalAlignment','center','fontsize',14)
end

if n3bd>10
 l3b=NaN*ii;
 l3b(i3bd)=ax(4);
 plot(ii,l3b,'-r','linewidth',8)
 text(mean(ax(1:2)),ax(3),['found ',int2str(n3bd),'% profiles 3 beam solutions'],...
 'VerticalAlignment','bottom','HorizontalAlignment','center','fontsize',14)
end

ylabel('range [m]')
xlabel('ensemble')
if existf(p,'name')
  streamer([p.name,' Figure 2']);
end
title(' W as function of bindepth and time')

if existf(d,'tsd_m')
 subplot(427)
 plot(d.tsd_m(1:length(d.izd),:),-d.zd)
 axis tight
 ax=axis;
 if existf(d,'tsu_m')
  hold on
  plot(d.tsu_m(1:length(d.izu),:),d.zu)
  axis tight
  ax=axis;
  plot(ax(1:2),ax(1:2)*0,'-k')
 end

 t=d.tsd_m;
 checkbeam(t,ax,1)

 if existf(d,'tsu_m')
  t=d.tsu_m;
  checkbeam(t,ax,0)
 end
 axis(ax)
 ylabel('distance [m]')
 xlabel('echo amplidute [dB]')
 title('Beam Performance')
end

if existf(d,'cmd_m')
 subplot(428)
 plot(d.cmd_m(1:length(d.izd),:),-d.zd)
 axis tight
 ax=axis;
 [dum,dum,dum,x,y]=makebars(-d.zd,sum(isfinite(d.weight(d.izd,:))'));
 hold on
 fill(-y/max(y)*10,x,'r')

 if existf(d,'cmu_m')
  hold on
  plot(d.cmu_m(1:length(d.izu),:),d.zu)
  plot(ax(1:2),ax(1:2)*0,'-k')
  [dum,dum,dum,x,y]=makebars(d.zu,sum(isfinite(d.weight(d.izu,:))'));
  fill(-y/max(y)*10,x,'g')
  axis tight
  ax=axis;
  for i=1:4
   [dum,ir]=min(abs(d.cmu_m(1:length(d.izu),i)-max(d.cmu_m(1,:))*0.3));
   text((0.12*i+0.27)*ax(2),ax(4),int2str(d.zu(ir)),'VerticalAlignment','top')
   p.up_range(i)=d.zu(ir);
  end
  text(0.01*ax(2),ax(4),['#',int2str(p.up_sn),...
         ' range:'],'VerticalAlignment','top')
 end
 ax(1)=-13; 

 for i=1:4
  [dum,ir]=min(abs(d.cmd_m(1:length(d.izd),i)-max(d.cmd_m(1,:))*0.3));
  p.dn_range(i)=d.zd(ir);
  text((0.12*i+0.27)*ax(2),ax(3),int2str(d.zd(ir)),'VerticalAlignment','bottom')
 end
 text(0.01*ax(2),ax(3),['#',int2str(p.down_sn),...
         ' range:'],'VerticalAlignment','bottom')

 axis(ax)
 ylabel('distance [m]')
 xlabel('correlation ')
 title('Range of good data')
  
end

if existf(d,'z')
 subplot(813)
 plot(d.z)
 ylabel('depth')
 ax=axis;
 ax(3)=min(d.z*1.05);
 ax(4)=0;
 ax(2)=length(d.z);
 axis(ax)
end

if existf(d,'tilt')
 subplot(814)
 plot(d.tilt(1,:))
 ylabel('tilt [deg]')
 ax=axis;
 ax(2)=length(d.z);
 ax(4)=30;
 axis(ax)
end

if existf(d,'hdg')
 subplot(815)
 plot(d.hdg(1,:))
 ylabel('heading [deg]')
 ax=axis;
 ax(4)=360;
 ax(2)=length(d.z);
 axis(ax)
 set(gca,'YTick',[0 90 180 270 360])
end

if existf(d,'xmv')
 subplot(816)
 plot(d.xmv(1,:))
 text(length(d.xmv)/2,mean(d.xmv(1,:)),[' mean: ',num2str(fix(p.xmv(1)*10)/10)])
 ylabel('X-mit volt down')
 xlabel('ensemble')
 ax=axis;
 ax(2)=length(d.z);
 axis(ax)
end

%=============================================
function checkbeam(t,ax,is_bottom)
% check beam performance

%%%bl=size(t,2); 
bl = length(t);	%%% BH fix

iend=fix(bl/2):bl;

tax=mean(ax(1:2));
if is_bottom
 tay=ax(3);
 tflag='bottom';
else
 tay=ax(4);
 tflag='top';
end

for i=1:4;
% first correct for source level
 t(:,i)=t(:,i)-mean(t(:,i));
% find noise level
 tn(i)=medianan(t(iend,i),2);
 s2n(i)=mean(t(1:2,i))-tn(i);
end

ifail=s2n<max(s2n)*0.5;
ibad=~ifail & s2n<max(s2n)*0.65;
iweak=~ifail & ~ibad & s2n<max(s2n)*0.8;
for i=1:4;
text(ax(1)+0.2*i*diff(ax(1:2)),tay,[int2str(s2n(i)./max(s2n)*100),'%'],...
 'VerticalAlignment',tflag)
end

if sum(ifail)>0
 it=find(ifail==1);
 text(tax,tay*0.5,[' beam ',int2str(it),' broken'],...
 'VerticalAlignment',tflag,'HorizontalAlignment','center','fontsize',18)
end

if sum(ibad)>0
 it=find(ibad==1);
 text(tax,tay*0.65,[' beam ',int2str(it),' bad'],...
 'VerticalAlignment',tflag,'HorizontalAlignment','center','fontsize',15)
end

if sum(iweak)>0
 it=find(iweak==1);
 text(tax,tay*0.8,[' beam ',int2str(it),' weak'],...
 'VerticalAlignment',tflag,'HorizontalAlignment','center','fontsize',12)
end

return
