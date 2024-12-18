function plotinv(dr,d,p,ps)
% function plotinv(dr,d,p,ps)
%                              
% - plot final velocity profile
%
% Martin Visbeck and Gerd Krahmann, LDEO, April-2000
%

%======================================================================
%                    P L O T I N V . M 
%                    doc: Fri Mar 18 09:37:33 2005
%                    dlm: Fri Mar 18 09:56:54 2005
%                    (c) 2005 A.M. Thurnherr
%                    uE-Info: 20 67 NIL 0 0 72 0 2 4 NIL ofnI
%======================================================================

% CHANGES BY ANT:
%	Mar 18, 2005: - in cases of bad range, [geterr.m] does not return
%					finite velocity errors; therefore the code was
%					changed so as not to bomb on figure(1) any more
%   Aug 15, 2022: EPA: disable LaTeX interpreter for software version
%   July 25, 2024: LC: altered legend of Figure 1, added the meaning of 
%                  colour letter to the left panel depth vs velocity.

if existf(dr,'range_do');
 zpmax=p.maxdepth+maxnan([0;dr.range_do]);
else
 zpmax=p.maxdepth+maxnan([0;dr.range]);
end

if existf(p,'zbottom')
 if isfinite(p.zbottom)
  zpmax=p.zbottom;
 end
end

if ~existf(p,'name')
 p.name=' ';
end


ua=dr.u;
va=dr.v;
z=dr.z;

ur=maxnan([maxnan(abs(ua)), maxnan(abs(va))])*120;
if (ur<=0 | isnan(ur)), ur=1; end
ps=setdefv(ps,'urange',ur);

% plot final absolute velocity
clf
axes('position',[0.1 0.23 0.4 0.7])
plot(ua*100,-z,'-r','linewidth',2.5)
grid
hold on
plot(va*100,-z,'--g','linewidth',2.5)
ct=[' U(r-) V(g--)']; % also list the colours in the legend
if existf(dr,'uerr')==1
 plot([va-dr.uerr]*100,-z,':g','linewidth',1.8)
 plot([va+dr.uerr]*100,-z,':g','linewidth',1.8)
 plot([ua-dr.uerr]*100,-z,':r','linewidth',1.8)
 plot([ua+dr.uerr]*100,-z,':r','linewidth',1.8)
end
if existf(dr,'u_do')
 plot((dr.u_do+dr.ubar)*100,-z,'-r','linewidth',0.5)
 plot((dr.u_up+dr.ubar)*100,-z,'-r','linewidth',0.5)
 plot((dr.v_do+dr.vbar)*100,-z,'--g','linewidth',0.5)
 plot((dr.v_up+dr.vbar)*100,-z,'--g','linewidth',0.5)
 iz=3:5:length(dr.u_do);
 plot((dr.u_do(iz)+dr.ubar)*100,-z(iz),'.b','markersize',6)
 plot((dr.v_do(iz)+dr.vbar)*100,-z(iz),'.b','markersize',6)
 ct=[ct,'; blue dots down cast'];
end
if existf(dr,'u_shear_method')
 plot((dr.u_shear_method+dr.ubar)*100,-dr.z,'-r','linewidth',0.9)
 plot((dr.v_shear_method+dr.vbar)*100,-dr.z,'--g','linewidth',0.9)
 iz=2:4:length(dr.u_shear_method);
 plot((dr.u_shear_method(iz)+dr.ubar)*100,-dr.z(iz),'.r','markersize',7)
 plot((dr.v_shear_method(iz)+dr.vbar)*100,-dr.z(iz),'.g','markersize',7)
 ct=[ct,'; dotted shear'];
end
if existf(dr,'u_sadcp')
  plot(dr.u_sadcp*100,-dr.z_sadcp,'-r','linewidth',1.3)
  plot(dr.v_sadcp*100,-dr.z_sadcp,'-g','linewidth',1.3)
  plot(dr.u_sadcp*100,-dr.z_sadcp,'pr','linewidth',0.9)
  plot(dr.v_sadcp*100,-dr.z_sadcp,'pg','linewidth',0.9)
  ct=[ct,'; pentagon SADCP'];
end
if (p.btrk_used>0 & existf(dr,'zbot') )
 plot(dr.ubot*100,-dr.zbot,'r-','linewidth',1.3)
 plot(dr.vbot*100,-dr.zbot,'g--','linewidth',1.3)
 plot(dr.ubot*100,-dr.zbot,'r^','linewidth',1.3)
 plot(dr.vbot*100,-dr.zbot,'g^','linewidth',1.3)
end


ax=axis;
ax(4)=0;
ax(3)=-zpmax;
ps=setdefv(ps,'zrange',ax(3:4));
ax(1:2)=[-1 1]*ps.urange;
ax(3:4)=ps.zrange;
if existf(dr,'onlyshear')
 if dr.onlyshear
  text(ax(1),ax(3)*0.95,'  SHEAR SOLUTION ONLY','color','b','fontsize',15)
 end
end
if ps.up_dn_looker==2
  text(ax(1),ax(3)*0.90,'  DOWN LOOKER ONLY','color','b','fontsize',13)
end
if ps.up_dn_looker==3
  text(ax(1),ax(3)*0.90,'    UP LOOKER ONLY','color','b','fontsize',13)
end
plot([0 0],ax(3:4),'-k')
axis(ax)

title(ct)
ylabel('depth [m]')

% plot bottom track data

if (p.btrk_used>0 & existf(dr,'zbot') )
 plot([-1 1]*ps.urange,[1 1]*zpmax,'-k','linewidth',2)
 set(gca,'XTickLabel',[]);
 zbprofr=-dr.zbot+p.zbottom;
 axes('position',[0.1 0.1 0.4 0.12])
 iz=find(zbprofr>30 & zbprofr <250);
 plot(dr.ubot(iz)*100,zbprofr(iz),'r-','linewidth',2.5)
 grid
 hold on
 plot(dr.vbot(iz)*100,zbprofr(iz),'g--','linewidth',2.5)
 if existf(dr,'uerrbot')
  plot([dr.ubot(iz)+dr.uerrbot(iz)]*100,zbprofr(iz),':r','linewidth',1.8)
  plot([dr.ubot(iz)-dr.uerrbot(iz)]*100,zbprofr(iz),':r','linewidth',1.8)
  plot([dr.vbot(iz)+dr.uerrbot(iz)]*100,zbprofr(iz),':g','linewidth',1.8)
  plot([dr.vbot(iz)-dr.uerrbot(iz)]*100,zbprofr(iz),':g','linewidth',1.8)
 end
 plot([0 0],[0 -maxnan(-[250, p.zbottom])],'-k')
 xlabel('velocity [cm/s]')
 ylabel('above bottom [m]')
 ax2=[[-1 1]*ps.urange 0 min([250, p.zbottom])];
 axis(ax2);
 if p.btrk_used>1 
  text(0.9*ax2(1),ax2(4)*0.92,'post processed bottom track')
 else
  text(0.9*ax2(1),ax2(4)*0.92,'RDI bottom track')
 end
end


% write some text information
subplot(322)

poss=p.poss;
pose=p.pose;

iy=1.2;
ix=-0.1;
idy=1/8;

iy=iy-idy;
[slat,slon] = pos2str(poss(1)+poss(2)/60,poss(3)+poss(4)/60);
text(ix,iy,[' Start:'])
text(ix+0.2,iy,[slat])
text(ix+0.7,iy,[slon])

date1=p.time_start;  
date1(1)=y2k(date1(1));
ds=datenum(date1(1),date1(2),date1(3),date1(4),date1(5),date1(6));
iy=iy-idy;
text(ix+0.2,iy,[datestr(ds,0)])

[slat,slon] = pos2str(pose(1)+pose(2)/60,pose(3)+pose(4)/60);
iy=iy-idy;
text(ix,iy,[' End:'])
text(ix+0.2,iy,[slat])
text(ix+0.7,iy,[slon])

date2=p.time_end;    
date2(1)=y2k(date2(1));
ds=datenum(date2(1),date2(2),date2(3),date2(4),date2(5),date2(6));
iy=iy-idy;
text(ix+0.2,iy,[datestr(ds,0)])

iy=iy-idy;
text(ix,iy,sprintf('u-mean: %3.0f [cm/s]    v-mean %3.0f [cm/s]',...
    meannan(ua)*100,meannan(va)*100))

iy=iy-idy;
if length(d)>0
 if length(d.zu>2), zu=d.zu; else, zu=[0 0]; end
 text(ix,iy,['binsize do: ',num2str(diff(d.zd([1,2]))),' [m]  binsize up:  ',...
            int2str(diff(zu([1,2]))),' [m]'])
end

idy=1/10;
iy=iy-idy;
te=['mag. deviation ',num3str(p.drot,4,1)];
if existf(p,'deviat_used')
 te=[te,'  dev. dn: ',int2str(p.deviat_used(1))];
 te=[te,'  up: ',int2str(p.deviat_used(2))];
end
text(ix,iy,te)

iy=iy-idy;
text(ix,iy,['wdiff: ',num2str(p.wlim),...
'  pglim: ',num2str(p.pglim),'  elim ',num2str(p.elim)])

iy=iy-idy;
dum=' ';
if ps.smoofac~=0
 dum=[dum,'smo:',num3str(ps.smoofac,3,2)];
end
if ps.dragfac~=0
 dum=[dum,' dra:',num3str(ps.dragfac,3,2)]; 
end
if sum(ps.smallfac(:,2))>0
 dum=[dum,' smal:',int2str(ps.smallfac(1,1)),'-',int2str(ps.smallfac(end,1))]; 
end
if ps.barofac~=0
 dum=[dum,' bar:',num3str(ps.barofac,3,1)];
end
if ps.botfac~=0
 dum=[dum,' bot:',num3str(ps.botfac,3,1)];
end
if ps.sadcpfac~=0
 dum=[dum,' sad:',num3str(ps.sadcpfac,3,1)]; 
end
text(ix,iy,dum)

iy=iy-idy;
text(ix,iy,['weightmin ',num3str(ps.weightmin,3,1),...
'  weightpower: ',num3str(ps.weightpower,3,1)])


iy=iy-idy;
ctext=['max depth: ',int2str(p.maxdepth),' [m]'];
if p.zbottom>-Inf, ctext=[ctext,'   bottom: ',int2str(p.zbottom),' [m]']; end
text(ix,iy,ctext)

streamer(['Station : ',p.name,'  Figure 1']);
axis off

% plot profiles of target strength
axes('position',[0.58 0.39 0.1 .25])
plot(dr.ts,-dr.z/1000,'b-','linewidth',1.5)
hold on
plot(dr.ts_out,(-max(d.zd)-dr.z)/1000,'k-')
if isfinite(maxnan(dr.ts))
 ax(1)=-0.9*maxnan(-dr.ts_out); 
 ax(2)=maxnan(dr.ts)*1.1; 
end
ax(3:4)=ax(3:4)/1000;
axis(ax)
ylabel('depth [km]')
xlabel('target strength [dB]','color','b')
set(gca,'fontsize',10);

% plot profiles of super ensemble data range
axes('position',[0.70 0.39 0.1 .25])

plot(sum(dr.range,2),-dr.z/1000,'r-','linewidth',1.5)
if existf(dr,'range_up')
 hold on
 plot(dr.range_up,-dr.z/1000,'-k')
 plot(dr.range_do,-dr.z/1000,'-k')
 iz=round(linspace(2,length(dr.z)-2,10));
 plot(dr.range_do(iz),-dr.z(iz)/1000,'.b')
end
if isfinite(maxnan(sum(dr.range,2)))
 ax(1)=0; ax(2)=1.1*maxnan(sum(dr.range,2)); 
end
axis(ax)
title('range of instuments [m]','color','r')
set(gca,'YtickLabel',[])
set(gca,'fontsize',10)


if existf(dr,'uerr')==1

% plot profiles of velocity error
axes('position',[0.82 0.39 0.1 .25])

ue=dr.uerr;
plot(ue,-dr.z/1000,'k-','linewidth',1.5)
hold on
ax(1)=0; ax(2)=3*medianan(ue);
if ~isfinite(ax(2)), ax(2) = 1; end;
axis(ax)
xlabel('vel error (-k) [m/s]')
set(gca,'fontsize',10)
set(gca,'YtickLabel',[])

end

% plot single ping accuary also
if existf(dr,'ensemble_vel_err')==1
 plot(dr.ensemble_vel_err,-dr.z/1000,'-b')
 ax(2)=max([2.5*medianan(ue),2.5*medianan(dr.ensemble_vel_err)]);
 axis(ax)
 title('single ping (-b)','color','b')
end

% plot position of CTD

subplot(326)

xctd=dr.xctd;
yctd=dr.yctd;
ii=fix(linspace(1,length(xctd),10));
[m,ib]=min(dr.zctd);
plot(xctd,yctd,'linewidth',2)
hold on
plot(xctd(ii),yctd(ii),'r.','markersize',10)
plot(xctd(ib),yctd(ib),'g+','markersize',9)
if existf(dr,'xship')
 plot(dr.xship,dr.yship,'-g',dr.xship(ii),dr.yship(ii),'k.','markersize',10)
 plot([xctd(ii);dr.xship(ii)],[yctd(ii); dr.yship(ii)],'-y','linewidth',0.5)
 xlabel('CTD-position (blue) and ship (green) east-west [m]')
else
 xlabel('CTD-position east-west [m]')
end
uship_a=p.uship+sqrt(-1)*p.vship;
if ~(abs(uship_a)==0 & p.lat==0 & p.lon==0)
 plot(p.xdisp,p.ydisp,'gp')
 text(p.xdisp,p.ydisp,' GPS-end')
end
text(xctd(ib),yctd(ib),' bottom')
text(xctd(end),yctd(end),' end')
axis equal
text(xctd(1),yctd(1),' start')
ylabel('north-south [m]')
grid
set(gca,'fontsize',10)

axes('position',[0.6 0.05 0.01 0.01])
text(0,0,p.software,'interpreter','none')
axis off


orient tall


%========================================================================
function d=y2k(d)
% fix date string
if d<80, d=2000+d; end
if d<100, d=1900+d; end

