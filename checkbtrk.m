function p=checkbtrk(d,di,de,dr,p);
% compare RDI to own bottom track
% 

%======================================================================
%                    C H E C K B T R K . M 
%                    doc: Sat Jul  3 23:57:25 2004
%                    dlm: Wed Sep  4 19:02:39 2019
%                    (c) 2004 ladcp@
%                    uE-Info: 24 0 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

% CHANGES BY ANT:
%  Jul  4, 2004: - adapted to [edit_data.m], i.e. plot unedited data
%		   in black & edited in blue; also use unedited data
%		   to calculate bias slope
%  Jul 16, 2004: - adapted to changes in [edit_data.m]
%  Mar  5, 2010: - debugging???
%  Oct  1, 2010: - BUG: code bombed when bottom was detected (i.e. p.zbottom
%			was set) but when all these bottom returns turned out
%			to be outside the valid range (i.e. almost certainly 
%			false)
%  Sep  4, 2019: - fixed plotting BUGs reported by GK

if ~existf(p,'zbottom'), return , end
if ~isfinite(p.zbottom), return, end
if ~existf(d,'bvel_own'), d.bvel_own=d.bvel*nan; d.hbot_own=d.hbot*nan; end
if ~existf(d,'bvel_rdi'), d.bvel_rdi=d.bvel*nan; d.hbot_rdi=d.hbot*nan; end

  ii=find(abs(d.z+p.zbottom)<max(p.btrk_range) & ...
          abs(d.z+p.zbottom)>min(p.btrk_range));

  if isempty(ii), p.zbottom = nan, return, end
 
  if existf(d,'bvel_rdi')
   bu_rdi=d.bvel_rdi(ii,1);
   bv_rdi=d.bvel_rdi(ii,2);
   bw_rdi=d.bvel_rdi(ii,3);
  else
   bu_rdi=d.bvel(ii,1)*nan;
   bv_rdi=d.bvel(ii,2)*nan;
   bw_rdi=d.bvel(ii,3)*nan;
  end

  if existf(d,'hbot_rdi')
   bh_rdi=d.hbot_rdi(ii);
  else
   bh_rdi=d.hbot(ii)*nan;
  end

  u_o=interp1(dr.z,dr.u,-d.z(ii)');
  u_c=-interp1(dr.tim',dr.uctd',d.time_jul(ii)');

  bu_own=d.bvel_own(ii,1);
  bv_own=d.bvel_own(ii,2);
  bw_own=d.bvel_own(ii,3);
  bh_own=d.hbot_own(ii);
  v_o=interp1(dr.z,dr.v,-d.z(ii)');
  v_c=-interp1(dr.tim',dr.vctd',d.time_jul(ii)');
  
  rw=d.rw(d.izd,ii);
  ij=find((d.izm(d.izd,ii)+p.zbottom)<d.zd(1));
  rw(ij)=nan;
  if existf(d,'wctd') ==-1
   w_c=d.wctd(ii)';
   disp(' use CTD pressure derived W')
  else
   w_c=medianan(rw,3)';
   disp(' use reference layer W')
  end
  w_c2=w_c-[0;diff(w_c)/2 ];

  bu_used=d.bvel(ii,1);
  bv_used=d.bvel(ii,2);
  bw_used=d.bvel(ii,3);
  bh_used=d.hbot(ii);

  blen=d.down.Cell_length/100;
 
% extract raw velocity data near bottom
rube = d.ru(d.izd,ii);
rvbe = d.rv(d.izd,ii);
rwbe = d.rw(d.izd,ii);

% the same, but only velocities that passed editing
ru = d.ru; ru(find(~isfinite(d.weight))) = NaN; ru = ru(d.izd,ii);
rv = d.rv; rv(find(~isfinite(d.weight))) = NaN; rv = rv(d.izd,ii);
rw = d.rw; rw(find(~isfinite(d.weight))) = NaN; rw = rw(d.izd,ii);

% correct velocities with solved U_ctd

rua=ru-repmat(u_c',[length(d.izd) 1]); ruabe=rube-repmat(u_c',[length(d.izd) 1]);
rva=rv-repmat(v_c',[length(d.izd) 1]); rvabe=rvbe-repmat(v_c',[length(d.izd) 1]);
rwa=rw-repmat(w_c',[length(d.izd) 1]); rwabe=rwbe-repmat(w_c',[length(d.izd) 1]);
rwan=abs(rw)./abs(repmat(w_c',[length(d.izd) 1]));
rwanbe=abs(rwbe)./abs(repmat(w_c',[length(d.izd) 1]));



bua_rdi=bu_rdi-u_c;
bva_rdi=bv_rdi-v_c;

% try to use w inbetween pings for real bottom track
bwa1_rdi=bw_rdi-w_c;
bwa2_rdi=bw_rdi-w_c2;
if stdnan(bwa1_rdi)>stdnan(bwa2_rdi)
 bwa_rdi=bwa2_rdi;
 disp(' use inbetween w-ref')
else
 bwa_rdi=bwa1_rdi;
end


bua_own=bu_own-u_c;
bva_own=bv_own-v_c;
bwa_own=bw_own-w_c;

bua_used=bu_used-u_c;
bva_used=bv_used-v_c;
bwa_used=bw_used-w_c;

% save bias and std of bottom track anomaly
p.btrk_u_bias=medianan(bua_used,6);
p.btrk_u_std=stdnan(bua_used);
p.btrk_v_bias=medianan(bva_used,6);
p.btrk_v_std=stdnan(bva_used);
p.btrk_w_bias=medianan(bwa_used,6);
p.btrk_w_std=stdnan(bwa_used);

disp(['CHECKBTRK: check bottom track against U_ctd solution '])
disp([' profiles within give acceptable range: ',int2str(length(ii))])
disp(['  U bias :',num3str(p.btrk_u_bias,6,3),' [m/s]  std: ',...
       num3str(p.btrk_u_std,5,3),' [m/s]'])
disp(['  V bias :',num3str(p.btrk_v_bias,6,3),' [m/s]  std: ',...
       num3str(p.btrk_v_std,5,3),' [m/s]'])
disp(['  W bias :',num3str(p.btrk_w_bias,6,3),' [m/s]  std: ',...
       num3str(p.btrk_w_std,5,3),' [m/s]'])

if abs(p.btrk_u_bias)>0.1
 warn=[' large U bottom track bias ',num2str(p.btrk_u_bias)];
 disp(warn)
 p.warnp(size(p.warnp,1)+1,1:length(warn))=warn;
end

if abs(p.btrk_v_bias)>0.1
 warn=[' large V bottom track bias ',num2str(p.btrk_v_bias)];
 disp(warn)
 p.warnp(size(p.warnp,1)+1,1:length(warn))=warn;
end

% correct nominal depth with detected bottom depth
izm=d.izm(d.izd,ii);
zcor=(-d.hbot(ii)+d.z(ii));
zcor_std=stdnan(zcor);
p.btrk_rough=zcor_std;
izm=izm-repmat(zcor,[length(d.izd) 1]);

% compute bias W slope
 dz=abs(diff(d.zd(1:2)));
 wbias=[];
 zbias=[];
 nbias=[];
 for i=1:6
  dzi=dz*i/2;
  ij=find(abs(izm+dzi)<(dz/4));
  nok=sum(isfinite(rwanbe(ij)));
  if nok>5
   wbias=[wbias, 1-medianan(rwanbe(ij))];
   zbias=[zbias, dzi];
   nbias=[nbias, nok];
  end
 end
 if exist('wbias')==1
  p.wbslope=polyfit(zbias,wbias,1);
 else
  wbias=nan;
  zbias=nan;
  p.wbslope=[nan nan];
 end

disp(['  W slope fact :',num3str(p.wbslope(1),6,4),' [1/m] lower W below bottom '])

% compute 'wave' action
p.btrk_wdiff=medianan(abs(diff(w_c)));
disp(['  W diff :',num3str(p.btrk_wdiff,6,4),' [m/s] ping to ping w difference '])



disp(['  H std :',num3str(zcor_std,6,1),' [m]  large means bottom is rough/sloped'])

p.btrk_tilt_mean=mean(d.tilt(ii));
p.btrk_tilt_std=std(d.tilt(ii));
disp(['  Tilt mean :',num3str(p.btrk_tilt_mean,3,1),' +/- ',...
           num3str(p.btrk_tilt_std,3,1),' [^o]  '])


% select bins used for own bottom track
iz=find(abs(izm+blen*p.btrk_below)<blen);


% plot results
tx2=0.25;
figure(13)
clf
  
 subplot(221); hold on;
 plot(ruabe,izm,'.k');			% all velocities
 plot(rua,izm,'.g'); 			% edited velocities
 plot(ruabe(iz),izm(iz),'.r');		% used-for-BT velocities
 axis tight;
 ax=axis;
 ax(1:2)=[-1 1]*0.8;
 ax(3)=-160;
 ax(4)=50;
 vh=floor(ax(1)*20)/20+0.0125 : 0.025: ceil(ax(2)*20)/20;
 plot(ax(1:2),[0 0],'-')
 plot([0 0],ax(3:4),'k-','linewidth',1.5)
 plot([0 0]+p.btrk_u_bias,ax(3:4),'r-','linewidth',1.5)
 text(tx2,10,' Bottom')

 if sum(isfinite(bua_rdi))>0
  [y,x]=hist(bua_rdi,vh);
  y([1,end])=0;
  ys=sum(y);
  [x,y]=stairs(x,y);
  fill(x,y/max(y)*30-50,'-g')
  plot([0 0]+medianan(bua_rdi,6),ax(3:4),'-g')
  text(tx2,-40,['RDI n: ',int2str(ys)])
  text(0.1,-10,[' bias ',num3str(medianan(bua_rdi,6),6,3)],'color','g')
 end

 [y,x]=hist(bua_own,vh);
 y([1,end])=0;
 ys=sum(y);
 [x,y]=stairs(x,y);
 fill(x,y/max(y)*30-100,'-b')
 plot([0 0]+medianan(bua_own,6),ax(3:4),'-b')
 text(tx2,-90,[' own n: ',int2str(ys)])
 text(0.1,-60,[' bias ',num3str(medianan(bua_own,6),6,3)],'color','b')

 [y,x]=hist(real(de.bvel-de.uctd(:,1).'),vh);
 y([1,end])=0;
 ys=sum(y);
 [x,y]=stairs(x,y);
 fill(x,y/max(y)*30-150,'-r')
 text(tx2,-140,['Sup. ens. n: ',int2str(ys)])
 text(0.1,-110,[' bias ',num3str(p.btrk_u_bias,6,3)],'color','r')
 title([' U bot-tr mean U_{ctd} : ',num3str(-meannan(real(de.uctd(:,1))),4,2),...
        ' [m/s] added'])
 xlabel(' U [m/s]')
 ylabel('depth [m]')
 axis(ax)
 text(-0.7,-154,['Bottom roughness: ',int2str(p.btrk_rough),' [m]'])

subplot(222); hold on;
 plot(rvabe,izm,'.k'); 			% all velocities
 plot(rva,izm,'.g'); 			% edited velocities
 plot(rvabe(iz),izm(iz),'.r')		% used-for-BT velocities
 axis tight
 ax=axis;
 ax(1:2)=[-1 1]*0.8;
 ax(3)=-160;
 ax(4)=50;
 plot(ax(1:2),[0 0],'-')
 plot([0 0],ax(3:4),'k-')
 plot([0 0]+p.btrk_v_bias,ax(3:4),'r-','linewidth',1.5)

 text(tx2,10,' Bottom')
 if sum(isfinite(bva_rdi))>0
  [y,x]=hist(bva_rdi,vh);
  y([1,end])=0;
  ys=sum(y);
  [x,y]=stairs(x,y);
  fill(x,y/max(y)*30-50,'-g')
  plot([0 0]+medianan(bva_rdi,6),ax(3:4),'-g')
  text(tx2,-40,[' RDI n: ',int2str(ys)])
  text(0.1,-10,[' bias ',num3str(medianan(bva_rdi,6),6,3)],'color','g')
 end

 [y,x]=hist(bva_own,vh);
 y([1,end])=0;
 ys=sum(y);
 [x,y]=stairs(x,y);
 fill(x,y/max(y)*30-100,'-b')
 plot([0 0]+medianan(bva_own,6),ax(3:4),'-b')
 text(tx2,-90,[' own n: ',int2str(ys)])
 text(0.1,-60,[' bias ',num3str(medianan(bva_own,6),6,3)],'color','b')

 [y,x]=hist(imag(de.bvel-de.uctd(:,1).'),vh);
 y([1,end])=0;
 ys=sum(y);
 [x,y]=stairs(x,y);
 fill(x,y/max(y)*30-150,'-r')
 text(tx2,-140,['Sup. ens. n: ',int2str(ys)])
 text(0.1,-110,[' bias ',num3str(p.btrk_v_bias,6,3)],'color','r')
 title([' V bot-tr mean V_{ctd} : ',num3str(-meannan(imag(de.uctd(:,1))),4,2),...
        ' [m/s] added'])
 xlabel(' V [m/s]')
 ylabel('depth [m]')
 text(-0.7,-154,['Tilt : ',num3str(p.btrk_tilt_mean,3,1),' +/- ',...
                 num3str(p.btrk_tilt_std,3,1),' [^o]'])
 axis(ax)

tx2=0.05;

subplot(223); hold on;
 plot(rwabe,izm,'.k');				% all velocities
 plot(rwa,izm,'.g');				% edited velocities
 plot(rwabe(iz),izm(iz),'.r');			% used-for-BT velocities
 axis tight
 ax=axis;
 ax(1:2)=[-1 1]*0.4;
 vh=floor(ax(1)*20)/20+0.005 : 0.01: ceil(ax(2)*20)/20;
 ax(3)=-160;
 ax(4)=50;
 plot(ax(1:2),[0 0],'-')
 plot([0 0],ax(3:4),'k-')
 plot([0 0]+p.btrk_w_bias,ax(3:4),'r-','linewidth',1.5)

 text(tx2,10,' Bottom')
 if sum(isfinite(bwa_rdi))>0
  [y,x]=hist(bwa_rdi,vh);
  y([1,end])=0;
  ys=sum(y);
  [x,y]=stairs(x,y);
  fill(x,y/max(y)*30-50,'-g')
  plot([0 0]+medianan(bwa_rdi,6),ax(3:4),'-g')
  text(tx2,-40,[' RDI n: ',int2str(ys)])
  text(0.1,-10,[' bias ',num3str(medianan(bva_rdi,6),6,3)],'color','g')
 end

 [y,x]=hist(bwa_own,vh);
 y([1,end])=0;
 ys=sum(y);
 [x,y]=stairs(x,y);
 fill(x,y/max(y)*30-100,'-b')
 plot([0 0]+medianan(bwa_own,6),ax(3:4),'-b')
 text(tx2,-90,[' own n: ',int2str(ys)])
 text(0.1,-60,[' bias ',num3str(medianan(bwa_own,6),6,3)],'color','b')

 [y,x]=hist(di.bvel(3,:)-medianan(di.rw),vh);
 y([1,end])=0;
 ys=sum(y);
 [x,y]=stairs(x,y);
 fill(x,y/max(y)*30-150,'-r')
 text(tx2,-140,['Sup. ens. n: ',int2str(ys)])
 text(0.1,-110,[' bias ',num3str(p.btrk_w_bias,6,3)],'color','r')
 text(-0.3,-154,['W diff: ',num3str(p.btrk_wdiff,6,3),' [m/s]'])

 title([' W bottom track [m/s] '])
 xlabel(' W [m/s]')
 ylabel('depth [m]')
 axis(ax)

subplot(224); hold on;
 plot(rwanbe,izm,'.k');				% all velocities
 plot(rwan,izm,'.g');				% edited velocities
 plot(rwan(iz),izm(iz),'.r');			% used-for-BT velocities
 axis tight
 ax=axis;
 ax(1:2)=[0.5 1.2];
 ax(3)=-60;
 ax(4)=60;
 title([' W bottom track factor '])
 xlabel(' abs(W)/abs(W_{ref}) ')
 ylabel('depth [m]')
 axis(ax)
 plot(ax(1:2),[0 0],'-')
 plot([1 1],ax(3:4),'k-','linewidth',1.5)
 text(tx2,10,' Bottom')
 plot(1-polyval(p.wbslope,[40 0]),[-40 0],'-b')
 plot(1-wbias,-zbias,'pb')
 grid
 text(0.54,-50,[' bias slope ',num3str(p.wbslope(1),6,3),' [1/m]'])
 text(0.54,-58,[' offset  ',num3str(p.wbslope(2),6,3) ])
 text(0.54,50,[' beam angle  ',int2str(d.down.Beam_angle)])
 text(0.54,42,[' bin length  ',int2str(d.down.Cell_length/100),' [m]'])
 text(0.54,34,[' Pings/Ens  ',int2str(d.down.Pings_per_Ensemble)])
 text(0.54,26,[' Frequency  ',int2str(d.down.Frequency),' [kHz]'])

 text(0.54,10,[' Zbottom  ',int2str(p.zbottom),' [m]'])

 
  streamer([p.name,' Figure 13']);

  orient tall

