function p=battery(p)
% function p=battery(p)
% try to find know battery calibration and issue warning if low

% battery level first low second critical
batlevel=[40 37];

%  WH 149
if p.instid(1)==102206758
 disp(' found CPU board of serial 149 ')
 p.battery=0.3*p.xmv(1);

%  WH 754
elseif p.instid(1)==2474849359
 disp(' found CPU board of serial 754 ')
 p.battery=0.37*p.xmv(1);
else
 disp(' do not know calibration of this instrument make a guess: ')
 p.battery=0.33*p.xmv(1);
end

if p.battery>batlevel(1)
 bc='g';
elseif p.battery>batlevel(2)
 bc='r';
else
  warn=([' Battery voltage is low : ',num2str(round(p.battery*10)/10),' V'])
    p.warn(size(p.warn,1)+1,1:length(warn))=warn;
 bc='r';
end
text(0,0,['Battery Voltage is ',num2str(round(p.battery*10)/10),' V'],'color',bc,...
       'fontsize',14,'fontweight','bold',...
       'Tag', 'Battery_do')

disp([' Battery Voltage is ',num2str(round(p.battery*10)/10),' V'])

