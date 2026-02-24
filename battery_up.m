function p=battery_up(p)
% function p=battery_up(p)
% try to find know battery calibration and issue warning if low

% battery level first low second critical
batlevel=[40 37];

%  WH 149
if p.instid(2)==102206758
 disp(' found CPU board of serial 149 ')
 p.battery_up=0.3*p.xmv(2);

%  WH 754
elseif p.instid(2)==2474849359
 disp(' found CPU board of serial 754 ')
 p.battery_up=0.37*p.xmv(2);
else
 disp(' do not know calibration of this instrument make a guess: ')
 p.battery_up=0.33*p.xmv(2);
end

if p.battery_up>batlevel(1)
 bc='g';
elseif p.battery_up>batlevel(2)
 bc='r';
else
  warn=([' Battery up voltage is low : ',num2str(round(p.battery_up*10)/10),' V'])
    p.warn(size(p.warn,1)+1,1:length(warn))=warn;
 bc='r';
end
hTxt = findobj(gcf,'Type','text','Tag','Battery_do');
hTxt.String = {
    hTxt.String
    ['Battery up Voltage is ', num2str(round(p.battery_up*10)/10), ' V'] };


disp([' Battery up Voltage is ',num2str(round(p.battery_up*10)/10),' V'])

