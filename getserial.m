%======================================================================
%                    G E T S E R I A L . M 
%                    doc: Wed Jan  7 16:30:11 2009
%                    dlm: Fri Mar  5 15:48:31 2010
%                    (c) 2009 A.M. Thurnherr
%                    uE-Info: 101 13 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

% CHANGES BY ANT:
%   Jan  7, 2009: - tightened use of exist()
%   Aug 15, 2022: - EPA: added support for separate folder for logs
%                   EPA: added ".TXT" as suffix for logs
%                   EPA: support ADCP serial numbers above 9999

function p=getadcpserial(f,p)
% try to read .LOG or RECOVE.LOG
% decode serial number
%

p=setdefv(p,'down_sn',nan);
p=setdefv(p,'up_sn',nan);

% check in local lookup table
if exist('adcpserial.tab','file')
 disp(' found local lookup table for serial # to instid ')
 load adcpserial.tab 
 if existf(p,'instid')
  [dum(1),ii]=min(abs(adcpserial(:,1)-p.instid(1)));
  if dum(1)==0
   p.down_sn=adcpserial(ii,2);
  end
  if length(p.instid)>1
   [dum(2),ii]=min(abs(adcpserial(:,1)-p.instid(2)));
   if dum(2)==0
    p.up_sn=adcpserial(ii,2);
   end
  end
 end
end

% could not find instrument in lookup table, try log file

if isfield(f,'ladcplogdo')
 flog=findlogwh(f.ladcplogdo);
else
 flog=findlogwh(f.ladcpdo);
end

if length(flog)>1
 [p.down_sn,p.down_rawlog]=readserialwh(flog);
else
 flog=findlogbb(f.ladcpdo);
 if length(flog)>3
  [p.down_sn,p.down_rawlog]=readserialbb(flog);
 end
end
disp([' down looking serial number is : ',int2str(p.down_sn)])

% make sure to put serial number in data base
if length(flog)>1
 if exist('adcpserial.tab','file')
  load adcpserial.tab 
  [dum,ii]=min(abs(adcpserial(:,2)-p.down_sn));
 else
  adcpserial=[0 0];
  dum=1;
 end
 if isfinite(p.down_sn)
  if dum~=0
   disp(' found new instrument serial number id pair save to ADCPSERIAL.TAB ') 
   adcpserial(end+1,:)=[p.instid(1) p.down_sn];
  else
   if p.instid(1)~=adcpserial(ii,1)
    disp(['found new CPU board for ',int2str(p.down_sn),' add to ADCPSERIAL.TAB '])
   end
  end
  fid = fopen('adcpserial.tab','w');
  fprintf(fid,'%16.0f  %5.0f\n',adcpserial');
  fclose(fid);
 else
  disp(['found new CPU board id ',int2str(p.instid(1)),...
        ' add it to ADCPSERIAL.TAB with serial number'])
 end
end
 

% up looker 
if isfield(f,'ladcplogup')
 flog=findlogwh(f.ladcplogup);
else
 flog=findlogwh(f.ladcpup);
end

if length(flog)>1
 [p.up_sn,p.up_rawlog]=readserialwh(flog);
else
 flog=findlogbb(f.ladcpup);
 if length(flog)>1
  [p.up_sn,p.up_rawlog]=readserialbb(flog);
 end
end
disp(['   up looking serial number is : ',int2str(p.up_sn)])

% make sure to put serial number in data base
if length(flog)>1
 if exist('adcpserial.tab','file')
  load adcpserial.tab
  [dum,ii]=min(abs(adcpserial(:,2)-p.up_sn));
 else
  adcpserial=[0 0];
  dum=1;
 end
 if isfinite(p.up_sn)
  if dum~=0
   disp(' found new instrument serial number id pair save to ADCPSERIAL.TAB ')
   adcpserial(end+1,:)=[p.instid(2) p.up_sn];
  else
   if p.instid(2)~=adcpserial(ii,1)
    disp(['found new CPU board for ',int2str(p.down_sn),' add to ADCPSERIAL.TAB '])
   end
  end

  fid = fopen('adcpserial.tab','w');
  fprintf(fid,'%16.0f  %5.0f\n',adcpserial');
  fclose(fid);
 else
  disp(['found new CPU board id ',int2str(p.instid(1)),...
        ' add it to ADCPSERIAL.TAB with serial number'])
 end
end

 


%====================================
function flog=findlogwh(fin)
% look for logfile in data directory

% ignore empty file
flog=' ';
if length(fin)<2, return, end

disp([' FINDLOG: search near file: ',fin])
% check if .LOG exists
dum=find(fin=='.');
if dum>0
 flog=[fin(1:max(dum)),'LOG'];
 if exist(flog,'file')
  return
 end
 flog=[fin(1:max(dum)),'log'];
 if exist(flog,'file')
  return
 end

 % check if .TXT exists
 flog=[fin(1:max(dum)),'TXT'];
 if exist(flog,'file')
  return
 end
 flog=[fin(1:max(dum)),'txt'];
 if exist(flog,'file')
  return
 end

% check if 000.LOG exists
 fin1=fin(1:(max(dum)-4));
 flog=[fin1,'000.LOG'];
 if exist(flog,'file')
  return
 end
 flog=[fin1,'000.log'];
 if exist(flog,'file')
  return
 end
end

% look in data directory for RECOVER.LOG
dum=find(fin=='/' | fin=='\') ;
if length(dum)>0
 i1=1:dum(end);
 fdir=fin(i1);
else
 fdir='';
end

flog=[fdir,'RECOVER.LOG'];
if exist(flog,'file')
 return
end

flog=[fdir,'recover.log'];
if exist(flog,'file')
 return
end

disp([' FINDLOG: give up'])
flog='';

return

%====================================
function flog=findlogbb(fin)
% look for logfile in data directory

% ignore empty file
flog=' ';
if length(fin)<2, return, end

disp([' FINDLOG: search near file: ',fin])
% check if .LOG exists
dum=find(fin=='/' | fin=='\') ;

if length(dum)<1, return, end

if length(dum)>1, i1=dum(end-1)+1; else i1=1; end

fname=fin(i1:[dum(end)-1]);
flog=[fin(1:max(dum)),lower(fname)];
if exist(flog,'file')
  return
end

flog=[fin(1:max(dum)),upper(fname)];
if exist(flog,'file')
  return
end

disp([' FINDLOG: give up'])
flog='';

return

%=========================================
function [s,A]=readserialwh(f)
% open file ans parse for serial number

 disp([' READSERIAL_WH: parse file ',f])
 id=fopen(f);
 A=fread(id);
 fclose(id);
 As=A;
 ii=find(As<32); As(ii)=[]; 
 ii=findstr(As','  Instrument S/N:');
 if length(ii)>0
  s=sscanf(char(As(ii(1)+17+[1:9])'),'%g')';
  if (s>0 & s<99999)
%  disp([' serial number is : ',int2str(s)])
  else
   s=nan;
  end
 else
  disp([' Can not find RDI serial number '])
  ii=findstr(A','Board Serial Number Data:');
  if length(ii)>0
     disp([' Check BOARD-ID: ',char(A(ii(1)+[26:208])')])
   s=sum(A(ii(1)+[26:208]));
   switch s
    case 8600
     s = 149;
    case 8560;
     s = 150.2;
    case 8547;
     s = 150.1;
    case 8477
     s = 299;
    case 8521
     s = 754;
    otherwise
     disp([' unknown BOARD-ID';])
     disp([' madeup serial number is : ',int2str(s)])
   end
  else
   s=nan;
  end
 end
 A=char(A)';
 

%=========================================
function [s,A]=readserialbb(f)
% open file ans parse for serial number

 disp([' READSERIAL_BB: parse file ',f])
 id=fopen(f);
 A=fread(id);
 fclose(id);
 As=A;
 ii=find(As<32); As(ii)=[]; 
 ii=findstr(As','Xducer Ser #:');
 s=sscanf(char(As(ii+13+[1:9])'),'%g')';
 if (s>0 & s<9999)
%  disp([' serial number is : ',int2str(s)])
 else
  s=nan;
 end
 A=char(A)';
 

