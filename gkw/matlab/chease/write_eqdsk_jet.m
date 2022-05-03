function  jet=write_eqdsk_jet(shot,t1,t2,eqm,usrnm,seq,pth,filename);
% Reads JET EFIT equilibria from PPF using MDSplus and writes a G-EQDSK file,
% takes a time average of all the equilibria in a given time interval
%
% write_eqdsk_jet(shot,t1,t2,[eqm],[usrnm],[seq],[pth],[file]);
%
% Where
%    shot:   integer  JET pulse number
%      t1:   real     start of the time interval to average
%      t2:   real     end of the time interval to average
%     eqm:   string   optional DDA name (default: 'efit')
%   usrnm:   string   optional PPF username (default: 'jetppf')
%     seq:   integer  optional PPF sequence number (default: 0)
%     pth:   string   optional directory to write in (default './')
%    file:   string   optional output filename (default g_eqdsk_[dda]_[shot]_[t1]_[t2]_[seq])
%                     use blank [] to supress file output
%
% Example
%
%   data = write_eqdsk_jet(85307,50.8,51.0,'eftm',[],0,'~/work/chease/eqdsk')
%
% Produces EQDSK which can be read by chease_90_10_2 with COCOS_IN = 3, 
% with the same signs as those produced directly by EFIT (using abs(q)).
% See http://dx.doi.org/10.1016/j.cpc.2012.09.010
%
% The data, metadata, sequence nubmer, and filename are also returned 
% to the matlab structure assigned to the function ('data', in the example)
%
% Requirements: * MDSplus is installed 
%               * You are at a Fusion lab that can access JET data by MDSplus
%
% F.J. Casson, 19/2/2014 
% with contributions from N. Hawkes and C. Angioni
%
% Updated versions of this script may appear at
% https://bitbucket.org/gkw/gkw/src/develop/matlab/chease/write_eqdsk_jet.m
% http://git.ccfe.ac.uk/fcasson/matlab/blob/master/gkw/chease/write_eqdsk_jet.m

if (~exist('eqm','var') || isempty(eqm))
   eqm='efit';
end

if (~exist('usrnm','var') || isempty(usrnm))
  usrnm='jetppf';
else  
  disp(usrnm)
end

if (~exist('seq','var') || isempty(seq))
  seq=0;
end

if (~exist('pth','var')  || isempty(pth))
  pth='./'
end

if (~exist('filename','var'))
  filename=1;
end

if (shot < 72124)
    error('JET shots prior to 72124 do not contain the PSI grid in the EFIT PPFs');
    % other PPFs, e.g. EFTM may be even later
end

mdsconnect('mdsplus.jet.efda.org');

%%%%%% 2D quantities

% get the psi grid - needs reshaping into 33 x 33 or 65 x 65 array
% not in JET data handbook, but can be found with PPFSUM
% (N Hawkes mentioned it by email)
tmp = mdsgdat(shot, 'ppf', [eqm '/psi/' num2str(seq)], usrnm);
inds=find(tmp.t>=t1 & tmp.t <= t2);
times = tmp.t(inds)

if (size(inds,1)<1)
    error('Selected PPF does not exist, or does not contain requested time interval');
end

tmin = min(tmp.t)
tmax = max(tmp.t)

if (size(inds,2)<1)
    %tmp.t
    error('Selected PPF does not contain requested time interval');
end

sizepsi=sqrt(size(tmp.r,1));
jet.psirz=reshape(mean(tmp.data(:,inds),2),sizepsi,sizepsi);

%%%%%% 1D quantities

% Get the RZ values of the flux grid

tmp = mdsgdat(shot, 'ppf', [eqm '/psir/' num2str(seq)], usrnm);
jet.psir = tmp.data;
tmp = mdsgdat(shot, 'ppf', [eqm '/psiz/' num2str(seq)], usrnm);
jet.psiz = tmp.data;

% Derived quantities needed for EQDSK
jet.rwid  = max(jet.psir)-min(jet.psir);      % Horizontal dimension of computaitonal box
jet.zhei  = max(jet.psiz)-min(jet.psiz);      % Vertical dimension of computational box
jet.rleft = min(jet.psir);                    % Minimum R of computational box
jet.zmid  = (max(jet.psiz)+min(jet.psiz))/2;  % Z of centre of computational box

% Pressure
tmp= mdsgdat(shot, 'ppf', [eqm '/p/' num2str(seq)], usrnm);
jet.pres = mean(tmp.data(:,inds),2);

% Pressure gradient
tmp = mdsgdat(shot, 'ppf', [eqm '/dpdp/' num2str(seq)], usrnm);
jet.pprime = mean(tmp.data(:,inds),2);

% F= R Bt
tmp = mdsgdat(shot, 'ppf', [eqm '/f/' num2str(seq)], usrnm);
jet.fpol = mean(tmp.data(:,inds),2);

% Gradient of F: (dF / dPsi) / Mu0
tmp = mdsgdat(shot, 'ppf', [eqm '/dfdp/' num2str(seq)], usrnm);
jet.ffprim = mean(tmp.data(:,inds),2)*1.2566e-6;

% safety factor q
tmp = mdsgdat(shot, 'ppf', [eqm '/q/' num2str(seq)], usrnm);
jet.qpsi = mean(tmp.data(:,inds),2);

% LCFS Boundary coords
tmp = mdsgdat(shot, 'ppf', [eqm '/zbnd/' num2str(seq)], usrnm);
jet.zbbbs = mean(tmp.data(1:end-1,inds),2);
tmp = mdsgdat(shot, 'ppf', [eqm '/rbnd/' num2str(seq)], usrnm);
jet.rbbbs = mean(tmp.data(1:end-1,inds),2);

% Limiter contour (hard coded for ILW 2013 below)
[jet.rlim jet.zlim] = limiter;

%%%%%% Scalar Quantities

% Plasma Current (Amperes)
tmp = mdsgdat(shot, 'ppf', [eqm '/xip/' num2str(seq)], usrnm);
jet.current = mean(tmp.data(inds));

% Magnetic axis coords
tmp = mdsgdat(shot, 'ppf', [eqm '/zmag/' num2str(seq)], usrnm);
jet.zmaxis = mean(tmp.data(inds));
tmp = mdsgdat(shot, 'ppf', [eqm '/rmag/' num2str(seq)], usrnm);
jet.rmaxis = mean(tmp.data(inds));

%poloidal flux at axis (Weber / 2 Pi)
tmp = mdsgdat(shot, 'ppf', [eqm '/faxs/' num2str(seq)], usrnm);
jet.simag = mean(tmp.data(inds));

% poloidal flux at boundary (Weber / 2 Pi)
tmp = mdsgdat(shot, 'ppf', [eqm '/fbnd/' num2str(seq)], usrnm);
jet.sibry = mean(tmp.data(inds));

% Vacuum magnetic field (Tesla) at given R
jet.rcentr = 2.96;
tmp = mdsgdat(shot, 'ppf', [eqm '/bvac/' num2str(seq)], usrnm);
jet.bcentr = mean(tmp.data(inds));

% Total magnetic field on axis (Tesla)
% Not strictly correct to commute the time averaging here
tmp = mdsgdat(shot, 'ppf', [eqm '/btax/' num2str(seq)], usrnm);
jet.baxis=interp1(tmp.r,mean(tmp.data(:,inds),2),jet.rmaxis)

%jet.bcentr = mean(tmp.data(inds));

%jet.simag = -jet.simag;
%jet.sibry = -jet.sibry;
%jet.psirz = -jet.psirz;
%jet.fpol = -jet.fpol;

%yyyo = jet;

%jet.psirz = fliplr(yyyo.psirz);
%jet.zbbbs = -yyyo.zbbbs;
%jet.zlim = -yyyo.zlim;
%jet.zmaxis = -yyyo.zmaxis;

%If possible get the metadata on the sequence number - don't yet know how to do this with mdsplus
try
  [dum1 dum2 dum3 dum4 dum5 dum6 dum7 tseq]=ppfread(squeeze(shot),upper(eqm),'ZMAG',seq,upper(usrnm));
  clear dum*
  jet.seq=tseq;
catch
  %u=mdsvalue('_sig=ppfuid("jetppf")')
  %tseq=mdsvalue('_sig=pdmseq("jetppf")')   
  %tseq=mdsvalue('_sig=ppfseq("jetppf")')   
  jet.seq=seq;
end    

% return the other metadata
jet.date=date;
jet.eqm=eqm;
jet.usrnm=usrnm;
jet.times=times;

if isempty(filename) return; end

% Write the data to eqdsk, for definition, see:
% https://fusion.gat.com/conferences/snowmass/working/mfe/physics/p3/equili
% bria/g_eqdsk_s.pdf

if (filename==1)
  if (size(inds,2) > 1)
     filename=['g_JET_' eqm '_' num2str(shot) '_t' num2str(t1,'%7.4f') '_' num2str(t2,'%7.4f') '_' num2str(jet.seq,'%i')]
  else
     filename=['g_JET_' eqm '_' num2str(shot) '_t' num2str(jet.times,'%7.4f') '_' num2str(jet.seq,'%i')]
  end
end

jet.pth=pth;
jet.filename=filename;

[pth '/' filename]
fptr = fopen([pth '/' filename], 'w');

if (size(inds,2) > 1)
  fprintf(fptr,['  JET  ' eqm '   # ' num2str(shot,'%06d') '  ' num2str(t1,'%4.1f') '  ' num2str(t2,'%4.1f') '  seq: ' num2str(jet.seq,'%04d') '     ']);
else
  fprintf(fptr,['  JET  ' eqm '   # ' num2str(shot,'%06d') '  ' num2str(jet.times,'%7.4f') '     seq: ' num2str(jet.seq,'%04d') '     ']);
end  
%fprintf(fptr,'  JET      02/24/2009    # 23529  3600ms          ')

fprintf(fptr, '%2i%4i%4i\n', 3, sizepsi, sizepsi);

xdum = 0.0;

zformat = '%16.9E%16.9E%16.9E%16.9E%16.9E\n';
%zformat = '%16.9E %16.9E %16.9E %16.9E %16.9E \n';

fprintf(fptr, zformat, jet.rwid, jet.zhei, jet.rcentr, jet.rleft, jet.zmid);
fprintf(fptr, zformat, jet.rmaxis, jet.zmaxis, jet.simag, jet.sibry, jet.bcentr);
fprintf(fptr, zformat, jet.current, jet.simag, xdum, jet.rmaxis, xdum);
fprintf(fptr, zformat, jet.zmaxis, xdum, jet.sibry, xdum, xdum);

% Profiles
fprintf(fptr, zformat, jet.fpol);
if (mod(size(jet.fpol),5)~=0); fprintf(fptr,'\n'); end

fprintf(fptr, zformat, jet.pres);
if (mod(size(jet.pres),5)~=0); fprintf(fptr,'\n'); end

fprintf(fptr, zformat, jet.ffprim);
if (mod(size(jet.ffprim),5)~=0); fprintf(fptr,'\n'); end

fprintf(fptr, zformat, jet.pprime);
if (mod(size(jet.pprime),5)~=0); fprintf(fptr,'\n'); end

% Psi grid
fprintf(fptr, zformat, jet.psirz);
if (mod(sizepsi*sizepsi,5)~=0); fprintf(fptr,'\n'); end

% q profile
fprintf(fptr, zformat, jet.qpsi);
if (mod(size(jet.qpsi),5)~=0); fprintf(fptr,'\n'); end

%fprintf(fptr,'\n');
fprintf(fptr, '%5i%5i\n',length(jet.rbbbs), length(jet.rlim));

clear aaa;
for ij=1:length(jet.rbbbs);
    ijn = 2*ij-1;
    aaa(ijn) = jet.rbbbs(ij);
    aaa(ijn+1) = jet.zbbbs(ij);
end
fprintf(fptr, zformat, aaa);
if (mod(size(aaa),5)~=0); fprintf(fptr,'\n'); end

clear aaa;
for ij=1:length(jet.rlim);
    ijn = 2*ij-1;
    aaa(ijn) = jet.rlim(ij);
    aaa(ijn+1) = jet.zlim(ij);
end

fprintf(fptr, zformat, aaa);
if (mod(size(aaa),5)~=0) fprintf(fptr,'\n'); end

fclose(fptr);

mdsdisconnect;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r,z] = limiter
% Returns the limiter contour for shot 85307 (ILW, 2013)

r = [3.2832,3.3119,3.3284,3.3524,3.3732,3.4233,3.4473,3.4943,3.5168,3.5591,3.5805,3.6193,3.639,3.6736,3.688,3.7223,3.735,3.7682,3.7758,3.8026,3.8108,3.8333,3.8396,3.8565,3.8621,3.8752,3.8784,3.8857,3.8879,3.8905,3.8914,3.8886,3.8884,3.8803,3.8789,3.8656,3.8635,3.8445,3.8409,3.8171,3.813,3.7837,3.779,3.7441,3.7382,3.6987,3.6918,3.6749,3.6674,3.6373,3.6421,3.6214,3.6208,3.5941,3.5935,3.5668,3.5662,3.5395,3.5389,3.5122,3.5116,3.4849,3.4844,3.4577,3.4571,3.4304,3.4298,3.4031,3.4025,3.3818,3.3315,3.2818,3.1863,3.1366,3.001,2.8689,2.8525,2.7756,2.7576,2.679,2.6589,2.5779,2.5606,2.4801,2.466,2.3884,2.3742,2.298,2.2875,2.1954,2.1824,2.1657,2.1653,2.149,2.1486,2.1323,2.1319,2.1156,2.1152,2.0989,2.0985,2.0822,2.0818,2.0682,2.0676,2.0547,2.0544,2.0415,2.0412,2.0283,2.028,2.0151,2.0148,2.0019,2.0016,1.9888,1.9885,1.9756,1.9753,1.9613,1.9299,1.927,1.9226,1.9425,1.9273,1.925,1.9124,1.9097,1.8984,1.8957,1.8858,1.883,1.8731,1.8704,1.8605,1.8584,1.8499,1.8482,1.8418,1.841,1.8372,1.837,1.8359,1.8364,1.8379,1.8391,1.8431,1.845,1.8517,1.8542,1.8634,1.8666,1.8785,1.8823,1.8967,1.9013,1.9182,1.9234,1.9427,1.9486,1.9706,1.9597,1.9618,2.0091,2.0204,2.0207,2.0345,2.0348,2.0486,2.049,2.0628,2.0631,2.0769,2.0772,2.091,2.0913,2.1051,2.1055,2.1193,2.1196,2.1334,2.1337,2.1475,2.1478,2.1616,2.1619,2.1757,2.1761,2.1899,2.1902,2.2015,2.1446,2.2936,2.2936,2.2954,2.3599,2.3962,2.4091,2.4122,2.4129,2.4129,2.4122,2.4076,2.398,2.4192,2.4212,2.4188,2.4163,2.4057,2.315,2.3535,2.3743,2.4274,2.4462,2.5237,2.5246,2.5591,2.553,2.5739,2.633,2.6337,2.6938,2.6943,2.7544,2.7552,2.8147,2.8147,2.8043,2.857,2.8785,2.9364,2.9573,2.987,2.8977,2.882,2.8816,2.9005,2.8905,2.8879,2.8859,2.8859,2.8895,2.9008,2.9133,2.9635,3.0097,3.06,3.194,3.2022,3.3063,3.2832];
z = [-1.1244,-1.0832,-1.0631,-1.0387,-1.0172,-0.95992,-0.93229,-0.87186,-0.84284,-0.78155,-0.75015,-0.68679,-0.6539,-0.58769,-0.56088,-0.48496,-0.4576,-0.37137,-0.35239,-0.2692,-0.24458,-0.15797,-0.13481,-0.05075,-0.02427,0.06852,0.08885,0.17625,0.1963,0.29045,0.31339,0.40756,0.42674,0.51993,0.53958,0.63228,0.64971,0.74347,0.76273,0.85325,0.87053,0.96125,0.97699,1.0672,1.0834,1.1705,1.1866,1.2188,1.2364,1.3339,1.4077,1.4266,1.4271,1.4514,1.4519,1.4762,1.4767,1.501,1.5016,1.5258,1.5264,1.5507,1.5512,1.5755,1.576,1.6003,1.6008,1.6251,1.6257,1.6445,1.7041,1.7387,1.8175,1.8521,1.8834,1.9396,1.9453,1.968,1.972,1.9829,1.984,1.9829,1.9811,1.9677,1.9642,1.9389,1.9331,1.8947,1.8884,1.8228,1.8237,1.7908,1.7901,1.7579,1.7572,1.7251,1.7243,1.6922,1.6915,1.6593,1.6586,1.6264,1.6257,1.5988,1.5982,1.5645,1.5637,1.53,1.5293,1.4956,1.4949,1.4612,1.4604,1.4267,1.426,1.3923,1.3915,1.3578,1.3571,1.3206,1.273,1.261,1.254,1.2346,1.1583,1.138,1.0607,1.04,0.96207,0.94009,0.86238,0.84067,0.76296,0.74168,0.66408,0.64287,0.56549,0.54523,0.4671,0.44616,0.36833,0.34737,0.27002,0.2485,0.17099,0.14967,0.07136,0.05067,-0.02749,-0.04782,-0.12484,-0.14554,-0.22313,-0.24315,-0.32025,-0.34019,-0.41681,-0.43655,-0.51204,-0.53185,-0.60693,-0.62652,-0.65756,-0.78399,-0.8113,-0.81204,-0.84537,-0.84612,-0.87945,-0.88019,-0.91352,-0.91426,-0.9476,-0.94834,-0.98167,-0.98241,-1.0157,-1.0167,-1.05,-1.0507,-1.0841,-1.0848,-1.1181,-1.1189,-1.1522,-1.153,-1.1863,-1.187,-1.2204,-1.2211,-1.2484,-1.2749,-1.3148,-1.3314,-1.3344,-1.3344,-1.3732,-1.4003,-1.422,-1.4315,-1.4685,-1.4768,-1.5044,-1.5164,-1.5922,-1.6102,-1.6428,-1.6561,-1.6897,-1.7387,-1.7387,-1.735,-1.7135,-1.7098,-1.7098,-1.7,-1.655,-1.638,-1.6018,-1.6171,-1.6199,-1.6355,-1.6382,-1.6548,-1.6566,-1.672,-1.7079,-1.7116,-1.7116,-1.716,-1.7414,-1.7459,-1.7459,-1.6823,-1.6228,-1.5916,-1.5104,-1.4984,-1.4892,-1.474,-1.4357,-1.4171,-1.3928,-1.3762,-1.3348,-1.3348,-1.2978,-1.214,-1.2089,-1.2089,-1.1244];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = mdsgdat(shot, strppf, strdata,usrid)
%
% function y = mdsgdat(shot, strppf, strdata, uid);
% basic function to read data via MDSplus
% if you want to include sequence number put it in the strdata
% namely e.g. strdata  = 'jps/ti/105'
% CLA from template by Maslov
%

if ~exist('usrid');
    usrid = 'jetppf';
end
if isempty(usrid);
    usrid = 'jetppf';
end
%

if ~strcmp(usrid, 'jetppf');
    str =['_uid=ppfuid("' usrid '")'];
    mdsvalue(str);
end

%

%mdsconnect('mdsplus.jet.efda.org');

str=['_sig=jet("' strppf '/' strdata '",',int2str(shot),')'];
[y.data,status]=mdsvalue(str);
y.t = transpose(mdsvalue('dim_of(_sig,1)'));
[a b] = size(y.data);
if a > 1;
    y.r = mdsvalue('dim_of(_sig,0)');
end

%
if ~strcmp(usrid, 'jetppf');
    str =['_uid=ppfuid("jetppf")'];
    mdsvalue(str);
end

%mdsdisconnect;

end