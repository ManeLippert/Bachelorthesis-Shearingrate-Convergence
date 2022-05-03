%This script defines a function which plots kyspec with a legend 
%The legend is taken from the krho file
%
%Run from a single run output folder
%Or if you are using a project kyspec folder first copy a krho file to it then run from there.
%
%Useage kyspec_label(krhomax,kyspec_file)
%
%Where kyspec_file is the kyspec filename (optional)
%
%You should provide krhomax since the values in file krho 
%have been normalised by kthnorm -a function of q and eps
%
%

function []=kyspec_label(krhomax,ky_file)

%Default values
if ~exist('krhomax')
	krhomax=1;
    disp('Used krhomax = 1')
end
if ~exist('ky_file')
    ky_file='kyspec';
end

set(0,'DefaultAxesLineStyleOrder','-|--|:')    
kyspec=load(ky_file);
plot(kyspec)
str=[ky_file '.krho'];

%load krho;
krho=load(str);
krho=krho*krhomax/krho(size(krho, 1));

for i=1:size(krho,1)
    labels(i)={num2str(krho(i,1))};
end

legend(labels)

end