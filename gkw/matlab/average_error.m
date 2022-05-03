function [av aver flmst] = fl_cal(in,nstart,nend) 
% function [av aver flmst] = fl_cal(in,nstart,nend) 
% 
% For a nonlinear run caculates the average of a signal 
% and determines the error on the average by deviding 
% the interval in several sub intervals over which the 
% average is performed. 
% 
% Input in      : the array of signal 
%       nstart  : starting point of the average 
%       nend    : end point of the average 
% Output av     : The averaged 
%        aver   : Error on the average
%        flstd  : Standard dev on the average


% Set the number of intervals 
ninter = 4; 

if (nend-nstart < ninter); 
    % The average 
    av = mean(in(nstart:nend));
    aver=NaN;
    flmst=NaN;    
    return; 
end

% The average 
av = mean(in(nstart:nend));

% The standard deviation 
flmst = std(in(nstart:nend));

% The real estimate of the error 
nstep = round((nend-nstart)/ninter);
nstart = round(nstart);

for i = 1:ninter-1 
  val(i) = mean(in(nstart+(i-1)*nstep:nstart+i*nstep)); 
end; 
val(ninter) = mean(in(nstart+(ninter-1)*nstep:nend)) ;
aver = min(std(val),flmst); 
