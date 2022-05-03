function [code reason] = gkw_why_stop(proj,sim_name)
%function [code reason] = gkw_why_stop(proj,sim_name)
% find the reason the code stopped from the GKW screen output file
%
% codes:
%   -3   Unstable: max_gr
%   -2   NaN in fluxes
%   -1   Abort
%    0   Converged linear run
%    1   Stable
%    2   NTIME complete
%    3   Walltime
%    4   External stop
%    5   Incomplete or running

[aborted str] = grep_for_string(proj,sim_name,'ABORT');
if (aborted)
   reason = str;   
   code = -1;
   return
end   

if grep_for_string(proj,sim_name,'Hit a NaN in fluxes: stop')
   reason = 'NaN in fluxes';
   code = -2;
   return
end  

if grep_for_string(proj,sim_name,'max_sec: stop')
   reason = 'walltime';
   code = 3;
   return
end    


if grep_for_string(proj,sim_name,'Internal: stop')
   if grep_for_string(proj,sim_name,'max_gr: Run unstable')
       reason = 'max_gr unstable';
       code = -3;   
       return              
   end   
   if grep_for_string(proj,sim_name,'Growth rate convergence reached: STOP')
       reason = 'Converged';
       code = 0;
       return
   end    
   if grep_for_string(proj,sim_name,'min_gr: Mode is stable: STOP') 
       reason = 'min_gr stable';
       code = 1;
       return
   end
end  

if grep_for_string(proj,sim_name,'External stop')
   reason = 'external';
   code = 4;
   return
end

if grep_for_string(proj,sim_name,'Run successfully completed')
   reason = 'NTIME complete ?';
   code = 2;
   return
end

% Otherwise, should be this 
reason = 'Incomplete or running';
code = 5;
return

end

