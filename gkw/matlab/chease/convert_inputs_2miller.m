function convert_inputs_hamada2miller(proj_in,proj_out,filename,hamada_path,eqm)
% convert_inputs_hamada2miller(proj_in,proj_out,filename*,eqm_path,eqm)
%
% Batch convert inputs from eqdsk or hamada to miller geom.
%
% Uses the hamada / eqdsk file specified in the gkw input, stripping "../../hamada"
% original chease gkw input file should not contain any of the miller parameters to avoid conflicts
% squareness and zmil are set to zero for now
%
% TO DO: Speedup by making hamada_miller only calculate a few surfaces of interest

if ~exist(eqm)
    eqm=='hamada'
end

if ~exist(hamada_path)
    hamada_path = gkwpath('hamada',proj_in)
end

files=dir([gkwpath('input',proj_in) filename]);

for i=1:length(files)
    
    if (files(i).isdir==1) continue; end
    
    disp(files(i).name);
    
    gkwin=read_gkwinput(files(i).name,proj_in);
    
    % get the miller parameters for the selected flux surface
    
    hamada_file=gkwin.GEOM.eqfile(14:end)  % This strips "../../hamada"
    %hamada_file=gkwin.GEOM.eqfile;
        
    if (eqm=='hamada')
      [dum mil]=hamada_miller(hamada_file,hamada_path,files(i).name,proj_in,'r',0,0,gkwin.GEOM.eps);
    elseif (eqm=='eqdsk')
      mil=eqdsk_miller(hamada_file,hamada_path,gkwin.GEOM.eps);    
      mil=mil.GEOM;
    else
      error('only hamada or eqdsk supported')  
    end    
        
    % modify the geom, remove conflicting params
    gkwin.GEOM.geom_type='miller';
    gkwin.GEOM=rmfield(gkwin.GEOM,'shat');
    gkwin.GEOM=rmfield(gkwin.GEOM,'q');
    gkwin.GEOM=rmfield(gkwin.GEOM,'eps');
    
    gkwin.GEOM=remove_miller(gkwin.GEOM)
    
    % add the miller params to the namelist
    gkwin.GEOM=mergestruct(gkwin.GEOM,mil);
    
    write_gkwinput(gkwin,files(i).name,proj_out,0)
    
end

end

function in=remove_miller(in)

in=rmfield(in,'kappa');
in=rmfield(in,'skappa');
in=rmfield(in,'delta');
in=rmfield(in,'sdelta');
%in=rmfield(in,'gradp');
in=rmfield(in,'gradp_type');
in=rmfield(in,'drmil');
in=rmfield(in,'dzmil');
in=rmfield(in,'zmil');
in=rmfield(in,'square');
in=rmfield(in,'ssquare');

end
