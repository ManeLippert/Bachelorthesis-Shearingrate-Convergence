% list the GKW scan file available for a given project and print their description
%	function []=list_gkwscan(proj)
% Inputs:
%	proj:	project name (optional)
%		path for the input files is obtained from the gkwpath function with "proj" in argument
%		If you give proj='ALL', all the scans of all the projects are listed

function []=list_gkwscan(proj)

% default 
if ~exist('proj')
	proj='ALL';
end

if strcmp(proj,'ALL')

	D=dir(gkwpath('top'));
	I=find([D.isdir]);
	for ii= 1:length(I)
		if ~(strcmp(D(I(ii)).name,'.')|strcmp(D(I(ii)).name,'..'))
		  fprintf('******************** Project %s*******************\n',D(I(ii)).name)
		  list_gkwscan(D(I(ii)).name)
		  fprintf('******************************************************\n\n',D(I(ii)).name)
		end
	end

else

	flpth=gkwpath('scan',proj);
	fl=dir(flpth);
	for ii=1:length(fl),
		if ~isempty(strfind(fl(ii).name,'.mat')),
			load([flpth fl(ii).name],'sinfo');
			sss=[fl(ii).name '\t' sinfo.date];
			for jj=1:sinfo.nbvar
				sss=[sss '\t' sinfo.var.name{jj}];	
				if isfield(sinfo.var,'coupled')
				  for kk=1:sinfo.var.coupled(jj)
					sss=[sss '+' sinfo.var.cpl(jj).name{kk}];
				  end
				end
			end
			sss=[sss '\n' sinfo.comments '\n\n'];
            
            
			fprintf(sss)
		end
    end
    
    dir(gkwpath('time',proj))

end