function convert_inputs(proj,filename,gs2_template)
% function convert_inputs(gkw_proj,gkw_filenames*,gs2_template)
%
% batch convert a set of gkw input files to gs2
% uses gkwpath.m for gkw
% gs2path is hard coded inside gkw2gs2_input (all files in same folder)

files=dir([gkwpath('input',proj) filename]);

for i=1:length(files)
    
  if (files(i).isdir==1) continue; end

  disp(files(i).name);
  gkw2gs2_input(files(i).name,proj,gs2_template);

  end

end
