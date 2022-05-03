#!/usr/bin/env octave

% This is a standalone GNU Octave script.


% You should add
%   addpath("~/gkw/octave");
% to your ~/.octaverc
% or whatever your GKW_HOME folder is, in order to use the
% files where often used functions are defined.
libgkwoct;

init()

set(0,'defaultfigurepaperposition',[0 0 660 450]*0.5)

%%%%%%%% LOAD DATA
pkg load hdf5oct

for dsetname = {'mat','mat_maxwll_background', 'mat_poisson','mat_field_diag'}
    rank = 0
    while(true)
      try
        dname=[dsetname{1},'_',num2str(rank,'%03d')]
        ## irs = reshape(irs,[],1);
        mat = complex(h5read("gkwdata.h5",["/diagnostic/diagnos_matrix/",dname,"_values_real"]),
	              h5read("gkwdata.h5",["/diagnostic/diagnos_matrix/",dname,"_values_imag"]));
        nmat=length(mat)
        jj = h5read("gkwdata.h5",["/diagnostic/diagnos_matrix/",dname,"_jj"]);
        length(jj)
        catch
          break
      end

        try
          irs = h5read("gkwdata.h5",["/diagnostic/diagnos_matrix/",dname,"_irs"]);

          nrows = length(irs)-1
          ## irs(end+1) = 0;
          [w,iw] = max(irs)

          row = [];
          col = [];
          val = [];

          for i = 1:nrows
            for j = irs(i):(irs(i+1)-1)
              %disp(['element',num2str(i),num2str(jj(j)),'is', num2str(mat(j))])
              try
                row(end+1) = i;
                col(end+1) = jj(j);
                val(end+1) = mat(j);
              catch
                lasterror
                exit
              end
            end
          end

        catch
          ii = h5read("gkwdata.h5",["/diagnostic/diagnos_matrix/",dname,"_ii"]);
          row = ii;
          col = jj;
          val = mat;
          lasterror
        end
        sparse_mat = sparse(row,col,val);

        %% display some information
        %spstats(sparse_mat)
        disp('Matrix type:')
        matrix_type(sparse_mat)


        ## spy(real(sparse_mat))
        ## print_all_formats(['sparse_',dname,'_real'])

        ## spy(imag(sparse_mat))
        ## print_all_formats(['sparse_',dname,'_imag'])

        spy(abs(sparse_mat))
        title(['rank ', num2str(rank)])
        filename = ['sparse_',dname,'_abs']
        print_all_formats(filename)
        mkdir('eps/')
        print(['eps/',filename,'.eps'], '-depsc','-color','-tight')

        rank++
      end
end
