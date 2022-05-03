#!/usr/bin/env octave

pkg load hdf5oct;

% define a function to be able to use the return statement to break out of the outer loop
function delete_all_3d_data()
  number_of_species = h5readatt('gkwdata.h5', '/', 'grid.number_of_species')
         
  fieldprefixes = {"diagnos_fields/Poten", "diagnos_fields/Spc3d", "diagnos_fields/Apc3d", "diagnos_fields/Apara", "diagnos_fields/Bpc3d", "diagnos_fields/Bpara"};

  for prefix = fieldprefixes
      disp(prefix{1})
      prefixes_flags.(prefix{1}) = true;
  end

  momentsprefixes = {};
  momentsprefixes_tmp = {"diagnos_moments/D3d", "diagnos_moments/E3d", "diagnos_moments/P3d"};
  
  for prefix = momentsprefixes_tmp
    for isp = 1:number_of_species
      momentsprefixes{end+1} = [prefix{1}, num2str(isp,"%02d")];
      disp(momentsprefixes{end})

      prefixes_flags.(momentsprefixes{end}) = true;
    end
  end

  i = input("Input the number to start from (default 1)")
  if(isempty(i) || ~isnumeric(i))
    i = 1
  endif
  
  while(true)
    for prefix = fieldprefixes
      try
        dsetname = ["/diagnostic/", prefix{1},num2str(i,"%08d")];
        h5delete("gkwdata.h5", dsetname)
        disp(["Deleted ", dsetname])

      catch
        disp(["Could not delete ", dsetname])
        prefixes_flags.(prefix{1}) = false
        if(~any(cell2mat((struct2cell(prefixes_flags)))))
          return
        end
      end
    end

    for prefix = momentsprefixes
      try
        dsetname = ["/diagnostic/", prefix{1},'_',num2str(i,"%06d")];
        h5delete("gkwdata.h5", dsetname)
        disp(["Deleted ", dsetname])
      catch
        disp(["Could not delete ", dsetname])
        prefixes_flags.(prefix{1}) = false
        if(~any(cell2mat((struct2cell(prefixes_flags)))))
          return
        end
      end
    end
    i = i + 1;
  end
end

delete_all_3d_data()

try
  disp("Calling h5repack to actually free diskspace:")
  commandline = ["h5repack gkwdata.h5 gkwdata.repacked.h5"];
  [status, commandoutput] = system (commandline);

  commandline = ["ls -lh gkwdata.h5 gkwdata.repacked.h5"];
  [status, commandoutput] = system (commandline);
  disp(commandoutput)

  if(yes_or_no("Do you want to replace gkwdata.h5 with the file gkwdata.repacked.h5 (WARNING: DESTROYED DATA CANNOT BE RECOVERED) ? "))
    rename("gkwdata.repacked.h5", "gkwdata.h5")
    disp("Done.")
  end
catch
  disp("Error when executing external programs")
  disp(lasterror.message)
end
  

