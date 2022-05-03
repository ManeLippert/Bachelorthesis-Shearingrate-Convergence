%This routine will extract and plot the kxspec information from the spc output files
%As with mkmovie you must first define nframes.
%Run from a run folder with spcXXXXX files
%

for j= 1: nframes 
  i = j-1;

s = int2str(i);
  if (i<10) 
    name = strcat('spc00000',s);
  else 
    if (i<100) 
      name = strcat('spc0000',s);
    else
      if (i<1000) 
        name = strcat('spc000',s);
      else
        if (i<10000) 
          name = strcat('spc00',s);
        else
          if (i<100000) 
            name = strcat('spc0',s);
          else
            name = strcat('spc',s);
          end
        end
      end
    end
  end
  
  spc=load(name);

  kxspec(j,:) = sum(spc,1);
    
  plot(kxspec)
  
end
  
  