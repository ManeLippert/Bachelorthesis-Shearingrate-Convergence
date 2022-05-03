% A function to read performance data from GKW perform.dat and make scaling graphs
% perform.dat is optional and is not required
%
% [struct rawdat] = read_perf(proj,filenames,file_ref,opttime)
%
% filenames is wildward string of file set
% file_ref is the name of the scaling reference file
% if opttime = 1, use perform instead of loop time for scaling

function [neat rawdat] = read_perf(proj,filename,file_ref,opttime)

files=dir([gkwpath('time',proj) filename]);
total = size(files)

if ~exist('opttime') opttime = 0; end

if (opttime==1)
  disp([gkwpath('perform',proj) file_ref])
  dum = importdata([gkwpath('perform',proj) file_ref]);
  ref_time = dum.data(1,2);
else
  ref_time=get_loop_time(proj,file_ref);
end

input=read_gkwinput(file_ref,proj,1);
ref_procs=input.GRIDSIZE.n_procs_x*input.GRIDSIZE.n_procs_sp*input.GRIDSIZE.n_procs_s*input.GRIDSIZE.n_procs_mu*input.GRIDSIZE.n_procs_vpar;

i=0

for ij = 1:total(1)
    
    if(files(ij).isdir==0 && files(ij).bytes > 0)
        
        i=i+1;
        
        input=read_gkwinput(files(ij).name,proj,1);
        procs = input.GRIDSIZE.n_procs_x*input.GRIDSIZE.n_procs_sp*input.GRIDSIZE.n_procs_s*input.GRIDSIZE.n_procs_mu*input.GRIDSIZE.n_procs_vpar;

        if (exist([gkwpath('perform',proj) files(ij).name],'file')==2)
          rawdat(i) = importdata([gkwpath('perform',proj) files(ij).name]);
        else
          continue
        end
              
        neat(i).procs=procs;
        if (opttime==1)        
         neat(i).time=rawdat(i).data(1,2);
        else
         neat(i).time=get_loop_time(proj,files(i).name);
        end
        neat(i).cputime=procs*neat(i).time;
        neat(i).speedup=ref_time*ref_procs/neat(i).time;
        neat(i).eff=neat(i).speedup/neat(i).procs;
        
              
        if (exist([gkwpath('perform',proj) files(ij).name],'file')==2)
            for j= 1:length(rawdat(i).data)
                field(j).name=regexprep(char(rawdat(i).textdata(j)),'[^\w'']','');
                field(j).name2=char(rawdat(i).textdata(j));
                neat(i).(field(j).name)=rawdat(i).data(j,2);
                neat(i).(['pc_' field(j).name])=rawdat(i).data(j,3);                
            end
        end
    end  
    
end

figure(101)
[sorted inds]=sort([neat(:).procs])
%[neat(inds).procs]
%[neat(inds).speedup]
%keyboard
loglog([neat(inds).procs],[neat(inds).speedup],'+-','DisplayName',proj)
hold all
xlabel('Cores')
ylabel('Speedup')
set(gca,'Xtick',unique([neat(:).procs]))
set(gca,'Ytick',unique([neat(:).procs]))
xlim([min([neat(:).procs]),max([neat(:).procs])+1000])
ylim([min([neat(:).procs]),max([neat(:).procs])+1000])
mm=max([neat(:).procs])
% perfect line
plot([neat(inds).procs],[neat(inds).procs],'k--');

figure(102)
semilogx([neat(inds).procs],[neat(inds).eff],'+-')
hold all
xlabel('Cores')
ylabel('Scaling efficiency')
set(gca,'Xtick',unique([neat(:).procs]))
xlim([min([neat(:).procs]),max([neat(:).procs])+1000])
ylim([0, 1.1])
set(gca,'Ytick',[0:0.1:1.0])
grid on

if (exist([gkwpath('perform',proj) files(1).name],'file')==2)
    figure    
    dum=rawdat(1).data(:,3);
    [sorted indx]=sort(dum,'descend');
    
    for j=1:length(rawdat(1).data)
      semilogx([neat(inds).procs],[neat(inds).(['pc_' field(indx(j)).name])],'+-','DisplayName',field(indx(j)).name2)       
      %doesn't work as indtended
      set(gca,'LineStyleOrder', '-+|-o|-^')  
      hold all
    end       
    
    set(gca,'Xtick',unique([neat(:).procs]))
    xlim([min([neat(:).procs]),max([neat(:).procs])+1000])
    
    xlabel('Cores')
    ylabel('Runtime percentage')  
    grid on    
end

end

function t=get_loop_time(proj,file)
    
    [err str]=system(['grep -i "(main loop time)" ' gkwpath('out',proj) file]);
    t=str2num(str(39:end));

end

