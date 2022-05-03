%Takes all the Phi3D files and recombines them to plot
%the whole tube, then outputs a hdf file for plotting
%in ViSit

%The number of s points must be changed here for correct
%Operation
%Interval is the number of steps between dumps
NumPs = 64;
Interval = 50;
NumFrames = 98;
maxstore = 0;
maxnum=0

for j=Interval:Interval:NumFrames*Interval
    for i=1:NumPs
        a{i}=load(sprintf('%s%4.4d%s%4.4d','PhiTD',i,'_',j));
        phy(:,:,i) = a{i};
    end
    b{j/Interval} = phy;
    maxn = max(max(max(phy)));
    if (maxn>maxnum)
        maxnum=maxn;
    end
    
   filename = ['Fluxtube',int2str(j/Interval),'.h5']
   hdf5write(filename,'/Phi',phy)
end
maxnum
maxdim = size(b{1});
for j=1:NumFrames
      j
      slice(b{j},0,0,[1,NumPs/8,NumPs/4,3*NumPs/8,NumPs/2,5*NumPs/8,3*NumPs/4,7*NumPs/8,NumPs]);shading interp;view(130,12);
      box on
      zlabel('Along field direction, s');
      xlabel('Radial direction');
      ylabel('Poloidal direction');colorbar('location','east');
      axis([1 maxdim(2) 1 maxdim(1) 1 maxdim(3) -1.1*maxnum 1.1*maxnum]);pbaspect([2 2 6]);
      M(j)=getframe;
      
end
movie2avi(M,'Fluxtube.avi','fps',10,'quality',100);
movie(M,2)
