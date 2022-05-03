function [aky, omavg, qheat, pflux, vflux, qmheat, pmflux] = read_gs2output(flnm, pth, opt)
%
% function [aky, omavg, qheat, pflux, vflux, qmheat, pmflux] = read_gs2output(flnm, [pth], opt)
%
% Input.
%       flnm :  filename
%       pth  :  path (optional: by default looks in user's directory ~/gs2run/output)
%       if pth == 1 looks in user's directory ~/gs2run/gs2archive/output)
% Output.
%       aky   : vector of k_theta values
%       omavg : matrix with complex frequency for each k_theta
%               omavg(:,1) : real frequency
%               omavg(:,2) : growth rate
%               same normalisations as in gs2 output !
%       qheat and pflux : matrix with heat and particle fluxes
%                         (same normalisations as in gs2)
%
% Calls read_gs2outtime to check time evolution if not converged if opt = 1 (default)
% Works for a single adiabatic species
%
% CLA 10.08.03
%

home = getenv('HOME')

if ~exist('pth')
    pth = [home '/gs2run/output/'];
end

if ~exist('opt'); opt = 1; end

totflnm = [pth flnm];
fptr = fopen(totflnm, 'r');
frewind(fptr);
%

iconv = 1;
yk = 0;
while yk == 0
    sss = fscanf(fptr,'%s', 1);
    if isempty(sss); break; end;
    if strcmp(sss, 'converged'), yk = 1; end;
end

if yk == 0
    disp([flnm ' not converged']);
    iconv = 0;
    if opt == 1
        ync = input('Do you want to check the output ? (0/1)  ');
        if ~ync
            aky = NaN;
            omavg = [NaN NaN];
            qheat = [NaN NaN];
            pflux = [NaN NaN];
            vflux = [NaN NaN];
            qmheat = [NaN NaN];
            pmflux = [NaN NaN];
            vmflux = [NaN NaN];
        else
            [aky, omavg, qheat, pflux, vflux, qmheat, pmflux, vmflux, ky] = read_gs2outtime(flnm, pth);
        end
    else
        aky = NaN;
        omavg = [NaN NaN];
        qheat = [NaN NaN];
        pflux = [NaN NaN];
        vflux = [NaN NaN];
        qmheat = [NaN NaN];
        pmflux = [NaN NaN];
        vmflux = [NaN NaN];
    end
end
%

if iconv
    ik = 0;
    yesik = 0;
    while yesik == 0
        % looks for locations of aky
        yk = 0;
        while yk == 0
            sss = fscanf(fptr,'%s', 1);
            if isempty(sss), yesik = 1; break; end;
            if strcmp(sss, 'aky='), yk = 1; end;
        end
        if yk == 1
            ik = ik+1;
            aky(ik) = fscanf(fptr, '%f', 1);
            zyk = 0;
            while zyk == 0
                sss = fscanf(fptr,'%s', 1);
                if strcmp(sss, 'omavg/(vt/a)='), zyk = 1; end;
            end
            omavg(ik,1) = fscanf(fptr, '%f', 1);
            omavg(ik,2) = fscanf(fptr, '%f', 1);
            zyk = 0;
            while zyk == 0
                sss = fscanf(fptr,'%s', 1);
                if strcmp(sss, 'qheat='), zyk = 1; end;
            end
            qheat(ik,1) = fscanf(fptr, '%f', 1);
            %qheat(ik,2) = fscanf(fptr, '%f', 1);
            aa = fscanf(fptr, '%f', 1);
            if ~isempty(aa); qheat(ik,2) = aa; end;
            aa = fscanf(fptr, '%f', 1);
            if ~isempty(aa); qheat(ik,3) = aa; end;
            zyk = 0;
            while zyk == 0
                sss = fscanf(fptr,'%s', 1);
                if strcmp(sss, 'pflux='), zyk = 1; end;
            end
            pflux(ik,1) = fscanf(fptr, '%f', 1);
            %pflux(ik,2) = fscanf(fptr, '%f', 1);
            aa = fscanf(fptr, '%f', 1);
            if ~isempty(aa); pflux(ik,2) = aa; end;
            aa = fscanf(fptr, '%f', 1);
            if ~isempty(aa); pflux(ik,3) = aa; end;
            %
            zyk = 0;ivflx=1; imgn=1;
            while zyk == 0
                sss = fscanf(fptr,'%s', 1);
                if isempty(sss); %no vflux in output
                    vflux(ik,1:length(qheat(ik,:)))=NaN;
                    ivflx=0;imgn=0; break;
                else;
                    if strcmp(sss, 'vflux='), zyk = 1;
                    elseif strcmp(sss, 'qmheat='), zyk = 1; ivflx=0;
                        vflux(ik,1:length(qheat(ik,:)))=NaN;end;
                end;end;
            if ivflx==1;
                vflux(ik,1) = fscanf(fptr, '%f', 1);
                %vflux(ik,2) = fscanf(fptr, '%f', 1);
                aa = fscanf(fptr, '%f', 1);
                if ~isempty(aa); vflux(ik,3) = aa; end;
            end;
            %
            if imgn==1;
                if ivflx ==1;
                    zyk = 0;
                    while zyk == 0;
                        sss = fscanf(fptr,'%s', 1);
                        if isempty(sss); %no qmheat in output
                            qmheat(ik,1:length(qheat(ik,:)))=NaN;
                            imgn=0; break;
                        else;
                            if strcmp(sss, 'qmheat='), zyk = 1;end;
                        end;end;end;end;
            if imgn==1;
                qmheat(ik,1) = fscanf(fptr, '%f', 1);
                %qmheat(ik,2) = fscanf(fptr, '%f', 1);
                aa = fscanf(fptr, '%f', 1);
                if ~isempty(aa); qmheat(ik,3) = aa; end;
            end;
            zyk = 0;
            while zyk == 0;
                sss = fscanf(fptr,'%s', 1);
                if isempty(sss); %no pmflux in output
                    pmflux(ik,1:length(qheat(ik,:)))=NaN;
                    imgn=0;break;
                else;
                    if strcmp(sss, 'pmflux='), zyk = 1;end;
                end;end;
            if imgn==1;
                pmflux(ik,1) = fscanf(fptr, '%f', 1);
                %pmflux(ik,2) = fscanf(fptr, '%f', 1);
                aa = fscanf(fptr, '%f', 1);
                if ~isempty(aa); pmflux(ik,3) = aa; end;
            else
                qmheat(ik,1:length(qheat(ik,:)))=NaN;
                pmflux(ik,1:length(qheat(ik,:)))=NaN;
            end;
            %
            %
        end %if yk
        %
    end % while yesik
    %
end % if iconv
fclose(fptr);
