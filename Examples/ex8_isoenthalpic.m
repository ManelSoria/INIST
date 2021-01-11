% INIST - 
% Interpolation of Nonideal Idiosyncratic Splendiferous Tables
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
% Example: Plot N2 isoenthalpic lines

clearvars
close all


pvec=logspace(log10(1),log10(200),10); % a set of pressure values


INIST_plotisobar('N2',pvec); 

hold on;

T0=linspace(70,250,10);

for j=1:length(T0)
    % plot isoenthalpic expansion
    hh=INIST('N2','h_pt',pvec(end),T0(j));
    is=[];
    iT=[];
    for i=1:length(pvec)
        if pvec(i)<INIST('N2','pcrit')
            hl=INIST('N2','hl_p',pvec(i));
            hv=INIST('N2','hv_p',pvec(i));
            sl=INIST('N2','sl_p',pvec(i));
            sv=INIST('N2','sv_p',pvec(i));
            xx=(hh-hl)/(hv-hl);
            if xx<=1 && xx>=0
                sat=1;
                is(i)=sl+xx*(sv-sl);
                iT(i)=INIST('N2','tsat_p',pvec(i));              
                continue;
            end
        end
        eq=@(Tx) INIST('N2','h_pt',pvec(i),Tx)-hh;
        tt=fsolve(eq,T0(j),optimset('Display','none'));
        iT(i)=tt;
        is(i)=INIST('N2','s_pt',pvec(i),tt);
    end
    plot(is,iT,'-r')
end

axis( [2 6 50 300])