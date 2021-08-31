% INIST - 
% Interpolation of Nonideal Idiosyncratic Splendiferous Tables
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
%
% Example:  Linde cycle for N2 liquefaction
% Multiple stage isentropic compression



error('THIS EXAMPLE HAS TO BE REVISED !!');

close all

plow=1; % bar
phigh=190; % bar

T1=300;
T3=164.5; % initial guess

mg=1 % kg/s


% point 1: cycle begins at ambient temperature and pressure

% Multiple stage compression
n=5; 

[~,qc,wc,Tvec,pvec,hvec,svec]=ns_isot_comp('N2',n,T1,plow,phigh);



T3=fsolve(@hastobezero,T3)

   


%%%%% NEED add_p
INIST_plotisobar(IND.N2,pvec);
%%%%%

% save('IND','IND');

% 
% plot(s1,T1,'ro');
% plot(s2,T2,'ro');
% plot(s3,T3,'ro');
% plot(s4,T4,'ro');
% plot(s5,T5,'ro');

T5

wc

ml


for i=2:length(pvec)
    plot(svec(i),Tvec(i),'ok')
end

% plot isoenthalpic expansion
ipres=linspace(p4,p3,40);
hh=h3;
is=[];
iT=[];
for i=1:length(ipres)
    if ipres(i)<INIST(IND.N2,'pcrit')
        hl=INIST('N2','hl_p',ipres(i));
        hv=INIST('N2','hv_p',ipres(i));
        sl=INIST('N2','sl_p',ipres(i));
        sv=INIST('N2','sv_p',ipres(i));
        xx=(hh-hl)/(hv-hl);
        is(i)=sl+xx*(sv-sl);
        iT(i)=INIST(IND.N2,'tsat_p',ipres(i));
    else
        eq=@(Tx) INIST('N2','h_pt',ipres(i),Tx)-hh;
        tt=fsolve(eq,T3,optimset('Display','none'));
        iT(i)=tt;
        is(i)=INIST('N2','s_pt',ipres(i),tt);
    end
end
plot(is,iT,'r')



 function e=hastobezero(T3)
        
        p1=plow;
        T1=Tvec(1);
        h1=hvec(1);
        s1=svec(1);

        % point 2: end of multiple stage compression
        p2=pvec(end); 
        T2=Tvec(end);
        h2=hvec(end);
        s2=svec(end);


        % point 3: end of regenerative cooling
        p3=p2;
        h3=INIST('N2','h_pt',p3,T3);
        s3=INIST('N2','s_pt',p3,T3);
        qi=mg*(h3-h2)


        % point 4: end of isoenthalpic expansion
        p4=plow;
        T4=INIST('N2','tsat_p',p4);
        h4=h3;
        hl=INIST('N2','hl_p',p4);
        hv=INIST('N2','hv_p',p4);
        sv=INIST('N2','sv_p',p4);
        sl=INIST('N2','sl_p',p4);

        x4=(h4-hl)/(hv-hl)
        s4=sl+x4*(sv-sl);
        if x4>1 
            error('Not wet steam, sorry ');
        end
        x4=(h4-hl)/(hv-hl)

        % point 5: end of separation
        T5=T4
        p5=plow;
        s5=sv;
        h5=hv;

        ml=mg*(1-x4)

        % heat that has to be released by 2->3 process
        abs(qi)

        % heat absorbed by 5->1 process
        (mg-ml)*(h1-h5)

        e=abs(qi)-(mg-ml)*(h1-h5);
    end

