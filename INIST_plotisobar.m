function INIST_plotisobar(dat,pv)
% INIST - 
% Interpolation of Nonideal Idiosyncratic Splendiferous Tables
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
% INIST_plotisobar:  plot isobar vector
% dat: data
% p: isobar
% Isobars must be in the database

global IND

try
    addpath('Database\')
catch
    error('Ups,... Database folder is not here pls download it')
end

try
    if isempty(IND) || ~isfield(IND,dat)  
        set = load(dat);
        IND.(dat) = set.(dat);
    end
catch
    error('%s not found',dat)
end



for j=1:length(pv)% plot isobar number j
   ok = 0;
   for  ii=1:length(IND.(dat).isoP)
       if IND.(dat).isoP{ii}.P==pv(j)
           plot(IND.(dat).isoP{ii}.s,IND.(dat).isoP{ii}.T);
           hold on
           ok = 1;
           break;
       end  
       if ii < length(IND.(dat).isoP) && IND.(dat).isoP{ii+1}.P > pv(j) && IND.(dat).isoP{ii}.P < pv(j)
           Snext = IND.(dat).isoP{ii+1}.s;
           Sprev = IND.(dat).isoP{ii}.s;
           Pnext = IND.(dat).isoP{ii+1}.P;
           Pprev = IND.(dat).isoP{ii}.P;
           S = (Snext - Sprev) ./ (Pnext -Pprev) .* (pv(j) - Pprev) + Sprev;
           Tprev = IND.(dat).isoP{ii}.T;
           Tnext = IND.(dat).isoP{ii+1}.T;
           T = (Tnext+Tprev)/2;
           
           plot(S,T);
           hold on
           ok = 1;
           break;
       end   
   end
   if ok == 0
      error('Isobar (%e) not found',pv(j)) 
   end
end

plot(IND.(dat).sl,IND.(dat).Tsat,'r');
plot(IND.(dat).sv,IND.(dat).Tsat,'r');


title(IND.(dat).name);
xlabel('s (kJ/kgK)');
ylabel('T (K)');

grid

return

end

