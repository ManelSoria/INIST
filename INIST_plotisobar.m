function INIST_plotisobar(dat,pv,varargin)
% INIST - 
% Interpolation of Nonideal Idiosyncratic Splendiferous Tables
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
% INIST_plotisobar(dat,p,color)
% INIST_plotisobar:  plot isobar vector
% dat: data
% p: isobar vector
% color (optional): isobar colors 
%
% examples:
% INIST_plotisobar('N2',[10,20],'k')

global IND

path = fileparts(which(mfilename));
addpath(genpath(path));

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
           if numel(varargin)>0
               plot(IND.(dat).isoP{ii}.s,IND.(dat).isoP{ii}.T,varargin{1});
           else
               plot(IND.(dat).isoP{ii}.s,IND.(dat).isoP{ii}.T);               
           end
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
           
           if numel(varargin)>0
               plot(S,T,varargin{1});
           else
               plot(S,T);
           end
           hold on
           ok = 1;
           break;
       end   
   end
   if ok == 0
      error('Isobar (%e) not found',pv(j)) 
   end
end

% plot saturation bell
plot(IND.(dat).sl,IND.(dat).Tsat,'r');
plot(IND.(dat).sv,IND.(dat).Tsat,'r');


title(IND.(dat).name);
xlabel('s (kJ/kgK)');
ylabel('T (K)');

grid

return

end

