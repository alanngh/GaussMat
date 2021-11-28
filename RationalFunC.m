function  fval = RationalFunC(x,E,F)
    
    k = length(x)/4;    
    %x = sort(x);
    p = x(1:k)       + x((2*k+1):3*k)*1i;
    q = x((k+1):2*k) + x((3*k+1):4*k)*1i;
    Nu = poly(p);
    De = poly(q);
    fval = max(abs(polyval(Nu,E))./abs(polyval(De,E))) ./ min(abs(polyval(Nu,F))./abs(polyval(De,F))) ;
    
%   Phi1 = ones(size(E));
%   for j=1:k 
%     Phi1 = Phi1.*((E-p(j))./(E-q(j)));
%   end
%   Num = max(abs(Phi1));
%   Phi2 = ones(size(F));
%   for j=1:k 
%     Phi2 = Phi2.*((F-p(j))./(F-q(j)));
%   end
%   Den = min(abs(Phi2));  
%   fval = Num/Den;
end 