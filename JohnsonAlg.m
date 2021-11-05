function varargout = JohnsonAlg(A,n,ax,T,C,opt)
%
% A : n by n matrix 
% n : # of tangfent lines
% ax : axis of the plot
% T = 'T' iindicates tangent lines and T = 'C' indicates Countour or 
% T = 'B' for both
% T = 'L'Countouur in log log
% T = 'N' no plot
%
% returns h1 (legend handle)
% could return W countour points if needed 
Col = 1;

if (nargin < 4) 
    T='T'
    Col = 0;
end

if (nargin < 6) 
    Col = 0;
end

t = -1e6:0.5e6:1e6;
Theta = 0:pi/(n-1):pi;

for k =1:n
    r = exp(1i*Theta(k));
    B = r*A;
    H = (B+B')/2;
    [X,D] = eig(H);
    [lamMax,IndMax] = max(diag(D));
    [lamMin,IndMin] = min(diag(D));    
    vMax = X(:,IndMax);
    vMin = X(:,IndMin);
    W(k) = vMax'*A*vMax/(vMax'*vMax);
    W(k+n) = vMin'*A*vMin/(vMin'*vMin);        
    if ( (T == 'T') || (T == 'B') )
        L1 = exp(-1i*Theta(k)).*(lamMax + t*1i);
        L2 = exp(-1i*Theta(k)).*(lamMin + t*1i);
        plot (real(L1),imag(L1),'--k',real(L2),imag(L2),'--k');
        hold on
    elseif (T == 'L')
        L1 = exp(-1i*Theta(k)).*(lamMax + t*1i);
        L2 = exp(-1i*Theta(k)).*(lamMin + t*1i);
        loglog(real(L1),imag(L1),'--k',real(L2),imag(L2),'--k');
    end
end



if(T=='C')
    if (Col == 0)
        h1 = plot (real(W),imag(W),C,'linewidth',2);        
       
    else
        h1 = plot (real(W),imag(W),C,'linewidth',2,'color',opt);
    end
    
elseif (T=='B')
    if (Col == 0)
        h1 = plot (real(W),imag(W),C,'linewidth',2);
    else
        h1 = plot (real(W),imag(W),C,'linewidth',2,'color',opt);
    end
elseif (T == 'L')
    if (Col == 0)
        h1 = loglog(real(W),imag(W),C,'linewidth',2);
    else
        h1 = loglog(real(W),imag(W),C,'linewidth',2,'color',opt);
    end
elseif (T == 'N')
    %% no plot   
    h1 = [];
end
axis(ax);


    for k = 1: max(nargout,1) 
        if k == 1
            varargout{1} = h1;
        elseif k ==2
            W = W';
            np = length(W);
            S = 0;
            for i = 2:np    
                S = S + norm([real(W(i)),imag(W(i))] - [real(W(i-1)),imag(W(i-1))]);    
            end
            S = S + norm([real(W(1)),imag(W(1))] - [real(W(np)),imag(W(np))]);
            varargout{2} = S;
        end
    end    
end
%axis equal

