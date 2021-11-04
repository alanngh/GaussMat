function [fval,C,R] = f(c,A,t,x0)

    [m,n] = size(A);
      
    s = 1;
    for i = 1:n
        for j = i:n             
            C(i,j) = c(s);
            s = s+1;
        end
    end    
    C = C*C.';
    C = (C+C.')/2;

    H = lyap(A',C);
    %eig(H)
    R = chol(H);
            
    RARi = R*(A/R);
    muH = max(real(eig(RARi+RARi')/2));
    
    %Ri = inv(R);
    %muH = -(1/2)*min(eig(R'*C*R));

    fval = cond(R)*exp(t*muH)*norm(x0);
    
    %display([fval, c'])
%    fprintf('fval = %10.5f; c = %10.5e, %10.5e, %10.5e\n', fval, c(1), c(2), c(3))
