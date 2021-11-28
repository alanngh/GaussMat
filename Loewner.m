function [Ar,Br,Cr,Er] = Loewner(A,B,C,E,mu,gamma,Le,Ri)

p = length(mu);
q = length(gamma);

[N ,N] = size(A);
[N1,m] = size(B);
[n,N2] = size(C);

if (N ~= N1 || N ~= N2)
    printf('\nDimension mismatch in A,B,C.\n')
    return
else
    
    H  = @(s) C*((s*E-A)\B);
    
    %Le = ones(m,p);   %% left interpolation directions
    %Ri = ones(n,q);   %% right interpolation directions

    H_mu = zeros(p,m);
    H_gamma = zeros(n,q);

    for i = 1:p
        H_mu(i,:)    = (Le(:,i)')*H(mu(i));
    end
    for i = 1:q
        H_gamma(:,i) = H(gamma(i))*Ri(:,i);
    end

    for i = 1:p
        for j = 1:q
            L(i,j) = (H_mu(i,:)*Ri(:,j) - Le(:,i)'*H_gamma(:,j))/(mu(i)-gamma(j));
            Ls(i,j)= (mu(i)*H_mu(i,:)*Ri(:,j) - gamma(j)*(Le(:,i)')*H_gamma(:,j))/(mu(i)-gamma(j));            
        end
    end

    Er  = -L;
    Ar  = -Ls;
    Br  = H_mu;
    Cr  = H_gamma;
end