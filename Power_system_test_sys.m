clc, close all, clear all
load('dynamicdata2.mat');


p = 500;   %% how much data form the case we take
p = length(M(:,1));
r = 100;

M = M(1:p,1:p);
L = L(1:p,1:p);
D = D(1:p,1:p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = round(L,5); 
for i = 1:p
    L(i,i) = 0;
    S = sum(L(i,:));
    L(i,i) = -S;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I = eye(p);
Z = zeros(p);

E = [I,Z;Z,M];
A = [Z,I;-L,-D];
B = [Z ; I];
%C = [I,Z];
C = [ones(1,p) , zeros(1,p)]/p;


[V,D,W] =  eig(full(A),full(E));
LL = diag(D);

% 
% a = 1;
% b = 100;
% %points  =  (linspace(b,a,2*r) + r/2 +.5)*1i - 0.2 ;  %% b
% LL = sort(LL);
% Indx = find(abs(LL) > 1e-10) ;
% points = -LL(Indx);
% points = points(1:2*r);
% abs(points)

a = 1;
b = 100;
%points  =  (linspace(b,a,2*r) + r/2 +.5)*1i - 0.2 ;  %% b
LL = sort(LL);

Indx =  find(imag(LL) > 0);
Pts = LL(Indx);
[Pts2, Indx]  = sort(abs(LL(Indx)));

Pts = Pts(Indx) ;
Pts = Pts(1:r)  ; %%  r complex points


M = max(abs(Pts))

Indx =  find(imag(LL) == 0);
Pts2 = LL(Indx);
Indx = find(abs(Pts2) > 1e-10 )
Pts2 = Pts2(Indx);
[Pts3, Indx]  = sort(abs(Pts2));
Pts2 = Pts2(Indx);
Indx = find( abs(Pts2) <=  M);
Pts2 = Pts2(Indx)       %%% some posible real points 

Nrp = length(Pts2);

if ( mod(Nrp,2) == 0 )
    Pts = Pts(1:(r-Nrp/2));
    points = - [Pts2 ; conj(Pts) ; flip(Pts) ];
else
    PtsA = Pts(1:(r-(Nrp-1)/2));
    PtsB = Pts(1:(r-(Nrp+1)/2));   %%% flip takes the conjugate 
    points = - [Pts2 ; conj(PtsB) ; flip(PtsA) ]
end

abs(points)

mu      = points(1:2:2*r) + 0.1;   
gamma   = points(2:2:2*r) + 0.1;


Espect =  figure('DefaultAxesFontSize',18)
plot(real(LL),imag(LL),'.k','markersize',15)
hold on
plot(real(mu),imag(mu),'.b','markersize',15)
plot(real(gamma),imag(gamma),'.r','markersize',15)
grid on
box on


[n,pp] = size(B)
[mm,n] = size(C)
Le = rand(mm,r);
Ri = rand(pp,r);
[Ar,Br,Cr,Er] = Loewner(A,B,C,E,mu,gamma,Le,Ri);
L = -Er;  

s1 = svd(L);
SVDPlot = figure('DefaultAxesFontSize',18)
semilogy(s1/s1(1),'linewidth',2)
grid on
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 2*p;

    R = zeros(N,r);
    O = zeros(r,N);

    for i = 1:r      
       R(:,i) = (gamma(i)*E-A)\B*Ri(:,i); 
       O(i,:) = (Le(:,i))'*C/(mu(i)*E-A);
    end

    S1_O = svd(O);
    S1_R = svd(R);

    
    SO_min = min(S1_O)
    SR_min = min(S1_R)
     
    Zk_Tol = 1e+80;
    TolS = 1e-50;

    
    
 Zk = 10000*ones(r,1);
    for i = 1:N
            for j = 1:r
                p_temp = LL(i);
                q_temp = mu(j);

                Nu = poly(p_temp);
                De = poly(q_temp);

                Zk_temp = max(abs(polyval(Nu,LL))./abs(polyval(De,LL))) / min(abs(polyval(Nu,mu))./abs(polyval(De,mu))) ;

                if (Zk_temp < Zk(1))
                    p(1) = p_temp;
                    q(1) = q_temp;
                    Zk(1)   = Zk_temp;
                end
            end
    end
        

        Zk(1) = min([1,Zk(1)]);
        
        p_temp = p(1);
        q_temp = q(1);
        Zk_0 = Zk(1);
        
        
 
                          
    for k = 2:r        
       clear p;
       clear q;
       clear PB;
       
       
       Zk(k) = Zk(k-1);
       p(1) = p_temp;
       q(1) = q_temp;
      
        for j = 2:k
            Nu = poly(p);
            De = poly(q);
            [M,Idx] = max(abs(polyval(Nu,LL))./abs(polyval(De,LL)));
            p(j) = LL(Idx);
            [m,Idx] = min(abs(polyval(Nu,mu))./abs(polyval(De,mu)));
            q(j) = mu(Idx);
            Nu = poly(p);
            De = poly(q);
            Zk_temp = max(abs(polyval(Nu,LL))./abs(polyval(De,LL))) / min(abs(polyval(Nu,mu))./abs(polyval(De,mu))) ;
            
            if (Zk_temp < Zk(k))        
                Zk(k) = Zk_temp;
            end
                           
        end
        
                
            x = [real(p),real(q),imag(p),imag(q)] ;            
            opt = optimset('MaxFunEvals',200000,'Display','off');         
            opt = optimset(opt, 'Algorithm','active-set');
            
            if ( RationalFunC(x,LL,mu) < Zk_Tol ) 
                xf = fminunc(@(x) RationalFunC(x,LL,mu),x,opt);
                nn = length(x);    
                rf = nn/4;
                px = xf(1:rf) + xf(2*rf+1:3*rf)*1i;
                qx = xf(rf+1:2*rf) + xf(3*rf+1:4*rf)*1i;

                 Nu_o = poly(px);
                 De_o = poly(qx);

                 Zk_opt = max(abs(polyval(Nu_o,LL))./abs(polyval(De_o,LL))) / min(abs(polyval(Nu_o,mu))./abs(polyval(De_o,mu)));

                if (Zk_opt < Zk(k))        
                    Zk(k) = Zk_opt;
                end    
            end

        
    end    
    Zk
    
    
    Zk2 = 10000*ones(r,1);
    for i = 1:N
            for j = 1:r
                p_temp = LL(i);
                q_temp = gamma(j);

                Nu = poly(p_temp);
                De = poly(q_temp);

                Zk_temp = max(abs(polyval(Nu,LL))./abs(polyval(De,LL))) / min(abs(polyval(Nu,gamma))./abs(polyval(De,gamma))) ;

                if (Zk_temp < Zk2(1))
                    p(1) = p_temp;
                    q(1) = q_temp;
                    Zk2(1)   = Zk_temp;
                end
            end
    end
        

        Zk2(1) = min([1,Zk2(1)]);
        
        p_temp = p(1);
        q_temp = q(1);
        Zk_0 = Zk2(1);
        
        
 
                          
    for k = 2:r        
       clear p;
       clear q;
       clear PB;
       
       
       Zk2(k) = Zk2(k-1);
       p(1) = p_temp;
       q(1) = q_temp;
      
        for j = 2:k
            Nu = poly(p);
            De = poly(q);
            [M,Idx] = max(abs(polyval(Nu,LL))./abs(polyval(De,LL)));
            p(j) = LL(Idx);
            [m,Idx] = min(abs(polyval(Nu,gamma))./abs(polyval(De,gamma)));
            q(j) = mu(Idx);
            Nu = poly(p);
            De = poly(q);
            Zk_temp = max(abs(polyval(Nu,LL))./abs(polyval(De,LL))) / min(abs(polyval(Nu,gamma))./abs(polyval(De,gamma))) ;
            
            if (Zk_temp < Zk(k))        
                Zk2(k) = Zk_temp;
            end
                           
        end
        
                
            x = [real(p),real(q),imag(p),imag(q)] ;            
            opt = optimset('MaxFunEvals',200000,'Display','off');         
            opt = optimset(opt, 'Algorithm','active-set');
            
            if ( RationalFunC(x,LL,gamma) < Zk_Tol ) 
                xf = fminunc(@(x) RationalFunC(x,LL,gamma),x,opt);
                nn = length(x);    
                rf = nn/4;
                px = xf(1:rf) + xf(2*rf+1:3*rf)*1i;
                qx = xf(rf+1:2*rf) + xf(3*rf+1:4*rf)*1i;

                 Nu_o = poly(px);
                 De_o = poly(qx);

                 Zk_opt = max(abs(polyval(Nu_o,LL))./abs(polyval(De_o,LL))) / min(abs(polyval(Nu_o,gamma))./abs(polyval(De_o,gamma)));

                if (Zk_opt < Zk2(k))        
                    Zk2(k) = Zk_opt;
                end    
            end

        
    end
    
    Zk2     
 
    
    UUP = 1e5*ones(1,r);   % inicial in pos          
    UP = 1e5*ones(1,r);   % inicial in pos       
    ConstV = cond(V)
    ConstW = cond(W) 
    s1(1)
                       
      for k= 1:r         
        for pos = (1+pp*k):r           
            if k == 1
                UP(pos) = 1*(S1_O(1))*(S1_R(1))*ConstV/s1(1);                
            else
                UP(pos) = Zk(k-1)*(S1_O(1))*(S1_R(1))*ConstV/s1(1);
            end
            for i = 0:max(find(S1_O > TolS))        %min(SwapIndx,length(S1_O))
                for j = 0:max(find(S1_R > TolS))    %min(SwapIndx,length(S1_R))
                    %[pos , pos+i+j]
                    if ((pos+i+j) <= r)
                        if k == 1
                            tempB = 1*(S1_O(1+i))*(S1_R(1+j))*ConstV/s1(1);
                        else
                            tempB = Zk(k-1)*(S1_O(1+i))*(S1_R(1+j))*ConstV/s1(1);
                        end                        
                        %tempB = Zk(k-1)*(S1_O(1+i))*(S1_R(1+j))*ConstV; 
                        if (tempB <= UUP(pos+i+j) ) %&& tempB > s1(pos+i+j)/s1(1))  %removing noisy data 
                            UUP(pos+i+j) = tempB;
                        end
                    end
                end
            end
        end
      end
      
    UUP2 = 1e5*ones(1,r);   % inicial in pos          
    UP2  = 1e5*ones(1,r);   % inicial in pos   
      
        for k= 1:r         
        for pos = (1+mm*k):r          
            if k == 1
                UP2(pos) = 1*(S1_O(1))*(S1_R(1))*ConstW/s1(1);                
            else
                UP2(pos) = Zk2(k-1)*(S1_O(1))*(S1_R(1))*ConstW/s1(1);                
            end
            for i = 0:max(find(S1_O > TolS))        %min(SwapIndx,length(S1_O))
                for j = 0:max(find(S1_R > TolS))    %min(SwapIndx,length(S1_R))
                    %[pos , pos+i+j]
                    if ((pos+i+j) <= r)
                        if k == 1
                            tempB = 1*(S1_O(1+i))*(S1_R(1+j))*ConstW/s1(1);                            
                        else
                            tempB = Zk2(k-1)*(S1_O(1+i))*(S1_R(1+j))*ConstW/s1(1);                            
                        end                        
                        %tempB = Zk(k-1)*(S1_O(1+i))*(S1_R(1+j))*ConstV; 
                        if (tempB <= UUP2(pos+i+j) ) %&& tempB > s1(pos+i+j)/s1(1))  %removing noisy data 
                            UUP2(pos+i+j) = tempB;
                        end
                    end
                end
            end
        end
      end

%%  plot UUP




hold on            
semilogy([1+pp*1:r],UUP(1+pp*1:r),'-r.','linewidth',2,'markersize',15)
semilogy([1+pp*1:r],UP(1+pp*1:r),'-g.','linewidth',2,'markersize',15)
semilogy([1+mm*1:r],UUP2(1+mm*1:r),':r.','linewidth',2,'markersize',15)
semilogy([1+mm*1:r],UP2(1+mm*1:r),':g.','linewidth',2,'markersize',15)
%semilogy([21:r],PBS,':b*','linewidth',2)   %% k + j  - 1  = 21 (since 1 extra comes from completing Z)
%legend('\sigma_j(L)/\sigma_1(L)','Aprox Bound','location','southwest')

         
figure(Espect)
F = gcf;
set(F,'PaperOrientation','landscape');
set(F, 'Position', get(0, 'Screensize'));
print(F,'~/Matlab/GaussMat/plots/Power/PS_spectraA','-dpdf','-fillpage')

figure(SVDPlot)
F = gcf;
set(F,'PaperOrientation','landscape');
set(F, 'Position', get(0, 'Screensize'));
print(F,'~/Matlab/GaussMat/plots/Power/PS_SVDA','-dpdf','-fillpage')
         
        


  AA1.Eig   = LL;
  AA1.Left  = mu;
  AA1.Right = gamma;
  AA1.IntPoints = points; 
  AA1.IniPos = 1+mm*1;
  AA1.FinPos  = r;
  AA1.CoroBnd = UUP2;
  AA1.ThmBnd  = UP2;
  AA1.SingVal = s1/s1(1);
  
  AA1.Ar = Ar;
  AA1.Br = Br;
  AA1.Cr = Cr;
  AA1.Er = Er;
  AA1.Le = Le;
  AA1.Ri = Ri;
  
  name = ['PowerSysTestDat_',num2str(N),'.mat'];
  save(name, '-struct', 'AA1');   
  
  
