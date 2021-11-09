clc, close all, clear all
load('dynamicdata2.mat');




M = M(2:end,2:end);
L = L(2:end,2:end);
D = D(2:end,2:end);

[m,n]=  size(M);

I = eye(n);
Z = zeros(n);

%% i'm kalling K to matrix D in my pdf...



E = [I,Z;Z,M];
A = [Z,I;-L,-D];
B = [Z ; I];
C = [I,Z];

[n,mm] = size(B)
[pp,n] = size(C)
r = 3;

Le = ones(pp,r);
Ri = ones(mm,r);

Idn = eye(pp);
for i =1:r
    Le(:,i) = Idn(:,mod(i,r)+1);
    Ri(:,i) = Idn(:,mod(i,r)+1);
end

Tol = 1e-10;

[m,n] = size(A)
%E  = eye(n,n);

a = -1;
b =  10;
Vr = zeros(m,r);
Wr = zeros(m,r);
sig = logspace(a,b,r);
sig = sig.';

%%
%% need directions r and l

fprintf("\n computing V and W ...")
for  i = 1:r
    % multiply by directions
    Vr(:,i) = sparse(((sig(i)*E-A)\B)*Ri(:,i));
    Wr(:,i) = sparse(((sig(i)*E'-A')\(C'))*Le(:,i));
end

s_old = ones(r,1);
res = 1;
it = 0;


h = figure('DefaultAxesFontSize',18);

L  = eig(full(A),full(E));
plot(real(L),imag(L),'.r','markersize',20);
hold on


%tStart = tic;   
fprintf("\n Iterating IRKA ...\n")
while (res > Tol && it < 500)
    s_old = sig;
    Ar = sparse((Wr'*A)*Vr);
    Er = sparse((Wr'*E)*Vr);    
    Eps = 1e-12;
      
    [V,D,FLAG] = eigs((Ar),(Er)+Eps*eye(r,r),r);    
    sig = -diag(D);    
    %plot(real(sig),imag(sig),'ob','markersize',15);
    for  i = 1:1:r              
        % multiply by directions
        Vr(:,i) = sparse(((sig(i)*E-A)\B)*Ri(:,i));
        Wr(:,i) = sparse(((sig(i)*E'-A')\(C'))*Le(:,i));
    end
    it = it + 1;   
    s = svd(full(Er)); 
    %semilogy(s/s(1),'.-','linewidth',2,'markersize',20)
    I = ones(size(sig));
    res = min(abs(kron(sig,I) - kron(I,s_old)));  % i need a beeter stoping criteria
    %fprintf(" it = %d \t res = %.8f \n",it,res);    
    fprintf("  %d \t  %.4e \t %d \t %.4e \n",it,res,r,s(end)/s(1));        
    if (mod(it,5) == 0)
        sv = s/s(1);
        if (sv(end) > 1e-10)            
            sig = [sig ; sig(end)+10*1i];             
            %Le(:,r+1) = rand(pp,1);
            %Ri(:,r+1) = rand(mm,1);
            Le(:,r+1) = Idn(:,mod(it,r)+1);
            Ri(:,r+1) = Idn(:,mod(it,r)+1);
            
            Vr(:,r+1) = sparse((sig(i)*E-A)\B*Ri(:,r+1));
            Wr(:,r+1) = sparse((sig(i)*E'-A')\(C')*Le(:,r+1));
            r = r+1;
            s = [s;1];
        end
    end 
    hold on
    grid on
    box on
    
end
%tEnd = toc(tStart)


hold on
plot(real(sig),imag(sig),'.b','markersize',20);
grid on
box on

r
sig
s = svd(full(Er)); 
s/s(1)


Ar = sparse(Wr'*A*Vr);
Er = sparse(Wr'*E*Vr);
Br = sparse(Wr'*B);
Cr = sparse(C*Vr);

H  = @(s) C*((s*E-A)\B);
Hr = @(s) Cr*((s*Er-Ar)\Br);

sys  = dss(full(A),full(B),full(C),[],full(E));
sysR = dss(full(Ar),full(Br),full(Cr),[],full(Er));

figure('DefaultAxesFontSize',18)
h = sigmaplot(sys,'-r',sysR,'-.b'); 
ax = gca;
h = findobj(gca,'Type','line')

NL = length(h)
figure('DefaultAxesFontSize',18)
SG = (NL-1)/2; 
newcolors = distinguishable_colors(2*SG,[1 1 1;1 0 0;0 0 1; 0 0 0]);
for i = 1:(NL-1)/2
    x1 = h(i+1).XData;
    x2 = h(i+1+SG).XData;
    if length(x1) > 2
        y1 = h(i+1).YData;
        y2 = h(i+1+SG).YData;
                
        semilogx(x2,y2,'-','color',newcolors(i+SG,:),'linewidth',2)            
        hold on
        semilogx(x1,y1,'--','color',newcolors(i,:),'linewidth',2)
                        
    end
end
axis([ax.XLim , ax.YLim ]);
grid on
box on

figure
F = gcf;
set(F,'PaperOrientation','landscape');
set(F, 'Position', get(0, 'Screensize'));
print(F,'~/Matlab/GaussMat/plots/IRKA/Aprox1','-dpdf','-fillpage')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%r = 50;



%Le = rand(pp,r);
%Ri = rand(mm,r);

Vr = zeros(m,r);
Wr = zeros(m,r);
sig = logspace(a,b,r);
sig = sig.';


for  i = 1:r
    % multiply by directions
    Vr(:,i) = (((sig(i)*E-A)\B)*Ri(:,i));
    Wr(:,i) = (((sig(i)*E'-A')\(C'))*Le(:,i));
end

s_old = ones(r,1);
res = 1;
it = 0;

%tStart2 = tic;   
while (res > Tol && it < 500)
    s_old = sig;
    Ar = sparse((Wr'*A)*Vr);
    Er = sparse((Wr'*E)*Vr);    
    Eps = 1e-12;
   
    [V,D,FLAG] = eigs((Ar),(Er)+Eps*eye(r,r),r);
    
    sig = -diag(D);
    %plot(real(sig),imag(sig),'ob','markersize',15);
    for  i = 1:1:r              
        % multiply by directions
        Vr(:,i) = (((sig(i)*E-A)\B)*Ri(:,i));
        Wr(:,i) = (((sig(i)*E'-A')\(C'))*Le(:,i));
    end
    it = it + 1;   
    s = svd(full(Er)); 
    %semilogy(s/s(1),'.-','linewidth',2,'markersize',20)
    I = ones(size(sig));
    res = min(abs(kron(sig,I) - kron(I,s_old)));  % i need a beeter stoping criteria
    %fprintf(" it = %d \t res = %.8f \n",it,res);    
    fprintf("  %d \t  %.4e \t %d \t %.4e \n",it,res,r,s(end)/s(1));        
%     if (mod(it,2) == 0)
%         sv = s/s(1);
%         if (sv(end) > 1e-10)            
%             sig = [sig ; sig(end)+10*1i];             
%             Le(:,r+1) = rand(pp,1);
%             Ri(:,r+1) = rand(mm,1);
%             Vr(:,r+1) = sparse((sig(i)*E-A)\B*Ri(:,r+1));
%             Wr(:,r+1) = sparse((sig(i)*E'-A')\(C')*Le(:,r+1));
%             r = r+1;
%             s = [s;1];
%         end
%     end 
    hold on
    grid on
    box on
end


s = svd(full(Er));
s = s/s(1);
indx = length(find(s > Tol))
Ar = sparse((Wr(:,1:indx)'*A)*Vr(:,1:indx));
Er = sparse((Wr(:,1:indx)'*E)*Vr(:,1:indx));    

%tEnd2 = toc(tStart2)

% 
% tEnd
% tEnd2


%Ar = sparse(Wr'*A*Vr);
%Er = sparse(Wr'*E*Vr);
Br = sparse(Wr(:,1:indx)'*B);
Cr = sparse(C*Vr(:,1:indx));

%H  = @(s) C*((s*E-A)\B);
Hr2 = @(s) Cr*((s*Er-Ar)\Br);

%sys  = dss(full(A),full(B),full(C),[],full(E));
sysR2 = dss(full(Ar),full(Br),full(Cr),[],full(Er));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure('DefaultAxesFontSize',18)
h3 = sigmaplot(sys,'-r',sysR2,'-.b'); 
legend
ax = gca;
h3 = findobj(gca,'Type','line')

NL2 = length(h3)
figure('DefaultAxesFontSize',18)
SG2 = (NL2-1)/2; 
newcolors = distinguishable_colors(2*SG+2*SG2,[1 1 1;1 0 0;0 0 1; 0 0 0]);
for i = 1:(NL-1)/2
    x1 = h(i+1).XData;
    x2 = h(i+1+SG).XData;
    
    x3 = h3(i+1).XData;
    x4 = h3(i+1+SG2).XData;
    if length(x1) > 2
        y1 = h(i+1).YData;
        y2 = h(i+1+SG).YData;
        
        y3 = h3(i+1).YData;
        y4 = h3(i+1+SG2).YData;
                
        semilogx(x2,y2,'-','color',newcolors(i+SG,:),'linewidth',2)            
        hold on
        semilogx(x1,y1,'--','color',newcolors(i,:),'linewidth',2)
        semilogx(x3,y3,'o','color',newcolors(i+2*SG,:),'linewidth',2)
        semilogx(x4,y4,':','color',newcolors(i+2*SG+SG2,:),'linewidth',2)
                        
    end
end
axis([ax.XLim , ax.YLim ]);
grid on
box on
%legend('Original','Rank Estimation','ROM','Original','Rank Estimation','ROM','Original','Rank Estimation','ROM','Original','Rank Estimation','ROM','Original','Rank Estimation','ROM')

%%%%%%%%%%%%%%%%%%
hold on 

figure
F = gcf;
set(F,'PaperOrientation','landscape');
set(F, 'Position', get(0, 'Screensize'));
print(F,'~/Matlab/GaussMat/plots/IRKA/approximations','-dpdf','-fillpage')


figure
sigma(sys-sysR,'-r',sys-sysR2,'-b')
legend('Rank update','full')

figure
F = gcf;
set(F,'PaperOrientation','landscape');
set(F, 'Position', get(0, 'Screensize'));
print(F,'~/Matlab/GaussMat/plots/IRKA/erros','-dpdf','-fillpage')

%legend

% 
% t = logspace(-6,2,100);
% nn = length(t);
% y  = zeros(nn,1);
% y2 = zeros(nn,1);
% for i=1:nn
%     y(i)  = 20*log10(abs(Hr(t(i))*1i));
%     y2(i) = 20*log10(abs(H(t(i))*1i));
% end
% semilogx(t,y2,'-b',t,y,'*r','linewidth',2)
% grid on
% legend('H','Hr')

