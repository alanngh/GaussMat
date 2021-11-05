function example_butter(n,a,b,step) 
clc,  close all

[A,B,C,D] = butter(n,1,'s');
X0 = ones(n,1);

fprintf('\n Computing eigen decomposition...')
[V,D] =eig(A);       
L = diag(D)
    
if (abs(prod(L)) <= 1e-20)
    fprintf('\n Computing jordan decomposition...')
    [V,J] = jordan(A);        
end

NC = cond(V)
alpha = max(real(L));           %% spectral abscissa 
omega = max(real(eig(A+A')/2)); %% numerical abscissa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting ||exp(At)||

t = a:0.01:b;
m = length(t);

Bnd= zeros(m,1);

fprintf('\n Computing spectral bounds...')

for k = 1:m    
    Bnd(k) = NC*exp(alpha*t(k))*norm(X0);
    Y     = expm(A*t(k))*X0;
    NY(k) = norm(Y);
end


T = a:step:b;
m = length(T)
NP = m;

newcolors = distinguishable_colors(NP,[1 1 1;1 0 0; 0 0 1; 0 0 0]);
colororder(newcolors)


h1 = figure('DefaultAxesFontSize',18);
semilogy(t,NY,'-b',t,Bnd,'--r','markersize',20,'linewidth',2)
grid on
axis([0 5 1e-1 1e2])
%plot(t,NY,'color',newcolors(1,:),'markersize',20,'linewidth',2)
xlabel('t')
title('\kappa(R)  exp(\omega_H(A) t) || X_0 ||')


R = zeros(size(A));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lab = cell(m+1,1);
Lab{1} = 'Exact';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h2 = figure('DefaultAxesFontSize',18);
title('W(RAR^-^1)')
hold on
h(1) = plot(real(L),imag(L),'.k','markersize',15);
[h(2),CntL] = JohnsonAlg(A,500,[-3 3 -3 3],'C','-b') ;

hL = '';

fprintf('\n Computing Godunov bounds...')


for p = 1:m
    clear C;
    clear H;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % optimize C at time t
    if p == 1
        R0 = eye(size(A));
    else
        R0 = R;        
    end
    
     fprintf('\n optimizing C at t = %d ...',T(p))
    ok = 0;
    [C,R,ok] = C_Opti(A,T(p),X0,R0,CntL)
       
    Lab{p+1} = ['t = ',num2str(T(p),'%.2f')];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                             
    fprintf('\n Computing H for previous C...')
    H = lyap(A',C);        
    R = chol(H);
            
    if (ok == 1)
        RARi = R*(A/R);
        muH(p) = max(real(eig(RARi+RARi')/2));
        CR(p) = cond(R);

        GBnd_t = CR(p).*exp(t.*muH(p))*norm(X0);   
        figure(h1);
        hold on
        semilogy(t,GBnd_t,'linewidth',2,'color',newcolors(p,:)) ;
        grid on
    end

    TopN = max([max(NY),max(Bnd),1e+1]);
    TopN = num2str(TopN);
    ind  = strfind(TopN,'+');  	
    if (length(ind) > 0)
        Ax1 = TopN((ind+1):end);
        Ax1 = str2num(['1e+',num2str(Ax1)]);
    else
    	TopN = num2str(TopN);
    	ind = length(TopN);
        Ax1 = str2num(['1e+',num2str(ind-1)]);
    end

    TopN = min([min(NY),min(Bnd),1e-1]);
    TopN = num2str(TopN);
    ind = strfind(TopN,'-');
    
    if (length(ind) > 0)
    	Ax2 = TopN((ind+1):end);
    	Ax2 = str2num(['1e-',num2str(str2num(Ax2))]);
    else
        ind = strfind(TopN,'.');
        Ax2 = length(TopN((ind+1):end));
        Ax2 = str2num(['1e-',num2str(Ax2)]);
    end
    axis([a b Ax2 Ax1])
    
    if (ok==1)
        figure(h2);
        hold on
        h(2+p) = JohnsonAlg(R*A*inv(R),500,[-3 1 -2 2],'C','-',newcolors(p,:))    ;        
    end
end

figure(h2);      
hold on
h(2+ NP +1) = line([0,0], ylim,'color','k','LineStyle','--', 'LineWidth', 1); % Draw line for Y axis.
grid on


figure('DefaultAxesFontSize',18); 
plot(T,CR,'-*b','linewidth',2)
title(['\kappa(R)  vs \kappa(V) = ',num2str(cond(V))])
xlabel('t')
grid on


figure('DefaultAxesFontSize',18);
plot(T,muH,'-*r','linewidth',2)
xlabel('t')
title(['\omega_H(A) vs \alpha(A) = ',num2str(alpha), ' vs \omega(A) = ',num2str(omega)])
grid on

fprintf('\n\n')

figure(2)
F = gcf;
set(F,'PaperOrientation','landscape');
set(F, 'Position', get(0, 'Screensize'));
print(F,['~/Matlab/GaussMat/plots/bounds',num2str(n)],'-dpdf','-fillpage')

figure(3)
F = gcf;
set(F,'PaperOrientation','landscape');
set(F, 'Position', get(0, 'Screensize'));
print(F,['~/Matlab/GaussMat/plots/spectra',num2str(n)],'-dpdf','-fillpage')

figure(4)
F = gcf;
set(F,'PaperOrientation','landscape');
set(F, 'Position', get(0, 'Screensize'));
print(F,['~/Matlab/GaussMat/plots/plot1-',num2str(n)],'-dpdf','-fillpage')

figure(5)
F = gcf;
set(F,'PaperOrientation','landscape');
set(F, 'Position', get(0, 'Screensize'));
print(F,['~/Matlab/GaussMat/plots/plot2-',num2str(n)],'-dpdf','-fillpage')
