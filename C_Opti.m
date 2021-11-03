function [Cop,Rop] = C_Opti(A,t,x0,C0)
    [m,n] = size(A);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% optimizing with an alias for A and t    
    c_opt = zeros(n*(n+1)/2,1);
    fmin = 1e10;
    
    s = 1;
    c0 = zeros(n*(n+1)/2,1);
    for i = 1:n
        for j = i:n             
            c0(s) = C0(i,j);
            s = s+1;
        end
    end    
    
    
    %opt = optimoptions('fmincon')
    opt = optimset('MaxFunEvals',10000,'Display','on'); %% TO DO: increase tol  (double check!)
    % Set OptimalityTolerance to 1e-3
    %opt = optimoptions( 'OptimalityTolerance', 1e-12); 
    %opt.StepTolerance = 1e-12;
    %opt = optimset(opt, 'Algorithm','active-set');
    opt = optimset(opt,'TolFun',1e-16,'TolX',1e-16,'Algorithm','active-set','MaxFunEvals',100000000,'Display','off'); %% TO DO: increase tol  
    fprintf('\n runinng fminunc')
    c_opt = fminunc(@(x) f(x,A,t,x0),c0,opt);              
    fprintf('\n getting matrices obtained')
    [fmin,Cop,Rop] = f(c_opt,A,t,x0);
    
    for i = 1:15        
        fprintf('\n Iterating solution %d for fmin = %e' ,i,fmin)
        %c_opt = fminunc(@(x) f(x,A,t,x0), c_opt+0.5*rand(size(c_opt)),opt);               
        c_opt = fminunc(@(x) f(x,A,t,x0), c_opt,opt);               
        [f0,C0,R0] = f(c_opt,A,t,x0);
        if (f0 < fmin)
            fprintf('\n better optimal founded fmin = %e ...',f0)
            fmin = f0;
            Cop = C0;
            Rop = R0;
        end
    end

   for i = 1:20
        fprintf('\n Iterating solution %d for fmin = %e' ,i,fmin)
        c_opt1 = fminunc(@(x) f(x,A,t,x0), c_opt+0.1*rand(size(c_opt)),opt);
        %c_opt = fminunc(@(x) f(x,A,t,x0), c_opt,opt);
        [f0,C0,R0] = f(c_opt1,A,t,x0);
        if (f0 < fmin)
            fprintf('\n better optimal founded fmin = %e ...',f0)
            fmin = f0;
            Cop = C0;
            Rop = R0;
        end
    end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% unpacking c
%     s = 1;
%     for i = 1:n
%         for j = i:n             
%             C(i,j) = c_opt(s);
%             s = s+1;
%         end
%     end
%     C = C*C';
    
    
