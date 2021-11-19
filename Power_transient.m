clc, clear all, close all


load('dynamicdata2.mat');

[m,n]=  size(M);

I = eye(n);
Z = zeros(n);

%% i'm calling K to matrix D in my pdf...



E = [I,Z;Z,M];
A = [Z,I;-L,-D];
B = [Z ; I];
C = [I,Z];

[n,mm] = size(B)
[pp,n] = size(C)

At = E\A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%t = 0:0.1:30;
t = 0:0.01:5;
N = length(t);
y = zeros(N,1);

for i = 1:N
    [i, t(i)]
    y(i) = norm(expm(At*t(i)));
end

  AA1.y = t;
  AA1.y = y;
  save('ExpoBnd2.mat', '-struct', 'AA1');  



