close all
clear
clc

pc = parcluster('local')
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

disp("This is mobius");

d = 1000;
upper = 7;
count = 7;
step = upper/count;
percent = 50;
x = zeros(percent,1);
for j = percent:-1:1
    x(j) = 1 - 0.01 * j;
end
delta = zeros(count,1);
for j = 1:(count)
    delta(j) = 2^(-step*j);
end
n = 1000;
normlist = zeros(percent,count);
A = rand(d);
A = A * A.';
B = rand(d) + 1i * rand(d);
B = (B - ctranspose(B))/2;
Z = A + B;
Z = Z/norm(Z,2);
parfor k = 1:percent
    disp(k);
    for j = 1:count
        Around = eye(d) + delta(j) * Z;
        S = lambert(d,n,x(k),Around);
        Snorm = norm(S);
        normlist(k,j) = Snorm;
    end
end

filename = 'mobius.mat';
save(filename);
figure(1)
hold on
plot(x,normlist(:,2),'LineWidth',1);
plot(x,normlist(:,3),'LineWidth',1);
plot(x,normlist(:,4),'LineWidth',1);
plot(x,normlist(:,5),'LineWidth',1);
plot(x,normlist(:,6),'LineWidth',1);
plot(x,normlist(:,7),'LineWidth',1);
leg2 = legend('$\delta=2^{-2}$','$\delta=2^{-3}$','$\delta=2^{-4}$','$\delta=2^{-5}$','$\delta=2^{-6}$','$\delta=2^{-7}$' ,'Fontsize',16,'Interpreter','latex');
set(leg2,'Interpreter','latex');
ylabel('norm','Fontsize',16,'Interpreter','latex')
xlabel('$x$','Fontsize',16,'Interpreter','latex')
hold off
print(figure(1),'mobius','-dpdf')

function mu = sng_mu(n)
    if n == 1, mu = 1; return; end
    a=factor(n);
    b=length(a);
    mu = (prod([0 a]-[a 0])~=0)*((2*(floor(b/2)==b/2))-1);
end
function S = lambert(d,n,x,A)
    [U,Sigma,V] = svd(A);
    sing = diag(Sigma);
    S = zeros(d);
    for k = 1:n
        w = k.^ sing;
        w2 = 1./ w;
        powerA = U * diag(w2) * V.';
        temp = (1-x)/(1-x^k) * k * x^k * sng_mu(k) * powerA;
        S = S + temp;
    end
end