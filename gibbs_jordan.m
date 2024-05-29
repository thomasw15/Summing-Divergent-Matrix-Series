close all
clear
clc

pc = parcluster('local');
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

dim = 1000;
div = 100;
parts = dim/div;
Freq = 200;
size = 800;

rate = rand(div,1);
signs = 2*randi(2,div,1)-3;
rate = rate.* signs;

one_vec = ones(1,parts-1);
super = diag(one_vec,1);
id = eye(dim);
idnorm = norm(id,"fro");

Xvals = linspace(-0.3*pi,0.3*pi,size);
regular_list=zeros(size, 1);
Cesaro_list=zeros(size, 1);
norm_list=zeros(size,1);

parfor i = 1:size
    disp(i);
    X = zeros(dim);
    for j = 1:10
        vals = rate(j)*Xvals(i)*ones(parts,1);
        Xtemp = diag(vals)+super;
        X((j-1)*parts+1:j*parts,(j-1)*parts+1:j*parts) = Xtemp;
    end
    for j = 11:div
        vals = rate(j)*Xvals(i)*ones(parts,1);
        Xtemp = diag(vals);
        X((j-1)*parts+1:j*parts,(j-1)*parts+1:j*parts) = Xtemp;
    end
    Xn = zeros(dim);
    sums = zeros(dim);
    Sigma = zeros(dim);
    c = zeros(dim);  %compensation for Kahan for sum
    d = zeros(dim);  %compensation for Kahan for Sigma

    for k = 1:Freq
        %Kahan summation for sum
        input = 2/(pi*k)*(1-(-1)^k)*funm(k*X,@sin);
        y = input - c;
        t = sums + y;
        c = (t - sums) - y;
        sums = t;
        
        %Kahan summation for Sima  
        z = sums - d;
        s = Sigma + sums;
        d = (s - Sigma) - z;
        Sigma = s;
    end

    Cesaro = Sigma/(Freq+1);

    regular_norm = norm(sums,"fro");
    regular_list(i)= regular_norm;
    Cesaro_norm = norm(Cesaro,"fro");
    Cesaro_list(i) = Cesaro_norm;
    norm_list(i) = idnorm;
end

filename = 'Jordan.mat';
save(filename)

hold on
fig = figure(1);
plot(Xvals,regular_list,'LineWidth',1)
plot(Xvals,Cesaro_list,'LineWidth',1)
plot(Xvals,norm_list,'LineWidth',1)
title('Norm of Approximations of $\mathrm{sign}(X)$','Interpreter','latex','Fontsize',16)
leg1 = legend('Truncated Fourier Series','Cesaro Summation');
set(leg1,'Interpreter','latex');
ylabel('norms','Fontsize',16,'Interpreter','latex')
xlabel('Eigenvalues','Fontsize',16,'Interpreter','latex')
hold off

print(fig,'Jordan','-dpdf')
