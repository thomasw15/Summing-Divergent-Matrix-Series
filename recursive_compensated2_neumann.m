close all
clear
clc

disp('rec100');

pc = parcluster('local')
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

n = 100;
cap = 5000;
lower = 1;
length = zeros(cap-lower,1);
regular_list = zeros(cap,1);
kahan_list = zeros(cap,1);
regular_list_backward = zeros(cap,1);
kahan_list_backward = zeros(cap,1);

R = rand(n);
R = 0.8 * R/norm(R,2);
OneMat = ones(n);
I = eye(n);

parfor k = lower:cap
    R = rand(n);
    R = 0.8 * R/norm(R,2);
    OneMat = ones(n);
    I = eye(n);
    length(k) = k;
    disp(k);    
    Target = (I - R^k) * inv(I-R);
    term = eye(n);
    sum = eye(n);
    for i = 1:k-1
        term = term * R;
        sum = sum + term;
    end
    term = eye(n);
    sum_kahan = eye(n);
    C = zeros(n);
    for i = 1: k-1
        term = term *R;
        Y = term - C;
        T = sum_kahan + Y;
        C = (T- sum_kahan) -Y;
        sum_kahan = T;
    end
    
    regular_error = norm(sum-Target,"fro");
    kahan_error = norm(sum_kahan-Target,"fro");
    regular_list(k) = regular_error;
    kahan_list(k) = kahan_error;
    regular_list_backward(k) = norm(sum*(I - R) - (I - R^k),"fro");
kahan_list_backward(k) = norm(sum_kahan*(I - R) - (I - R^k),"fro");
end

filename = 'rec_comp2_neumann_100.mat';
save(filename)

hold on
fig = figure(1);
plot(length,regular_list,'LineWidth',1)
plot(length,kahan_list,'LineWidth',1)
%title('Errors from summation methods','Interpreter','latex','Fontsize',16)
leg1 = legend('recursive', 'compensated','Fontsize',16,'Interpreter','latex','location','northwest');
set(leg1,'Interpreter','latex');
ylabel('error','Fontsize',16,'Interpreter','latex')
xlabel('length of series','Fontsize',16,'Interpreter','latex')
hold off

fig2 = figure(2);
hold on
plot(length,regular_list_backward,'LineWidth',1)
plot(length,kahan_list_backward,'LineWidth',1)
%title('Errors from summation methods','Interpreter','latex','Fontsize',16)
leg1 = legend('recursive', 'compensated','Fontsize',16,'Interpreter','latex','location','northwest');
set(leg1,'Interpreter','latex');
ylabel('backward error','Fontsize',16,'Interpreter','latex')
xlabel('length of series','Fontsize',16,'Interpreter','latex')
hold off

print(fig,'recursive_compensated2_neumann_100','-dpdf')