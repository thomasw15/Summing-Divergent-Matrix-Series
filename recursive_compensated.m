close all
clear
clc

pc = parcluster('local')
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

k = 1000;
cap = 1000;
regular_list = zeros(cap,1);
kahan_list = zeros(cap,1);
regular_list_backward = zeros(cap,1);
kahan_list_backward = zeros(cap,1);

parfor n = 1:cap
    disp(n);
    R = randn(n);
    R = R/(5*norm(R,2));
    term = eye(n);
    sum_kahan = eye(n);
    C = zeros(n);
    for i = 1: k-1
        term = term * R;
        Y = term - C;
        T = sum_kahan + Y;
        C = (T- sum_kahan) -Y;
        sum_kahan = T;
    end

    term = eye(n);
    sum = eye(n);
    for i = 1:k-1
        term = term * R;
        sum = sum + term;
    end
    I = eye(n);
    inverse = (I - R^k) * inv(I -R);
    regular_list(n) = norm(sum-inverse,'fro');
    kahan_list(n) = norm(sum_kahan - inverse, 'fro');

    regular_list_backward(n) = norm(sum*(I -R) - (I - R^k),'fro');
    kahan_list_backward(n) = norm(sum_kahan*(I -R) - (I - R^k), 'fro');
end

filename = 'rec_comp_neumann.mat';
save(filename)

hold on
fig = figure(1);
plot(regular_list,'LineWidth',1)
plot(kahan_list,'LineWidth',1)
leg1 = legend('recursive', 'compensated','Fontsize',16,'Interpreter','latex','location','southeast');
set(leg1,'Interpreter','latex');
ylabel('forward error','Fontsize',16,'Interpreter','latex')
xlabel('matrix dimension','Fontsize',16,'Interpreter','latex')
hold off

fig2 = figure(2);
hold on
plot(regular_list_backward,'LineWidth',1)
plot(kahan_list_backward,'LineWidth',1)
leg1 = legend('recursive', 'compensated','Fontsize',16,'Interpreter','latex','location','southeast');
set(leg1,'Interpreter','latex');
ylabel('backward error','Fontsize',16,'Interpreter','latex')
xlabel('matrix dimension','Fontsize',16,'Interpreter','latex')
hold off

print(fig,'recursive_compensated_forward','-dpdf')
print(fig,'recursive_compensated_backward','-dpdf')