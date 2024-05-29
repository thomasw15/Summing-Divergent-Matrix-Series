close all
clear
clc

pc = parcluster('local')
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

n = 1000;
cap = 5000;
regular_list = zeros(cap,1);
kahan_list = zeros(cap,1);

R = rand(n);
OneMat = ones(n);

parfor k = 4950:cap
    disp(k);    
    Target = (OneMat - R.^k)./(OneMat-R);
    
    term = ones(n);
    sum = ones(n);
    
    for i = 1:k-1
        term = term .* R;
        sum = sum + term;
    end
    
    term = ones(n);
    sum_kahan = ones(n);
    C = zeros(n);
    for i = 1: k-1
        term = term .*R;
        Y = term - C;
        T = sum_kahan + Y;
        C = (T- sum_kahan) -Y;
        sum_kahan = T;
    end
    
    regular_error = (norm(sum-Target,"fro"));
    kahan_error = (norm(sum_kahan-Target,"fro"));
    regular_list(k) = regular_error;
    kahan_list(k) = kahan_error;
end

filename = 'rec_comp2.mat';
save(filename)

hold on
fig = figure(1);
plot(regular_list,'LineWidth',1)
plot(kahan_list,'LineWidth',1)
title('Errors from summation methods','Interpreter','latex','Fontsize',16)
leg1 = legend('recursive', 'compensated','Fontsize',16,'Interpreter','latex','location','northwest');
set(leg1,'Interpreter','latex');
ylabel('error','Fontsize',16,'Interpreter','latex')
xlabel('length of series','Fontsize',16,'Interpreter','latex')
hold off

print(fig,'recursive_compensated2','-dpdf')