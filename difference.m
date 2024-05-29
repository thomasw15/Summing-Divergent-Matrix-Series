function [regular_error, kahan_error] = difference(n,k)

    R = rand(n);
    OneMat = ones(n);
    
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
end    