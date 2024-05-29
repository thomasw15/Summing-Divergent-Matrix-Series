function [sum_Amq] = eqsum_simplified(m_max, q, Z)
    
    Z = double(Z);
    n = size(Z);
    II = double(eye(n));
    
    iter = 0;
    ratio = double(q/(1+q));
    rec = double(1/(1+q));
    Zq = double(ratio * II + Z * rec);
    Amq = double(II);
    sum_Amq = Amq;

    c = zeros(n);  %compensation for Kahan summation
    
    while iter < m_max

        iter = iter + 1;
        Amq = double(Amq * Zq);
        y = double(Amq - c);
        t = double(sum_Amq + y);
        %c = (t - sum_Amq) - y; %update compensation for Kahan summation
        c = double((t - sum_Amq) - y); 
        sum_Amq = double(t);

    end
    sum_Amq = double(sum_Amq* rec);
end