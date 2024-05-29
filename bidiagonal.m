function  [V,Vinv] = bidiagonal(n,p)
    a = 0.8*rand(n,1);
    for i = 1:n
        coin = rand;
        if coin > p
            a(i) = -a(i);
        end
    end
    b = 2*rand(n-1,1)-1;
    c = zeros(n-1,1);
    
    V = (diag(a)+diag(b,1)+diag(c,-1));

    a = 1-a;
    b = -b;
    Theta = zeros(n+1,1);
    Theta(1) = 1; Theta(2) = a(1);
    Phi = zeros(n+1,1);
    Phi(n+1) =1; Phi(n) = a(n);
    
    Vinv = zeros(n);
    for i = 2:n
        Theta(i+1) = a(i)*Theta(i) - b(i-1)*c(i-1)*Theta(i-1);
        j = n+1-i;
        Phi(j) = a(j)*Phi(j+1) - b(j)*c(j)*Phi(j+2);
    end
    
    for j = 1:n
        for i = 1:j-1
            Vinv(i,j) = (-1)^(i+j)* prod(b(i:j-1))* Theta(i)*Phi(j+1)/Theta(n+1);
        end
        Vinv(j,j) = Theta(j)*Phi(j+1)/Theta(n+1);
        for i = j+1:n
            Vinv(i,j) = (-1)^(i+j)* prod(c(j:i-1))* Theta(j)*Phi(i+1)/Theta(n+1);
        end
    end