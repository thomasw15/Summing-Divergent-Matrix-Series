function F = SchurParlett(A,delta,m_max, q)
%This implementation follows that of Xiaobo Liu (2024). 
% mp-spalg (https://github.com/xiaobo-liu/mp-spalg), GitHub. 
% Retrieved November 15, 2024.
if nargin < 3 || isempty(delta)
    delta = 0.1;
end
ord = [];
[m,n] = size(A);
if  ~isfloat(A) || ~ismatrix(A) || m ~= n
   error(message('MATLAB:funm:InputDim'));
end

if isequal(A,triu(A))
   T = A; U = eye(n);
   diagT = diag(T);
else
   [U,T] = schur(A,'complex');
   diagT = diag(T);
end
if isequal(T,diag(diagT)) 
   F = U*diag(fun(diagT))*U';
   return
end
if isequal(delta,'norm')
    max_tii = norm(diagT,Inf);
    delta = 0.1*max_tii;
end
% Determine reordering of Schur form into block form.
if isempty(ord), ord = blocking(T,delta); end
[ord, ind] = swapping(ord);  % Gives the blocking.
ord = max(ord)-ord+1;        % Since ORDSCHUR puts highest index top left.
[U,T] = ordschur(U,T,ord);
m = length(ind);
%disp("the new T is");
%disp(T);


% Calculate F(T); this is where we change a bit
F = zeros(n);
for col=1:m
   j = ind{col};
   % F(j,j) = trim_diagpertub(T(j,j),fun); we replace this line with the
   % approximation by Euler
   F(j,j) = eqsum_simplified(m_max, q, T(j,j));
   for row=col-1:-1:1
      i = ind{row};
      if length(i) == 1 && length(j) == 1
         % Scalar case.
         k = i+1:j-1;
         temp = T(i,j)*(F(i,i) - F(j,j)) + F(i,k)*T(k,j) - T(i,k)*F(k,j);
         F(i,j) = temp/(T(i,i)-T(j,j));
      else
         k = cat(2,ind{row+1:col-1});
         rhs = F(i,i)*T(i,j) - T(i,j)*F(j,j) + F(i,k)*T(k,j) - T(i,k)*F(k,j);
         F(i,j) = sylv_tri(T(i,i),-T(j,j),rhs);
      end
   end
end
F = U*F*U';
if isreal(A) && norm(imag(F),1) <= 10*n*eps*norm(F,1)
   F = real(F);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = blocking(A,delta)
a = diag(A); n = length(a);
m = zeros(1,n); maxM = 0;
if nargin < 2 || isempty(delta), delta = 0.1; end
for i = 1:n
    if m(i) == 0
        m(i) = maxM + 1; % If a(i) hasn`t been assigned to a set
        maxM = maxM + 1; % then make a new set and assign a(i) to it.
    end
    for j = i+1:n
        if m(i) ~= m(j)    % If a(i) and a(j) are not in same set.
            if abs(a(i)-a(j)) <= delta
                if m(j) == 0
                    m(j) = m(i); % If a(j) hasn`t been assigned to a
                                 % set, assign it to the same set as a(i).
                else
                    p = max(m(i),m(j)); q = min(m(i),m(j));
                    m(m==p) = q; % If a(j) has been assigned to a set
                                 % place all the elements in the set
                                 % containing a(j) into the set
                                 % containing a(i) (or vice versa).
                    m(m>p) = m(m>p) -1;
                    maxM = maxM - 1;
                                 % Tidying up. As we have deleted set
                                 % p we reduce the index of the sets
                                 % > p by 1.
                end
            end
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mm,ind] = swapping(m)
mmax = max(m); mm = zeros(size(m));
g = zeros(1,mmax); h = zeros(1,mmax);
for i = 1:mmax
    p = find(m==i);
    h(i) = length(p);
    g(i) = sum(p)/h(i);
end
[~,y] = sort(g);
h = [0 cumsum(h(y))];
ind = cell(mmax,1);
for i = 1:mmax
    mm(m==y(i)) = i;
    ind{i} = h(i)+1:h(i+1);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = sylv_tri(T,U,B)
m = length(T);
n = length(U);
X = zeros(m,n);
opts.UT = true;
% Forward substitution.
for i = 1:n
    X(:,i) = linsolve(T + U(i,i)*eye(m), B(:,i) - X(:,1:i-1)*U(1:i-1,i), opts);
end
end
