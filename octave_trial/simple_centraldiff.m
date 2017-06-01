dr = 0.05; dr2 = dr*dr;
R1 = 1;
N = R1/dr;

A = sparse(N, N);
for i = 2:N-1
	A(i,i) = -2;
	A(i,i-1) = 1;
	A(i,i+1) = 1;
end
A(1,1) = -2;
A(1,2) = 1;
A(N,N) = -2;
A(N,N-1) = 1;

A = A/dr2;

b = ones(N,1);

plot(A\b)
