function symvec = symmat2vec(S)
%  converts symmetric matrix to a vector over the lower triangle
n = size(S,1);
symvec = zeros(n*(n+1)/2,1);
m = 0;
for i=1:n
    for j=1:i
        m = m + 1;
        symvec(m) = S(i,j);
    end
end