function symmat = symvec2symmat(Svec)
%  converts symmetric matrix to a vector over the lower triangle
Svec = Svec(:);
m = size(Svec,1);
n = (-1 + sqrt(1 + 8*m))/2;
symmat = zeros(n);
m = 0;
for i=1:n
    for j=1:i
        m = m + 1;
        symmat(i,j) = Svec(m);
        symmat(j,i) = Svec(m);
    end
end