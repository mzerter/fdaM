function symvec = symmat2symvec(S)
%  converts symmetric matrix to a vector over the lower triangle
% if any(any(S ~= S'))
%     disp(S)
%     error('S is not symmetric');
% end
n = size(S,1);
symvec = zeros(n*(n+1)/2,1);
m = 0;
for i=1:n
    for j=1:i
        m = m + 1;
        symvec(m) = S(i,j);
    end
end