function norms = W2norms(W,groups);
norms = zeros(length(groups),1);
for i=1:length(groups)
    norms(i) = norm(W(groups{i}));
end