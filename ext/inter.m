function [interactions] = inter(nFact)
%getints Get interactions matrix for anovan
%   Produces matrix of all combinations between anovan factors

set=[1:nFact];

% interactions=false(factorial(nFact),nFact);

r=1;
for k=1:nFact;
    clear combs
    combs = combnk(set,k);
    for c=1:size(combs,1)
        interactions(r,combs(c,:)) = 1;
        r=r+1;
    end
end

end

