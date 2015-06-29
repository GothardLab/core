function x = makerow(x)
%N Killian 110927
tf = size(x,2)==1;
if tf
    x = x';
end