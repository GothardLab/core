function x = makecol(x)
%N Killian 110927
tf = size(x,1)==1;
if tf
    x = x';
end