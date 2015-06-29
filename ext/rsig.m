%round to significant figures past decimal
function out = rsig(in,nsigfig)
out = round(in*10^nsigfig)/10^nsigfig;