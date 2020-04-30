% threshoding
function st = thresh(s,tau)
if abs(s) < tau
    st = 0;
elseif s <= -tau
    st = s + tau;
else
    st= s - tau;
end
