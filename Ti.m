%this is the funtion file to caluclae Temp for given P 

function T = Ti(A,B,C,P)
    T = (B/(A-log(P)))-C;
end