function r = Gauss3Kronrod7
% Gauss-Kronrod (3,7) pair.

%   Copyright 2011-2020 The MathWorks, Inc.

r.Nodes = [ ...
    -0.9604912687080202,-0.7745966692414834,-0.4342437493468026, ...
    0, ...
    0.4342437493468026,0.7745966692414834,0.9604912687080202];
r.LowWeights = [0, 5/9, 0, 8/9, 0, 5/9, 0];
r.HighWeights = [ ...
    0.1046562260264672,0.2684880898683334,0.4013974147759622, ...
    0.4509165386584744, ...
    0.4013974147759622,0.2684880898683334,0.1046562260264672];
