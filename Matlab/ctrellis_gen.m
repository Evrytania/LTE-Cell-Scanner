%Writed by jxj.All rights reserved.
function trellis=ctrellis_gen(L,G,fbG)
if nargin == 2
    trellis = poly2trellis(L, G);
elseif nargin == 3
    if isempty(fbG)
        trellis = poly2trellis(L, G);
    else
        trellis = poly2trellis(L, G, fbG);
    end
else
    disp('Incorrect number of input arguments!!!');
    trellis = -inf;
    return;
end
% trellis=poly2trellis(ConstraintLength,CodeGenerator,FeedbackGenerator);
trellis.outputs=oct2dec(trellis.outputs);
for s=1:trellis.numStates
    [tmp,ttmp]=find(trellis.nextStates==s-1);
    tttmp(1:2:2*length(tmp))=tmp-1;
    tttmp(2:2:2*length(tmp))=ttmp-1;
    trellis.priorStates(s,1:2*length(tmp))=tttmp;
end;
trellis.L1=length(trellis.priorStates(1,:))./2;
for out=1:trellis.numOutputSymbols
    [tmp,ttmp]=find(trellis.outputs==out-1);
    tttmp(1:2:2*length(tmp))=tmp-1;
    tttmp(2:2:2*length(tmp))=ttmp-1;
    trellis.m1ivo(out,1:2*length(tmp))=tttmp;
end;
trellis.L2=length(trellis.m1ivo(1,:))./2;
