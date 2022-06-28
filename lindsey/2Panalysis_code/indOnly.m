function out = indOnly(A,ind);
    %extracts indices only from array
    indA=sub2ind(size(A),(1:size(A,1))',ind);
    out=A(indA);
end