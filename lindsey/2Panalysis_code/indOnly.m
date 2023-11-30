function out = indOnly(A,ind);
    %extracts indices only from array
    %A is a 2d vector, e.g. cells x direction
    %ind is list of indices to extract, e.g. direction for each cell
    %out is the value for those indices, e.g. for each cell
    indA=sub2ind(size(A),(1:size(A,1))',ind);
    out=A(indA);
end