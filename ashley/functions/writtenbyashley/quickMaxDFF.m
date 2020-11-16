function maxDFF = quickMaxDFF(img_stack)

F = mean(img_stack,3);
dF = bsxfun(@minus, img_stack,F);
dFF = bsxfun(@rdivide,dF,F);

maxDFF = max(dFF,[],3);
end