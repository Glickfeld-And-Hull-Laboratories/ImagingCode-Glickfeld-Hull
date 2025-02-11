shallowToDeepGood=0;
%type in shallow to deep coordinates- column one is unit, column two is
%depth
shallowToDeepGood(:,3) = 2364.4-shallowToDeepGood(:,2);
shallowToDeepGood(:,4) = shallowToDeepGood(:,3)*sind(75);
shallowToDeepGood(:,4) = round(shallowToDeepGood(:,4));