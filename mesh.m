function [edof,q] = mesh(a,nElements, elementLength)

% Allozieren
edof = zeros(nElements ,2);
q = zeros(nElements + 1, 1);
q(1:size(q, 1)) = a*ones(1, nElements+1)+(0:nElements)*elementLength;

% Zuweisung Knotenfreiheitsgrade <-> Element    
for i = 1:nElements
        edof(i,:) = [i i+1]; 
end




