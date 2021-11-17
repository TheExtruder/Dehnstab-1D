function [Ke,Fe] = element(elementLength,integrationParameter,C)
% Elementroutine
%
%   Eingang: elementData - Zeile der Elementgeometriedatentabelle
%            elementData = [ x1 x2 x3 ... xn];
%            integrationParameter - Parameter Integrationsstuetzpunkte
%            C  - Materialmatrix



% Anzahl der Integrationspunkte
nGausspoints = integrationParameter;

% Fallunterscheidung - Stuetzpunkte, Gewichte
if nGausspoints == 1
    gausspoint1 = 0.0; 
    weight1 = 2.0;
    gausspoint = gausspoint1; 
    weight  = weight1;
    
elseif nGausspoints == 2
    gausspoint1 = 0.577350269189626; 
    weight1 = 1;
    gausspoint = [-gausspoint1; gausspoint1]; 
    weight  = [ weight1; weight1];    
        
elseif nGausspoints == 3
    gausspoint1 = 0.774596669241483; 
    gausspoint2 = 0;
    weight1 = 0.555555555555555; 
    weight2 = 0.888888888888888;
    gausspoint = [-gausspoint1; gausspoint2; gausspoint1];
    weight  = [ weight1; weight2; weight1];
   
else
    disp('Used number of intergratins poits not implemented')
    return
end

xi = gausspoint(:);

%% Lineare Formfunktionen
%   N = [   N1(xi1)     N2(xi1)
%           N1(xi2)     N2(xi2) ]
N(:,1)=(1-xi)/2;  
N(:,2)=(1+xi)/2;   

%% Elementsteifigkeitmatrizen/-lastvektoren aufbauen
% Allozieren
Ke = zeros(2, 2);
Fe = zeros(2, 1);
B = [-1/elementLength; 1/elementLength]; %   B = [   dN1dx  dN2dx ]
J = elementLength/2;
    detJ = det(J);

for i = 1:nGausspoints
    % Elementsteifigkeitsmatrix
    Ke(1,1) =  Ke(1,1) + B(1)'*C*B(1)*detJ*weight(i);
    Ke(2,2) =  Ke(2,2) + B(2)'*C*B(2)*detJ*weight(i);
    Ke(1,2) =  Ke(1,2) + B(1)'*C*B(2)*detJ*weight(i);
    Ke(2,1) =  Ke(1,2);
end

   

