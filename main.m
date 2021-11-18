close all; clear all; clc;

disp('-------------------------------------------------------------------')
disp('-                           Dehnstab                              -') 
disp('-------------------------------------------------------------------')

%% Konstanten / Vorgaben

% Anzahl Elemente
disp('Geometry')
nElements = input('Enter number of elements: ');             % in x-Richtung
integrationParameter              = 1;           

% Geometrie
a   = 0;  b  = 4; 

%   (horizontaler) Lastvektor
load = 25; 
    
% Materialdaten
emod    = 300;                 % E-Modul
% Materialmatrix
C = emod;

%% Netzgenerator
% Netz generieren (Elementfreiheitsgradzuordnungstabelle, Geometrievektor)

elementLength =  b / nElements;
disp(elementLength);
[edof,q]       = mesh(a, nElements, elementLength);

% Geometriedatentabelle bzgl. Gebiet (1D Elemente)
[elementData]           = extract(edof,q);

% bzgl. Randpunkt
BoundaryEdofRight = edof(end , 2);
[BoundaryEdRight] = extract(BoundaryEdofRight,q);

%% Assembly
%   Allozieren / Speicher reservieren
nDof = size(q,1);
K    = sparse(nDof,nDof);
F    = sparse(nDof,1);

% Elementschleife
for i = 1:nElements
    [Ke,Fe] = element(elementLength,integrationParameter,C);
    [K]     = assem(K,Ke,edof(i,:));
    [F]     = assem(F,Fe,edof(i,:));
end

%% Randbedingungen einbauen
% Dirichlet-Rand
boundaryCondition1    = 0;
dirichletBoundary      = ones(1,2);
dirichletBoundary(1,2) = boundaryCondition1;
% Neumann-Rand
boundaryCondition2    = load;

% Randelement
    [Fe]    = boundaryCondition2;
    [F]     = assem(F,Fe(1),BoundaryEdofRight);


%% Gleichungsloeser
d = solveq(K,F,dirichletBoundary);

%% Darstellung Ausgangssituation
y = zeros(1, nElements + 1);
plot(q,y,'r-o','LineWidth',3);
hold on

%% Darstellung Gleichgewichtslage
q = q + d;                      % Update, Konfiguration Gleichgewichtslage
plot(q,y,'b-o','LineWidth',2);
hold off

[newElementData] = extract(edof,q); 

disp('Ausgangslage:'  )
disp(elementData);
disp('Verschiebungsvektor:'  )
disp(d);
disp('Gleichgewichtslage:'  )
disp(newElementData);


