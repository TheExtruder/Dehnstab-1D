function q = solveq(K,F,dirichletBoundary)
%% Loesen des lineare Gleichungssystems unter Beruecksichtiung Dirichlet-Randbedingungen
%
% q = solveq(K,F,dirichletBoundarycondition)
%
% q: Loesungsvektor
% K: Steifigkeitsmatrix
% F: Lastvektor
% dirichletBoundary: Dirichlet-Freiheitsgrade
%
% dirichletBoundary =  [idxDiri1 valueDiri1
%                       dxDiri2 valueDiri2
%                               :              
%                       dxDiriN valueDiriN]

if ~isempty(dirichletBoundary)
    % Eingaben einlesen
    dofDiri = dirichletBoundary(:,1);   % 1
    qDiri   = dirichletBoundary(:,2);   % 0

    % Freiheitsgrade mit Dirichlet-Raender identifizieren
    nDof = numel(F);
    DoSolve = true(nDof,1);         % logisches Array
    DoSolve(dofDiri) = false;       % erster Eintrag = 0

    % Nach Freiheitsgraden, die nicht Dirichlet-Rand sind loesen, bei
    % anderen Werten Dirichletrand einsetzen
    q = zeros(size(K,1),1);         
    q(DoSolve) = K(DoSolve,DoSolve)\(F(DoSolve) - K(DoSolve,~DoSolve)*qDiri);
    q(~DoSolve) = qDiri;
else
    % Nach allen Werten loesen
    q = K\F;
end

end