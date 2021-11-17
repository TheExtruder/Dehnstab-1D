function elementData = extract(edof,q)
%% Extrahieren der Elementdatentabelle aus der Elementfreiheitsgradzuordnungsdatentabelle und dem Geometrievektor q
%
% elementData = extract(edof,q)
% 
% elementData:   Elementdatentabelle
% edof: Elementfreiheitsgradzuordnungsdatentabelle
% q:    Geometrievektor

elementData = q(edof')';

end