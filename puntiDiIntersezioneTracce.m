% Carica i dati dal file CSV
percorsoFile = 'C:\Users\alexf\OneDrive\Desktop\Progetto_PCS_2024\Debug\tracce.csv'; % Cambia questo percorso con il percorso corretto
dati = readtable(percorsoFile);

% Estrai i dati delle tracce
X1 = dati.X1;
Y1 = dati.Y1;
Z1 = dati.Z1;
X2 = dati.X2;
Y2 = dati.Y2;
Z2 = dati.Z2;

% Plot delle tracce in 3D
figure;
hold on;
for i = 1:height(dati)
    plot3([X1(i) X2(i)], [Y1(i) Y2(i)], [Z1(i) Z2(i)], '-o');
end
hold off;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Visualizzazione Tracce');
view(3);
grid on;
axis equal;