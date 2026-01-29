%% ============================
%  CONFIGURACIÓN Y PARÁMETROS
% ============================
archivos = {
    '11062025.txt'
    '16062025.txt'
    '17062025.txt'
    '18062025.txt'
    '19062025.txt'
    '23062025.txt'
    '24062025.txt'
    '25062025.txt'
    '26062025V2.txt'
    '27062025.txt'
};

LAT  = 1; LON = 2; ALT = 3; VEL = 4; TIME = 5;

convert_to_kmh   = false;
umbral_stop_kmh  = 3;
umbral_acel_mps2 = 0.2;

% Parámetros físicos del vehículo de combustión
rho   = 0.907;
Cd    = 0.35;
A     = 2.9;
Cr    = 0.0115; 
m     = 2100;
g     = 9.81; 
delta = 1.05;
R     = 6371000;

%% ============================
%  PROCESAMIENTO DE TODAS LAS RUTAS
% ============================
n = numel(archivos);
v_list = cell(n,1);
t_list = cell(n,1);
Eacum_list = cell(n,1);
fecha_str = strings(n,1);

v_max_global = -Inf;
E_max_global = -Inf;
t_max_global = -Inf;

for k = 1:n
    file = archivos{k};
    D = readmatrix(file);
    
    v_kmh = D(:,VEL);
    if convert_to_kmh, v_kmh = v_kmh * 3.6; end
    v_ms = v_kmh / 3.6;

    t_raw = D(:,TIME);
    t_corr = double(t_raw);
    for i = 2:length(t_corr)
        if t_corr(i) < t_corr(i-1)
            t_corr(i:end) = t_corr(i:end) + 2^32;
        end
    end
    t_s = (t_corr - t_corr(1)) * 1e-6;
    t_min = t_s / 60;

    dt = diff(t_s); valido = dt > 0;
    v_seg = (v_kmh(1:end-1) + v_kmh(2:end)) / 2; v_seg = v_seg(valido);
    v_seg_ms = v_seg / 3.6;
    a_seg = diff(v_ms) ./ diff(t_s); a_seg = a_seg(valido);
    dt = dt(valido);

    lat = D(:,LAT); lon = D(:,LON);
    dphi = deg2rad(diff(lat)); dphi = dphi(valido);
    dl   = deg2rad(diff(lon)); dl = dl(valido);
    phi1 = deg2rad(lat(1:end-1)); phi1 = phi1(valido);
    phi2 = deg2rad(lat(2:end));   phi2 = phi2(valido);
    a = sin(dphi/2).^2 + cos(phi1).*cos(phi2).*sin(dl/2).^2;
    c = 2*atan2(sqrt(a), sqrt(1-a));
    dist_seg_m = R * c;

    % Energía total (sin regeneración)
    F_aero = 0.5 * rho * Cd * A .* v_seg_ms.^2 .* (v_seg_ms > 10/3.6);
    F_roll = Cr * m * g;
    F_iner = delta * m .* a_seg;
    F_pos = F_aero + F_roll + max(F_iner, 0);
    E_J = F_pos .* dist_seg_m;

    Eacum_kWh = cumsum(E_J) / 3.6e6;
    t_seg = t_min(2:end); t_seg = t_seg(valido);

    v_list{k} = v_seg;
    t_list{k} = t_seg;
    Eacum_list{k} = Eacum_kWh;

    v_max_global = max(v_max_global, max(v_seg));
    E_max_global = max(E_max_global, max(Eacum_kWh));
    t_max_global = max(t_max_global, max(t_seg));

    [~, nameNoExt] = fileparts(file);
    if strlength(nameNoExt) >= 8
        dia  = extractBetween(nameNoExt,1,2);
        mes  = extractBetween(nameNoExt,3,4);
        anio = extractBetween(nameNoExt,5,8);
        fecha_str(k) = sprintf('%s/%s/%s', char(dia), char(mes), char(anio));
    else
        fecha_str(k) = nameNoExt;
    end
end

%% ============================
%  GRAFICAR CICLOS Y ENERGÍA
% ============================
figure('Units','centimeters','Position',[2 2 22 4*n]);
tiledlayout(n,1,'Padding','compact','TileSpacing','compact');

for k = 1:n
    nexttile;

    yyaxis left
    plot(t_list{k}, v_list{k}, 'k-', 'LineWidth', 1.2);
    ylabel('Speed (km/h)', 'FontSize', 10);
    ylim([0, ceil(v_max_global/10)*10]);

    [v_max, idx_vmax] = max(v_list{k});
    t_vmax = t_list{k}(idx_vmax);
    text(t_vmax, v_max + 1, sprintf('Max V = %.1f km/h', v_max), ...
        'FontSize', 9, 'Color', 'k', 'HorizontalAlignment', 'left');

    yyaxis right
    plot(t_list{k}, Eacum_list{k}, 'r-', 'LineWidth', 1.2);
    ylabel('Energy (kWh)', 'FontSize', 10);
    ylim([0, ceil(E_max_global*10)/10]);

    [E_max, idx_emax] = max(Eacum_list{k});
    t_emax = t_list{k}(idx_emax);
    text(t_emax, E_max + 0.1, sprintf('Max E = %.2f kWh', E_max), ...
        'FontSize', 9, 'Color', 'r', 'HorizontalAlignment', 'left');

    xlim([0, ceil(t_max_global)]);
    grid on; box on;
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
    
    text(0.01*t_max_global, 0.85*v_max_global, ...
        fecha_str(k), 'FontWeight', 'bold', 'FontSize', 10);

    if k == n
        xlabel('Time (min)', 'FontSize', 10);
    else
        set(gca, 'XTickLabel', []);
    end
end

set(gcf,'Color','w');
exportgraphics(gcf,'CiclosVelocidad_Energia_Combustion.pdf','ContentType','vector');
disp('Gráfico exportado como "CiclosVelocidad_Energia_Combustion.pdf"');
