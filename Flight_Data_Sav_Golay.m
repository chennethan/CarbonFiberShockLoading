function Flight_Data_Sav_Golay()
% Savitzky–Golay smoothing of flight data (robust to '#' commented header)

fname = 'for_ethan.csv';

%% --- Parse the commented header to get real column names ---
fid = fopen(fname,'r'); assert(fid>0, 'Could not open %s', fname);
hdr = '';
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if startsWith(strtrim(ln), 'Time (s)')
        hdr = strtrim(ln);
        break;
    end
end
fclose(fid);
assert(~isempty(hdr), 'Header line not found (expected one starting with "# Time (s)").');

rawNames   = strtrim(strsplit(erase(hdr,'#'), ','));         % e.g., 'Time (s)', 'Altitude (m)', ...
validNames = matlab.lang.makeValidName(rawNames,'ReplacementStyle','delete');

%% --- Read data rows (skip '#' comments), then apply the header we parsed ---
opts = detectImportOptions(fname, 'CommentStyle','#', 'ReadVariableNames', false);
T    = readtable(fname, opts);

% Align widths (just in case) and apply names
ncol = min(width(T), numel(validNames));
T    = T(:, 1:ncol);
T.Properties.VariableNames = validNames(1:ncol);

% Helper to find a column by original (raw) header text
colIdx = @(pat) find(contains(rawNames(1:ncol), pat, 'IgnoreCase', true), 1, 'first');
idxT   = colIdx('Time (s)');
idxAlt = colIdx('Altitude');
idxAcc = colIdx('Vertical acceleration');

assert(~isempty(idxT) && ~isempty(idxAlt) && ~isempty(idxAcc), ...
    'Couldn''t find expected columns. Found: %s', strjoin(rawNames(1:ncol), ', '));

t     = T{:, idxT};
alt   = T{:, idxAlt};
aVert = T{:, idxAcc};

%% --- Savitzky–Golay settings ---
order   = 3;       % polynomial order (2–4 typical)
win_sec = 0.50;    % window length (seconds)

% Convert seconds -> odd frame length in samples (SG requires odd frame)
dt    = median(diff(t));
frame = max(order+2, round(win_sec/dt));
frame = 2*floor(frame/2) + 1;

% Smooth series
alt_s = sgolayfilt(alt,  order, frame);
a_s   = sgolayfilt(aVert, order, frame);

%% --- Plots: raw vs smoothed ---
figure('Color','w'); tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile;
plot(t, alt,  '-', 'DisplayName','Altitude (raw)'); hold on;
plot(t, alt_s,'-', 'LineWidth',1.5, 'DisplayName','Altitude (SG)');
xlabel('Time (s)'); ylabel('Altitude (m)'); grid on; legend('Location','best');
title(sprintf('Savitzky–Golay — frame %d (%.2fs), order %d', frame, frame*dt, order));

nexttile;
plot(t, aVert, '-', 'DisplayName','Vert accel (raw)'); hold on;
plot(t, a_s,   '-', 'LineWidth',1.5, 'DisplayName','Vert accel (SG)');
xlabel('Time (s)'); ylabel('Vertical accel (m/s^2)'); grid on; legend('Location','best');

%% --- Optional: zoom near descent decel peak ---
[~, apIdx] = max(alt);
decel      = max(0, -aVert(apIdx:end));
[~, kpk]   = max(decel);
iPeak      = apIdx + kpk - 1;
twin       = 3;  % seconds around the peak

mask = t >= (t(iPeak)-twin) & t <= (t(iPeak)+twin);
figure('Color','w');
plot(t(mask), aVert(mask), '-', 'DisplayName','Raw'); hold on;
plot(t(mask), a_s(mask),   '-', 'LineWidth',1.5, 'DisplayName','SG');
xlabel('Time (s)'); ylabel('Vertical accel (m/s^2)'); grid on; legend('Location','best');
title('Zoom: descent peak (raw vs. SG)');

end