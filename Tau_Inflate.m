%% tau_inflate from time/altitude CSV (5–95% impulse duration)
% Set your filename:
fname = 'for_ethan.csv';

%% 1) Read CSV with '#' comments and preserve header if it's commented
try
    % First pass: try readtable with comments skipped
    T = readtable(fname, 'FileType','text', 'Delimiter',',', ...
        'CommentStyle','#', 'VariableNamingRule','preserve');
    varNames = string(T.Properties.VariableNames);

    % If expected names are missing (commented header got skipped), do a header recovery
    if ~all(ismember(["Time (s)","Altitude (m)"], varNames))
        % Scan file to find the commented header line
        fid = fopen(fname,'r'); assert(fid>0, 'Cannot open file.');
        hdrLine = ""; hdrIdx = 0; ln = 0;
        while true
            L = fgetl(fid); if ~ischar(L), break; end; ln = ln + 1;
            if startsWith(strtrim(L), '#') && contains(L, 'Time') && contains(L, 'Altitude')
                hdrLine = string(erase(strtrim(L), '#')); hdrIdx = ln;
                break
            end
        end
        fclose(fid);
        assert(hdrIdx>0, 'Could not locate a header line with Time/Altitude.');

        % Parse variable names
        V = strtrim(strsplit(hdrLine, ','));
        V = regexprep(V, '^\s*#\s*', '');  % remove leading '# ' if any

        % Build import options to start right after header
        opts = detectImportOptions(fname, 'FileType','text', 'Delimiter',',');
        opts.VariableNames = V;
        opts.VariableNamingRule = 'preserve';
        opts.DataLines = [hdrIdx+1, inf];   % data starts after header
        T = readtable(fname, opts);
        varNames = string(T.Properties.VariableNames);
    end
catch
    error('Failed to read %s', fname);
end

%% 2) Extract Time (s) and Altitude (m) columns (with fallbacks)
if all(ismember(["Time (s)","Altitude (m)"], varNames))
    t   = T.("Time (s)");
    alt = T.("Altitude (m)");
else
    % Fallback: pick likely columns by name, else first two columns
    names = lower(varNames);
    it = find(contains(names,"time") | names=="t", 1, 'first');
    ia = find(contains(names,"alt")  | contains(names,"height"), 1, 'first');
    if isempty(it) || isempty(ia)
        assert(width(T) >= 2, 'Need at least two columns for time and altitude.');
        it = 1; ia = 2;
    end
    t   = T{:, it};
    alt = T{:, ia};
end

% Column vectors, drop NaNs, sort by time
t   = t(:);  alt = alt(:);
keep = isfinite(t) & isfinite(alt);
t = t(keep); alt = alt(keep);
[ t, order ] = sort(t, 'ascend');
alt = alt(order);

dt = diff(t(:));
dt_med = median(dt);
fprintf('dt_med = %.6f s  (~%.1f Hz)\n', dt_med, 1/dt_med);

assert(numel(t) >= 5, 'Not enough samples after cleaning.');


%% 3) Smooth altitude, then differentiate to v and a
% Diagnostic
dt = diff(t(:));
dt_med = median(dt);
fprintf('dt_med = %.6f s  (~%.1f Hz)\n', dt_med, 1/dt_med);

% 1) Resample to a uniform grid (kills 1 ms spikes)
fs_target = 100;                        % Hz (adjust 50–200 as needed)
tq   = (t(1):1/fs_target:t(end)).';     % uniform time vector [s]
altq = interp1(t, alt, tq, 'pchip', 'extrap');

% 2) Smooth altitude before differentiating
span = max(11, 2*round(0.30*fs_target)+1);  % ~0.30 s window, odd
try
    alt_s = sgolayfilt(altq, 3, span);
catch
    alt_s = smoothdata(altq, 'sgolay', span);
end

% 3) Derivatives on the uniform grid (up is +)
v = gradient(alt_s, 1/fs_target);       % m/s
a = gradient(v,      1/fs_target);      % m/s^2

% 4) Use DESCENT only (after apogee)
[~, iApo] = max(alt_s);
t_d = tq(iApo:end);
a_d = a(iApo:end);

% 5) Deceleration magnitude during opening (upward accel while descending)
decel = max(0, a_d);
assert(any(decel>0), 'No deceleration bump detected after apogee.');

% 6) Robust event window: threshold + guards
[~, kpk] = max(decel);
peak     = decel(kpk);

a_floor     = 0.3;                       % m/s^2 floor vs noise (tune 0.2–0.5)
th          = max(0.05*peak, a_floor);
idx         = find(decel > th);
edges       = [0; find(diff(idx)~=1); numel(idx)];

min_dur_s   = 0.30;                      % ignore events shorter than this
min_samples = max(10, ceil(min_dur_s*fs_target));

best = struct('J',-inf,'i0',[],'i1',[]);
for b = 1:numel(edges)-1
    seg = idx(edges(b)+1 : edges(b+1));
    if numel(seg) < min_samples, continue; end
    if (t_d(seg(end)) - t_d(seg(1))) < min_dur_s, continue; end
    % Impulse per unit mass on this segment
    cum = [0; cumsum(0.5*(decel(seg(1:end-1))+decel(seg(2:end))).*diff(t_d(seg)))];
    Jtot = cum(end);
    if Jtot > best.J
        best.J = Jtot; best.i0 = seg(1); best.i1 = seg(end);
    end
end
assert(isfinite(best.J), 'No segment passed duration/sample guards.');

% 7) 5–95%% impulse duration and tau_inflate
seg   = best.i0:best.i1;
t_win = t_d(seg);
a_win = decel(seg);

cum   = [0; cumsum(0.5*(a_win(1:end-1)+a_win(2:end)).*diff(t_win))];
Jtot  = cum(end);
t05   = interp1(cum, t_win, 0.05*Jtot, 'linear', 'extrap');
t95   = interp1(cum, t_win, 0.95*Jtot, 'linear', 'extrap');

dt_05_95   = t95 - t05;
tau_inflate = dt_05_95 / 0.713;          % half-sine duration [s]
a_peak      = max(a_win);

%% Line-stretch speed just before the decel pulse
v_d = v(iApo:end);                 % velocity after apogee [m/s], up is +
pad = 0.20;                        % average over the 0.2 s before t05
t0  = max(tq(1), t05 - pad);
sel = (tq >= t0) & (tq <= t05);

% Two equivalent conventions:
V0_mag = mean(abs(v(sel)));       % speed magnitude [m/s]
V0_down = max(0, -mean(v(sel)));  % downward-positive speed [m/s]

DeltaV       = Jtot;
a_peak_tau   = (pi/2) * (DeltaV / tau_inflate);

fprintf('V0 ≈ %.2f m/s (|v|),  or  %.2f m/s (downward +)\n', V0_mag, V0_down)

fprintf('Δt_05–95 = %.3f s | tau_inflate = %.3f s | a_peak = %.3f m/s^2\n', ...
        dt_05_95, tau_inflate, a_peak);
fprintf('Window: [%.3f, %.3f] s | samples = %d | fs_target = %d Hz\n', ...
        t05, t95, numel(t_win), fs_target);

save('tau_outputs.mat', 'tau_inflate','DeltaV','a_peak','a_peak_tau','t05','t95','dt_05_95');