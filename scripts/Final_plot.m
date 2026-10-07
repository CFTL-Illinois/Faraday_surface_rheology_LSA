% ============================================================
% JFM_mechanism_plot.m
%
% Plot-only script for the JFM mechanism figures.
% No stability calculation is performed here.
%
% Required MAT file in the current folder:
%   JFM_final_data.mat
%
% Output: 600 dpi PNG figures only.
% ============================================================

clear; clc; close all;

%% Repository setup
this_file = mfilename('fullpath');
script_dir = fileparts(this_file);
repo_root = fileparts(script_dir);

data_dir   = fullfile(repo_root,'data');
figure_dir = fullfile(repo_root,'figures');

if ~exist(figure_dir,'dir')
    mkdir(figure_dir);
end

DATA_FILE = fullfile(data_dir,'JFM_final_data.mat');

if ~isfile(DATA_FILE)
    error(['Could not find JFM_final_data.mat. Run scripts/Final_data.m first, ' ...
           'or place the supplied JFM_final_data.mat in the data folder.']);
end


set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% ---------------- USER SETTINGS ----------------
SAVE_FIGURES = true;
MARK_EVERY   = 4;
LINE_W       = 1.35;
MARKER_SIZE  = 4.2;
FONT_SIZE    = 10.5;
AXIS_W       = 0.9;

%% ---------------- LOAD UNIFIED DATA ----------------
S = load(DATA_FILE);

THRESHOLD = S.THRESHOLD;
FW_CM     = S.FW_CM;
FW_CS     = S.FW_CS;
EIG_CM_P  = S.EIG_CM;
EIG_CS_P  = S.EIG_CS;
META      = S.META;

CM = get_results_cell(FW_CM);
CS = get_results_cell(FW_CS);
EIG_CM = get_eigen_cases(EIG_CM_P);
EIG_CS = get_eigen_cases(EIG_CS_P);

Omega = FW_CM.metadata.Omega;

ncase = min(numel(CM),numel(CS));
color_order = lines(max(ncase,7));
marker_list = {'o','s','d','^','v','>','<','p','h'};

%% ============================================================
% FIGURE 1: DENSE FARADAY ONSET THRESHOLD
% cols: C_M, C_S
%% ============================================================
fig0 = figure('Color','w','Units','inches','Position',[0.8 0.8 7.1 3.2]);
tl0 = tiledlayout(fig0,1,2,'TileSpacing','compact','Padding','compact');
ax01 = nexttile(tl0,1); hold(ax01,'on');
ax02 = nexttile(tl0,2); hold(ax02,'on');

for jc = 1:ncase
    Rc = THRESHOLD.CM_cases{jc};
    Rs = THRESHOLD.CS_cases{jc};
    c  = color_order(jc,:);
    mk = marker_list{mod(jc-1,numel(marker_list))+1};

    vm = isfinite(Rc.CM) & Rc.CM>0 & isfinite(Rc.ac);
    vs = isfinite(Rs.CS) & Rs.CS>0 & isfinite(Rs.ac);

    plot_with_sparse_markers(ax01,Rc.CM(vm),Rc.ac(vm),c,mk, ...
        LINE_W,MARKER_SIZE,MARK_EVERY,eta_label(Rc.mu2));
    plot_with_sparse_markers(ax02,Rs.CS(vs),Rs.ac(vs),c,mk, ...
        LINE_W,MARKER_SIZE,MARK_EVERY,'');

    % Exact case-dependent physical C_M limit at beta=beta_limit.
    if isfinite(Rc.CM_limit)
        xline(ax01,Rc.CM_limit,'--','Color',[1 0 0],'LineWidth',1.1,'HandleVisibility','off');
    end
end

for a=[ax01 ax02]
    format_axis(a,FONT_SIZE,AXIS_W);
    set(a,'XScale','log');
end

ylabel(ax01,'$a_c/g$');
xlabel(ax01,'$C_M(k_c)\,[\mathrm{Pa\,s\,m^{-1}}]$');
xlabel(ax02,'$C_S(k_c)\,[\mathrm{Pa\,s\,m^{-1}}]$');
panel_label(ax01,'(a)');
panel_label(ax02,'(b)');
lg0=legend(ax01,'Location','northoutside','Orientation','horizontal','NumColumns',min(ncase,3));
lg0.Layout.Tile='north';

if SAVE_FIGURES
    exportgraphics(fig0,fullfile(figure_dir,'JFM_threshold_CM_CS.png'),'Resolution',600);
end

%% ============================================================
% FIGURE 2: FREE WAVE AT FIXED CLEAN k_c
% rows: damping, dimensionless detuning
% cols: C_M, C_S
%% ============================================================
fig1 = figure('Color','w','Units','inches','Position',[0.8 0.8 7.1 5.0]);
tl = tiledlayout(fig1,2,2,'TileSpacing','compact','Padding','compact');

ax11 = nexttile(tl,1); hold(ax11,'on');
ax12 = nexttile(tl,2); hold(ax12,'on');
ax13 = nexttile(tl,3); hold(ax13,'on');
ax14 = nexttile(tl,4); hold(ax14,'on');

for jc = 1:ncase
    Rm = CM{jc};
    Rs = CS{jc};
    c  = color_order(jc,:);
    mk = marker_list{mod(jc-1,numel(marker_list))+1};

    vm = Rm.CM_fixed>0 & isfinite(Rm.gamma_fixed) & isfinite(Rm.detune_fixed);
    vs = Rs.CS_fixed>0 & isfinite(Rs.gamma_fixed) & isfinite(Rs.detune_fixed);

    plot_with_sparse_markers(ax11,Rm.CM_fixed(vm),Rm.gamma_fixed(vm)/(Omega/2),c,mk, ...
        LINE_W,MARKER_SIZE,MARK_EVERY,eta_label(Rm.mu2));
    plot_with_sparse_markers(ax13,Rm.CM_fixed(vm),Rm.detune_fixed(vm)/(Omega/2),c,mk, ...
        LINE_W,MARKER_SIZE,MARK_EVERY,'');

    plot_with_sparse_markers(ax12,Rs.CS_fixed(vs),Rs.gamma_fixed(vs)/(Omega/2),c,mk, ...
        LINE_W,MARKER_SIZE,MARK_EVERY,'');
    plot_with_sparse_markers(ax14,Rs.CS_fixed(vs),Rs.detune_fixed(vs)/(Omega/2),c,mk, ...
        LINE_W,MARKER_SIZE,MARK_EVERY,'');
end

% Case-dependent C_M limits for the fixed-k free-wave panels.
for jc=1:ncase
    Rm=CM{jc}; c=color_order(jc,:);
    if isfield(Rm,'CM_limit_fixed') && isfinite(Rm.CM_limit_fixed)
        xline(ax11,Rm.CM_limit_fixed,'--','Color',[1 0 0],'LineWidth',1.1,'HandleVisibility','off');
        xline(ax13,Rm.CM_limit_fixed,'--','Color',[1 0 0],'LineWidth',1.1,'HandleVisibility','off');
    end
end

format_axis(ax11,FONT_SIZE,AXIS_W); format_axis(ax12,FONT_SIZE,AXIS_W);
format_axis(ax13,FONT_SIZE,AXIS_W); format_axis(ax14,FONT_SIZE,AXIS_W);

set([ax11 ax12 ax13 ax14],'XScale','log');
ylabel(ax11,'$\gamma_d/(\Omega/2)$');
ylabel(ax13,'$|\omega_0-\Omega/2|/(\Omega/2)$');
xlabel(ax13,'$C_M(k_{c,0})\,[\mathrm{Pa\,s\,m^{-1}}]$');
xlabel(ax14,'$C_S(k_{c,0})\,[\mathrm{Pa\,s\,m^{-1}}]$');

panel_label(ax11,'(a)'); panel_label(ax12,'(b)');
panel_label(ax13,'(c)'); panel_label(ax14,'(d)');

lg = legend(ax11,'Location','northoutside','Orientation','horizontal','NumColumns',min(ncase,3));
lg.Layout.Tile = 'north';

if SAVE_FIGURES
    exportgraphics(fig1,fullfile(figure_dir,'JFM_freewave_fixed_k.png'),'Resolution',600);
end

%% ============================================================
% FIGURE 3: FREE WAVE ALONG ACTUAL FARADAY CRITICAL PATH
% rows: damping, dimensionless detuning
% cols: C_M, C_S
%% ============================================================

fig2 = figure('Color','w','Units','inches','Position',[0.8 0.8 7.1 5.0]);

tl2 = tiledlayout(fig2,2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

ax21 = nexttile(tl2,1); hold(ax21,'on');
ax22 = nexttile(tl2,2); hold(ax22,'on');
ax23 = nexttile(tl2,3); hold(ax23,'on');
ax24 = nexttile(tl2,4); hold(ax24,'on');

for jc = 1:ncase

    Rm = CM{jc};
    Rs = CS{jc};

    c  = color_order(jc,:);
    mk = marker_list{mod(jc-1,numel(marker_list))+1};

    vm = Rm.CM_path > 0 & ...
         isfinite(Rm.gamma_path) & ...
         isfinite(Rm.detune_path);

    vs = Rs.CS_path > 0 & ...
         isfinite(Rs.gamma_path) & ...
         isfinite(Rs.detune_path);

    % ---------------- C_M ----------------
    plot_with_sparse_markers( ...
        ax21, ...
        Rm.CM_path(vm), ...
        Rm.gamma_path(vm)/(Omega/2), ...
        c,mk,LINE_W,MARKER_SIZE,MARK_EVERY, ...
        eta_label(Rm.mu2));

    plot_with_sparse_markers( ...
        ax23, ...
        Rm.CM_path(vm), ...
        Rm.detune_path(vm)/(Omega/2), ...
        c,mk,LINE_W,MARKER_SIZE,MARK_EVERY,'');

    % ---------------- C_S ----------------
    plot_with_sparse_markers( ...
        ax22, ...
        Rs.CS_path(vs), ...
        Rs.gamma_path(vs)/(Omega/2), ...
        c,mk,LINE_W,MARKER_SIZE,MARK_EVERY,'');

    plot_with_sparse_markers( ...
        ax24, ...
        Rs.CS_path(vs), ...
        Rs.detune_path(vs)/(Omega/2), ...
        c,mk,LINE_W,MARKER_SIZE,MARK_EVERY,'');

end

% Case-dependent C_M limits along the actual critical branch.
for jc=1:ncase
    Rm=CM{jc}; c=color_order(jc,:);
    if isfield(Rm,'CM_limit_path') && isfinite(Rm.CM_limit_path)
        xline(ax21,Rm.CM_limit_path,'--','Color',[1 0 0],'LineWidth',1.1,'HandleVisibility','off');
        xline(ax23,Rm.CM_limit_path,'--','Color',[1 0 0],'LineWidth',1.1,'HandleVisibility','off');
    end
end

%% ---------------- AXIS FORMAT ----------------

axs = [ax21 ax22 ax23 ax24];

for a = axs
    format_axis(a,FONT_SIZE,AXIS_W);
    set(a,'XScale','log');
end

%% ---------------- LABELS ----------------

ylabel(ax21,'$\gamma_d/(\Omega/2)$');
ylabel(ax23,'$|\omega_0-\Omega/2|/(\Omega/2)$');

xlabel(ax23,'$C_M(k_c)\,[\mathrm{Pa\,s\,m^{-1}}]$');
xlabel(ax24,'$C_S(k_c)\,[\mathrm{Pa\,s\,m^{-1}}]$');

%% ---------------- PANEL LABELS ----------------

panel_label(ax21,'(a)');
panel_label(ax22,'(b)');
panel_label(ax23,'(c)');
panel_label(ax24,'(d)');

%% ---------------- LEGEND ----------------

lg2 = legend(ax21, ...
    'Location','northoutside', ...
    'Orientation','horizontal', ...
    'NumColumns',min(ncase,3));

lg2.Layout.Tile = 'north';

%% ---------------- EXPORT ----------------

if SAVE_FIGURES

    exportgraphics( ...
        fig2, ...
        fullfile(figure_dir,'JFM_freewave_critical_path.png'), ...
        'Resolution',600);

end

%% ============================================================
% FIGURE 4: COMBINED C_M / C_S EIGENVECTOR DIAGNOSTICS
% Row 1: Gamma', tau_M, u_s versus C_M
% Row 2: blank,  tau_S, u_s versus C_S
% One shared legend for all panels
%% ============================================================
fig3 = figure('Color','w','Units','inches','Position',[0.6 0.6 7.1 5.0]);

tl3 = tiledlayout(fig3,2,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

ax31 = nexttile(tl3,1); hold(ax31,'on');
ax32 = nexttile(tl3,2); hold(ax32,'on');
ax33 = nexttile(tl3,3); hold(ax33,'on');

ax34 = nexttile(tl3,4);
ax35 = nexttile(tl3,5); hold(ax35,'on');
ax36 = nexttile(tl3,6); hold(ax36,'on');


%% ============================================================
% ASYMPTOTIC SCALING TEST
%
% Expected:
%   Gamma' ~ C_M^(-1)
%   u_s    ~ C_M^(-1)
%   tau_M  ~ C_M^( 0)
%
%   u_s    ~ C_S^(-1)
%   tau_S  ~ C_S^( 0)
%
% Fit is performed over the largest decade:
%       C >= C_max/10
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('ASYMPTOTIC SCALING TEST\n');
fprintf('Fit range: largest decade of C_M or C_S\n');
fprintf('============================================================\n');

%% ---------------- C_M scaling ----------------

fprintf('\nC_M SWEEP\n');
fprintf(['eta2 [Pa s]       ', ...
         'p_Gamma      R2_Gamma      ', ...
         'p_us         R2_us         ', ...
         'p_tauM       R2_tauM\n']);

for jc = 1:numel(EIG_CM)

    R = EIG_CM(jc);

    if isfield(R,'run_CM')
        R = R.run_CM;
    end

    if ~isfield(R,'CM_actual')
        continue;
    end

    x = R.CM_actual(:);

    [pG,r2G] = fit_power_law_tail( ...
        x,R.Gamma_RMS_norm(:));

    [pu,r2u] = fit_power_law_tail( ...
        x,R.us_RMS_norm(:));

    [pt,r2t] = fit_power_law_tail( ...
        x,R.tauM_RMS_norm(:));

    mu2 = EIG_CM(jc).mu2;

    fprintf('% .3e      % .4f      %.4f        % .4f      %.4f        % .4f      %.4f\n', ...
        mu2,pG,r2G,pu,r2u,pt,r2t);

end


%% ---------------- C_S scaling ----------------

fprintf('\nC_S SWEEP\n');
fprintf(['eta2 [Pa s]       ', ...
         'p_us         R2_us         ', ...
         'p_tauS       R2_tauS\n']);

for jc = 1:numel(EIG_CS)

    R = EIG_CS(jc);

    if isfield(R,'run_CS')
        R = R.run_CS;
    end

    if ~isfield(R,'CS_actual')
        continue;
    end

    x = R.CS_actual(:);

    [pu,r2u] = fit_power_law_tail( ...
        x,R.us_RMS_norm(:));

    [pt,r2t] = fit_power_law_tail( ...
        x,R.tauS_RMS_norm(:));

    mu2 = EIG_CS(jc).mu2;

    fprintf('% .3e      % .4f      %.4f        % .4f      %.4f\n', ...
        mu2,pu,r2u,pt,r2t);

end


%% ---------------- C_M ROW ----------------
neig_cm = min(numel(EIG_CM),size(color_order,1));

for jc = 1:neig_cm

    R = EIG_CM(jc);
    if isfield(R,'run_CM'), R = R.run_CM; end
    if ~isfield(R,'CM_actual'), continue; end

    c  = color_order(jc,:);
    mk = marker_list{mod(jc-1,numel(marker_list))+1};
    x  = R.CM_actual(:);

    vg = isfinite(x) & x>0 & isfinite(R.Gamma_RMS_norm(:));
    vt = isfinite(x) & x>0 & isfinite(R.tauM_RMS_norm(:));
    vu = isfinite(x) & x>0 & isfinite(R.us_RMS_norm(:));

    lab = '';
    if isfield(EIG_CM(jc),'mu2')
        lab = eta_label(EIG_CM(jc).mu2);
    end

    plot_with_sparse_markers(ax31,x(vg),R.Gamma_RMS_norm(vg), ...
        c,mk,LINE_W,MARKER_SIZE,MARK_EVERY,lab);

    plot_with_sparse_markers(ax32,x(vt),R.tauM_RMS_norm(vt), ...
        c,mk,LINE_W,MARKER_SIZE,MARK_EVERY,'');

    plot_with_sparse_markers(ax33,x(vu),R.us_RMS_norm(vu), ...
        c,mk,LINE_W,MARKER_SIZE,MARK_EVERY,'');
end

% Case-dependent C_M limits for the eigenvector panels.
for jc=1:neig_cm
    R=EIG_CM(jc); c=color_order(jc,:);
    if isfield(R,'run_CM'), R=R.run_CM; end
    if isfield(R,'CM_limit') && isfinite(R.CM_limit)
        xline(ax31,R.CM_limit,'--','Color',[1 0 0],'LineWidth',1.1,'HandleVisibility','off');
        xline(ax32,R.CM_limit,'--','Color',[1 0 0],'LineWidth',1.1,'HandleVisibility','off');
        xline(ax33,R.CM_limit,'--','Color',[1 0 0],'LineWidth',1.1,'HandleVisibility','off');
    end
end

%% ---------------- C_S ROW ----------------
neig_cs = min(numel(EIG_CS),size(color_order,1));

for jc = 1:neig_cs

    R = EIG_CS(jc);
    if isfield(R,'run_CS'), R = R.run_CS; end
    if ~isfield(R,'CS_actual'), continue; end

    c  = color_order(jc,:);
    mk = marker_list{mod(jc-1,numel(marker_list))+1};
    x  = R.CS_actual(:);

    vt = isfinite(x) & x>0 & isfinite(R.tauS_RMS_norm(:));
    vu = isfinite(x) & x>0 & isfinite(R.us_RMS_norm(:));

    plot_with_sparse_markers(ax35,x(vt),R.tauS_RMS_norm(vt), ...
        c,mk,LINE_W,MARKER_SIZE,MARK_EVERY,'');

    plot_with_sparse_markers(ax36,x(vu),R.us_RMS_norm(vu), ...
        c,mk,LINE_W,MARKER_SIZE,MARK_EVERY,'');
end

%% ---------------- AXIS FORMAT ----------------
for a = [ax31 ax32 ax33 ax35 ax36]
    format_axis(a,FONT_SIZE,AXIS_W);
    set(a,'XScale','log');
end

xlabel(ax31,'$C_M(k_c)\,[\mathrm{Pa\,s\,m^{-1}}]$');
xlabel(ax32,'$C_M(k_c)\,[\mathrm{Pa\,s\,m^{-1}}]$');
xlabel(ax33,'$C_M(k_c)\,[\mathrm{Pa\,s\,m^{-1}}]$');

ylabel(ax31,'$\|\Gamma^\prime\|/\|\zeta\|$');
ylabel(ax32,'$\|\tau_M\|/\|\zeta\|$');
ylabel(ax33,'$\|u_s\|/\|\zeta\|$');

xlabel(ax35,'$C_S(k_c)\,[\mathrm{Pa\,s\,m^{-1}}]$');
xlabel(ax36,'$C_S(k_c)\,[\mathrm{Pa\,s\,m^{-1}}]$');

ylabel(ax35,'$\|\tau_S\|/\|\zeta\|$');
ylabel(ax36,'$\|u_s\|/\|\zeta\|$');

%% ---------------- BLANK LOWER-LEFT TILE ----------------
axis(ax34,'off');

text(ax34,0.50,0.56,'$C_S$ sweep', ...
    'Units','normalized', ...
    'Interpreter','latex', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontWeight','bold', ...
    'FontSize',FONT_SIZE+0.5);

text(ax34,0.50,0.40,'$\beta_s=0$', ...
    'Units','normalized', ...
    'Interpreter','latex', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',FONT_SIZE);

%% ---------------- PANEL LABELS ----------------
panel_label(ax31,'(a)');
panel_label(ax32,'(b)');
panel_label(ax33,'(c)');
panel_label(ax34,'(d)');
panel_label(ax35,'(e)');
panel_label(ax36,'(f)');

%% ---------------- ONE SHARED LEGEND ----------------
if neig_cm > 0
    lg3 = legend(ax31, ...
        'Location','northoutside', ...
        'Orientation','horizontal', ...
        'NumColumns',min(neig_cm,3));
    lg3.Layout.Tile = 'north';
end

%% ---------------- EXPORT ----------------
if SAVE_FIGURES
    exportgraphics(fig3,fullfile(figure_dir,'JFM_eigenvector_diagnostics_combined.png'),'Resolution',600);
end

fprintf('All final JFM figures regenerated successfully.\n');

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================
function results = get_results_cell(S)
if isfield(S,'results')
    results=S.results;
else
    f=fieldnames(S);
    results=S.(f{1});
end
if isstruct(results) && ~iscell(results)
    results=num2cell(results);
end
end

function cases = get_eigen_cases(S)
if isfield(S,'results')
    R=S.results;
else
    f=fieldnames(S); R=S.(f{1});
end
if isstruct(R) && isfield(R,'cases')
    cases=R.cases;
elseif isstruct(R)
    cases=R;
else
    error('Could not identify eigenvector diagnostic cases in %s.',inputname(1));
end
end

function Omega = get_omega(S1,S2)
Omega=nan;
if isfield(S1,'metadata') && isfield(S1.metadata,'Omega'), Omega=S1.metadata.Omega; end
if ~isfinite(Omega) && isfield(S2,'metadata') && isfield(S2.metadata,'Omega'), Omega=S2.metadata.Omega; end
if ~isfinite(Omega), Omega=2*pi*100; end
end

function plot_with_sparse_markers(ax,x,y,c,mk,lw,ms,markEvery,displayName)
x=x(:); y=y(:);
if isempty(x), return; end
idx=1:markEvery:numel(x);
if isempty(displayName)
    plot(ax,x,y,'-','Color',c,'LineWidth',lw,'HandleVisibility','off');
    plot(ax,x(idx),y(idx),mk,'Color',c,'MarkerSize',ms,'LineWidth',0.8,'HandleVisibility','off');
else
    plot(ax,x,y,'-','Color',c,'LineWidth',lw,'DisplayName',displayName);
    plot(ax,x(idx),y(idx),mk,'Color',c,'MarkerSize',ms,'LineWidth',0.8,'HandleVisibility','off');
end
end

function format_axis(ax,fs,lw)
set(ax,'FontSize',fs,'LineWidth',lw,'TickDir','out','Layer','top');
box(ax,'on');
grid(ax,'off');
end

function panel_label(ax,s)
text(ax,0.03,0.94,s,'Units','normalized','FontWeight','bold','FontSize',10.5, ...
    'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex');
end

function s = eta_label(mu)
e=floor(log10(mu)); c=mu/10^e;
s=sprintf('$\\eta_2=%.3g\\times10^{%d}\\;\\mathrm{Pa\\,s}$',c,e);
end

function [p,R2] = fit_power_law_tail(x,y)

x = x(:);
y = y(:);

% keep only valid positive values
valid = isfinite(x) & isfinite(y) & x>0 & y>0;

x = x(valid);
y = y(valid);

if numel(x) < 4
    p  = NaN;
    R2 = NaN;
    return;
end

% largest decade
xmax = max(x);
xmin_fit = xmax/10;

tail = x >= xmin_fit;

xf = x(tail);
yf = y(tail);

if numel(xf) < 4
    p  = NaN;
    R2 = NaN;
    return;
end

X = log10(xf);
Y = log10(yf);

coef = polyfit(X,Y,1);

p = coef(1);

Yfit = polyval(coef,X);

SSres = sum((Y-Yfit).^2);
SStot = sum((Y-mean(Y)).^2);

if SStot > 0
    R2 = 1 - SSres/SStot;
else
    R2 = NaN;
end

end
