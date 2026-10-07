% ============================================================
% JFM_mechanism_data.m
% Compute and SAVE all data needed for final JFM mechanism figures.
% No figures are generated here.
%
% Outputs:
%   JFM_freewave_CM_data.mat
%   JFM_freewave_CS_data.mat
%   JFM_CM_eigenvector_diagnostics_data.mat
%   JFM_CS_eigenvector_diagnostics_data.mat
%
% Required external solver functions:
%   compute_eigenvalues.m
%   matrix_An.m
%   banded_matrix_B.m
%   fhs_general.m, fhs_general_deep.m
%   fhs_special.m, fhs_special_deep.m
% ============================================================
clear; clc; close all;

%% Repository setup
this_file = mfilename('fullpath');
script_dir = fileparts(this_file);
repo_root = fileparts(script_dir);

src_dir  = fullfile(repo_root,'src');
data_dir = fullfile(repo_root,'data');

addpath(src_dir);

if ~exist(data_dir,'dir')
    mkdir(data_dir);
end

required_files = { ...
    'compute_eigenvalues.m', ...
    'matrix_An.m', ...
    'banded_matrix_B.m', ...
    'fhs_general.m', ...
    'fhs_general_deep.m', ...
    'fhs_special.m', ...
    'fhs_special_deep.m'};

for ireq = 1:numel(required_files)
    f = fullfile(src_dir,required_files{ireq});
    if ~isfile(f)
        error('Missing required source file: %s',f);
    end
end


fprintf('\n=== 1/5: Dense onset-threshold data ===\n');
THRESHOLD = compute_dense_threshold_data();

fprintf('\n=== 2/5: C_M free-wave data ===\n');
FW_CM = compute_fwcm_data();

fprintf('\n=== 3/5: C_S free-wave data ===\n');
FW_CS = compute_fwcs_data();

fprintf('\n=== 4/5: C_M eigenvector diagnostics ===\n');
EIG_CM = compute_eigdiag_data();

fprintf('\n=== 5/5: C_S eigenvector diagnostics ===\n');
EIG_CS = compute_eigdiag_CS_data();

META = struct();
META.description = 'Unified JFM onset/free-wave/eigenvector data';
META.CM_units = 'Pa s m^{-1}';
META.CS_units = 'Pa s m^{-1}';
META.A_limit = 2.36;
META.Gamma0 = 0.2;
META.beta_base = 1.0;
META.beta_limit = META.A_limit*(1-META.Gamma0)/META.Gamma0;

save(fullfile(data_dir,'JFM_final_data.mat'), ...
    'THRESHOLD','FW_CM','FW_CS','EIG_CM','EIG_CS','META','-v7.3');

fprintf('\nSaved unified data file:\n%s\n',fullfile(data_dir,'JFM_final_data.mat'));

function PACK = compute_dense_threshold_data()
%% Dense onset threshold for final a_c/g figure
N=32; g_dim=9.8066; Omega_dim=2*pi*100;
rho_dim=[998,1.2]; mu1_dim=1e-3; sigma_dim=0.06;
lc=sqrt(sigma_dim/((rho_dim(1)-rho_dim(2))*g_dim)); h_dim=[5*lc,5*lc];
Gamma0_dim=0.2; Ds_dim=2.5e-9;
mu2_list=[1.82e-5,1e-4,1e-3,5e-3,1e-2];

% Much denser than the mechanism sweeps.
nBetaDense=201; nCSDense=201;
beta_base=1.0; A_limit=2.36;
beta_limit=A_limit*(1-Gamma0_dim)/Gamma0_dim; % 9.44
beta_dense=unique(sort([0,logspace(-4,4,nBetaDense-1),beta_base,beta_limit]));
CS_target_dense=[0,logspace(-2,5,nCSDense-1)];

w_dim=pi; cols=3:12; tol_im=1e-6;
m_min_coarse=400; m_max_coarse=1800; m_step_coarse=20;
refine_halfwidth=40; continuation_halfwidth=80; jump_tol_rel=0.30;
continuous_halfwidth=2.0; continuous_tol_k=1e-4;

ncase=numel(mu2_list); CM_cases=cell(ncase,1); CS_cases=cell(ncase,1);
for jc=1:ncase
    mu2=mu2_list(jc); mu=[mu1_dim,mu2]; nu=mu./rho_dim;
    fprintf('\n[DENSE THRESHOLD] eta_2 = %.8e Pa s\n',mu2);

    [kc_CM,ac_CM,winner_CM]=fwcm_compute_faraday_CM_path( ...
        beta_dense,N,Omega_dim,nu,mu,rho_dim,g_dim,sigma_dim,h_dim, ...
        Gamma0_dim,Ds_dim,w_dim,cols,tol_im, ...
        m_min_coarse,m_max_coarse,m_step_coarse, ...
        refine_halfwidth,continuation_halfwidth,jump_tol_rel, ...
        continuous_halfwidth,continuous_tol_k);
    CM_axis=fwcm_compute_CM_reference(beta_dense,kc_CM,sigma_dim,Gamma0_dim,Ds_dim,Omega_dim);
    [~,ibase]=min(abs(beta_dense-beta_base)); [~,ilim]=min(abs(beta_dense-beta_limit));
    Rc=struct(); Rc.mu2=mu2; Rc.beta=beta_dense(:); Rc.kc=kc_CM(:); Rc.ac=ac_CM(:); Rc.winner=winner_CM(:); Rc.CM=CM_axis(:);
    Rc.beta_base=beta_dense(ibase); Rc.beta_limit=beta_dense(ilim); Rc.kc_base=kc_CM(ibase); Rc.kc_limit=kc_CM(ilim); Rc.CM_base=CM_axis(ibase); Rc.CM_limit=CM_axis(ilim);
    CM_cases{jc}=Rc;

    i0=find(beta_dense==0,1); k_clean=kc_CM(i0); muS_dense=CS_target_dense/k_clean^2;
    [kc_CS,ac_CS,winner_CS]=fwcs_compute_faraday_CS_path( ...
        muS_dense,N,Omega_dim,nu,mu,rho_dim,g_dim,sigma_dim,h_dim,0, ...
        w_dim,cols,tol_im,m_min_coarse,m_max_coarse,m_step_coarse, ...
        refine_halfwidth,continuation_halfwidth,jump_tol_rel, ...
        continuous_halfwidth,continuous_tol_k);
    CS_axis=muS_dense(:).*kc_CS(:).^2;
    Rs=struct(); Rs.mu2=mu2; Rs.muS=muS_dense(:); Rs.kc=kc_CS(:); Rs.ac=ac_CS(:); Rs.winner=winner_CS(:); Rs.CS=CS_axis(:);
    CS_cases{jc}=Rs;
end
metadata=struct('rho',rho_dim,'mu1',mu1_dim,'sigma',sigma_dim,'g',g_dim,'Omega',Omega_dim,'Gamma0',Gamma0_dim,'Ds',Ds_dim,'h',h_dim,'mu2_list',mu2_list,'beta_base',beta_base,'beta_limit',beta_limit,'A_limit',A_limit,'CM_units','Pa s m^{-1}','CS_units','Pa s m^{-1}');
PACK=struct('CM_cases',{CM_cases},'CS_cases',{CS_cases},'metadata',metadata);
end


function PACK = compute_fwcm_data()
%% USER SETTINGS
N         = 32;
g_dim     = 9.8066;
Omega_dim = 2*pi*100;

rho_dim   = [998,1.2];
mu1_dim   = 1e-3;
sigma_dim = 0.06;

lc = sqrt(sigma_dim/((rho_dim(1)-rho_dim(2))*g_dim));
h_dim = [5*lc,5*lc];

Gamma0_dim = 0.2;
Ds_dim     = 2.5e-9;

mu2_list = [1.82e-5,1e-4,1e-3,5e-3,1e-2];

% Dense C_M free-wave sweep.
nb = 121;
beta_base  = 1.0;
A_limit    = 2.36;
beta_limit = A_limit*(1-Gamma0_dim)/Gamma0_dim;  % 9.44 for Gamma0=0.2
beta_list = unique(sort([0,logspace(-4,4,nb-1),beta_base,beta_limit]));

% Faraday threshold search
w_dim  = pi;     % with this convention k=m
cols   = 3:12;
tol_im = 1e-6;

m_min_coarse  = 400;
m_max_coarse  = 1800;
m_step_coarse = 20;
refine_halfwidth      = 40;
continuation_halfwidth = 80;
jump_tol_rel = 0.30;  % retained for compatibility; path search below is global

% JFM-quality k refinement: global coarse scan at every beta, followed by
% integer local refinement and continuous sub-grid minimization.
continuous_halfwidth = 2.0;   % 1/m
continuous_tol_k     = 1e-4;  % 1/m

% nonlinear free-wave root
ROOT_TOL     = 1e-9;
ROOT_MAXITER = 80;

RUN_ALL_MU2 = true;
DEBUG_CASE_INDEX = 1;


%% SELECT CASES
if RUN_ALL_MU2
    case_indices = 1:numel(mu2_list);
else
    case_indices = DEBUG_CASE_INDEX;
end

ncase = numel(case_indices);
% MATLAB struct arrays must have identical fields.
% Store each case in a cell array to avoid MATLAB struct-field mismatch issues.
results = cell(ncase,1);

%% MAIN LOOP
for jc = 1:ncase
    ic = case_indices(jc);

    mu2 = mu2_list(ic);
    mu  = [mu1_dim,mu2];
    nu  = mu./rho_dim;

    fprintf('\n============================================================\n');
    fprintf('FREE-WAVE TEST: eta_2 = %.8e Pa s\n',mu2);
    fprintf('============================================================\n');

    % 1) Actual Faraday critical path
    [kc_path,ac_path,winner_path] = fwcm_compute_faraday_CM_path( ...
        beta_list,N,Omega_dim,nu,mu,rho_dim,g_dim,sigma_dim,h_dim, ...
        Gamma0_dim,Ds_dim,w_dim,cols,tol_im, ...
        m_min_coarse,m_max_coarse,m_step_coarse, ...
        refine_halfwidth,continuation_halfwidth,jump_tol_rel, ...
        continuous_halfwidth,continuous_tol_k);

    if ~isfinite(kc_path(1))
        error('Could not determine the clean critical wavenumber.');
    end

    k_clean = kc_path(1);

    fprintf('\nClean critical k_c = %.10g 1/m\n',k_clean);
    fprintf('Clean a_c/g        = %.10g\n',ac_path(1));

    % 2) Fixed-k free mode
    lambda_fixed = complex(nan(size(beta_list)));
    gamma_fixed  = nan(size(beta_list));
    omega_fixed  = nan(size(beta_list));
    detune_fixed = nan(size(beta_list));
    res_fixed    = nan(size(beta_list));
    ok_fixed     = false(size(beta_list));

    lambda_guess = fwcm_initial_free_wave_guess( ...
        k_clean,rho_dim,mu,g_dim,sigma_dim,h_dim);

    for ib = 1:numel(beta_list)
        beta = beta_list(ib);

        [lam,ok,resnorm] = fwcm_solve_free_gravity_capillary_mode( ...
            lambda_guess,k_clean,nu,mu,rho_dim,g_dim,sigma_dim,h_dim, ...
            beta,Gamma0_dim,Ds_dim,ROOT_TOL,ROOT_MAXITER);

        if ~ok
            guess2 = fwcm_initial_free_wave_guess( ...
                k_clean,rho_dim,mu,g_dim,sigma_dim,h_dim);
            [lam2,ok2,res2] = fwcm_solve_free_gravity_capillary_mode( ...
                guess2,k_clean,nu,mu,rho_dim,g_dim,sigma_dim,h_dim, ...
                beta,Gamma0_dim,Ds_dim,ROOT_TOL,ROOT_MAXITER);
            if ok2 || res2 < resnorm
                lam = lam2; ok = ok2; resnorm = res2;
            end
        end

        lambda_fixed(ib) = lam;
        res_fixed(ib) = resnorm;
        ok_fixed(ib) = ok;

        if isfinite(real(lam)) && isfinite(imag(lam))
            gamma_fixed(ib)  = -real(lam);
            omega_fixed(ib)  = abs(imag(lam));
            detune_fixed(ib) = abs(omega_fixed(ib)-Omega_dim/2);
            lambda_guess = lam;
        end

        fprintf(['[fixed k] beta=% .4e | lambda=% .6e %+ .6ei | ' ...
                 'gamma=% .6e | omega0=% .6e | res=% .3e | ok=%d\n'], ...
            beta,real(lam),imag(lam),gamma_fixed(ib), ...
            omega_fixed(ib),resnorm,ok);
    end

    % 3) Free mode along k_c(beta)
    lambda_path = complex(nan(size(beta_list)));
    gamma_path  = nan(size(beta_list));
    omega_path  = nan(size(beta_list));
    detune_path = nan(size(beta_list));
    res_path    = nan(size(beta_list));
    ok_path     = false(size(beta_list));

    lambda_guess = lambda_fixed(1);

    for ib = 1:numel(beta_list)
        beta = beta_list(ib);
        k_now = kc_path(ib);

        if ~isfinite(k_now), continue; end

        if ib>1 && isfinite(kc_path(ib-1))
            if abs(k_now-kc_path(ib-1))/kc_path(ib-1) > 0.15
                lambda_guess = fwcm_initial_free_wave_guess( ...
                    k_now,rho_dim,mu,g_dim,sigma_dim,h_dim);
            end
        end

        [lam,ok,resnorm] = fwcm_solve_free_gravity_capillary_mode( ...
            lambda_guess,k_now,nu,mu,rho_dim,g_dim,sigma_dim,h_dim, ...
            beta,Gamma0_dim,Ds_dim,ROOT_TOL,ROOT_MAXITER);

        if ~ok
            guess2 = fwcm_initial_free_wave_guess( ...
                k_now,rho_dim,mu,g_dim,sigma_dim,h_dim);
            [lam2,ok2,res2] = fwcm_solve_free_gravity_capillary_mode( ...
                guess2,k_now,nu,mu,rho_dim,g_dim,sigma_dim,h_dim, ...
                beta,Gamma0_dim,Ds_dim,ROOT_TOL,ROOT_MAXITER);
            if ok2 || res2 < resnorm
                lam = lam2; ok = ok2; resnorm = res2;
            end
        end

        lambda_path(ib) = lam;
        res_path(ib) = resnorm;
        ok_path(ib) = ok;

        if isfinite(real(lam)) && isfinite(imag(lam))
            gamma_path(ib)  = -real(lam);
            omega_path(ib)  = abs(imag(lam));
            detune_path(ib) = abs(omega_path(ib)-Omega_dim/2);
            lambda_guess = lam;
        end

        fprintf(['[critical k] beta=% .4e | k=% .6g | ' ...
                 'lambda=% .6e %+ .6ei | gamma=% .6e | ' ...
                 'omega0=% .6e | res=% .3e | ok=%d\n'], ...
            beta,k_now,real(lam),imag(lam),gamma_path(ib), ...
            omega_path(ib),resnorm,ok);
    end

    % 4) map beta -> Faraday-reference C_M
    CM_fixed = fwcm_compute_CM_reference( ...
        beta_list,k_clean*ones(size(beta_list)), ...
        sigma_dim,Gamma0_dim,Ds_dim,Omega_dim);

    CM_path = fwcm_compute_CM_reference( ...
        beta_list,kc_path,sigma_dim,Gamma0_dim,Ds_dim,Omega_dim);

    R = struct();
    R.mu2 = mu2;
    R.mu  = mu;
    R.nu  = nu;
    R.beta = beta_list(:);
    R.k_clean = k_clean;
    R.kc_path = kc_path(:);
    R.ac_path = ac_path(:);
    R.winner_path = winner_path(:);
    R.CM_fixed = CM_fixed(:);
    R.CM_path  = CM_path(:);
    R.lambda_fixed = lambda_fixed(:);
    R.gamma_fixed  = gamma_fixed(:);
    R.omega_fixed  = omega_fixed(:);
    R.detune_fixed = detune_fixed(:);
    R.res_fixed    = res_fixed(:);
    R.ok_fixed     = ok_fixed(:);
    R.lambda_path = lambda_path(:);
    R.gamma_path  = gamma_path(:);
    R.omega_path  = omega_path(:);
    R.detune_path = detune_path(:);
    R.res_path    = res_path(:);
    R.ok_path     = ok_path(:);

    [~,ibase] = min(abs(R.beta-1.0));
    beta_limit = 2.36*(1-Gamma0_dim)/Gamma0_dim;
    [~,ilim] = min(abs(R.beta-beta_limit));
    R.beta_base = R.beta(ibase);
    R.beta_limit = R.beta(ilim);
    R.kc_base = R.kc_path(ibase);
    R.kc_limit = R.kc_path(ilim);
    R.CM_base_fixed = R.CM_fixed(ibase);
    R.CM_limit_fixed = R.CM_fixed(ilim);
    R.CM_base_path = R.CM_path(ibase);
    R.CM_limit_path = R.CM_path(ilim);

    results{jc} = R;
end
metadata = struct();
metadata.rho = rho_dim;
metadata.mu1 = mu1_dim;
metadata.sigma = sigma_dim;
metadata.g = g_dim;
metadata.Omega = Omega_dim;
metadata.Gamma0 = Gamma0_dim;
metadata.Ds = Ds_dim;
metadata.h = h_dim;
metadata.beta_list = beta_list;
metadata.definition_CM = ...
    'Faraday-reference CM evaluated with i*Omega/2 in denominator';

PACK = struct('results',{results},'metadata',metadata);
end

function PACK = compute_fwcs_data()
%% parameters
N         = 32;
g_dim     = 9.8066;
Omega_dim = 2*pi*100;

rho_dim   = [998,1.2];
mu1_dim   = 1e-3;
sigma_dim = 0.06;


lc    = sqrt(sigma_dim/((rho_dim(1)-rho_dim(2))*g_dim));
h_dim = [5*lc,5*lc];

Gamma0_dim = 0.2; 
Ds_dim     = 2.5e-9; 
beta_dim   = 0;      % C_M = 0 exactly

mu2_list = [1.82e-5,1e-4,1e-3,5e-3,1e-2];

% rho_dim   = [519.933,1];
% mu1_dim   = 3.908e-5;
% sigma_dim = 2.181e-6;
% 
% lc    = sqrt(sigma_dim/((rho_dim(1)-rho_dim(2))*g_dim));
% h_dim = [5*lc,5*lc];
% 
% Gamma0_dim = 0.2;
% Ds_dim     = 2.5e-9;
% beta_dim   = 0;
% 
% mu2_list = [2e-5,3e-5,3.908e-5,5e-5,6e-5];



% Desired fixed-k C_S range.  muS_tot will be generated from this
% after k_c,clean is known for each viscosity case.
nCS = 121;
CS_target = [0,logspace(-2,5,nCS-1)];

% threshold search
w_dim = pi;
cols   = 3:12;
tol_im = 1e-6;

% m_min_coarse  = 6000;
% m_max_coarse  = 20000;
% m_step_coarse = 50;

m_min_coarse  = 400;
m_max_coarse  = 1800;
m_step_coarse = 20;

refine_halfwidth       = 40;
continuation_halfwidth = 80;
jump_tol_rel           = 0.30;  % retained for compatibility; path search below is global

% JFM-quality k refinement: global coarse scan at every parameter value,
% integer local refinement, then continuous sub-grid minimization.
continuous_halfwidth    = 2.0;   % 1/m, bracket around each integer branch minimum
continuous_tol_k        = 1e-4;  % 1/m, fminbnd tolerance

% nonlinear root
ROOT_TOL = 1e-9;
ROOT_MAXITER = 80;

RUN_ALL_MU2 = true;
DEBUG_CASE_INDEX = 1;


if RUN_ALL_MU2
    case_indices = 1:numel(mu2_list);
else
    case_indices = DEBUG_CASE_INDEX;
end

ncase = numel(case_indices);
results = cell(ncase,1);

%% ============================================================
% MAIN LOOP
%% ============================================================
for jc = 1:ncase

    ic = case_indices(jc);
    mu2 = mu2_list(ic);
    mu = [mu1_dim,mu2];
    nu = mu./rho_dim;

    fprintf('\n============================================================\n');
    fprintf('FREE-WAVE C_S TEST: eta_2 = %.8e Pa s\n',mu2);
    fprintf('C_M = 0 (beta = 0)\n');
    fprintf('============================================================\n');

    % --------------------------------------------------------
    % 1) clean critical point
    % --------------------------------------------------------
    R0 = fwcs_search_k_list_CS( ...
        m_min_coarse:m_step_coarse:m_max_coarse, ...
        w_dim,cols,tol_im,N,Omega_dim,nu,mu,rho_dim,g_dim, ...
        sigma_dim,h_dim,beta_dim,0);

    if ~isfinite(R0.kc)
        error('Could not determine clean critical wavenumber.');
    end

    m0 = max(1,round(R0.kc*w_dim/pi));
    R1 = fwcs_search_k_list_CS( ...
        max(1,m0-refine_halfwidth):(m0+refine_halfwidth), ...
        w_dim,cols,tol_im,N,Omega_dim,nu,mu,rho_dim,g_dim, ...
        sigma_dim,h_dim,beta_dim,0);

    if isfinite(R1.kc), R0 = R1; end

    % Continuous sub-grid refinement of the clean critical point.
    R0 = fwcs_refine_k_continuous_CS(R0,continuous_halfwidth,continuous_tol_k, ...
        N,Omega_dim,nu,mu,rho_dim,g_dim,sigma_dim,h_dim,beta_dim,0,cols,tol_im);

    k_clean = R0.kc;
    ac_clean = R0.ac;

    fprintf('clean kc = %.8g 1/m, ac/g = %.8g\n',k_clean,ac_clean);

    % physical total surface viscosity [N s / m]
    muS_list = CS_target/k_clean^2;

    % --------------------------------------------------------
    % 2) actual forced Faraday critical path
    % --------------------------------------------------------
    [kc_path,ac_path,winner_path] = fwcs_compute_faraday_CS_path( ...
        muS_list,N,Omega_dim,nu,mu,rho_dim,g_dim,sigma_dim,h_dim, ...
        beta_dim,w_dim,cols,tol_im,m_min_coarse,m_max_coarse, ...
        m_step_coarse,refine_halfwidth,continuation_halfwidth,jump_tol_rel, ...
        continuous_halfwidth,continuous_tol_k);

    % --------------------------------------------------------
    % 3) free wave at fixed k_clean
    % --------------------------------------------------------
    lambda_fixed = complex(nan(size(muS_list)));
    gamma_fixed  = nan(size(muS_list));
    omega_fixed  = nan(size(muS_list));
    detune_fixed = nan(size(muS_list));
    res_fixed    = nan(size(muS_list));
    ok_fixed     = false(size(muS_list));

    lambda_guess = fwcs_initial_free_wave_guess(k_clean,rho_dim,mu,g_dim,sigma_dim,h_dim);

    for is = 1:numel(muS_list)
        muS_tot = muS_list(is);

        [lam,ok,resnorm] = fwcs_solve_free_mode_CS( ...
            lambda_guess,k_clean,nu,mu,rho_dim,g_dim,sigma_dim,h_dim, ...
            muS_tot,ROOT_TOL,ROOT_MAXITER);

        if ~ok
            guess2 = fwcs_initial_free_wave_guess(k_clean,rho_dim,mu,g_dim,sigma_dim,h_dim);
            [lam2,ok2,res2] = fwcs_solve_free_mode_CS( ...
                guess2,k_clean,nu,mu,rho_dim,g_dim,sigma_dim,h_dim, ...
                muS_tot,ROOT_TOL,ROOT_MAXITER);
            if ok2 || res2 < resnorm
                lam=lam2; ok=ok2; resnorm=res2;
            end
        end

        lambda_fixed(is)=lam;
        res_fixed(is)=resnorm;
        ok_fixed(is)=ok;

        if isfinite(real(lam)) && isfinite(imag(lam))
            gamma_fixed(is)  = -real(lam);
            omega_fixed(is)  = abs(imag(lam));
            detune_fixed(is) = abs(omega_fixed(is)-Omega_dim/2);
            lambda_guess = lam;
        end

        fprintf(['[fixed k] CS=% .4e | lambda=% .6e %+ .6ei | ' ...
                 'gamma=% .6e | omega0=% .6e | res=% .3e | ok=%d\n'], ...
            muS_tot*k_clean^2,real(lam),imag(lam), ...
            gamma_fixed(is),omega_fixed(is),resnorm,ok);
    end

    % --------------------------------------------------------
    % 4) free wave along k = kc(CS)
    % --------------------------------------------------------
    lambda_path = complex(nan(size(muS_list)));
    gamma_path  = nan(size(muS_list));
    omega_path  = nan(size(muS_list));
    detune_path = nan(size(muS_list));
    res_path    = nan(size(muS_list));
    ok_path     = false(size(muS_list));

    lambda_guess = lambda_fixed(1);

    for is = 1:numel(muS_list)
        k_now = kc_path(is);
        muS_tot = muS_list(is);
        if ~isfinite(k_now), continue; end

        if is>1 && isfinite(kc_path(is-1)) && ...
                abs(k_now-kc_path(is-1))/kc_path(is-1) > 0.15
            lambda_guess = fwcs_initial_free_wave_guess(k_now,rho_dim,mu,g_dim,sigma_dim,h_dim);
        end

        [lam,ok,resnorm] = fwcs_solve_free_mode_CS( ...
            lambda_guess,k_now,nu,mu,rho_dim,g_dim,sigma_dim,h_dim, ...
            muS_tot,ROOT_TOL,ROOT_MAXITER);

        if ~ok
            guess2 = fwcs_initial_free_wave_guess(k_now,rho_dim,mu,g_dim,sigma_dim,h_dim);
            [lam2,ok2,res2] = fwcs_solve_free_mode_CS( ...
                guess2,k_now,nu,mu,rho_dim,g_dim,sigma_dim,h_dim, ...
                muS_tot,ROOT_TOL,ROOT_MAXITER);
            if ok2 || res2 < resnorm
                lam=lam2; ok=ok2; resnorm=res2;
            end
        end

        lambda_path(is)=lam;
        res_path(is)=resnorm;
        ok_path(is)=ok;

        if isfinite(real(lam)) && isfinite(imag(lam))
            gamma_path(is)  = -real(lam);
            omega_path(is)  = abs(imag(lam));
            detune_path(is) = abs(omega_path(is)-Omega_dim/2);
            lambda_guess = lam;
        end

        fprintf(['[critical k] CS=% .4e | k=% .6g | ' ...
                 'lambda=% .6e %+ .6ei | gamma=% .6e | omega0=% .6e | ' ...
                 'res=% .3e | ok=%d\n'], ...
            muS_tot*k_now^2,k_now,real(lam),imag(lam), ...
            gamma_path(is),omega_path(is),resnorm,ok);
    end

    CS_fixed = muS_list*k_clean^2;
    CS_path  = muS_list.*kc_path.^2;

    R = struct();
    R.mu2 = mu2;
    R.mu = mu;
    R.nu = nu;
    R.muS_list = muS_list(:);
    R.k_clean = k_clean;
    R.ac_clean = ac_clean;
    R.CS_fixed = CS_fixed(:);
    R.CS_path = CS_path(:);
    R.kc_path = kc_path(:);
    R.ac_path = ac_path(:);
    R.winner_path = winner_path(:);
    R.lambda_fixed = lambda_fixed(:);
    R.gamma_fixed = gamma_fixed(:);
    R.omega_fixed = omega_fixed(:);
    R.detune_fixed = detune_fixed(:);
    R.res_fixed = res_fixed(:);
    R.ok_fixed = ok_fixed(:);
    R.lambda_path = lambda_path(:);
    R.gamma_path = gamma_path(:);
    R.omega_path = omega_path(:);
    R.detune_path = detune_path(:);
    R.res_path = res_path(:);
    R.ok_path = ok_path(:);
    results{jc}=R;
end

metadata = struct('rho',rho_dim,'mu1',mu1_dim,'sigma',sigma_dim,'g',g_dim, ...
    'Omega',Omega_dim,'h',h_dim,'beta',0,'CM',0, ...
    'definition_CS','C_S(k)=(kappa_0^S+eta_0^S)k^2');
PACK = struct('results',{results},'metadata',metadata);

end

function PACK = compute_eigdiag_data()
%% ============================================================
% USER SETTINGS
%% ============================================================

N         = 32;
g_dim     = 9.8066;
omega_dim = 2*pi*100;

rho_dim   = [998,1.2];
mu1_dim   = 1e-3;
gamma_dim = 0.06;

lc    = sqrt(gamma_dim/(rho_dim(1)-rho_dim(2))/g_dim);
h_dim = [5*lc,5*lc];

Gamma0_dim   = 0.2;
Ds_dim       = 2.5e-9;
sigma_values = 0;

mu2_list = [1.82e-5,1e-4,1e-3,5e-3,1e-2];

nb = 121;
beta_base  = 1.0;
A_limit    = 2.36;
beta_limit = A_limit*(1-Gamma0_dim)/Gamma0_dim;  % 9.44 for Gamma0=0.2
beta_list = unique(sort([logspace(log10(1e-4),log10(1e4),nb), ...
    beta_base,beta_limit]));

w_dim  = pi;
cols   = 3:12;
tol_im = 1e-6;

USE_REFINE       = true;
USE_CONTINUATION = true;
ALLOW_FALLBACK   = true;

m_max_coarse     = 1800;
refine_halfwidth = 120;
jump_tol_rel     = 0.35;
edge_tol_m       = 5;

verbose_print = true;
%% ============================================================
% RUN C_M SWEEPS
%% ============================================================

ncase = numel(mu2_list);

cases = repmat(struct( ...
    'mu2',[],'mu_dim',[],'nu_dim',[],'run_CM',[]),ncase,1);

for ic = 1:ncase

    mu2_dim = mu2_list(ic);
    mu_dim  = [mu1_dim,mu2_dim];
    nu_dim  = mu_dim./rho_dim;

    fprintf('\n============================================================\n');
    fprintf('C_M CASE %d / %d : eta_2 = %.12g Pa s\n',ic,ncase,mu2_dim);
    fprintf('============================================================\n');

    run_CM = eig_run_one_CM_sweep( ...
        beta_list, ...
        N,g_dim,omega_dim,rho_dim,mu_dim,nu_dim,gamma_dim,h_dim, ...
        Gamma0_dim,Ds_dim,sigma_values, ...
        w_dim,cols,tol_im, ...
        USE_REFINE,USE_CONTINUATION,ALLOW_FALLBACK, ...
        m_max_coarse,refine_halfwidth,jump_tol_rel,edge_tol_m, ...
        verbose_print);

    run_CM.CM_actual = eig_compute_CM_axis( ...
        run_CM.beta_list,run_CM.kc_list, ...
        gamma_dim,Gamma0_dim,Ds_dim,omega_dim);

    [~,ibase] = min(abs(run_CM.beta_list-1.0));
    beta_limit = 2.36*(1-Gamma0_dim)/Gamma0_dim;
    [~,ilim] = min(abs(run_CM.beta_list-beta_limit));
    run_CM.beta_base = run_CM.beta_list(ibase);
    run_CM.beta_limit = run_CM.beta_list(ilim);
    run_CM.kc_base = run_CM.kc_list(ibase);
    run_CM.kc_limit = run_CM.kc_list(ilim);
    run_CM.CM_base = run_CM.CM_actual(ibase);
    run_CM.CM_limit = run_CM.CM_actual(ilim);

    cases(ic).mu2    = mu2_dim;
    cases(ic).mu_dim = mu_dim;
    cases(ic).nu_dim = nu_dim;
    cases(ic).run_CM = run_CM;
end

metadata = struct();
metadata.rho = rho_dim;
metadata.mu1 = mu1_dim;
metadata.sigma = gamma_dim;
metadata.g = g_dim;
metadata.Omega = omega_dim;
metadata.Gamma0 = Gamma0_dim;
metadata.Ds = Ds_dim;
metadata.h = h_dim;
metadata.normalization = 'Each critical eigenvector is rescaled so that ||zeta||_RMS = 1';
results = struct();
results.mu2_list = mu2_list;
results.cases = cases;
PACK = struct('results',results,'metadata',metadata);
end

function PACK = compute_eigdiag_CS_data()
%% ============================================================
% C_S EIGENVECTOR DIAGNOSTICS
% Computes critical-mode u_s and actual surface-viscous traction tau_S
% along the Faraday critical branch for a pure surface-viscosity sweep.
%% ============================================================

N         = 32;
g_dim     = 9.8066;
omega_dim = 2*pi*100;

rho_dim   = [998,1.2];
mu1_dim   = 1e-3;
gamma_dim = 0.06;

lc    = sqrt(gamma_dim/(rho_dim(1)-rho_dim(2))/g_dim);
h_dim = [5*lc,5*lc];

beta_dim   = 0;
Gamma0_dim = 0.2;
Ds_dim     = 2.5e-9;

mu2_list = [1.82e-5,1e-4,1e-3,5e-3,1e-2];

nCS = 121;
CS_target = logspace(-2,5,nCS);

w_dim  = pi;
cols   = 3:12;
tol_im = 1e-6;

USE_REFINE       = true;
USE_CONTINUATION = true;
ALLOW_FALLBACK   = true;

m_max_coarse     = 1800;
refine_halfwidth = 120;
jump_tol_rel     = 0.35;
edge_tol_m       = 5;

verbose_print = true;

ncase = numel(mu2_list);
cases = repmat(struct('mu2',[],'mu_dim',[],'nu_dim',[],'run_CS',[]),ncase,1);

for ic = 1:ncase
    mu2_dim = mu2_list(ic);
    mu_dim  = [mu1_dim,mu2_dim];
    nu_dim  = mu_dim./rho_dim;

    fprintf('\n============================================================\n');
    fprintf('C_S CASE %d / %d : eta_2 = %.12g Pa s\n',ic,ncase,mu2_dim);
    fprintf('============================================================\n');

    R0 = eig_find_critical_point( ...
        1,nan, ...
        N,g_dim,omega_dim,rho_dim,mu_dim,nu_dim,gamma_dim,h_dim, ...
        beta_dim,Gamma0_dim,Ds_dim,0,0, ...
        w_dim,cols,tol_im, ...
        USE_REFINE,false,ALLOW_FALLBACK, ...
        m_max_coarse,refine_halfwidth,jump_tol_rel,edge_tol_m);

    if ~isfinite(R0.kc) || ~isfinite(R0.ac)
        error('Could not determine clean critical point for eta_2 = %.6e.',mu2_dim);
    end

    k_clean = R0.kc;
    muS_list = CS_target/k_clean^2;

    run_CS = eig_run_one_CS_sweep( ...
        muS_list, ...
        N,g_dim,omega_dim,rho_dim,mu_dim,nu_dim,gamma_dim,h_dim, ...
        beta_dim,Gamma0_dim,Ds_dim, ...
        w_dim,cols,tol_im, ...
        USE_REFINE,USE_CONTINUATION,ALLOW_FALLBACK, ...
        m_max_coarse,refine_halfwidth,jump_tol_rel,edge_tol_m, ...
        verbose_print);

    run_CS.CS_actual = run_CS.muS_list(:).*run_CS.kc_list(:).^2;
    run_CS.CS_reference = run_CS.muS_list(:)*k_clean^2;
    run_CS.k_clean = k_clean;
    run_CS.ac_clean = R0.ac;

    cases(ic).mu2     = mu2_dim;
    cases(ic).mu_dim  = mu_dim;
    cases(ic).nu_dim  = nu_dim;
    cases(ic).run_CS  = run_CS;
end

metadata = struct();
metadata.rho = rho_dim;
metadata.mu1 = mu1_dim;
metadata.sigma = gamma_dim;
metadata.g = g_dim;
metadata.Omega = omega_dim;
metadata.h = h_dim;
metadata.beta = beta_dim;
metadata.Gamma0 = Gamma0_dim;
metadata.Ds = Ds_dim;
metadata.normalization = 'Each critical eigenvector is rescaled so that ||zeta||_RMS = 1';
metadata.definition_CS = 'C_S(k_c)=(kappa_0^S+eta_0^S)k_c^2';
metadata.definition_tauS = 'tau_S,n = -C_S u_s,n';

results = struct();
results.mu2_list = mu2_list;
results.cases = cases;

PACK = struct('results',results,'metadata',metadata);
end

% ======================== C_M FREE-WAVE HELPERS ========================
function [kc,ac,winner] = fwcm_compute_faraday_CM_path( ...
    beta_list,N,Omega,nu,mu,rho,g,sigma0,h,Gamma0,Ds, ...
    w_dim,cols,tol_im,m_min_coarse,m_max_coarse,m_step_coarse, ...
    refine_halfwidth,continuation_halfwidth,jump_tol_rel, ...
    continuous_halfwidth,continuous_tol_k)
% JFM-quality critical-wavenumber search.
% Every beta is searched independently to avoid continuation-window locking.
% Harmonic and subharmonic minima are refined separately in continuous k.
%
% continuation_halfwidth and jump_tol_rel remain in the signature only for
% backward compatibility.

nb = numel(beta_list);
kc = nan(size(beta_list)); ac = nan(size(beta_list)); winner = nan(size(beta_list));

for ib = 1:nb
    beta = beta_list(ib); surface = 0;
    fprintf('\nFaraday threshold beta = %.5e\n',beta);

    % ---- Stage 1: GLOBAL coarse search at every beta ----
    m_coarse = m_min_coarse:m_step_coarse:m_max_coarse;
    R = fwcm_search_k_list(m_coarse,w_dim,cols,tol_im,N,Omega,nu,mu, ...
        rho,g,sigma0,h,beta,Gamma0,Ds,surface);

    % ---- Stage 2: integer-k local refinement ----
    if isfinite(R.kc)
        m0=max(1,round(R.kc*w_dim/pi));
        m_ref=max(1,m0-refine_halfwidth):(m0+refine_halfwidth);
        R2=fwcm_search_k_list(m_ref,w_dim,cols,tol_im,N,Omega,nu,mu, ...
            rho,g,sigma0,h,beta,Gamma0,Ds,surface);
        if isfinite(R2.kc), R=R2; end
    end

    % ---- Stage 3: continuous sub-grid refinement ----
    R=fwcm_refine_k_continuous_CM(R,continuous_halfwidth,continuous_tol_k, ...
        N,Omega,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds,surface,cols,tol_im);

    kc(ib)=R.kc; ac(ib)=R.ac; winner(ib)=R.winner;
    fprintf('  kc = %.8f, ac/g = %.10g, winner = %g\n',kc(ib),ac(ib),winner(ib));
end
end

function R = fwcm_search_k_list(m_list,w_dim,cols,tol_im,N,Omega,nu,mu, ...
    rho,g,sigma0,h,beta,Gamma0,Ds,surface)

R = struct('kc',nan,'ac',nan,'winner',nan, ...
           'kh',nan,'ah',inf,'ks',nan,'as',inf);
best_h = inf; best_s = inf; kh = nan; ks = nan;

for ii = 1:numel(m_list)
    k = m_list(ii)*pi/w_dim;

    out_h = compute_eigenvalues(0,0,true,k,N,Omega,nu,mu,rho,g, ...
        sigma0,h,beta,Gamma0,Ds,surface);
    out_s = compute_eigenvalues(0,Omega/2,false,k,N,Omega,nu,mu,rho,g, ...
        sigma0,h,beta,Gamma0,Ds,surface);

    vh = fwcm_valid_positive_thresholds(out_h,cols,tol_im);
    vs = fwcm_valid_positive_thresholds(out_s,cols,tol_im);

    if ~isempty(vh)
        ah=min(vh); if ah<best_h, best_h=ah; kh=k; end
    end
    if ~isempty(vs)
        as=min(vs); if as<best_s, best_s=as; ks=k; end
    end
end

R.kh=kh; R.ah=best_h; R.ks=ks; R.as=best_s;
if isfinite(best_h) && isfinite(best_s)
    if best_h<=best_s, R.kc=kh; R.ac=best_h; R.winner=1;
    else, R.kc=ks; R.ac=best_s; R.winner=2; end
elseif isfinite(best_h)
    R.kc=kh; R.ac=best_h; R.winner=1;
elseif isfinite(best_s)
    R.kc=ks; R.ac=best_s; R.winner=2;
end
end

function R = fwcm_refine_k_continuous_CM(R,halfwidth,tol_k,N,Omega,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds,surface,cols,tol_im)
% Continuous refinement removes the integer-k staircase in k_c(beta).
opts=optimset('Display','off','TolX',tol_k,'MaxIter',100,'MaxFunEvals',250);
kh=R.kh; ah=R.ah; ks=R.ks; as=R.as;

if isfinite(kh)
    a=max(1e-12,kh-halfwidth); b=kh+halfwidth;
    obj=@(k) fwcm_threshold_at_k_CM(k,1,cols,tol_im,N,Omega,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds,surface);
    [ktry,atry]=fminbnd(obj,a,b,opts);
    if isfinite(atry) && atry<1e90, kh=ktry; ah=atry; end
end
if isfinite(ks)
    a=max(1e-12,ks-halfwidth); b=ks+halfwidth;
    obj=@(k) fwcm_threshold_at_k_CM(k,2,cols,tol_im,N,Omega,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds,surface);
    [ktry,atry]=fminbnd(obj,a,b,opts);
    if isfinite(atry) && atry<1e90, ks=ktry; as=atry; end
end

R.kh=kh; R.ah=ah; R.ks=ks; R.as=as;
if isfinite(ah) && isfinite(as)
    if ah<=as, R.kc=kh; R.ac=ah; R.winner=1; else, R.kc=ks; R.ac=as; R.winner=2; end
elseif isfinite(ah), R.kc=kh; R.ac=ah; R.winner=1;
elseif isfinite(as), R.kc=ks; R.ac=as; R.winner=2;
else, R.kc=nan; R.ac=nan; R.winner=nan;
end
end

function acrit = fwcm_threshold_at_k_CM(k,branch,cols,tol_im,N,Omega,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds,surface)
if ~(isfinite(k) && k>0), acrit=1e100; return; end
if branch==1
    out=compute_eigenvalues(0,0,true,k,N,Omega,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds,surface);
else
    out=compute_eigenvalues(0,Omega/2,false,k,N,Omega,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds,surface);
end
v=fwcm_valid_positive_thresholds(out,cols,tol_im);
if isempty(v), acrit=1e100; else, acrit=min(v); end
end

function vals = fwcm_valid_positive_thresholds(out,cols,tol_im)
if isempty(out), vals=[]; return; end
cols_use = cols(cols<=numel(out));
if isempty(cols_use), vals=[]; return; end
v = out(cols_use);
good = isfinite(v) & real(v)>0 & abs(imag(v))<tol_im;
vals = real(v(good));
end

function [lambda,ok,resnorm] = fwcm_solve_free_gravity_capillary_mode( ...
    lambda0,k,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds,tol,maxiter)

x = [real(lambda0);abs(imag(lambda0))];
if x(2)<=0 || ~isfinite(x(2)), x(2)=1; end
ok = false;

for iter = 1:maxiter
    F = fwcm_free_wave_residual_real(x,k,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds);
    fn = norm(F,2);
    if fn<tol, ok=true; break; end

    J = zeros(2,2);
    for j=1:2
        hj = 1e-6*max(1,abs(x(j)));
        xp=x; xm=x;
        xp(j)=xp(j)+hj; xm(j)=xm(j)-hj;
        Fp = fwcm_free_wave_residual_real(xp,k,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds);
        Fm = fwcm_free_wave_residual_real(xm,k,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds);
        J(:,j)=(Fp-Fm)/(2*hj);
    end

    if ~all(isfinite(J(:))) || rcond(J)<1e-13, break; end
    dx = -J\F;

    accepted=false; step=1;
    for ls=1:12
        xt=x+step*dx;
        if xt(2)<=0, xt(2)=abs(xt(2)); end
        Ft = fwcm_free_wave_residual_real(xt,k,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds);
        if all(isfinite(Ft)) && norm(Ft,2)<fn
            x=xt; accepted=true; break;
        end
        step=step/2;
    end
    if ~accepted, break; end
end

F = fwcm_free_wave_residual_real(x,k,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds);
if norm(F,2)>=tol
    obj = @(xx) sum(fwcm_free_wave_residual_real([xx(1);abs(xx(2))], ...
        k,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds).^2);
    opts = optimset('Display','off','TolX',1e-10,'TolFun',1e-14, ...
        'MaxIter',1000,'MaxFunEvals',3000);
    xtry = fminsearch(obj,x,opts);
    xtry(2)=abs(xtry(2));
    Ftry = fwcm_free_wave_residual_real(xtry,k,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds);
    if norm(Ftry,2)<norm(F,2), x=xtry; F=Ftry; end
end

lambda = x(1)+1i*x(2);
resnorm = norm(F,2);
ok = resnorm < max(tol,1e-7);
end

function F = fwcm_free_wave_residual_real(x,k,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds)
lambda = x(1)+1i*x(2);
R = fwcm_free_wave_residual_complex(lambda,k,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds);
rho_jump = rho(2)-rho(1);
Krest = sigma0*k^4 - rho_jump*g*k^2;
Kscale = max(abs(Krest),1);
F = [real(R);imag(R)]/Kscale;
if ~all(isfinite(F)), F=[1e20;1e20]; end
end

function R = fwcm_free_wave_residual_complex(lambda,k,nu,mu,rho,g,sigma0,h,beta,Gamma0,Ds)
D = fwcm_interface_derivatives_per_unit_zeta_lambda(lambda,k,nu,mu,sigma0,h,beta,Gamma0,Ds);

Ntotal = ...
    ((rho(2)*lambda+3*k^2*mu(2))*D.dz2 - mu(2)*D.dzzz2) ...
  - ((rho(1)*lambda+3*k^2*mu(1))*D.dz1 - mu(1)*D.dzzz1);

rho_jump = rho(2)-rho(1);
Krest = sigma0*k^4 - rho_jump*g*k^2;
R = Ntotal - Krest;
end

function D = fwcm_interface_derivatives_per_unit_zeta_lambda( ...
    lambda,k,nu,mu,sigma0,h,beta,Gamma0,Ds)

if beta~=0 && Gamma0~=0
    Mn_M = beta*Gamma0*sigma0*k^2 / ...
        ((lambda+Ds*k^2)*(1-Gamma0));
else
    Mn_M = 0;
end
Mn = Mn_M;

q1 = sqrt(k^2+lambda/nu(1));
q2 = sqrt(k^2+lambda/nu(2));
if real(q1)<0, q1=-q1; end
if real(q2)<0, q2=-q2; end

KH_DEEP=8.0;
USE_DEEP=min(abs(k*h(1)),abs(k*h(2)))>=KH_DEEP;

if abs(lambda)>1e-14
    if ~USE_DEEP
        A=fhs_general(k,mu,h,q1,q2,Mn);
        b=zeros(8,1); b(1)=lambda;
        c=fwcm_scaled_solve(A,b);
        dz1=[k,-k,q1,-q1]*c(1:4);
        dz2=[k,-k,q2,-q2]*c(5:8);
        dzzz1=[k^3,-k^3,q1^3,-q1^3]*c(1:4);
        dzzz2=[k^3,-k^3,q2^3,-q2^3]*c(5:8);
    else
        A=fhs_general_deep(k,mu,q1,q2,Mn);
        b=zeros(4,1); b(1)=lambda;
        c=fwcm_scaled_solve(A,b);
        dz1=[k,q1]*c(1:2);
        dz2=[-k,-q2]*c(3:4);
        dzzz1=[k^3,q1^3]*c(1:2);
        dzzz2=[-k^3,-q2^3]*c(3:4);
    end
else
    if ~USE_DEEP
        A=fhs_special(k,mu,h,Mn);
        b=zeros(8,1); b(1)=lambda;
        c=fwcm_scaled_solve(A,b);
        dz1=[k,-k,1,1]*c(1:4);
        dz2=[k,-k,1,1]*c(5:8);
        dzzz1=[k^3,-k^3,3*k^2,3*k^2]*c(1:4);
        dzzz2=[k^3,-k^3,3*k^2,3*k^2]*c(5:8);
    else
        A=fhs_special_deep(k,mu,Mn);
        b=zeros(4,1); b(1)=lambda;
        c=fwcm_scaled_solve(A,b);
        dz1=[k,1]*c(1:2);
        dz2=[-k,1]*c(3:4);
        dzzz1=[k^3,3*k^2]*c(1:2);
        dzzz2=[-k^3,3*k^2]*c(3:4);
    end
end

D=struct('dz1',dz1,'dz2',dz2,'dzzz1',dzzz1,'dzzz2',dzzz2, ...
    'q1',q1,'q2',q2,'Mn_M',Mn_M,'USE_DEEP',USE_DEEP);
end

function x = fwcm_scaled_solve(A,b)
if ~all(isfinite(A(:))) || ~all(isfinite(b(:)))
    error('FreeWave:NonFiniteHydroMatrix','Hydrodynamic matrix contains NaN/Inf.');
end
rn=sqrt(sum(abs(A),2)); cn=sqrt(sum(abs(A'),2));
rn(rn==0)=1; cn(cn==0)=1;
D1=diag(1./rn); D2=diag(1./cn);
xs=(D1*A*D2)\(D1*b);
x=D2*xs;
end

function lambda0 = fwcm_initial_free_wave_guess(k,rho,mu,g,sigma0,h)
den = rho(1)*coth(max(k*h(1),1e-12)) + ...
      rho(2)*coth(max(k*h(2),1e-12));
omega2 = ((rho(1)-rho(2))*g*k + sigma0*k^3)/den;
omega0 = sqrt(max(real(omega2),1e-16));
nu_eff = (mu(1)+mu(2))/(rho(1)+rho(2));
gamma_guess = 2*nu_eff*k^2;
lambda0 = -gamma_guess + 1i*omega0;
end

function CM = fwcm_compute_CM_reference(beta,k,sigma0,Gamma0,Ds,Omega)
beta=beta(:); k=k(:);
if isscalar(k), k=repmat(k,size(beta)); end
CM=nan(size(beta));
for i=1:numel(beta)
    if ~isfinite(k(i)) || k(i)<=0, continue; end
    Mref = sigma0*beta(i)*Gamma0*k(i)^2 / ...
        ((1-Gamma0)*(1i*Omega/2+Ds*k(i)^2));
    CM(i)=abs(Mref);
end
end


% ======================== C_S FREE-WAVE HELPERS ========================
function [kc,ac,winner] = fwcs_compute_faraday_CS_path(muS_list,N,Omega,nu,mu,rho,g,sigma0,h,beta,w_dim,cols,tol_im,m_min_coarse,m_max_coarse,m_step_coarse,refine_halfwidth,continuation_halfwidth,jump_tol_rel,continuous_halfwidth,continuous_tol_k)
% JFM-quality critical-wavenumber search.
% Each parameter value is searched independently so the result cannot be
% trapped by the previous continuation window.  The search has three stages:
%   1) global coarse scan,
%   2) integer-k local refinement,
%   3) continuous sub-grid minimization for harmonic/subharmonic branches.
%
% continuation_halfwidth and jump_tol_rel are retained in the function
% signature only for backward compatibility with older calls.

ns=numel(muS_list); kc=nan(size(muS_list)); ac=kc; winner=kc;

for is=1:ns
    surface=muS_list(is);

    % ---- Stage 1: GLOBAL coarse search at every C_S ----
    mcoarse = m_min_coarse:m_step_coarse:m_max_coarse;
    R=fwcs_search_k_list_CS(mcoarse,w_dim,cols,tol_im,N,Omega,nu,mu,rho,g,sigma0,h,beta,surface);

    % ---- Stage 2: integer-k refinement around coarse winner ----
    if isfinite(R.kc)
        m0=max(1,round(R.kc*w_dim/pi));
        mref=max(1,m0-refine_halfwidth):(m0+refine_halfwidth);
        R2=fwcs_search_k_list_CS(mref,w_dim,cols,tol_im,N,Omega,nu,mu,rho,g,sigma0,h,beta,surface);
        if isfinite(R2.kc), R=R2; end
    end

    % ---- Stage 3: continuous sub-grid refinement ----
    R=fwcs_refine_k_continuous_CS(R,continuous_halfwidth,continuous_tol_k, ...
        N,Omega,nu,mu,rho,g,sigma0,h,beta,surface,cols,tol_im);

    kc(is)=R.kc; ac(is)=R.ac; winner(is)=R.winner;
    fprintf('  surface=% .4e | CS(actual)=% .4e | kc=%.8f | ac/g=%.10g | winner=%g\n', ...
        surface,surface*kc(is)^2,kc(is),ac(is),winner(is));
end
end

function R = fwcs_search_k_list_CS(m_list,w_dim,cols,tol_im,N,Omega,nu,mu,rho,g,sigma0,h,beta,surface)
R=struct('kc',nan,'ac',nan,'winner',nan, ...
         'kh',nan,'ah',inf,'ks',nan,'as',inf);
best_h=inf; best_s=inf; kh=nan; ks=nan;
for ii=1:numel(m_list)
    k=m_list(ii)*pi/w_dim;
    oh=compute_eigenvalues(0,0,true,k,N,Omega,nu,mu,rho,g,sigma0,h,beta,0,0,surface);
    os=compute_eigenvalues(0,Omega/2,false,k,N,Omega,nu,mu,rho,g,sigma0,h,beta,0,0,surface);
    vh=fwcs_valid_positive_thresholds(oh,cols,tol_im); vs=fwcs_valid_positive_thresholds(os,cols,tol_im);
    if ~isempty(vh)
        ah=min(vh);
        if ah<best_h, best_h=ah; kh=k; end
    end
    if ~isempty(vs)
        as=min(vs);
        if as<best_s, best_s=as; ks=k; end
    end
end
R.kh=kh; R.ah=best_h; R.ks=ks; R.as=best_s;
if isfinite(best_h) && isfinite(best_s)
    if best_h<=best_s, R.kc=kh; R.ac=best_h; R.winner=1; else, R.kc=ks; R.ac=best_s; R.winner=2; end
elseif isfinite(best_h), R.kc=kh; R.ac=best_h; R.winner=1;
elseif isfinite(best_s), R.kc=ks; R.ac=best_s; R.winner=2;
end
end

function R = fwcs_refine_k_continuous_CS(R,halfwidth,tol_k,N,Omega,nu,mu,rho,g,sigma0,h,beta,surface,cols,tol_im)
% Refine harmonic and subharmonic branch minima independently in continuous k.
% This removes the integer-k staircase that generated the visible jaggedness.
opts=optimset('Display','off','TolX',tol_k,'MaxIter',100,'MaxFunEvals',250);

kh=R.kh; ah=R.ah; ks=R.ks; as=R.as;

if isfinite(kh)
    a=max(1e-12,kh-halfwidth); b=kh+halfwidth;
    obj=@(k) fwcs_threshold_at_k_CS(k,1,cols,tol_im,N,Omega,nu,mu,rho,g,sigma0,h,beta,surface);
    [ktry,atry]=fminbnd(obj,a,b,opts);
    if isfinite(atry) && atry<1e90, kh=ktry; ah=atry; end
end
if isfinite(ks)
    a=max(1e-12,ks-halfwidth); b=ks+halfwidth;
    obj=@(k) fwcs_threshold_at_k_CS(k,2,cols,tol_im,N,Omega,nu,mu,rho,g,sigma0,h,beta,surface);
    [ktry,atry]=fminbnd(obj,a,b,opts);
    if isfinite(atry) && atry<1e90, ks=ktry; as=atry; end
end

R.kh=kh; R.ah=ah; R.ks=ks; R.as=as;
if isfinite(ah) && isfinite(as)
    if ah<=as, R.kc=kh; R.ac=ah; R.winner=1; else, R.kc=ks; R.ac=as; R.winner=2; end
elseif isfinite(ah), R.kc=kh; R.ac=ah; R.winner=1;
elseif isfinite(as), R.kc=ks; R.ac=as; R.winner=2;
else, R.kc=nan; R.ac=nan; R.winner=nan;
end
end

function acrit = fwcs_threshold_at_k_CS(k,branch,cols,tol_im,N,Omega,nu,mu,rho,g,sigma0,h,beta,surface)
% branch = 1 harmonic, branch = 2 subharmonic
if ~(isfinite(k) && k>0), acrit=1e100; return; end
if branch==1
    out=compute_eigenvalues(0,0,true,k,N,Omega,nu,mu,rho,g,sigma0,h,beta,0,0,surface);
else
    out=compute_eigenvalues(0,Omega/2,false,k,N,Omega,nu,mu,rho,g,sigma0,h,beta,0,0,surface);
end
v=fwcs_valid_positive_thresholds(out,cols,tol_im);
if isempty(v), acrit=1e100; else, acrit=min(v); end
end

function vals = fwcs_valid_positive_thresholds(out,cols,tol_im)
if isempty(out), vals=[]; return; end
cu=cols(cols<=numel(out)); if isempty(cu), vals=[]; return; end
v=out(cu); good=isfinite(v)&real(v)>0&abs(imag(v))<tol_im; vals=real(v(good));
end

function [lambda,ok,resnorm] = fwcs_solve_free_mode_CS(lambda0,k,nu,mu,rho,g,sigma0,h,muS_tot,tol,maxiter)
x=[real(lambda0);abs(imag(lambda0))]; if x(2)<=0 || ~isfinite(x(2)), x(2)=1; end
ok=false;
for iter=1:maxiter
    F=fwcs_free_wave_residual_real_CS(x,k,nu,mu,rho,g,sigma0,h,muS_tot); fn=norm(F,2);
    if fn<tol, ok=true; break; end
    J=zeros(2,2);
    for j=1:2
        hj=1e-6*max(1,abs(x(j))); xp=x; xm=x; xp(j)=xp(j)+hj; xm(j)=xm(j)-hj;
        Fp=fwcs_free_wave_residual_real_CS(xp,k,nu,mu,rho,g,sigma0,h,muS_tot);
        Fm=fwcs_free_wave_residual_real_CS(xm,k,nu,mu,rho,g,sigma0,h,muS_tot);
        J(:,j)=(Fp-Fm)/(2*hj);
    end
    if ~all(isfinite(J(:))) || rcond(J)<1e-13, break; end
    dx=-J\F; accepted=false; step=1;
    for ls=1:12
        xt=x+step*dx; if xt(2)<=0, xt(2)=abs(xt(2)); end
        Ft=fwcs_free_wave_residual_real_CS(xt,k,nu,mu,rho,g,sigma0,h,muS_tot);
        if all(isfinite(Ft)) && norm(Ft,2)<fn, x=xt; accepted=true; break; end
        step=step/2;
    end
    if ~accepted, break; end
end
F=fwcs_free_wave_residual_real_CS(x,k,nu,mu,rho,g,sigma0,h,muS_tot);
if norm(F,2)>=tol
    obj=@(xx) sum(fwcs_free_wave_residual_real_CS([xx(1);abs(xx(2))],k,nu,mu,rho,g,sigma0,h,muS_tot).^2);
    opts=optimset('Display','off','TolX',1e-10,'TolFun',1e-14,'MaxIter',1000,'MaxFunEvals',3000);
    xtry=fminsearch(obj,x,opts); xtry(2)=abs(xtry(2)); Ftry=fwcs_free_wave_residual_real_CS(xtry,k,nu,mu,rho,g,sigma0,h,muS_tot);
    if norm(Ftry,2)<norm(F,2), x=xtry; F=Ftry; end
end
lambda=x(1)+1i*x(2); resnorm=norm(F,2); ok=resnorm<max(tol,1e-7);
end

function F = fwcs_free_wave_residual_real_CS(x,k,nu,mu,rho,g,sigma0,h,muS_tot)
lambda=x(1)+1i*x(2); R=fwcs_free_wave_residual_complex_CS(lambda,k,nu,mu,rho,g,sigma0,h,muS_tot);
rho_jump=rho(2)-rho(1); Krest=sigma0*k^4-rho_jump*g*k^2; Kscale=max(abs(Krest),1);
F=[real(R);imag(R)]/Kscale; if ~all(isfinite(F)), F=[1e20;1e20]; end
end

function R = fwcs_free_wave_residual_complex_CS(lambda,k,nu,mu,rho,g,sigma0,h,muS_tot)
D=fwcs_interface_derivatives_lambda_CS(lambda,k,nu,mu,h,muS_tot);
Ntotal=((rho(2)*lambda+3*k^2*mu(2))*D.dz2-mu(2)*D.dzzz2)-((rho(1)*lambda+3*k^2*mu(1))*D.dz1-mu(1)*D.dzzz1);
rho_jump=rho(2)-rho(1); Krest=sigma0*k^4-rho_jump*g*k^2; R=Ntotal-Krest;
end

function D = fwcs_interface_derivatives_lambda_CS(lambda,k,nu,mu,h,muS_tot)
Mn=k^2*muS_tot; q1=sqrt(k^2+lambda/nu(1)); q2=sqrt(k^2+lambda/nu(2));
if real(q1)<0, q1=-q1; end; if real(q2)<0, q2=-q2; end
KH_DEEP=8; USE_DEEP=min(abs(k*h(1)),abs(k*h(2)))>=KH_DEEP;
if abs(lambda)>1e-14
    if ~USE_DEEP
        A=fhs_general(k,mu,h,q1,q2,Mn); b=zeros(8,1); b(1)=lambda; c=fwcs_scaled_solve(A,b);
        dz1=[k,-k,q1,-q1]*c(1:4); dz2=[k,-k,q2,-q2]*c(5:8);
        dzzz1=[k^3,-k^3,q1^3,-q1^3]*c(1:4); dzzz2=[k^3,-k^3,q2^3,-q2^3]*c(5:8);
    else
        A=fhs_general_deep(k,mu,q1,q2,Mn); b=zeros(4,1); b(1)=lambda; c=fwcs_scaled_solve(A,b);
        dz1=[k,q1]*c(1:2); dz2=[-k,-q2]*c(3:4);
        dzzz1=[k^3,q1^3]*c(1:2); dzzz2=[-k^3,-q2^3]*c(3:4);
    end
else
    if ~USE_DEEP
        A=fhs_special(k,mu,h,Mn); b=zeros(8,1); b(1)=lambda; c=fwcs_scaled_solve(A,b);
        dz1=[k,-k,1,1]*c(1:4); dz2=[k,-k,1,1]*c(5:8);
        dzzz1=[k^3,-k^3,3*k^2,3*k^2]*c(1:4); dzzz2=[k^3,-k^3,3*k^2,3*k^2]*c(5:8);
    else
        A=fhs_special_deep(k,mu,Mn); b=zeros(4,1); b(1)=lambda; c=fwcs_scaled_solve(A,b);
        dz1=[k,1]*c(1:2); dz2=[-k,1]*c(3:4);
        dzzz1=[k^3,3*k^2]*c(1:2); dzzz2=[-k^3,3*k^2]*c(3:4);
    end
end
D=struct('dz1',dz1,'dz2',dz2,'dzzz1',dzzz1,'dzzz2',dzzz2,'q1',q1,'q2',q2,'Mn_S',Mn);
end

function x = fwcs_scaled_solve(A,b)
if ~all(isfinite(A(:))) || ~all(isfinite(b(:))), error('FreeWaveCS:NonFiniteHydroMatrix','Hydrodynamic matrix contains NaN/Inf.'); end
rn=sqrt(sum(abs(A),2)); cn=sqrt(sum(abs(A'),2)); rn(rn==0)=1; cn(cn==0)=1;
D1=diag(1./rn); D2=diag(1./cn); xs=(D1*A*D2)\(D1*b); x=D2*xs;
end

function lambda0 = fwcs_initial_free_wave_guess(k,rho,mu,g,sigma0,h)
den=rho(1)*coth(max(k*h(1),1e-12))+rho(2)*coth(max(k*h(2),1e-12));
omega2=((rho(1)-rho(2))*g*k+sigma0*k^3)/den; omega0=sqrt(max(real(omega2),1e-16));
nu_eff=(mu(1)+mu(2))/(rho(1)+rho(2)); gamma_guess=2*nu_eff*k^2; lambda0=-gamma_guess+1i*omega0;
end

% ======================== EIGENVECTOR DIAGNOSTIC HELPERS ========================
function out = eig_run_one_CM_sweep( ...
    beta_list, ...
    N,g_dim,omega_dim,rho_dim,mu_dim,nu_dim,gamma_dim,h_dim, ...
    Gamma0_dim,Ds_dim,sigma_values, ...
    w_dim,cols,tol_im, ...
    USE_REFINE,USE_CONTINUATION,ALLOW_FALLBACK, ...
    m_max_coarse,refine_halfwidth,jump_tol_rel,edge_tol_m, ...
    verbose_print)

nx = numel(beta_list);

kc_list = nan(nx,1);
ac_list = nan(nx,1);
winner  = nan(nx,1);

Gamma_RMS_norm = nan(nx,1);
tauM_RMS_norm  = nan(nx,1);
us_RMS_norm    = nan(nx,1);
dz_RMS_norm    = nan(nx,1);
zeta_RMS_check = nan(nx,1);

kc_prev = nan;

for ix = 1:nx

    beta = beta_list(ix);
    surface = 0;

    R = eig_find_critical_point( ...
        ix,kc_prev, ...
        N,g_dim,omega_dim,rho_dim,mu_dim,nu_dim,gamma_dim,h_dim, ...
        beta,Gamma0_dim,Ds_dim,surface,sigma_values, ...
        w_dim,cols,tol_im, ...
        USE_REFINE,USE_CONTINUATION,ALLOW_FALLBACK, ...
        m_max_coarse,refine_halfwidth,jump_tol_rel,edge_tol_m);

    if isfinite(R.kc) && isfinite(R.ac)

        kc_list(ix)=R.kc;
        ac_list(ix)=R.ac;
        winner(ix)=R.winner_flag;
        kc_prev=R.kc;

        if R.winner_flag==1
            alpha=0;
            is_harmonic=true;
        else
            alpha=omega_dim/2;
            is_harmonic=false;
        end

        P = eig_critical_CM_diagnostics( ...
            R.kc,R.ac,alpha,is_harmonic, ...
            N,omega_dim,nu_dim,mu_dim,rho_dim,g_dim,gamma_dim,h_dim, ...
            beta,Gamma0_dim,Ds_dim);

        Gamma_RMS_norm(ix)=P.Gamma_RMS;
        tauM_RMS_norm(ix)=P.tauM_RMS;
        us_RMS_norm(ix)=P.us_RMS;
        dz_RMS_norm(ix)=P.dz_RMS;
        zeta_RMS_check(ix)=P.zeta_RMS;

    else
        kc_prev=nan;
    end

    if verbose_print
        fprintf([ ...
            '[CM] beta=% .5e | kc=% .7g | ac/g=% .7g | winner=%g | ' ...
            'Gamma=% .5e | tauM=% .5e | us=% .5e | zetaRMS=% .5e\n'], ...
            beta,kc_list(ix),ac_list(ix),winner(ix), ...
            Gamma_RMS_norm(ix),tauM_RMS_norm(ix), ...
            us_RMS_norm(ix),zeta_RMS_check(ix));
    end
end

out = struct();

out.beta_list = beta_list(:);
out.kc_list   = kc_list;
out.ac_list   = ac_list;
out.winner    = winner;

out.Gamma_RMS_norm = Gamma_RMS_norm;
out.tauM_RMS_norm  = tauM_RMS_norm;
out.us_RMS_norm    = us_RMS_norm;
out.dz_RMS_norm    = dz_RMS_norm;
out.zeta_RMS_check = zeta_RMS_check;

end

function out = eig_run_one_CS_sweep( ...
    muS_list, ...
    N,g_dim,omega_dim,rho_dim,mu_dim,nu_dim,gamma_dim,h_dim, ...
    beta_dim,Gamma0_dim,Ds_dim, ...
    w_dim,cols,tol_im, ...
    USE_REFINE,USE_CONTINUATION,ALLOW_FALLBACK, ...
    m_max_coarse,refine_halfwidth,jump_tol_rel,edge_tol_m, ...
    verbose_print)

nx = numel(muS_list);
kc_list = nan(nx,1); ac_list = nan(nx,1); winner = nan(nx,1);
us_RMS_norm = nan(nx,1); tauS_RMS_norm = nan(nx,1);
dz_RMS_norm = nan(nx,1); zeta_RMS_check = nan(nx,1);
kc_prev = nan;

for ix = 1:nx
    surface = muS_list(ix);

    R = eig_find_critical_point( ...
        ix,kc_prev, ...
        N,g_dim,omega_dim,rho_dim,mu_dim,nu_dim,gamma_dim,h_dim, ...
        beta_dim,Gamma0_dim,Ds_dim,surface,0, ...
        w_dim,cols,tol_im, ...
        USE_REFINE,USE_CONTINUATION,ALLOW_FALLBACK, ...
        m_max_coarse,refine_halfwidth,jump_tol_rel,edge_tol_m);

    if isfinite(R.kc) && isfinite(R.ac)
        kc_list(ix)=R.kc; ac_list(ix)=R.ac; winner(ix)=R.winner_flag; kc_prev=R.kc;
        if R.winner_flag==1
            alpha=0; is_harmonic=true;
        else
            alpha=omega_dim/2; is_harmonic=false;
        end

        P = eig_critical_CS_diagnostics( ...
            R.kc,R.ac,alpha,is_harmonic, ...
            N,omega_dim,nu_dim,mu_dim,rho_dim,g_dim,gamma_dim,h_dim,surface);

        us_RMS_norm(ix)=P.us_RMS;
        tauS_RMS_norm(ix)=P.tauS_RMS;
        dz_RMS_norm(ix)=P.dz_RMS;
        zeta_RMS_check(ix)=P.zeta_RMS;
    else
        kc_prev=nan;
    end

    if verbose_print
        CS_now = surface*kc_list(ix)^2;
        fprintf(['[CS] muS=% .5e | CS=% .5e | kc=% .7g | ac/g=% .7g | winner=%g | ' ...
            'tauS=% .5e | us=% .5e | zetaRMS=% .5e\n'], ...
            surface,CS_now,kc_list(ix),ac_list(ix),winner(ix), ...
            tauS_RMS_norm(ix),us_RMS_norm(ix),zeta_RMS_check(ix));
    end
end

out = struct();
out.muS_list = muS_list(:);
out.kc_list = kc_list;
out.ac_list = ac_list;
out.winner = winner;
out.us_RMS_norm = us_RMS_norm;
out.tauS_RMS_norm = tauS_RMS_norm;
out.dz_RMS_norm = dz_RMS_norm;
out.zeta_RMS_check = zeta_RMS_check;
end

function P = eig_critical_CS_diagnostics( ...
    k,ac_over_g,alpha,is_harmonic, ...
    N,omega,nu,mu,rho,g,gamma,h,surface)

E = eig_get_critical_eigenvector( ...
    k,ac_over_g,alpha,is_harmonic, ...
    N,omega,nu,mu,rho,g,gamma,h,0,0,0,surface);

nkeep = N+1;
dz = complex(zeros(nkeep,1));
us = complex(zeros(nkeep,1));
tauS = complex(zeros(nkeep,1));

for ij = 1:nkeep
    n = ij-1;
    Dn = eig_interface_derivatives_per_unit_zeta( ...
        0,alpha,n,omega,k,nu,mu,rho,g,gamma,h,0,0,0,surface);
    zn = E.zeta(ij);
    dz1 = Dn.dz1_over_zeta*zn;
    dz2 = Dn.dz2_over_zeta*zn;
    dz(ij)=0.5*(dz1+dz2);
    us(ij)=1i*dz(ij)/k;
    tauS(ij)=-Dn.Mn_S*us(ij);
end

P = struct();
P.zeta_RMS = eig_floquet_norm(E.zeta,E.weights);
P.dz_RMS   = eig_floquet_norm(dz,E.weights);
P.us_RMS   = eig_floquet_norm(us,E.weights);
P.tauS_RMS = eig_floquet_norm(tauS,E.weights);
end

function P = eig_critical_CM_diagnostics( ...
    k,ac_over_g,alpha,is_harmonic, ...
    N,omega,nu,mu,rho,g,gamma,h,beta,Gamma0,Ds)

E = eig_get_critical_eigenvector( ...
    k,ac_over_g,alpha,is_harmonic, ...
    N,omega,nu,mu,rho,g,gamma,h,beta,Gamma0,Ds,0);

nkeep = N+1;

dz    = complex(zeros(nkeep,1));
us    = complex(zeros(nkeep,1));
Gamma = complex(zeros(nkeep,1));
tauM  = complex(zeros(nkeep,1));

for ij = 1:nkeep

    n = ij-1;

    Dn = eig_interface_derivatives_per_unit_zeta( ...
        0,alpha,n,omega,k,nu,mu,rho,g,gamma,h, ...
        beta,Gamma0,Ds,0);

    zn = E.zeta(ij);

    dz1 = Dn.dz1_over_zeta*zn;
    dz2 = Dn.dz2_over_zeta*zn;

    dz(ij)=0.5*(dz1+dz2);

    us(ij)=1i*dz(ij)/k;

    denom = Dn.fexp + Ds*k^2;

    if Gamma0~=0 && abs(denom)>1e-30
        Gamma(ij)=Gamma0*dz(ij)/denom;
    else
        Gamma(ij)=0;
    end

    % physical Marangoni tangential traction
    tauM(ij) = -1i*k*gamma*beta/(1-Gamma0)*Gamma(ij);
end

P = struct();

P.zeta_RMS = eig_floquet_norm(E.zeta,E.weights);
P.dz_RMS   = eig_floquet_norm(dz,E.weights);
P.us_RMS   = eig_floquet_norm(us,E.weights);
P.Gamma_RMS = eig_floquet_norm(Gamma,E.weights);
P.tauM_RMS  = eig_floquet_norm(tauM,E.weights);

end

function E = eig_get_critical_eigenvector( ...
    k,ac_over_g,alpha,is_harmonic, ...
    N,omega,nu,mu,rho,g,gamma,h,beta,Gamma0,Ds,surface)

sigma=0;
nkeep=N+1;

An_block=cell(nkeep,1);

for ij=1:nkeep

    n=ij-1;

    An=matrix_An( ...
        sigma,alpha,n,omega,k,nu,mu,rho,g,gamma,h, ...
        beta,Gamma0,Ds,surface);

    An_block{ij}=[real(An),-imag(An); ...
                  imag(An), real(An)];
end

A=blkdiag(An_block{:});
B=banded_matrix_B(N,is_harmonic);

Gmat=0.5*(rho(2)-rho(1))*k^2*g*full(B);

[V,D]=eig(A,Gmat,'vector');

valid=isfinite(D) & abs(imag(D))<1e-8;

if ~any(valid)
    error('No valid real generalized eigenvalue found.');
end

Dvalid=real(D(valid));
Vvalid=V(:,valid);

[~,jj]=min(abs(abs(Dvalid)-ac_over_g));
x=Vvalid(:,jj);

if numel(x)<2*nkeep
    error('Unexpected generalized eigenvector length.');
end

zeta=complex(x(1:2:2*nkeep),x(2:2:2*nkeep));

if is_harmonic
    weights=[1;2*ones(N,1)];
else
    weights=2*ones(nkeep,1);
end

zeta_RMS_raw=eig_floquet_norm(zeta,weights);

if ~isfinite(zeta_RMS_raw) || zeta_RMS_raw<1e-30
    error('Critical eigenvector has invalid zeta norm.');
end

% consistent normalization
zeta=zeta/zeta_RMS_raw;

% fix arbitrary global phase for reproducible visualization
if abs(zeta(1))>1e-30
    zeta=zeta*exp(-1i*angle(zeta(1)));
end

E=struct();

E.zeta=zeta;
E.weights=weights;
E.zeta_RMS=eig_floquet_norm(zeta,weights);
E.zeta_RMS_raw=zeta_RMS_raw;
E.alpha=alpha;
E.is_harmonic=is_harmonic;
E.k=k;
E.ac_over_g=ac_over_g;

end

function val = eig_floquet_norm(q,weights)

q=q(:);
weights=weights(:);

val=sqrt(sum(weights.*abs(q).^2));

end

function D = eig_interface_derivatives_per_unit_zeta( ...
    sigma,alpha,n,omega,k,nu,mu,rho,g,gamma,h,beta,Gamma0,Ds,surface)

fexp=sigma+1i*(alpha+n*omega);
sigma0=gamma;

Mn_M=0;

if beta~=0 && Gamma0~=0
    Mn_M=beta*Gamma0*sigma0*k^2 / ...
        ((fexp+Ds*k^2)*(1-Gamma0));
end

Mn_S=k^2*surface;
Mn=Mn_M+Mn_S;

q1=sqrt(k^2+fexp/nu(1));
q2=sqrt(k^2+fexp/nu(2));

KH_DEEP=8.0;
USE_DEEP=min(abs(k*h(1)),abs(k*h(2)))>=KH_DEEP;

if abs(fexp)>0

    if ~USE_DEEP

        lhs=fhs_general(k,mu,h,q1,q2,Mn);

        rhs=zeros(8,1);
        rhs(1)=fexp;

        c=eig_scaled_solve(lhs,rhs);

        dz1=[k,-k,q1,-q1]*c(1:4);
        dz2=[k,-k,q2,-q2]*c(5:8);

    else

        lhs=fhs_general_deep(k,mu,q1,q2,Mn);

        rhs=zeros(4,1);
        rhs(1)=fexp;

        c=eig_scaled_solve(lhs,rhs);

        dz1=[k,q1]*c(1:2);
        dz2=[-k,-q2]*c(3:4);
    end

else

    if ~USE_DEEP

        lhs=fhs_special(k,mu,h,Mn);

        rhs=zeros(8,1);
        rhs(1)=fexp;

        c=eig_scaled_solve(lhs,rhs);

        dz1=[k,-k,1,1]*c(1:4);
        dz2=[k,-k,1,1]*c(5:8);

    else

        lhs=fhs_special_deep(k,mu,Mn);

        rhs=zeros(4,1);
        rhs(1)=fexp;

        c=eig_scaled_solve(lhs,rhs);

        dz1=[k,1]*c(1:2);
        dz2=[-k,1]*c(3:4);
    end
end

D=struct();

D.dz1_over_zeta=dz1;
D.dz2_over_zeta=dz2;
D.fexp=fexp;
D.Mn_M=Mn_M;
D.Mn_S=Mn_S;
D.Mn=Mn;

end

function x = eig_scaled_solve(A,b)

if ~all(isfinite(A(:))) || ~all(isfinite(b(:)))
    error('CMViz:NonFiniteHydroMatrix', ...
        'Hydrodynamic matrix contains NaN/Inf.');
end

row_norm=sqrt(sum(abs(A),2));
col_norm=sqrt(sum(abs(A'),2));

row_norm(row_norm==0)=1;
col_norm(col_norm==0)=1;

D1=diag(1./row_norm);
D2=diag(1./col_norm);

As=D1*A*D2;
bs=D1*b;

xs=As\bs;
x=D2*xs;

end

function R = eig_find_critical_point( ...
    ix,kc_prev, ...
    N,g_dim,omega_dim,rho_dim,mu_dim,nu_dim,gamma,h_dim, ...
    beta,Gamma0,Ds,surface,sigma_values, ...
    w_dim,cols,tol_im, ...
    USE_REFINE,USE_CONTINUATION,ALLOW_FALLBACK, ...
    m_max_coarse,refine_halfwidth,jump_tol_rel,edge_tol_m)

R=struct('kc',nan,'ac',nan,'winner_flag',nan);

use_local_first= ...
    USE_CONTINUATION && ix>1 && isfinite(kc_prev);

result_ok=false;

if use_local_first

    m0=max(1,round(kc_prev*w_dim/pi));
    m_lo=max(1,m0-refine_halfwidth);
    m_hi=m0+refine_halfwidth;

    [a_h,a_s]=eig_sweep_over_m_list( ...
        m_lo:m_hi,sigma_values,w_dim,N,omega_dim, ...
        nu_dim,mu_dim,rho_dim,g_dim,gamma,h_dim, ...
        beta,Gamma0,Ds,surface);

    R=eig_extract_global_critical_from_maps(a_h,a_s,cols,tol_im);

    result_ok=isfinite(R.kc) && isfinite(R.ac);

    need_fallback=~result_ok;

    if result_ok

        rel_jump=abs(R.kc-kc_prev)/max(abs(kc_prev),1e-12);

        if rel_jump>jump_tol_rel
            need_fallback=true;
        end

        mcrit=round(R.kc*w_dim/pi);

        if (mcrit-m_lo<=edge_tol_m) || (m_hi-mcrit<=edge_tol_m)
            need_fallback=true;
        end
    end

    if need_fallback && ALLOW_FALLBACK

        [a_h,a_s]=eig_sweep_over_m_list( ...
            100:100:m_max_coarse,sigma_values,w_dim,N,omega_dim, ...
            nu_dim,mu_dim,rho_dim,g_dim,gamma,h_dim, ...
            beta,Gamma0,Ds,surface);

        R=eig_extract_global_critical_from_maps(a_h,a_s,cols,tol_im);

        result_ok=isfinite(R.kc) && isfinite(R.ac);
    end

else

    [a_h,a_s]=eig_sweep_over_m_list( ...
        100:100:m_max_coarse,sigma_values,w_dim,N,omega_dim, ...
        nu_dim,mu_dim,rho_dim,g_dim,gamma,h_dim, ...
        beta,Gamma0,Ds,surface);

    R=eig_extract_global_critical_from_maps(a_h,a_s,cols,tol_im);

    result_ok=isfinite(R.kc) && isfinite(R.ac);
end

if result_ok && USE_REFINE

    m0=max(1,round(R.kc*w_dim/pi));
    m_lo=max(1,m0-refine_halfwidth);
    m_hi=m0+refine_halfwidth;

    [a_h,a_s]=eig_sweep_over_m_list( ...
        m_lo:m_hi,sigma_values,w_dim,N,omega_dim, ...
        nu_dim,mu_dim,rho_dim,g_dim,gamma,h_dim, ...
        beta,Gamma0,Ds,surface);

    R_ref=eig_extract_global_critical_from_maps(a_h,a_s,cols,tol_im);

    if isfinite(R_ref.kc) && isfinite(R_ref.ac)
        R=R_ref;
    end
end

end

function [a_harmonic,a_subharmonic] = eig_sweep_over_m_list( ...
    m_list,sigma_values,w_dim,N,omega,nu,mu,rho,g,gamma,h, ...
    beta,Gamma0,Ds,surface)

nm=numel(m_list);

k0=m_list(1)*pi/w_dim;

out_h0=compute_eigenvalues( ...
    sigma_values,0,true, ...
    k0,N,omega,nu,mu,rho,g,gamma,h,beta,Gamma0,Ds,surface);

out_s0=compute_eigenvalues( ...
    sigma_values,omega/2,false, ...
    k0,N,omega,nu,mu,rho,g,gamma,h,beta,Gamma0,Ds,surface);

nh=size(out_h0,2);
ns=size(out_s0,2);

a_harmonic=nan(nm,nh);
a_subharmonic=nan(nm,ns);

for i=1:nm

    k=m_list(i)*pi/w_dim;

    out_h=compute_eigenvalues( ...
        sigma_values,0,true, ...
        k,N,omega,nu,mu,rho,g,gamma,h,beta,Gamma0,Ds,surface);

    out_s=compute_eigenvalues( ...
        sigma_values,omega/2,false, ...
        k,N,omega,nu,mu,rho,g,gamma,h,beta,Gamma0,Ds,surface);

    a_harmonic(i,1:numel(out_h))=out_h;
    a_subharmonic(i,1:numel(out_s))=out_s;

    a_harmonic(i,1)=k;
    a_subharmonic(i,1)=k;
end

end

function R = eig_extract_global_critical_from_maps( ...
    a_harmonic,a_subharmonic,cols,tol_im)

R=struct( ...
    'kc',nan,'ac',nan, ...
    'ah_min',inf,'as_min',inf, ...
    'kh_min',nan,'ks_min',nan, ...
    'winner_flag',nan);

if isempty(a_harmonic) || isempty(a_subharmonic)
    return;
end

k=a_harmonic(:,1);
k=k(:);

cols_h=cols(cols<=size(a_harmonic,2));
cols_s=cols(cols<=size(a_subharmonic,2));

if isempty(cols_h) || isempty(cols_s)
    return;
end

Yh=a_harmonic(:,cols_h);
Ys=a_subharmonic(:,cols_s);

Yh(~isfinite(Yh) | real(Yh)<=0 | abs(imag(Yh))>tol_im)=NaN;
Ys(~isfinite(Ys) | real(Ys)<=0 | abs(imag(Ys))>tol_im)=NaN;

Yh_vec=real(Yh(:));
Ys_vec=real(Ys(:));

Kh_vec=repmat(k,numel(cols_h),1);
Ks_vec=repmat(k,numel(cols_s),1);

vh=isfinite(Yh_vec);
vs=isfinite(Ys_vec);

Yh_vec=Yh_vec(vh);
Ys_vec=Ys_vec(vs);

Kh_vec=Kh_vec(vh);
Ks_vec=Ks_vec(vs);

if ~isempty(Yh_vec)
    [R.ah_min,ih]=min(Yh_vec);
    R.kh_min=Kh_vec(ih);
end

if ~isempty(Ys_vec)
    [R.as_min,isub]=min(Ys_vec);
    R.ks_min=Ks_vec(isub);
end

if isfinite(R.ah_min) && isfinite(R.as_min)

    if R.ah_min<=R.as_min
        R.ac=R.ah_min;
        R.kc=R.kh_min;
        R.winner_flag=1;
    else
        R.ac=R.as_min;
        R.kc=R.ks_min;
        R.winner_flag=2;
    end

elseif isfinite(R.ah_min)

    R.ac=R.ah_min;
    R.kc=R.kh_min;
    R.winner_flag=1;

elseif isfinite(R.as_min)

    R.ac=R.as_min;
    R.kc=R.ks_min;
    R.winner_flag=2;
end

end

function CM = eig_compute_CM_axis( ...
    beta_list,kc_list,gamma,Gamma0,Ds,omega)

beta_list=beta_list(:);
kc_list=kc_list(:);

CM=nan(size(beta_list));

for i=1:numel(beta_list)

    k=kc_list(i);

    if ~isfinite(k) || k<=0
        continue;
    end

    Pe=omega/(Ds*k^2);

    prefactor=gamma*Gamma0/(1-Gamma0) * ...
        k^2/(omega*sqrt((1/Pe)^2+(1/2)^2));

    CM(i)=prefactor*beta_list(i);
end

end
