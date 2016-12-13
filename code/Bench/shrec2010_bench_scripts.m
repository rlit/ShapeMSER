%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  vertex weight  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% vertex weight - t = 2048
params = GetDefaultBenchParams();
params.is_ew   = 0;
params.vw_fun    = 'heat_kernel';
params.removeZeroEval = 1;
params.timeVals	 = 2048;
params.showPlots = 1;
params.preCalcMsers = 1;
params.load_msers = 0;
params.benchName = ['vw+ t=' int2str(params.timeVals)];
try shrec2010_bench(params);catch me;me = struct(me);save([params.benchName '.FAIL.mat'],'me');end

%% vertex weight - CT
params = GetDefaultBenchParams();
params.is_ew   = 0;
params.vw_fun  = 'CT';
params.removeZeroEval = 1;
params.showPlots = 1;
params.preCalcMsers = 1;
params.load_msers = 0;
params.benchName = ['vw+ ' params.vw_fun];
try shrec2010_bench(params);catch me;me = struct(me);save([params.benchName '.FAIL.mat'],'me');end 

%% vertex weight - siHks_at0
params = GetDefaultBenchParams();
params.is_ew   = 0;
params.vw_fun  = 'siHks_at0';
params.removeZeroEval = 1;
params.showPlots = 1;
params.preCalcMsers = 1;
params.load_msers = 0;
params.benchName = ['vw+ ' params.vw_fun];
try shrec2010_bench(params);catch me;me = struct(me);save([params.benchName '.FAIL.mat'],'me');end 

%% vertex weight - siHksNorm
% params = GetDefaultBenchParams();
% params.is_ew   = 0;
% params.vw_fun  = 'siHksNorm';
% params.removeZeroEval = 1;
% params.showPlots = 1;
% params.preCalcMsers = 1;
% params.load_msers = 0;
% params.benchName = ['vw+ ' params.vw_fun];
% try shrec2010_bench(params);catch me;me = struct(me);save([params.benchName '.FAIL.mat'],'me');end 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  edge weight  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% edge weight - geodesic
params = GetDefaultBenchParams();
params.is_ew   = 1;
params.ew_distFun  = 'geodesic';
params.removeZeroEval = 1;
params.showPlots = 1;
params.preCalcMsers = 1;
params.benchName = ['ew ' params.ew_clusterType ' ' params.ew_distFun ];
try shrec2010_bench(params);catch me;me = struct(me);save([params.benchName '.FAIL.mat'],'me');end 

%% edge weight - |h(x,x)-h(y,y)| , t = 2048  scalar_HK
params = GetDefaultBenchParams();
params.is_ew   = 1;
params.timeVals	 = 2048;
params.mser_filters.MaxScore = 10^5.4;
params.ew_distFun  = 'scalar_HK';
params.removeZeroEval = 1;
params.showPlots = 1;
params.preCalcMsers = 1;
params.benchName = ['ew ' params.ew_clusterType ' ' params.ew_distFun ' t=' int2str(params.timeVals)];
try shrec2010_bench(params);catch me;me = struct(me);save([params.benchName '.FAIL.mat'],'me');end 

%% edge weight - |c(x,x)-c(y,y)|  scalar_CT
params = GetDefaultBenchParams();
params.is_ew   = 1;
params.timeVals	 = 2048;
params.mser_filters.MaxScore = 10^2;
params.ew_distFun  = 'scalar_CT';
params.removeZeroEval = 1;
params.showPlots = 1;
params.preCalcMsers = 1;
params.benchName = ['ew ' params.ew_clusterType ' ' params.ew_distFun ];
try shrec2010_bench(params);catch me;me = struct(me);save([params.benchName '.FAIL.mat'],'me');end 

%% edge weight - |h(x,*)-h(y,*)| , t = 2048  diffusion_distance
% params = GetDefaultBenchParams();
% params.is_ew   = 1;
% params.timeVals	 = 2048;
% params.mser_filters.MaxScore = 10^6.75;
% params.ew_distFun  = 'diffusion_distance';
% params.removeZeroEval = 1;
% params.showPlots = 1;
% params.preCalcMsers = 1;
% params.benchName = ['ew ' params.ew_clusterType ' ' params.ew_distFun ' t=' int2str(params.timeVals)];
% try shrec2010_bench(params);catch me;me = struct(me);save([params.benchName '.FAIL.mat'],'me');end 
%% edge weight - time_integral
params = GetDefaultBenchParams();
params.is_ew   = 1;
params.timeVals	 = [128 32768];
params.mser_filters.MaxScore = 10^7.2;
params.ew_distFun  = 'time_integral';
params.removeZeroEval = 1;
params.showPlots = 1;
params.preCalcMsers = 1;
params.benchName = ['ew ' params.ew_clusterType ' ' params.ew_distFun ' t=' int2str(params.timeVals)];
try shrec2010_bench(params);catch me;me = struct(me);save([params.benchName '.FAIL.mat'],'me');end 


%% edge weight - inverse_HKS, t = 2048
params = GetDefaultBenchParams();
params.is_ew   = 1;
params.timeVals	 = 2048;
params.ew_distFun  = 'inverse_HKS';
params.removeZeroEval = 1;
params.showPlots = 1;
params.preCalcMsers = 1;
params.benchName = ['ew ' params.ew_clusterType ' ' params.ew_distFun ' t=' int2str(params.timeVals)];
try shrec2010_bench(params);catch me;me = struct(me);save([params.benchName '.FAIL.mat'],'me');end 

%% inverse_CT
params = GetDefaultBenchParams();
params.is_ew   = 1;
params.ew_distFun  = 'inverse_CT';
params.mser_filters.MaxScore = 10^0;
params.removeZeroEval = 1;
params.showPlots = 1;
params.preCalcMsers = 1;
params.benchName = ['ew ' params.ew_clusterType ' ' params.ew_distFun ];
try shrec2010_bench(params);catch me;me = struct(me);save([params.benchName '.FAIL.mat'],'me');end 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  edge weight - SCALE INVARIANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% edge weight - inverse_SIHKS_diff, t = [...]
% params = GetDefaultBenchParams();
% params.is_ew   = 1;
% params.timeVals	 = 2.^(1: 1/16 : 25-1/16);
% params.ew_distFun  = 'inverse_SIHKS_diff';
% params.removeZeroEval = 0;
% params.showPlots = 1;
% params.preCalcMsers = 1;
% params.load_msers = 0;
% params.benchName = ['ew ' params.ew_clusterType ' ' params.ew_distFun ];
% try shrec2010_bench(params);catch me;me = struct(me);save([params.benchName '.FAIL.mat'],'me');end 

%% edge weight - inverse_SIHKS, t = [...]
params = GetDefaultBenchParams();
params.is_ew   = 1;
params.timeVals	 = 2.^(1: 1/16 : 25-1/16);
params.ew_distFun  = 'inverse_SIHKS';
params.mser_filters.MaxScore = 10^2;
params.removeZeroEval = 0;
params.showPlots = 1;
params.preCalcMsers = 1;
params.load_msers = 0;
params.benchName = ['ew ' params.ew_clusterType ' ' params.ew_distFun ];
try shrec2010_bench(params);catch me;me = struct(me);save([params.benchName '.FAIL.mat'],'me');end 

%% edge weight - SIHKS_at0
% params = GetDefaultBenchParams();
% params.is_ew   = 1;
% params.timeVals	 = 2.^(1: 1/16 : 25-1/16);
% params.ew_distFun  = 'SIHKS_at0';
% params.showPlots = 1;
% params.preCalcMsers = 1;
% params.load_msers = 0;
% params.benchName = ['ew ' params.ew_clusterType ' ' params.ew_distFun ];
% try shrec2010_bench(params);catch me;me = struct(me);save([params.benchName '.FAIL.mat'],'me');end 

%% edge weight - siHksNorm, all freq
% params = GetDefaultBenchParams();
% params.is_ew   = 1;
% params.timeVals	 = 2.^(1: 1/16 : 25-1/16);
% params.ew_distFun  = 'siHksNorm';
% params.showPlots = 1;
% params.preCalcMsers = 1;
% params.load_msers = 0;
% params.benchName = ['ew ' params.ew_clusterType ' ' params.ew_distFun ];
% try shrec2010_bench(params);catch me;me = struct(me);save([params.benchName '.FAIL.mat'],'me');end 

