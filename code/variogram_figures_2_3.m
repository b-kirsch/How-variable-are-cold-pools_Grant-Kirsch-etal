%% variogram_figures_2_3
% script to read in variogram data from various simulations
% and create figures 2-3 of GRL manuscript, "How variable are cold pools?"
%
% written by: Leah Grant (leah.grant@colostate.edu), Sept 2023
%
% tested with MATLAB_R2018a
% 
% needs access to some extra functions:  boundedline, multipan, 
% set_figsize, and colormap arrays redblu, orangepurple, and parula

clear

% set up some plotting defaults in case there is no default startup file
% ------
% default axes font size and line width and text font size
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultLineLineWidth',4)
set(0,'DefaultTextFontSize',16)

% figure background white
set(0,'DefaultFigureColor','w')

% add path to extra functions needed for this script
addpath(genpath('./matlab_functions'))


for v=1:38 % loop over variogram data files

    % input directory containing nc files, relative to current directory
    nc_dir = '../data/';
    % input directory containing nc files for Jennie's ensemble
    ensemble_dir = [nc_dir,'Jennie_ensemble/'];
    
    if v==1 % Jogi obs case
        name{v} = 'OBSJogi'; c=name{v};
        ttl.(c) = 'OBS-Jogi';
        nc_fname = 'variogram_obs_20210629.nc';
        dx.(c) = 0.; % not relevant for obs
        bin_lim.(c) = 15.; % km
        
    elseif v==2 % ICON 75m run
        name{v} = 'CSJogi75'; c=name{v};
        ttl.(c) = 'CS-Jogi-75m';
        nc_fname = 'variogram_icon_75m_20210629.nc';
        dx.(c) = 75.;
        bin_lim.(c) = 10.; % km
        
    elseif v==3 % ICON 156m run
        name{v} = 'CSJogi156'; c=name{v};
        ttl.(c) = 'CS-Jogi-156m';
        nc_fname = 'variogram_icon_156m_20210629.nc';
        dx.(c) = 156.; 
        bin_lim.(c) = 15.;
        
    elseif v==4 % ICON 312m run
        name{v} = 'CSJogi312'; c=name{v};
        ttl.(c) = 'CS-Jogi-312m';
        nc_fname = 'variogram_icon_312m_20210629.nc';
        dx.(c) = 312.; 
        bin_lim.(c) = 15.;
        
    elseif v==5 % ICON 625m run
        name{v} = 'CSJogi625'; c=name{v};
        ttl.(c) = 'CS-Jogi-625m';
        nc_fname = 'variogram_icon_625m_20210629.nc';
        dx.(c) = 312.; 
        bin_lim.(c) = 15.;
        
    elseif v==6 % RAMS idealized dry BL 50 m no shear case (Leah's simulation)
        name{v} = 'IDEALdryBL50'; c=name{v};
        ttl.(c) = 'IDEAL-DryBL-50m';
        nc_fname = 'variograms_IDEAL-DryBL-50.nc';
        dx.(c) = 50.; 
        bin_lim.(c) = 15.;
        
    elseif v==7 % RAMS idealized dry BL 100 m no shear case (Leah's simulation)
        name{v} = 'IDEALdryBL100'; c=name{v};
        ttl.(c) = 'IDEAL-DryBL-100m';
        nc_fname = 'variograms_IDEAL-DryBL-100.nc';
        dx.(c) = 100.; 
        bin_lim.(c) = 15.;
        
    elseif v==8 % RAMS idealized dry BL 100 m w/ shear case, downshear network placement (Christine's simulation)
        name{v} = 'IDEALdownshear100'; c=name{v};
        ttl.(c) = 'IDEAL-DownShear-100m';
        nc_fname = 'IDEAL-Downshear-100m-variogram_data_v04.nc';
        dx.(c) = 100.; 
        bin_lim.(c) = 15.;
        
    elseif v==9 % RAMS idealized dry BL 250 m w/ shear case, downshear network placement (Christine's simulation)
        name{v} = 'IDEALdownshear250'; c=name{v};
        ttl.(c) = 'IDEAL-DownShear-250m';
        nc_fname = 'IDEAL-Downshear-250m-variogram_data_v04.nc';
        dx.(c) = 250.; 
        bin_lim.(c) = 15.;
        
    elseif v==37 % RAMS idealized dry BL 100 m w/ shear case, upshear network placement (Christine's simulation)
        name{v} = 'IDEALupshear100'; c=name{v};
        ttl.(c) = 'IDEAL-UpShear-100m';
        nc_fname = 'IDEAL-Upshear-100m-variogram_data_v04.nc';
        dx.(c) = 100.; 
        bin_lim.(c) = 15.;
        
    elseif v==38 % RAMS idealized dry BL 250 m w/ shear case, upshear network placement (Christine's simulation)
        name{v} = 'IDEALupshear250'; c=name{v};
        ttl.(c) = 'IDEAL-UpShear-250m';
        nc_fname = 'IDEAL-Upshear-250m-variogram_data_v04.nc';
        dx.(c) = 250.; 
        bin_lim.(c) = 15.;
        
    elseif v==10 % RAMS tropical oceanic case study 100 m (Nick's simulation)
        name{v} = 'CStropOce100'; c=name{v};
        ttl.(c) = 'CS-TropOce-100m';
        nc_fname = 'CS-TropOce-100_variogram-data_v02.nc';
        dx.(c) = 100.; 
        bin_lim.(c) = 15.;
        
    elseif v==11 % RAMS tropical oceanic case study 100 m (Nick's simulation)
        name{v} = 'CStropOce1000'; c=name{v};
        ttl.(c) = 'CS-TropOce-1km';
        nc_fname = 'CS-TropOce-1000_variogram-data_v02.nc';
        dx.(c) = 1000.; 
        bin_lim.(c) = 15.;
        
    elseif v==12 % ICON full obs dataset
        name{v} = 'OBS'; c=name{v};
        ttl.(c) = 'OBS';
        nc_fname = 'variogram_obs_20210517-20210827.nc';
        dx.(c) = 0.; % not relevant for obs
        bin_lim.(c) = 15.;
        
    elseif v>=13 & v<=36 % ensemble of Jennie's simulations
        % set up directories
        ensref = 13; % reference for ensemble numbering - should match v above
        ensemble_files = dir([ensemble_dir,'*nc']);
        ifile = v-ensref+1; % file counter
        nc_fname = ensemble_files(ifile).name;
        name{v} = ['IDEALhaboob',nc_fname(1:end-3)]; c=name{v};
        ttl.(c) = ['IDEAL-Haboob-',nc_fname(1:end-3)];
            ttl.(c)(ttl.(c)=='_')='-';
            % rename the initial temperature perturbation ensemble
            if strcmp(ttl.(c)(end-6:end),'CP-day1'); ttl.(c)='IDEAL-Haboob-20KDay-150m';
            elseif strcmp(ttl.(c)(end-6:end),'CP-day2'); ttl.(c)='IDEAL-Haboob-17KDay-150m';
            elseif strcmp(ttl.(c)(end-6:end),'CP-day3'); ttl.(c)='IDEAL-Haboob-13KDay-150m';
            elseif strcmp(ttl.(c)(end-6:end),'CP-day4'); ttl.(c)='IDEAL-Haboob-10KDay-150m';
            elseif strcmp(ttl.(c)(end-8:end),'CP-night1'); ttl.(c)='IDEAL-Haboob-20KNight-150m';
            elseif strcmp(ttl.(c)(end-8:end),'CP-night2'); ttl.(c)='IDEAL-Haboob-17KNight-150m';
            elseif strcmp(ttl.(c)(end-8:end),'CP-night3'); ttl.(c)='IDEAL-Haboob-13KNight-150m';
            elseif strcmp(ttl.(c)(end-8:end),'CP-night4'); ttl.(c)='IDEAL-Haboob-10KNight-150m';
            end
        dx.(c) = 150.; 
        bin_lim.(c) = 15.;
        % re-write nc_dir for ensemble
        nc_dir = ensemble_dir;
        
    end


    %% read variograms from nc file
    
    nc_file = [nc_dir,nc_fname];
    c = name{v}; % short character for making structure arrays
    
    time.(c) = ncread( nc_file, 'time' );
    % convert time from integer seconds since Jan 1 1970 0Z to datetime array
    time.(c) = datetime( time.(c), 'convertfrom', 'posixtime' );
    
    bin_mids.(c) = ncread( nc_file, 'distance' );
    bin_units.(c) = ncreadatt( nc_file, 'distance', 'units' );
    if ~strcmp(bin_units.(c),'m'); 
        error(['check bin_units for case ',c]); 
    end
    % change bin_mids to km
    bin_mids.(c) = bin_mids.(c)/1000.; bin_units.(c) = 'km';
    
    variogram.(c) = ncread( nc_file, 'variogram' );
    variogram_units.(c) = ncreadatt( nc_file, 'variogram', 'units' );
    % add ^ so 2 becomes a superscript in the string label
    variogram_units.(c) = [variogram_units.(c)(1:end-1),'^',variogram_units.(c)(end)];
    
    cp_flag.(c) = ncread( nc_file, 'cp_flag' );
    
    
    % for the ICON 312 and 625 m simulations, remove the data off the hour
    % because those are NaNs (only hourly output for those two simulations)
    if any(strcmp(c,{'CSJogi312';'CSJogi625'}))
        itime_keep = time.(c).Minute==0;
        time.(c) = time.(c)(itime_keep);
        variogram.(c) = variogram.(c)(:,itime_keep);
        cp_flag.(c) = cp_flag.(c)(itime_keep);
    end
    
    % for IDEALhaboobSFC_night1, remove the last two entries, they are erroneous
    % can probably remove this code when Jennie re-processes the variograms
    if strcmp(c,'IDEALhaboobSFC_night1')
        time.(c) = time.(c)(1:end-2);
        variogram.(c) = variogram.(c)(:,1:end-2);
        cp_flag.(c) = cp_flag.(c)(1:end-2);
    end

        
    
    
    %% variogram processing for plots
    
    % variogram time series, averaged over (relevant) bins based on
    % specified bin_lim above
    ibin_lim.(c) = find( bin_mids.(c) <= bin_lim.(c), 1, 'last' );
    variogram_ts.(c) = nanmedian( variogram.(c)(1:ibin_lim.(c),:), 1 ); % median over dimension 1 (bins)
    
    % parse out cold pool periods. Everything should be relative to the
    % cold pool start time.
    
    cpstart = find(cp_flag.(c)==1);
    % most cases should just have one cold pool start time, but in case
    % there is more than one (e.g. OBS), loop over number of CP starts
    nCPs.(c) = length(cpstart); % save number of cold pools present in this run
    
    % for OBS, also make a pre-CP flag
    if strcmp(c,'OBS')
        precp_flag.(c) = zeros(size(cp_flag.(c)));
    end
    
    for n=1:nCPs.(c)
        cn = [c,num2str(n)]; % new character array for this cold pool
        
        % save mean variogram time series, average variogram over cold
        % pool, and time array relative to cold pool start
        time_cprel.(cn) = time.(c) - time.(c)(cpstart(n));
        % save time array from -180 min to + 240 min (-3h to +4h)
        it1 = find( time_cprel.(cn) <= minutes(-180), 1, 'last' );
        if isempty(it1); it1=1; end % check
        it2 = find( time_cprel.(cn) >= minutes(+240), 1, 'first' );
        if isempty(it2); it2=length(time_cprel.(cn)); end
        time_cprel.(cn) = time_cprel.(cn)(it1:it2);
        
        % save the preCP flag for OBS case
        if strcmp(c,'OBS')
            precp_flag.(c)(it1:cpstart(n)-1) = 1;
            % set to 0 if there are any cold pools during this time period
            precp_flag.(c)( cp_flag.(c)(it1:cpstart(n)-1)>=1 ) = 0;
        end
        
        % save the variogram relative time series for this cold pool
        variogram_ts_cprel.(cn) = variogram_ts.(c)(it1:it2);
        
        % also save cp flag relative time series
        cp_flag_cprel.(cn) = cp_flag.(c)(it1:it2);
        
        % save 2D variogram information over CP relative time period
        variogram_cprel.(cn) = variogram.(c)(:,it1:it2);
        
        % calculate mean/median variogram over the CP time period
        % first calculate index of cold pool start and end in the relative framework
        cp1 = find( time_cprel.(cn) == 0 ); % already defined to be time 0
        % use diff to find first instance where CP flag changes from 2 to 0 again
        % need to only look after cp1 time index, then add this in again at the end
        cp2 = find( diff(cp_flag_cprel.(cn)(cp1:end))==-2,1,'first')+cp1-1; 
        % check in case cold pool goes to the end of the record
        if isempty(cp2); cp2=length(time_cprel.(cn)); end
        
        % check that cold pool end time is not more than 4.5 h or less than 0
        if time_cprel.(cn)(cp2) > minutes(270)
            disp(['CP end time is more than 4.5 h after CP start for ',cn])
            disp('Setting CP end time to 4.5 h')
            cp2 = find(time_cprel.(cn) >= minutes(270),1,'first');
        elseif time_cprel.(cn)(cp2) < minutes(0)
            error(['CP end time is before CP start time for ',cn])
        end
        
        % for Jennie's ensemble of daytime haboobs, set CP duration to 1 h
        % NOTE: need Jennie to change this in the nc files
        if strcmp(c(1:min(length(c),11)),'IDEALhaboob')
            if strcmp(c(end-3:end-1),'day')
                cp2 = find(time_cprel.(cn) >= minutes(60),1,'first');
            end
        end
        
        
        % save variogram stats (mean, median, IQR)
        variogram_cpmean.(cn) = nanmean( variogram_cprel.(cn)(:,cp1:cp2),2 ) ; % mean over dimension 2 (time)
        % also do median, 25, and 75 percentiles
        variogram_cpmed.(cn) = nanmedian( variogram_cprel.(cn)(:,cp1:cp2),2 );
        variogram_cp_pctiles.(cn) = prctile( variogram_cprel.(cn)(:,cp1:cp2), [25 50 75], 2 );
        
        % also save pre-cold pool variograms
        preCP1 = find( time_cprel.(cn) >= minutes(-180),1,'first' );
        if isempty(preCP1); error('preCP1 empty'); end
        preCP2 = cp1-1; % index before CP1 time
        if isempty(preCP2); error('preCP2 empty'); end
        if preCP1>preCP2; error('preCP1 > preCP2'); end
        
        variogram_precpmean.(cn) = nanmean( variogram_cprel.(cn)(:,preCP1:preCP2),2 ); % mean over dimension 2 (time)
        variogram_precpmed.(cn) = nanmedian( variogram_cprel.(cn)(:,preCP1:preCP2),2 );
        variogram_precp_pctiles.(cn) = prctile( variogram_cprel.(cn)(:,preCP1:preCP2), [25 50 75], 2 );
        
        % calculate nighttime variograms for OBS-Jogi, ICON simulations,
        % and Haboob ensemble
        if strcmp(c,'OBSJogi') || strcmp(c(1:min(length(c),6)),'CSJogi')
            
            NightHrLT = [22 05]; NightHrZ = NightHrLT-2; % NOTE should match below
            iNight = time.(c).Hour>=NightHrZ(1) | time.(c).Hour<NightHrZ(2);
            variogram_night_pctiles.(cn) = prctile( variogram.(c)(:,iNight), [25 50 75], 2 );
            
        elseif strcmp(c(1:min(length(c),11)),'IDEALhaboob')
            if strcmp(c(end-5:end-1),'night')
                
                % for these cases, just use the 1-h pre-cold pool period defined above
                variogram_night_pctiles.(cn) = prctile( variogram_cprel.(cn)(:,preCP1:preCP2), [25 50 75], 2 );
            end
        end
        
        
    end % loop over # cold pools in this case
    
end % loop over variogram data files


%% composites and percentiles of observed cold pools

% for panel showing variogram for all cold pool and non-cold pool data
% (separate into day, night, and transition time. 11-18 for day, 22-05 for night)
DayHrLT = [11 18]; DayHrZ = DayHrLT-2;
NightHrLT = [22 05]; NightHrZ = NightHrLT-2; % NOTE should match above
% save cold pool times
iCP = cp_flag.OBS>=1;
CPtimesLT = time.OBS( iCP ) + hours(2); % times with cold pools; change to local time
% save no cold pool times
inoCP = cp_flag.OBS==0;
noCPtimesLT = time.OBS( inoCP ) + hours(2); % times without cold pools; change to local time

% index of no cold pool daytimes
inoCPday = ( inoCP & time.OBS.Hour>=DayHrZ(1) & time.OBS.Hour<DayHrZ(2) );
noCPtimesDayLT = time.OBS( inoCPday ) + hours(2); % day times without cold pools; change to local time

% index of no cold pool nighttimes
inoCPnight = ( inoCP & (time.OBS.Hour>=NightHrZ(1) | time.OBS.Hour<NightHrZ(2)) );
noCPtimesNightLT = time.OBS( inoCPnight ) + hours(2); % night times without cold pools; change to local time

% index of no cold pool transition times
inoCPtrans = ~iCP & ~inoCPday & ~inoCPnight;
noCPtimesTransLT = time.OBS( inoCPtrans ) + hours(2);

% calculate variogram percentiles for CP, non-CP, non-CP day only, non-CP night only
variogramOBS_CP_pctiles = prctile( variogram.OBS(:,iCP), [25 50 75], 2 );
variogramOBS_noCP_pctiles = prctile( variogram.OBS(:,inoCP), [25 50 75], 2 );
variogramOBS_noCPday_pctiles = prctile( variogram.OBS(:,inoCPday), [25 50 75], 2 );
variogramOBS_noCPnight_pctiles = prctile( variogram.OBS(:,inoCPnight), [25 50 75], 2 );
variogramOBS_noCPtrans_pctiles = prctile( variogram.OBS(:,inoCPtrans), [25 50 75], 2 );

% further separate CPs into daytime, nighttime, and transition
iCPday = ( iCP & time.OBS.Hour>=DayHrZ(1) & time.OBS.Hour<DayHrZ(2) );
CPdaytimesLT = time.OBS( iCPday ) + hours(2);
iCPnight = ( iCP & (time.OBS.Hour>=NightHrZ(1) | time.OBS.Hour<NightHrZ(2)) );
CPnighttimesLT = time.OBS( iCPnight ) + hours(2);
iCPtrans = iCP & ~iCPday & ~iCPnight;
CPtranstimesLT = time.OBS( iCPtrans ) + hours(2);

% calculate variogram percentiles for CP-day, CP-night, and CP-transition
variogramOBS_CPday_pctiles = prctile( variogram.OBS(:,iCPday), [25 50 75], 2 );
variogramOBS_CPnight_pctiles = prctile( variogram.OBS(:,iCPnight), [25 50 75], 2 );
variogramOBS_CPtrans_pctiles = prctile( variogram.OBS(:,iCPtrans), [25 50 75], 2 );

% calculate pre-CP variogram percentiles
ipreCP = precp_flag.OBS==1;
variogramOBS_preCP_pctiles = prctile( variogram.OBS(:,ipreCP), [25 50 75], 2 );



%% make figure 2 - compare observations overall and observations of Jogi case to ICON simulations


load('orangepurple'); load('redblu'); load('parula'); % load in colormaps needed for plotting
% set up colors for lines
colsorg = orangepurple( round(linspace(42,64,4)),: ); % oranges, for ICON runs
colspur = orangepurple( round(linspace(1,15,2)),: ); % purple, for RAMS DryBL runs
colsred = redblu( round(linspace(46,64,2)),: ); % red, for RAMS TropOce case study runs
colsgrn(1,:) = [0 .35 0]; % dark green, for RAMS Shear runs
colsgrn(2,:) = [0 .7 0];

clear lgdstr

% handle for fig2
ifig2 = figure; set_figsize(gcf,[15 12]);
% margins
om = [.07 .08 .03 .05]; im = [.06 .07];

% panel (a) - observations - median, with 25 and 75 percentiles shaded, for CP and Non-CP
multipan(ifig2,2,2,1,'om',om,'im',im);
hold all
l=0;
for iln=[3 1]; l=l+1; % plot Dyatime non-cold pool vs cold pool (indices 3 and 1). The rest below are kept for checks / flexibility
    if iln==1
        var2plot = variogramOBS_CP_pctiles;
        col = parula(10,:); % blue
        %col = 'k'; % black
        lgdstr{l} = 'OBS-All Cold Pool';
    elseif iln==2
        var2plot = variogramOBS_noCP_pctiles; 
        col = 'k'; % black
        lgdstr{l} = 'Non-Cold Pool';
    elseif iln==3
        var2plot = variogramOBS_noCPday_pctiles; 
        col = redblu(58,:); % red
        lgdstr{l} = 'OBS-All Daytime Non-Cold Pool';
    elseif iln==4
        var2plot = variogramOBS_noCPnight_pctiles; 
        %col = orangepurple(1,:); % purple
        col = 'k'; % black
        lgdstr{l} = 'Non-Cold Pool, Night';
    elseif iln==5
        var2plot = variogramOBS_noCPtrans_pctiles;
        col = colsgrn(1,:); % dark green
        lgdstr{l} = 'Non-Cold Pool, Transition';
    elseif iln==6
        var2plot = variogramOBS_CPday_pctiles;
        col = orangepurple(64,:); % orange
        lgdstr{l} = 'Cold Pool, Day';
    elseif iln==7
        var2plot = variogramOBS_CPnight_pctiles;
        col = parula(10,:); % blue
        lgdstr{l} = 'Cold Pool, Night';
    elseif iln==8
        var2plot = variogramOBS_CPtrans_pctiles;
        col = colsgrn(1,:); % dark green
        lgdstr{l} = 'Cold Pool, Transition';
    elseif iln==9 % pre-CP
        var2plot = variogramOBS_preCP_pctiles;
        col = redblu(58,:); % red
        lgdstr{l} = 'Pre-Cold Pool';
    end
    pct50 = var2plot(:,2); pct25 = var2plot(:,1); pct75 = var2plot(:,3);
    bds = [ pct50-pct25 pct75-pct50 ];
    boundedline( bin_mids.OBS, pct50, bds, '-', ...
        'alpha', 'transparency', .2, 'orientation', 'vert', 'nan', 'remove', 'cmap', col );
end
% bring all lines to the foreground
ln=findall(gca,'type','line'); uistack(ln,'top'); 
lgd = legend(lgdstr,'location','northwest');
lgd.Box = 'off';

box on
%xlabel(['Distance (',bin_units.OBS,')'])
xticklabels([]);
ylabel(['Temperature Variogram (',variogram_units.OBS,')'])
ylim([0 2]);
title('(a)  All FESSTVaL OBS Variograms')


% panel (b) - median variograms over cold pool time period for OBS-All, Jogi,
% and ICON simulation suite
clear lgdstr
multipan(ifig2,2,2,2,'om',om,'im',im);
hold all
% first overplot shading for cold pool obs
var2plot = variogramOBS_CP_pctiles;
lgdstr{1} = 'OBS-All Cold Pool';
pct50 = var2plot(:,2); pct25 = var2plot(:,1); pct75 = var2plot(:,3);
bds = [ pct50-pct25 pct75-pct50 ];
[ln,ptch] = boundedline( bin_mids.OBS, pct50, bds, '-', ...
    'alpha', 'transparency', .2, 'orientation', 'vert', 'nan', 'remove', 'cmap', parula(10,:) );
% remove the median line and add a legend entry to the shading
delete(ln); ptch.Annotation.LegendInformation.IconDisplayStyle='on';
l=1; % set to 1 since I already have a legend entry
for p=1:5; l=l+1;
    if p==1; c = 'OBSJogi'; col = 'k';
    elseif p==2; c = 'CSJogi75'; col = colsorg(4,:); % darkest
    elseif p==3; c = 'CSJogi156'; col = colsorg(3,:); 
    elseif p==4; c = 'CSJogi312'; col = colsorg(2,:); 
    elseif p==5; c = 'CSJogi625'; col = colsorg(1,:); 
    end
    lstyle='-'; if 8==5; lstyle=':'; end % dashed line to separate upshear case from downshear

    cn=[c,'1']; % cold pool #
    if p==1 % include shading for OBSJogi
        var2plot = variogram_cp_pctiles.(cn);
        pct50 = var2plot(:,2); pct25 = var2plot(:,1); pct75 = var2plot(:,3);
        bds = [ pct50-pct25 pct75-pct50 ];
        [ln,ptch] = boundedline( bin_mids.(c)(1:ibin_lim.(c)), pct50(1:ibin_lim.(c)), bds(1:ibin_lim.(c),:), ...
            lstyle, 'alpha', 'transparency', .2, 'orientation', 'vert', 'nan', 'remove', 'cmap', col );
    else
        plot( bin_mids.(c)(1:ibin_lim.(c)), variogram_cpmed.(cn)(1:ibin_lim.(c)), lstyle, 'color', col )
    end
    lgdstr{l} = ttl.(c);
end

lgd = legend( lgdstr, 'location', 'northwest' );
lgd.Box='off';

%xlabel(['Distance (',bin_units.(c),')'])
%ylabel(['Temperature Variogram (',variogram_units.(c),')'])
xticklabels([]);
title('(b)  Jogi Cold Pool Variograms')
ylim([0 8.5])
box on


% panel (c) - pre-cold pool comparisons between Jogi OBS, and pre-Jogi ICON runs
clear lgdstr
multipan(ifig2,2,2,3,'om',om,'im',im);
hold all
l=0;
for p=1:5; l=l+1;
    if p==1; c = 'OBSJogi'; col = 'k';
    elseif p==2; c = 'CSJogi75'; col = colsorg(4,:); % darkest
    elseif p==3; c = 'CSJogi156'; col = colsorg(3,:);
    elseif p==4; c = 'CSJogi312'; col = colsorg(2,:);
    elseif p==5; c = 'CSJogi625'; col = colsorg(1,:);
    end
    lstyle='-'; if p==8; lstyle=':'; end % dashed line to separate upshear case from downshear

    cn=[c,'1']; % cold pool #
    if p==1  % include shading for OBSJogi
        var2plot = variogram_precp_pctiles.(cn);
        pct50 = var2plot(:,2); pct25 = var2plot(:,1); pct75 = var2plot(:,3);
        bds = [ pct50-pct25 pct75-pct50 ];
        [ln,ptch] = boundedline( bin_mids.(c)(1:ibin_lim.(c)), pct50(1:ibin_lim.(c)), bds(1:ibin_lim.(c),:), ...
            lstyle, 'alpha', 'transparency', .2, 'orientation', 'vert', 'nan', 'remove', 'cmap', col );
    else
        plot( bin_mids.(c)(1:ibin_lim.(c)), variogram_precpmed.(cn)(1:ibin_lim.(c)), lstyle, 'color', col )
    end
    lgdstr{l} = ttl.(c);
end

lgd = legend( lgdstr, 'location', 'northwest' );
lgd.Box='off';

xlabel(['Distance (',bin_units.(c),')'])
ylabel(['Temperature Variogram (',variogram_units.(c),')'])
title('(c)  Jogi Pre-Cold Pool Variograms')
ylim([0 2]);
box on


% panel (d) - nighttime comparisons, include nighttime JOGI-OBS and ICON runs
clear lgdstr
multipan(ifig2,2,2,4,'om',om,'im',im);
hold all
l=0;
for p=1:5; l=l+1;
    if p==1; c = 'OBSJogi'; col = 'k';
    elseif p==2; c = 'CSJogi75'; col = colsorg(4,:); % darkest
    elseif p==3; c = 'CSJogi156'; col = colsorg(3,:);
    elseif p==4; c = 'CSJogi312'; col = colsorg(2,:);
    elseif p==5; c = 'CSJogi625'; col = colsorg(1,:);
    end

    cn=[c,'1']; % cold pool #
    if p==1  % include shading for OBSJogi
        var2plot = variogram_night_pctiles.(cn);
        pct50 = var2plot(:,2); pct25 = var2plot(:,1); pct75 = var2plot(:,3);
        bds = [ pct50-pct25 pct75-pct50 ];
        [ln,ptch] = boundedline( bin_mids.(c)(1:ibin_lim.(c)), pct50(1:ibin_lim.(c)), bds(1:ibin_lim.(c),:), ...
            lstyle, 'alpha', 'transparency', .2, 'orientation', 'vert', 'nan', 'remove', 'cmap', col );
    else
        plot( bin_mids.(c)(1:ibin_lim.(c)), variogram_night_pctiles.(cn)(1:ibin_lim.(c),2), '-', 'color', col ) % 50th percentile
    end
    lgdstr{l} = ttl.(c);
end

lgd = legend( lgdstr, 'location', 'northwest' );
lgd.Box='off';

xlabel(['Distance (',bin_units.(c),')'])
%ylabel(['Temperature Variogram (',variogram_units.(c),')'])
title('(d)  Jogi Nighttime Variograms')
ylim([0 2]);
box on

print('-dpng','Fig2_OBS-All_OBS-Jogi_CS-Jogi')
saveas(gcf,'Fig2_OBS-All_OBS-Jogi_CS-Jogi','fig')


%% figure 3 - compare variograms among observations, resolution tests, different environments

ifig3 = figure; set_figsize(gcf,[15 12]);

% specify margins just for the top two panels
% use multipanel to make it look like a 1-by-2 panel for these plots, since
% panel (c) needs to fill the space differently
om = [.07 .08 .03 .08]; im = [.06 .09];
omoffset = .43;


% panel (a) - median variograms over cold pool time period among resolution tests
clear lgdstr
multipan( ifig3, 1,2,1,'om', om+[0 omoffset 0 0], 'im', im );
hold all
l=0;
for p=2:13; l=l+1;
    if p==1; c = 'OBSJogi'; col = 'k';
    elseif p==2; c = 'CSJogi75'; col = colsorg(4,:); % darkest
    elseif p==3; c = 'CSJogi156'; col = colsorg(3,:);
    elseif p==4; c = 'CSJogi312'; col = colsorg(2,:);
    elseif p==5; c = 'CSJogi625'; col = colsorg(1,:);
    elseif p==6; c = 'IDEALdryBL50'; col = colspur(1,:); % darkest
    elseif p==7; c = 'IDEALdryBL100'; col = colspur(2,:);
    elseif p==8; c = 'IDEALdownshear100'; col = colsgrn(1,:); % darkest
    elseif p==9; c = 'IDEALdownshear250'; col = colsgrn(2,:);
    elseif p==10; c = 'IDEALupshear100'; col = colsgrn(1,:); % darkest
    elseif p==11; c = 'IDEALupshear250'; col = colsgrn(2,:);
    elseif p==12; c = 'CStropOce100'; col = colsred(2,:); % darkest
    elseif p==13; c = 'CStropOce1000'; col = colsred(1,:);
    end
    lstyle='-'; if p==10||p==11; lstyle=':'; end % dashed line to separate upshear case from downshear

    cn=[c,'1']; % cold pool #
    ln = plot( bin_mids.(c)(1:ibin_lim.(c)), variogram_cpmed.(cn)(1:ibin_lim.(c)), lstyle, 'color', col );
    % only include legend information for one of each simulation...
    % otherwise legend is too large!
    if any(p==[2 6 8 10 12])
        lgdstr{l} = ttl.(c);
        lgdstr{l} = lgdstr{l}(1:find(lgdstr{l}=='-',1,'last')-1); % remove the grid spacing information
    else
        l=l-1;
        ln.Annotation.LegendInformation.IconDisplayStyle='off';
    end
end

xlabel(['Distance (',bin_units.(c),')'])
ylabel(['Temperature Variogram (',variogram_units.(c),')'])
title({'(a)  Cold Pool Variograms:';'Environment and Resolution Sensitivity'})
ylim([0 6])

lgd = legend( lgdstr, 'location', 'northwest' ); lgd.Box='off';

box on


% panel (b) - Jennie's ensemble
multipan( ifig3, 1,2,2,'om', om+[0 omoffset 0 0], 'im', im );

% ifig=1 here (sensitivity to initial CP temp, day and night) is panel b
% additional code kept for checks / flexibility
for ifig=1:1;%3; % 3 sets of simulations, each with day and night
    %figure; set_figsize(gcf,[10 7]); 
    hold all
    clear lgdstr
    colsday = parula(round(linspace(46,57,4)),:); colsday=flipdim(colsday,1);
    colsnight = parula(round(linspace(1,22,4)),:);
    cols = [colsday; colsnight];
    if ifig==1
        istart = ensref;
        ttlstr = {'(b)  IDEAL-Haboob-150m Ensemble:';'Sensitivity to Initial Temp. and Time of Day'};
    elseif ifig==2
        istart = ensref+8;
        ttlstr = 'IDEAL-Haboob Ensemble: Sensitivity to Land Sfc Type';
    elseif ifig==3
        istart = ensref + 16;
        ttlstr = 'IDEAL-Haboob Ensemble: Sensitivity to Soil Moisture';
    end
    for v=istart:istart+7 % loop through simulations
        c=name{v}; cn=[c,'1'];
        plot( bin_mids.(c), variogram_cpmed.(cn),'-','color',cols(v-istart+1,:) )
        lgdstr{v-istart+1} = ttl.(c)(14:end-5);
    end
    lgd = legend( lgdstr, 'location','northwest' ); lgd.Box='off';
    box on
    xlabel(['Distance (',bin_units.(c),')'])
    %ylabel('Temperature Variogram (K^2)')
    ylim([0 6])
    title(ttlstr)
end




% panel (c) - time series
clear lgdstr
% this panel covers the entire lower part of the figure
% modify margins a little
multipan( ifig3, 1,1,1,'om', om+[.02 0 -.01 1-omoffset-.03] ); 
hold all
l=0;
for p=1:7; l=l+1;
    if p==1; c = 'OBSJogi'; col = 'k';
    elseif p==2; c = 'CSJogi156'; col = colsorg(4,:); % darkest
    elseif p==3; c = 'IDEALdryBL50'; col = colspur(1,:); % darkest
    elseif p==4; c = 'IDEALdownshear100'; col = colsgrn(1,:); % darkest
    elseif p==5; c = 'IDEALupshear100'; col = colsgrn(1,:); % darkest
    elseif p==6; c = 'CStropOce100'; col = colsred(2,:); % darkest
    elseif p==7; c = 'IDEALhaboobCP_day1'; col = colsday(1,:); % same yellow as in panel (b)
    end
    lstyle='-'; if p==5; lstyle=':'; end % dashed line to separate upshear case from downshear

    cn=[c,'1']; % cold pool #
    plot( time_cprel.(cn), variogram_ts_cprel.(cn), lstyle, 'color', col )
    lgdstr{l} = ttl.(c);
end

xlim(duration([minutes(-30) minutes(120)]))
xtickformat('hh:mm')

set(gca,'yscale','log'); ylim([5e-3 50])

xlabel(['Time (hh:mm since CP onset)'])
ylabel(['Temperature Variogram (',variogram_units.(c),')'])
title('(c)  Variogram Time Series')

lgd = legend( lgdstr, 'location', 'eastoutside' );

box on

print('-dpng','Fig3_Environment_Resolution_Sensitivity')
saveas(gcf,'Fig3_Environment_Resolution_Sensitivity','fig')


