
%% SETUP
% clear, make some variables
close all
warning('off')
processed = [];
myplots = [];

% read in the Herzberg data
herzbergData = readtable([pwd, '\supplemental\Herzberg_etal_2010_Supplemental.xls']);
% age data
herzbergAgeIndex = contains(table2cell(herzbergData(:, "Var1")), 'age (Ma)');
herzbergAgeData = 4.543 - (table2array(herzbergData(herzbergAgeIndex, 3:end)) / 1e3);
% temperature data
herzbergMgOIndex = contains(table2cell(herzbergData(:, "Var1")), 'MgO');
herzbergMgOData = table2array(herzbergData(herzbergMgOIndex, 3:end));
herzbergTempData = 1463 + 12.74.*herzbergMgOData - 2924./herzbergMgOData;

% read in the Condie data
condieData = readtable([pwd, '\supplemental\Condie_etal_2016_Supplemental.xlsx'], 'Sheet', "Komatiite");
% age data
condieAgeData = 4.543 - (table2array(condieData(:, 'Age_Ma_')) / 1e3);
% temperature data
condieMgOData = table2array(condieData(:, 'MgO'));
condieTempData = 1463 + 12.74.*condieMgOData - 2924./condieMgOData;
% average komatiites by age into bins so we can do statistics
condieUniqueAges = unique(condieAgeData(condieAgeData < 4));
condieTempAvgs = NaN(size(condieUniqueAges));
condieTempMinDiffs = NaN(size(condieUniqueAges));
condieTempMaxDiffs = NaN(size(condieUniqueAges));
for i = 1:length(condieUniqueAges)
    temps = condieTempData(condieAgeData == condieUniqueAges(i));
    condieTempAvgs(i) = mean(temps,'omitnan');
    condieTempMinDiffs(i) = condieTempAvgs(i) - min(temps, [], 'omitnan');
    condieTempMaxDiffs(i) = max(temps, [], 'omitnan') - condieTempAvgs(i);
end

% add them to the geochem structure
geochem.herzbergAgeData = herzbergAgeData;
geochem.herzbergTempData = herzbergTempData;
geochem.condieUniqueAges = condieUniqueAges;
geochem.condieTempAvgs = condieTempAvgs;
geochem.condieTempMinDiffs = condieTempMinDiffs;
geochem.condieTempMaxDiffs = condieTempMaxDiffs;

% clear everything that's not processed, myplots, geochem
clearvars -except processed myplots geochem t0conc t0fug

% present-day uncertainties
errors.T = [1350, 80, 95];
errors.Xum = [100, 50, 100];
errors.OM = [1, 0.1, 0.1];
errors.mantleOM = [4.625, 4.375, 4.375];
errors.eta_um = [1e21, 9e20, 9e21];
errors.u = [4.5, 0.5, 0.5];
errors.Qs = [32, 2, 2];
errors.RminusD_Korenaga = [3.5e11, 0.8e11, 1.5e11];
errors.RminusD = [4.5e11, 1.8e11, 6.7e11];
errors.UreyRatio = [0.29, 0.17, 0.20];
errors.fH2Oum = [31.825, 31.475, 31.475];

%% LOAD IN THE FWD BASE TIMESERIES DATA
% file paths
SectionBreak()
disp('NOW LOADING: FWD BASE TS')
processed.settings.baseDirInfo = dir([pwd, '\models\fwd\base\', '*TS.mat']);

% loading those cases
for j = 1:length(processed.settings.baseDirInfo)
    % display that we're loading things
    disp(['--- ', processed.settings.baseDirInfo(j).name])
    if contains(processed.settings.baseDirInfo(j).name, 'FUG')
        processed.fwd.fug.base = load([pwd, '\models\fwd\base\', processed.settings.baseDirInfo(j).name]);
    else
        processed.fwd.conc.base = load([pwd, '\models\fwd\base\', processed.settings.baseDirInfo(j).name]);
    end
end
disp(['--- COMPLETE CONC BASELINE REALIZATIONS: ', num2str(processed.fwd.conc.base.ECs.solved.nFullSims)])
disp(['--- COMPLETE FUG BASELINE REALIZATIONS: ', num2str(processed.fwd.fug.base.ECs.solved.nFullSims)])

myplots.TS.avgFig = figure;
myplots.TS.labels = ('a':'z');
myplots.TS.params = {'T', 'Xum', 'u', 'Qs', 'oceanmass', 'mantleH2Omass', 'eta_um', 'UreyRatio', 'RminusD', 'fH2Oum'};
myplots.TS.timeInterval = linspace(processed.fwd.fug.base.ECs.range.t(1) .* processed.fwd.fug.base.ECs.conversion.sec2yr ./ 1e9, ...
                                   processed.fwd.fug.base.ECs.range.t(end) .* processed.fwd.fug.base.ECs.conversion.sec2yr ./ 1e9, ...
                                   500);
myplots.steadystate.timeInterval = myplots.TS.timeInterval;
for i = 1:length(myplots.TS.params)
    % interpolating the good CONC data
    myplots.TS.interpConc = NaN(length(processed.fwd.conc.base.ECs.solved.TS), length(myplots.TS.timeInterval));
    for j = 1:length(processed.fwd.conc.base.ECs.solved.TS)
        myplots.TS.interpConc(j, :) = interp1(processed.fwd.conc.base.ECs.solved.TS(j).t, ...
                                              processed.fwd.conc.base.ECs.solved.TS(j).(myplots.TS.params{i}), ...
                                              myplots.TS.timeInterval);
    end
    % interpolating the good FUG data
    myplots.TS.interpFug = NaN(length(processed.fwd.fug.base.ECs.solved.TS), length(myplots.TS.timeInterval));
    for j = 1:length(processed.fwd.fug.base.ECs.solved.TS)
        myplots.TS.interpFug(j, :) = interp1(processed.fwd.fug.base.ECs.solved.TS(j).t, ...
                                             processed.fwd.fug.base.ECs.solved.TS(j).(myplots.TS.params{i}), ...
                                             myplots.TS.timeInterval);
    end

    % median, 25/75 percentiles, min/max
    myplots.TS.avgConc = median(myplots.TS.interpConc, 'omitnan');
    myplots.TS.prctile25Conc = prctile(myplots.TS.interpConc, 25);
    myplots.TS.prctile75Conc = prctile(myplots.TS.interpConc, 75);
    myplots.TS.minConc = min(myplots.TS.interpConc, [], 'omitnan');
    myplots.TS.maxConc = max(myplots.TS.interpConc, [], 'omitnan');
    myplots.TS.avgFug = median(myplots.TS.interpFug, 'omitnan');
    myplots.TS.prctile25Fug = prctile(myplots.TS.interpFug, 25);
    myplots.TS.prctile75Fug = prctile(myplots.TS.interpFug, 75);
    myplots.TS.minFug = min(myplots.TS.interpFug, [], 'omitnan');
    myplots.TS.maxFug = max(myplots.TS.interpFug, [], 'omitnan');

    % setting up uncertainty bands
    myplots.TS.timeBand = [myplots.TS.timeInterval, myplots.TS.timeInterval(end:-1:1)];
    myplots.TS.concBand = [myplots.TS.prctile75Conc, myplots.TS.prctile25Conc(end:-1:1)];
    myplots.TS.concExtremes = [myplots.TS.maxConc, myplots.TS.minConc(end:-1:1)];
    myplots.TS.fugBand = [myplots.TS.prctile75Fug, myplots.TS.prctile25Fug(end:-1:1)];
    myplots.TS.fugExtremes = [myplots.TS.maxFug, myplots.TS.minFug(end:-1:1)];

    % conditional by parameter
    if strcmp('T', myplots.TS.params{i})
        subplot(3, 4, [1,2]);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel(['T [',char(176),'C]']); ylim([800, 1850]);
        errorbar(processed.fwd.fug.base.ECs.range.t(end) .* processed.fwd.fug.base.ECs.conversion.sec2yr ./ 1e9, errors.T(1), errors.T(2), errors.T(3), 'k-', 'linewidth', 1)
        myplots.magni.conc.avgT = myplots.TS.avgConc;
        myplots.magni.conc.bandT = myplots.TS.concBand;
        myplots.magni.conc.extremeT = myplots.TS.concExtremes;
        myplots.magni.fug.avgT = myplots.TS.avgFug;
        myplots.magni.fug.bandT = myplots.TS.fugBand;
        myplots.magni.fug.extremeT = myplots.TS.fugExtremes;
    elseif strcmp('Xum', myplots.TS.params{i})
        subplot(3, 4, [3,4]);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('\chi_{um} [ppm]'); set(gca, 'YScale', 'log'); ylim([0.1 5000]); yticks([.1 10 1000]);
        errorbar(processed.fwd.fug.base.ECs.range.t(end) .* processed.fwd.fug.base.ECs.conversion.sec2yr ./ 1e9, errors.Xum(1), errors.Xum(2), errors.Xum(3), 'k-', 'linewidth', 1)
        myplots.TS.concBand(myplots.TS.concBand < 0) = 1e-10;
        myplots.TS.fugBand(myplots.TS.fugBand < 0) = 1e-10;
    elseif strcmp('u', myplots.TS.params{i})
        subplot(3, 4, 5);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('u [cm yr^{-1}]'); set(gca, 'YScale', 'log'); ylim([0.1 150]); yticks([1 10 100]);
        myplots.TS.concBand(myplots.TS.concBand < eps) = eps;
        myplots.TS.fugBand(myplots.TS.fugBand < eps) = eps;
        errorbar(processed.fwd.fug.base.ECs.range.t(end) .* processed.fwd.fug.base.ECs.conversion.sec2yr ./ 1e9, errors.u(1), errors.u(2), errors.u(3), 'k-', 'linewidth', 1)
        myplots.magni.conc.avgu = myplots.TS.avgConc;
        myplots.magni.conc.bandu = myplots.TS.concBand;
        myplots.magni.conc.extremeu = myplots.TS.concExtremes;
        myplots.magni.fug.avgu = myplots.TS.avgFug;
        myplots.magni.fug.bandu = myplots.TS.fugBand;
        myplots.magni.fug.extremeu = myplots.TS.fugExtremes;
    elseif strcmp('Qs', myplots.TS.params{i})
        subplot(3, 4, 6);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('Q_s [TW]'); set(gca, 'YScale', 'log'); ylim([5 150]); yticks([10 100]);
        errorbar(processed.fwd.fug.base.ECs.range.t(end) .* processed.fwd.fug.base.ECs.conversion.sec2yr ./ 1e9, errors.Qs(1), errors.Qs(2), errors.Qs(3), 'k-', 'linewidth', 1)
        myplots.TS.concBand(myplots.TS.concBand < 0) = 1e-10;
        myplots.TS.fugBand(myplots.TS.fugBand < 0) = 1e-10;
    elseif strcmp('oceanmass', myplots.TS.params{i})
        subplot(3, 4, 7);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('Surface Water [OM]'); ylim([0 11]);
        errorbar(processed.fwd.fug.base.ECs.range.t(end) .* processed.fwd.fug.base.ECs.conversion.sec2yr ./ 1e9, errors.OM(1), errors.OM(2), errors.OM(3), 'k-', 'linewidth', 1)
    elseif strcmp('mantleH2Omass', myplots.TS.params{i})
        subplot(3, 4, 8);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('Mantle Water [OM]'); ylim([0 11]);
        errorbar(processed.fwd.fug.base.ECs.range.t(end) .* processed.fwd.fug.base.ECs.conversion.sec2yr ./ 1e9, errors.mantleOM(1), errors.mantleOM(2), errors.mantleOM(3), 'k-', 'linewidth', 1)
    elseif strcmp('eta_um', myplots.TS.params{i})
        subplot(3, 4, 9);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('\eta_{um} [Pa s]'); set(gca, 'YScale', 'log'); ylim([1e19 1e24]); yticks([1e19 1e21 1e23])
        errorbar(processed.fwd.fug.base.ECs.range.t(end) .* processed.fwd.fug.base.ECs.conversion.sec2yr ./ 1e9, errors.eta_um(1), errors.eta_um(2), errors.eta_um(3), 'k-', 'linewidth', 1)
    elseif strcmp('UreyRatio', myplots.TS.params{i})
        subplot(3, 4, 10);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('Urey Ratio [-]'); ylim([0 1.75])
        errorbar(processed.fwd.fug.base.ECs.range.t(end) .* processed.fwd.fug.base.ECs.conversion.sec2yr ./ 1e9, errors.UreyRatio(1), errors.UreyRatio(2), errors.UreyRatio(3), 'k-', 'linewidth', 1)
    elseif strcmp('RminusD', myplots.TS.params{i})
        subplot(3, 4, 11);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('R - D [kg yr^{-1}]'); ylim([-1e12 4e12])
        errorbar(processed.fwd.fug.base.ECs.range.t(end) .* processed.fwd.fug.base.ECs.conversion.sec2yr ./ 1e9, errors.RminusD_Korenaga(1), errors.RminusD_Korenaga(2), errors.RminusD_Korenaga(3), 'k-', 'linewidth', 1)
    elseif strcmp('fH2Oum', myplots.TS.params{i})
        subplot(3, 4, 12);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('f_{H_2O, um} [%]'); ylim([0 65]);
        errorbar(processed.fwd.fug.base.ECs.range.t(end) .* processed.fwd.fug.base.ECs.conversion.sec2yr ./ 1e9, errors.fH2Oum(1), errors.fH2Oum(2), errors.fH2Oum(3), 'k-', 'linewidth', 1)
    end

    % plotting
    xlabel('Time [Gyr]'); xlim([myplots.TS.timeInterval(1) myplots.TS.timeInterval(end)]);
    title(['(',myplots.TS.labels(i),')'])
    if strcmp('T', myplots.TS.params{i})
        scatter(geochem.herzbergAgeData, geochem.herzbergTempData, 40, [46/255, 201/255, 123/255], 'filled')
        errorbar(geochem.condieUniqueAges, geochem.condieTempAvgs, geochem.condieTempMinDiffs, geochem.condieTempMaxDiffs, ...
                 "Marker", 'square', "MarkerSize", 6.5, ...
                 "MarkerEdgeColor", [219/255 186/255 39/255], ...
                 "MarkerFaceColor", [219/255 186/255 39/255], ...
                 'CapSize', 0, "Color", [219/255 186/255 39/255], "LineStyle", "none")
    end
    fill(myplots.TS.timeBand, myplots.TS.concExtremes, [.85 .325 .098], 'FaceAlpha', 0.15, 'EdgeColor', 'none')
    fill(myplots.TS.timeBand, myplots.TS.fugExtremes, [0 .447 .741], 'FaceAlpha', 0.15, 'EdgeColor', 'none')
    fill(myplots.TS.timeBand, myplots.TS.concBand, [.85 .325 .098], 'FaceAlpha', 0.35, 'EdgeColor', 'none')
    fill(myplots.TS.timeBand, myplots.TS.fugBand, [0 .447 .741], 'FaceAlpha', 0.35, 'EdgeColor', 'none')
    plot(myplots.TS.timeInterval, myplots.TS.avgConc, 'r-')
    plot(myplots.TS.timeInterval, myplots.TS.avgFug, 'b-')
end
% show plot
figure(myplots.TS.avgFig)

SectionBreak()
disp('NOW APPLYING: MAGNI FEEDBACK')
myplots.magni.viscDep = {'conc', 'fug'};
myplots.magni.params = {'u', 'delta', 'T'};
myplots.magni.labels = {'Varying [u]', 'Varying [a]', 'Varying [T]', 'Varying [const]'};
myplots.magni.timeInterval = myplots.TS.timeInterval;
myplots.magni.W = [1.06 0.14 -0.023 17 2.32];
myplots.magni.varyParam = [0.5 1 2];
myplots.magni.linestyles = {'-.', '-', '--'};

myplots.magni.sec2year = processed.fwd.fug.base.ECs.conversion.sec2yr;
myplots.magni.kappa = processed.fwd.fug.base.ECs.avg.kappa;
myplots.magni.L = processed.fwd.fug.base.ECs.avg.L;
myplots.magni.conc.Xps = NaN(length(processed.fwd.conc.base.ECs.solved.TS), 1);
myplots.magni.fug.Xps = NaN(length(processed.fwd.fug.base.ECs.solved.TS), 1);
for j = 1:length(processed.fwd.fug.base.ECs.solved.TS)
    myplots.magni.conc.Xps(j) = processed.fwd.conc.base.ECs.solved.TS(j).Xp(end);
    myplots.magni.fug.Xps(j) = processed.fwd.fug.base.ECs.solved.TS(j).Xp(end);
end
myplots.magni.conc.Xp = median(myplots.magni.conc.Xps, 'omitnan') * processed.fwd.conc.base.ECs.conversion.ppm2massfrac;
myplots.magni.fug.Xp = median(myplots.magni.fug.Xps, 'omitnan') * processed.fwd.fug.base.ECs.conversion.ppm2massfrac;
myplots.magni.rhop = 3315;
myplots.magni.A_O = processed.fwd.fug.base.ECs.avg.A_O;
myplots.magni.trench = 47650e3;
myplots.magni.fRscale = 1;
processed = [];

% plotting just the averages from the forward time direction as Magni parameters
myplots.magni.TSs = figure;
tiledlayout('flow','Padding','compact','TileSpacing','compact');

    % plot stuff 
    nexttile;
    hold on; grid on;
    fontsize(gca, scale=1.33);
    xlabel('Time [Gyr]'); xlim([myplots.magni.timeInterval(1) myplots.magni.timeInterval(end)]);
    ylabel('f_R(t) [%]'); set(gca, 'YScale', 'log'); ylim([1e-3 1e2]);
    title(['(',myplots.TS.labels(1),')'])

    % calculating average f_R(t)
    for i = 1:length(myplots.magni.viscDep)
        % calculate W
        myplots.magni.(myplots.magni.viscDep{i}).avgdelta = 2 .* sqrt((myplots.magni.kappa .* myplots.magni.L) ./ (myplots.magni.(myplots.magni.viscDep{i}).avgu .* myplots.magni.sec2year ./ 100));
        myplots.magni.(myplots.magni.viscDep{i}).avga = ((myplots.magni.(myplots.magni.viscDep{i}).avgdelta ./ myplots.magni.W(5)).^2) ./ myplots.magni.kappa .* myplots.magni.sec2year .* 1e-6;
        myplots.magni.(myplots.magni.viscDep{i}).avgW = 1e5 .* ((myplots.magni.W(1) .* myplots.magni.(myplots.magni.viscDep{i}).avgu) + ...
                                                                (myplots.magni.W(2) .* myplots.magni.(myplots.magni.viscDep{i}).avga) + ...
                                                                (myplots.magni.W(3) .* myplots.magni.(myplots.magni.viscDep{i}).avgT) + ...
                                                                myplots.magni.W(4));
        myplots.magni.(myplots.magni.viscDep{i}).banddelta = 2 .* sqrt((myplots.magni.kappa .* myplots.magni.L) ./ (myplots.magni.(myplots.magni.viscDep{i}).bandu .* myplots.magni.sec2year ./ 100));
        myplots.magni.(myplots.magni.viscDep{i}).banda = ((myplots.magni.(myplots.magni.viscDep{i}).banddelta ./ myplots.magni.W(5)).^2) ./ myplots.magni.kappa .* myplots.magni.sec2year .* 1e-6;
        myplots.magni.(myplots.magni.viscDep{i}).bandW = 1e5 .* ((myplots.magni.W(1) .* myplots.magni.(myplots.magni.viscDep{i}).bandu) + ...
                                                                 (myplots.magni.W(2) .* myplots.magni.(myplots.magni.viscDep{i}).banda) + ...
                                                                 (myplots.magni.W(3) .* myplots.magni.(myplots.magni.viscDep{i}).bandT) + ...
                                                                 myplots.magni.W(4));
        myplots.magni.(myplots.magni.viscDep{i}).extremedelta = 2 .* sqrt((myplots.magni.kappa .* myplots.magni.L) ./ (myplots.magni.(myplots.magni.viscDep{i}).extremeu .* myplots.magni.sec2year ./ 100));
        myplots.magni.(myplots.magni.viscDep{i}).extremea = ((myplots.magni.(myplots.magni.viscDep{i}).extremedelta ./ myplots.magni.W(5)).^2) ./ myplots.magni.kappa .* myplots.magni.sec2year .* 1e-6;
        myplots.magni.(myplots.magni.viscDep{i}).extremeW = 1e5 .* ((myplots.magni.W(1) .* myplots.magni.(myplots.magni.viscDep{i}).extremeu) + ...
                                                                 (myplots.magni.W(2) .* myplots.magni.(myplots.magni.viscDep{i}).extremea) + ...
                                                                 (myplots.magni.W(3) .* myplots.magni.(myplots.magni.viscDep{i}).extremeT) + ...
                                                                 myplots.magni.W(4));


        % calculate fR(t), bound from 0 to 1
        myplots.magni.(myplots.magni.viscDep{i}).avgfR_t = myplots.magni.fRscale .* (myplots.magni.(myplots.magni.viscDep{i}).avgW ./ (myplots.magni.(myplots.magni.viscDep{i}).Xp.*myplots.magni.rhop.*(myplots.magni.A_O./myplots.magni.trench)));
        myplots.magni.(myplots.magni.viscDep{i}).tempfR = myplots.magni.(myplots.magni.viscDep{i}).avgfR_t;
        myplots.magni.(myplots.magni.viscDep{i}).tempfR(myplots.magni.(myplots.magni.viscDep{i}).tempfR < eps) = eps; myplots.magni.(myplots.magni.viscDep{i}).tempfR(myplots.magni.(myplots.magni.viscDep{i}).tempfR >= 1) = 1-eps;
        myplots.magni.(myplots.magni.viscDep{i}).avgfR_t = myplots.magni.(myplots.magni.viscDep{i}).tempfR;
        myplots.magni.(myplots.magni.viscDep{i}).bandfR_t = myplots.magni.fRscale .* (myplots.magni.(myplots.magni.viscDep{i}).bandW ./ (myplots.magni.(myplots.magni.viscDep{i}).Xp.*myplots.magni.rhop.*(myplots.magni.A_O./myplots.magni.trench)));
        myplots.magni.(myplots.magni.viscDep{i}).tempfR = myplots.magni.(myplots.magni.viscDep{i}).bandfR_t;
        myplots.magni.(myplots.magni.viscDep{i}).tempfR(myplots.magni.(myplots.magni.viscDep{i}).tempfR < eps) = eps; myplots.magni.(myplots.magni.viscDep{i}).tempfR(myplots.magni.(myplots.magni.viscDep{i}).tempfR >= 1) = 1-eps;
        myplots.magni.(myplots.magni.viscDep{i}).bandfR_t = myplots.magni.(myplots.magni.viscDep{i}).tempfR;
        myplots.magni.(myplots.magni.viscDep{i}).extremefR_t = myplots.magni.fRscale .* (myplots.magni.(myplots.magni.viscDep{i}).extremeW ./ (myplots.magni.(myplots.magni.viscDep{i}).Xp.*myplots.magni.rhop.*(myplots.magni.A_O./myplots.magni.trench)));
        myplots.magni.(myplots.magni.viscDep{i}).tempfR = myplots.magni.(myplots.magni.viscDep{i}).extremefR_t;
        myplots.magni.(myplots.magni.viscDep{i}).tempfR(myplots.magni.(myplots.magni.viscDep{i}).tempfR < eps) = eps; myplots.magni.(myplots.magni.viscDep{i}).tempfR(myplots.magni.(myplots.magni.viscDep{i}).tempfR >= 1) = 1-eps;
        myplots.magni.(myplots.magni.viscDep{i}).extremefR_t = myplots.magni.(myplots.magni.viscDep{i}).tempfR;

        % actually plot
        if strcmp('conc',myplots.magni.viscDep{i})
            %fill(myplots.TS.timeBand, 100 .* myplots.magni.(myplots.magni.viscDep{i}).extremefR_t, [.85 .325 .098], 'FaceAlpha', 0.15, 'EdgeColor', 'none')
            %fill(myplots.TS.timeBand, 100 .* myplots.magni.(myplots.magni.viscDep{i}).bandfR_t, [.85 .325 .098], 'FaceAlpha', 0.125, 'EdgeColor', 'none')
            plot(myplots.magni.timeInterval, 100 .* myplots.magni.(myplots.magni.viscDep{i}).avgfR_t, 'r-')
        else
            %fill(myplots.TS.timeBand, 100 .* myplots.magni.(myplots.magni.viscDep{i}).extremefR_t, [0 .447 .741], 'FaceAlpha', 0.15, 'EdgeColor', 'none')
            %fill(myplots.TS.timeBand, 100 .* myplots.magni.(myplots.magni.viscDep{i}).bandfR_t, [0 .447 .741], 'FaceAlpha', 0.125, 'EdgeColor', 'none')
            plot(myplots.magni.timeInterval, 100 .* myplots.magni.(myplots.magni.viscDep{i}).avgfR_t, 'b-')
        end
    end

figure(myplots.magni.TSs)

% checking each parameter to see which one we should vary
myplots.magni.magniFig = figure;
for j = 1:length(myplots.magni.labels)
    % basic subplot stuff
    nexttile
    hold on; grid on;
    fontsize(gca, scale=1.33)
    xlabel('Time [Gyr]'); xlim([myplots.magni.timeInterval(1) myplots.magni.timeInterval(end)]);
    title(myplots.magni.labels{j})
    ylabel('f_R(t) [%]'); set(gca, 'YScale', 'log'); ylim([1e-3 1e2])
    % calculating, plotting fR(t) using tweaked Magni parameters
    for i = 1:length(myplots.magni.viscDep)
        % calculate W
        for k = 1:length(myplots.magni.varyParam)
            if strcmp('Varying [u]', myplots.magni.labels{j})
                % calculate W, fR
                myplots.magni.(myplots.magni.viscDep{i}).W(k,:) = 1e5 .* ((myplots.magni.W(1) .* myplots.magni.varyParam(k) .* myplots.magni.(myplots.magni.viscDep{i}).avgu) + ...
                                                                          (myplots.magni.W(2) .* myplots.magni.(myplots.magni.viscDep{i}).avga) + ...
                                                                          (myplots.magni.W(3) .* myplots.magni.(myplots.magni.viscDep{i}).avgT) + ...
                                                                         myplots.magni.W(4));
            elseif strcmp('Varying [a]', myplots.magni.labels{j})
                myplots.magni.(myplots.magni.viscDep{i}).W(k,:) = 1e5 .* ((myplots.magni.W(1) .* myplots.magni.(myplots.magni.viscDep{i}).avgu) + ...
                                                                          (myplots.magni.W(2) .* myplots.magni.varyParam(k) .* myplots.magni.(myplots.magni.viscDep{i}).avga) + ...
                                                                          (myplots.magni.W(3) .* myplots.magni.(myplots.magni.viscDep{i}).avgT) + ...
                                                                         myplots.magni.W(4));
            elseif strcmp('Varying [T]', myplots.magni.labels{j})
                myplots.magni.(myplots.magni.viscDep{i}).W(k,:) = 1e5 .* ((myplots.magni.W(1) .* myplots.magni.(myplots.magni.viscDep{i}).avgu) + ...
                                                                          (myplots.magni.W(2) .* myplots.magni.(myplots.magni.viscDep{i}).avga) + ...
                                                                          (myplots.magni.W(3) .* myplots.magni.varyParam(k) .* myplots.magni.(myplots.magni.viscDep{i}).avgT) + ...
                                                                         myplots.magni.W(4));
            elseif strcmp('Varying [const]', myplots.magni.labels{j})
                myplots.magni.(myplots.magni.viscDep{i}).W(k,:) = 1e5 .* ((myplots.magni.W(1) .* myplots.magni.(myplots.magni.viscDep{i}).avgu) + ...
                                                                          (myplots.magni.W(2) .* myplots.magni.(myplots.magni.viscDep{i}).avga) + ...
                                                                          (myplots.magni.W(3) .* myplots.magni.(myplots.magni.viscDep{i}).avgT) + ...
                                                                         myplots.magni.W(4) .* myplots.magni.varyParam(k));
            end
        end
        % calculate fR_t
        for k = 1:length(myplots.magni.varyParam)
            myplots.magni.(myplots.magni.viscDep{i}).fR_t(k,:) = myplots.magni.fRscale .* (myplots.magni.(myplots.magni.viscDep{i}).W(k,:) ./ (myplots.magni.(myplots.magni.viscDep{i}).Xp*myplots.magni.rhop*(myplots.magni.A_O./myplots.magni.trench)));
            % bound to [0, 1]
            myplots.magni.tempfR = myplots.magni.(myplots.magni.viscDep{i}).fR_t(k, :);
            myplots.magni.tempfR(myplots.magni.tempfR < eps) = eps; myplots.magni.tempfR(myplots.magni.tempfR >= 1) = 1-eps;
            myplots.magni.(myplots.magni.viscDep{i}).fR_t(k, :) = myplots.magni.tempfR;
        end

        % plot
        plot(NaN, NaN, 'k-.'); plot(NaN, NaN, 'k-'); plot(NaN, NaN, 'k--');
        if strcmp('conc',myplots.magni.viscDep{i})
            for k = 1:length(myplots.magni.varyParam)
                plot(myplots.magni.timeInterval, 100 .* myplots.magni.(myplots.magni.viscDep{i}).fR_t(k,:), 'Color', 'r', 'LineStyle', myplots.magni.linestyles{k})
            end
        else
            for k = 1:length(myplots.magni.varyParam)
                plot(myplots.magni.timeInterval, 100 .* myplots.magni.(myplots.magni.viscDep{i}).fR_t(k,:), 'Color', 'b', 'LineStyle', myplots.magni.linestyles{k})
            end
        end
        if j == length(myplots.magni.labels)
            legend('0.5*param', 'param', '2*param', 'Location', 'southeast')
        end
    end
end
figure(myplots.magni.magniFig)

% meshgrid of plate speeds, temperatures to find feasible range for Magni (re)scalings
myplots.magni.fRtContour = figure;
    hold on; grid on;
    fontsize(gca, scale=1.33);
    myplots.magni.scalings.Trange = linspace(1000, 2000, 1001);
    myplots.magni.scalings.urange = logspace(-1, 2, 1001);
    [myplots.magni.scalings.U, myplots.magni.scalings.T] = meshgrid(myplots.magni.scalings.urange, myplots.magni.scalings.Trange);
    myplots.magni.scalings.DELTA = 2 .* sqrt((myplots.magni.kappa .* myplots.magni.L) ./ (myplots.magni.scalings.U .* myplots.magni.sec2year ./ 100));
    myplots.magni.scalings.A = ((myplots.magni.scalings.DELTA ./ myplots.magni.W(5)).^2) ./ myplots.magni.kappa .* myplots.magni.sec2year .* 1e-6;
    myplots.magni.scalings.W = 1e5 .* (((myplots.magni.W(1) * 1e0) .* myplots.magni.scalings.U) + ...
                                       (myplots.magni.W(2) .* myplots.magni.scalings.A) + ...
                                       (myplots.magni.W(3) .* myplots.magni.scalings.T) + ...
                                       myplots.magni.W(4) * 6e-1);
    myplots.magni.scalings.fRt = 5e2 * (100 .* myplots.magni.scalings.W ./ (myplots.magni.fug.Xp.*myplots.magni.rhop.*(myplots.magni.A_O./myplots.magni.trench)));
    myplots.magni.scalings.fRt(myplots.magni.scalings.fRt < 0) = NaN;
    myplots.magni.scalings.fRt(myplots.magni.scalings.fRt > 100) = 100;
    contourf(myplots.magni.scalings.U, myplots.magni.scalings.T, myplots.magni.scalings.fRt)
    xlabel('u [cm/yr]'); ylabel(['T [',char(176),'C]']);
    set(gca, 'XScale', 'log')
    c = colorbar;
    set(gca, 'ColorScale', 'log')
    c.Label.String = 'f_R(t) [%]';
    myplots.magni.scalings.Earthregion = polyshape([4 1270; 4 1445; 5 1445; 5 1270]);
    plot(myplots.magni.scalings.Earthregion, 'FaceColor', 'None', 'LineWidth', 2)
    plot(myplots.magni.conc.avgu, myplots.magni.conc.avgT, 'r-')
    plot(myplots.magni.fug.avgu, myplots.magni.fug.avgT, 'b-')
    scatter(4.5, 1350, 100, 'k', 'filled', 'pentagram')
figure(myplots.magni.fRtContour)
myplots = rmfield(myplots, {'TS', 'magni'});

%% FIND THE RANGES OF PARAMETERS TO TARGET
SectionBreak()
disp('NOW LOADING: fR(t) PARAM RANGES')

myplots.MCsearch.conc = [];
myplots.MCsearch.fug = [];
processed.settings.timeDirections = {'fwd'};
for i = 1:length(processed.settings.timeDirections)
    % setting the file path to the baseline cases
    processed.settings.dirPath = [pwd,'\models\',processed.settings.timeDirections{i},'\'];

    % investigate the regular or targeted MC data
    % processed.settings.baseDirInfo = dir([processed.settings.dirPath,'var\', '*set*VarfR_MC.mat']);
    processed.settings.baseDirInfo = dir([processed.settings.dirPath,'var\', '*set*VarfR_VarphiRum_MC.mat']);

    % loading those cases
    for j = 1:length(processed.settings.baseDirInfo)
        % display that we're loading things
        if contains(processed.settings.baseDirInfo(j).name, "zeroedMagni")
            continue
        end
        disp(['--- ', processed.settings.baseDirInfo(j).name])

        % checking for fwd, rev
        if strcmp('fwd', processed.settings.timeDirections{i})
            % checking for fug, conc
            if contains(processed.settings.baseDirInfo(j).name, 'FUG')
                processed.fwd.fug.fRt = load([processed.settings.dirPath, 'var\', processed.settings.baseDirInfo(j).name]);
                myplots.MCsearch.fug = [myplots.MCsearch.fug, processed.fwd.fug.fRt.ECs.solved.MC];
            else
                processed.fwd.conc.fRt = load([processed.settings.dirPath, 'var\', processed.settings.baseDirInfo(j).name]);
                myplots.MCsearch.conc = [myplots.MCsearch.conc, processed.fwd.conc.fRt.ECs.solved.MC];
            end
        else
            if contains(processed.settings.baseDirInfo(j).name, 'FUG')
                processed.rev.fug.fRt = load([processed.settings.dirPath, 'var\', processed.settings.baseDirInfo(j).name]);
            else
                processed.rev.conc.fRt = load([processed.settings.dirPath, 'var\', processed.settings.baseDirInfo(j).name]);
            end
        end
    end
end

% inputs, outputs
% myplots.histo.fRt.inputs = {'T0', 'totalH2Omass', 'H2Opartition', 'Xppd', 'fRGinMantle', ...
%                             'C_U235', 'C_U238', 'C_Th232','C_K40', ...
%                             'E', 'eta_pd', 'fR0', 'MagniCScale', 'fum'};
myplots.histo.fRt.inputs = {'T0', 'totalH2Omass', 'H2Opartition', 'Xppd', 'fRGinMantle', ...
                            'C_U235', 'C_U238', 'C_Th232','C_K40', ...
                            'E', 'eta_pd', 'fR0', 'MagniCScale', 'phiRum'};

for i = 1:length(myplots.histo.fRt.inputs)
    % get the base data
    if isempty(myplots.MCsearch.conc)
        myplots.histo.success.concData = NaN;
    else
        myplots.histo.success.concData = [myplots.MCsearch.conc.(myplots.histo.fRt.inputs{i})];
        if ismember(myplots.histo.fRt.inputs{i}, {'C_U235', 'C_U238', 'C_Th232','C_K40'})
            myplots.histo.success.concData = myplots.histo.success.concData ./ ([myplots.MCsearch.conc.fRGinMantle]./100);
        end
    end    
    myplots.histo.success.fugData = [myplots.MCsearch.fug.(myplots.histo.fRt.inputs{i})];
    if ismember(myplots.histo.fRt.inputs{i}, {'C_U235', 'C_U238', 'C_Th232','C_K40'})
        myplots.histo.success.fugData = myplots.histo.success.fugData ./ ([myplots.MCsearch.fug.fRGinMantle]./100);
    end

    % print out appropriate ranges to use
    disp(['Explore range for successful ',myplots.histo.fRt.inputs{i},': [',num2str(0.9*min([myplots.histo.success.concData myplots.histo.success.fugData], [], 'omitnan')),' - ', ...
                                                                            num2str(1.1*max([myplots.histo.success.concData myplots.histo.success.fugData], [], 'omitnan')),']'])
end
processed = [];
myplots = [];

%% TARGETED OUTPUT ASSESSMENT USING JUST BEST-KNOWN CONSTRAINTS
SectionBreak()
disp('NOW LOADING: TARGETED fR(t) MC')

processed.fwd.conc.fRt.successes = [];
processed.fwd.fug.fRt.successes = [];
processed.fwd.conc.fRt.MC = [];
processed.fwd.fug.fRt.MC = [];
processed.settings.timeDirections = {'fwd'};
for i = 1:length(processed.settings.timeDirections)
    % setting the file path to the baseline cases
    processed.settings.dirPath = [pwd,'\models\',processed.settings.timeDirections{i},'\'];

    % investigate the regular or targeted MC data
    %processed.settings.baseDirInfo = dir([processed.settings.dirPath,'var\', '*set*VarfR_EssentialTargetedMC.mat']);
    %processed.settings.baseDirInfo = dir([processed.settings.dirPath,'var\', '*set*VarfR_VarphiRum_EssentialTargetedMC.mat']);
    %processed.settings.baseDirInfo = dir([processed.settings.dirPath,'var\', '*set*VarfR_FullTargetedMC.mat']);
    processed.settings.baseDirInfo = dir([processed.settings.dirPath,'var\', '*set*VarfR_VarphiRum_FullTargetedMC.mat']);

    % loading those cases
    for j = 1:length(processed.settings.baseDirInfo)
        % display that we're loading things
        disp(['--- ', processed.settings.baseDirInfo(j).name])

        % checking for fwd, rev
        if strcmp('fwd', processed.settings.timeDirections{i})
            % checking for fug, conc
            if contains(processed.settings.baseDirInfo(j).name, 'FUG')
                processed.temp = load([processed.settings.dirPath, 'var\', processed.settings.baseDirInfo(j).name]);
                processed.fwd.fug.fRt.successes = [processed.fwd.fug.fRt.successes, processed.temp.ECs.solved.successes];
                processed.fwd.fug.fRt.MC = [processed.fwd.fug.fRt.MC, processed.temp.ECs.solved.MC];
            else
                processed.temp = load([processed.settings.dirPath, 'var\', processed.settings.baseDirInfo(j).name]);
                processed.fwd.conc.fRt.successes = [processed.fwd.conc.fRt.successes, processed.temp.ECs.solved.successes];
                processed.fwd.conc.fRt.MC = [processed.fwd.conc.fRt.MC, processed.temp.ECs.solved.MC];
            end
        else
            if contains(processed.settings.baseDirInfo(j).name, 'FUG')
                processed.rev.fug.fRt = load([processed.settings.dirPath, 'var\', processed.settings.baseDirInfo(j).name]);
            else
                processed.rev.conc.fRt = load([processed.settings.dirPath, 'var\', processed.settings.baseDirInfo(j).name]);
            end
        end
    end
end
disp(['TOTAL SUCCESSES: CONC -- ', num2str(sum(processed.fwd.conc.fRt.successes)),', FUG -- ', num2str(sum(processed.fwd.fug.fRt.successes))])
myplots.histo.bins = 25;
myplots.histo.labels = ('a':'z');

%%%%%%%% SUCCESS OUTPUT FACTOR MAPPING %%%%%%%%
% actually plotting outputs
%myplots.histo.fRt.inputs = {'T0', 'totalH2Omass', 'H2Opartition', 'Xppd', 'fRGinMantle', ...
%                            'C_U235', 'C_U238', 'C_Th232','C_K40', 'E', 'eta_pd', 'fR0', 'MagniCScale', 'fum'};
myplots.histo.fRt.inputs = {'T0', 'totalH2Omass', 'H2Opartition', 'Xppd', 'fRGinMantle', ...
                           'C_U235', 'C_U238', 'C_Th232','C_K40', 'E', 'eta_pd', 'fR0', 'MagniCScale', 'phiRum'};
myplots.histo.fRt.outputs = {'Xumpd', 'fRpd', 'RminusDpd', 'Ureypd', 'tRminusD0'};

myplots.design.conc.outputParamsMat = NaN(sum(processed.fwd.conc.fRt.successes), length(myplots.histo.fRt.outputs));
for p = 1:length(myplots.histo.fRt.outputs)
    if strcmp('RminusDpd', myplots.histo.fRt.outputs{p})
        myplots.design.conc.outputData = log10([processed.fwd.conc.fRt.MC.(myplots.histo.fRt.outputs{p})]);
    else
        myplots.design.conc.outputData = [processed.fwd.conc.fRt.MC.(myplots.histo.fRt.outputs{p})];
    end
    myplots.design.conc.outputParamsMat(:, p) = myplots.design.conc.outputData;
end
myplots.design.fug.outputParamsMat = NaN(sum(processed.fwd.fug.fRt.successes), length(myplots.histo.fRt.outputs));
for p = 1:length(myplots.histo.fRt.outputs)
    if strcmp('RminusDpd', myplots.histo.fRt.outputs{p})
        myplots.design.fug.outputData = log10([processed.fwd.fug.fRt.MC.(myplots.histo.fRt.outputs{p})]);
    else
        myplots.design.fug.outputData = [processed.fwd.fug.fRt.MC.(myplots.histo.fRt.outputs{p})];
    end
    myplots.design.fug.outputParamsMat(:, p) = myplots.design.fug.outputData;
end
for u = 1:length(myplots.histo.fRt.outputs)
    if strcmp('Xumpd', myplots.histo.fRt.outputs{u})
        myplots.design.unc.(myplots.histo.fRt.outputs{u}) = [50 200];
    elseif strcmp('RminusDpd', myplots.histo.fRt.outputs{u})
        myplots.design.unc.(myplots.histo.fRt.outputs{u}) = log10([2.7e11 5.0e11]);
    elseif strcmp('Ureypd', myplots.histo.fRt.outputs{u})
        myplots.design.unc.(myplots.histo.fRt.outputs{u}) = [0.12 0.49];
    elseif strcmp('tRminusD0', myplots.histo.fRt.outputs{u})
        myplots.design.unc.(myplots.histo.fRt.outputs{u}) = [2 3];
    elseif strcmp('fRpd', myplots.histo.fRt.outputs{u})
        myplots.design.unc.(myplots.histo.fRt.outputs{u}) = [10 60];
    end
end
myplots.design.successOutputFig = figure;
for m = 1:(length(myplots.histo.fRt.outputs)-1)
    mm = m + 1;
    % actual plot info
    for n = 1:(length(myplots.histo.fRt.outputs)-1)
        spIdx = (m-1)*(length(myplots.histo.fRt.outputs)-1) + n;
        h = subplot(length(myplots.histo.fRt.outputs)-1, length(myplots.histo.fRt.outputs)-1, spIdx);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        % label info
        if mm == length(myplots.histo.fRt.outputs)
            if strcmp('Xumpd', myplots.histo.fRt.outputs{n})
                xlabel('\chi_{um} [ppm]');
            elseif strcmp('RminusDpd', myplots.histo.fRt.outputs{n})
                xlabel('log_{10}(R-D) [kg/yr]');
            elseif strcmp('fRpd', myplots.histo.fRt.outputs{n})
                xlabel('f_R [%]');
            elseif strcmp('Ureypd', myplots.histo.fRt.outputs{n})
                xlabel('Urey [-]');
            elseif strcmp('tRminusD0', myplots.histo.fRt.outputs{n})
                xlabel('t_{R=D} [Gyr]');
            end
        end
        if n == 1
            if strcmp('Xumpd', myplots.histo.fRt.outputs{mm})
                ylabel('\chi_{um} [ppm]');
            elseif strcmp('RminusDpd', myplots.histo.fRt.outputs{mm})
                ylabel('log_{10}(R-D) [kg/yr]');
            elseif strcmp('fRpd', myplots.histo.fRt.outputs{mm})
                ylabel('f_R [%]');
            elseif strcmp('Ureypd', myplots.histo.fRt.outputs{mm})
                ylabel('Urey [-]');
            elseif strcmp('tRminusD0', myplots.histo.fRt.outputs{mm})
                ylabel('t_{R=D} [Gyr]');
            end
        end
        if n >= mm
            set(h,'Visible','off');
            set(get(h,'Children'),'Visible','off');
        else
            scatter(myplots.design.conc.outputParamsMat(:, n), myplots.design.conc.outputParamsMat(:, mm), 10, [.85 .325 .098], 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
            scatter(myplots.design.fug.outputParamsMat(:, n), myplots.design.fug.outputParamsMat(:, mm), 10, [0 .447 .741], 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
            myplots.design.xunc = myplots.design.unc.(myplots.histo.fRt.outputs{n});
            myplots.design.yunc = myplots.design.unc.(myplots.histo.fRt.outputs{mm});
            myplots.design.uncregion = polyshape([myplots.design.xunc(1) myplots.design.yunc(1); ...
                                                  myplots.design.xunc(1) myplots.design.yunc(end); ...
                                                  myplots.design.xunc(end) myplots.design.yunc(end); ...
                                                  myplots.design.xunc(end) myplots.design.yunc(1)]);
            plot(myplots.design.uncregion, 'FaceColor', 'None', 'LineWidth', 2)
        end
    end
end
figure(myplots.design.successOutputFig)

% print out a table of input param ranges that also hit Xum, R-D, Urey present-day constraints
myplots.histo.fRt.vsuccessConc = ([processed.fwd.conc.fRt.MC.Xumpd] >= myplots.design.unc.Xumpd(1) & [processed.fwd.conc.fRt.MC.Xumpd] <= myplots.design.unc.Xumpd(end)) & ...
                                 ([processed.fwd.conc.fRt.MC.RminusDpd] >= 10.^(myplots.design.unc.RminusDpd(1)) & [processed.fwd.conc.fRt.MC.RminusDpd] <= 10.^(myplots.design.unc.RminusDpd(end))) & ...
                                 ([processed.fwd.conc.fRt.MC.Ureypd] >= myplots.design.unc.Ureypd(1) & [processed.fwd.conc.fRt.MC.Ureypd] <= myplots.design.unc.Ureypd(end)) & ...
                                 ([processed.fwd.conc.fRt.MC.fRpd] >= myplots.design.unc.fRpd(1) & [processed.fwd.conc.fRt.MC.fRpd] <= myplots.design.unc.fRpd(end));
myplots.histo.fRt.vsuccessFug = ([processed.fwd.fug.fRt.MC.Xumpd] >= myplots.design.unc.Xumpd(1) & [processed.fwd.fug.fRt.MC.Xumpd] <= myplots.design.unc.Xumpd(end)) & ...
                                ([processed.fwd.fug.fRt.MC.RminusDpd] >= 10.^(myplots.design.unc.RminusDpd(1)) & [processed.fwd.fug.fRt.MC.RminusDpd] <= 10.^(myplots.design.unc.RminusDpd(end))) & ...
                                ([processed.fwd.fug.fRt.MC.Ureypd] >= myplots.design.unc.Ureypd(1) & [processed.fwd.fug.fRt.MC.Ureypd] <= myplots.design.unc.Ureypd(end)) & ...
                                ([processed.fwd.fug.fRt.MC.fRpd] >= myplots.design.unc.fRpd(1) & [processed.fwd.fug.fRt.MC.fRpd] <= myplots.design.unc.fRpd(end));
for i = 1:length(myplots.histo.fRt.inputs)
    % get the base data
    myplots.histo.vsuccess.concData = [processed.fwd.conc.fRt.MC.(myplots.histo.fRt.inputs{i})];
    myplots.histo.vsuccess.fugData = [processed.fwd.fug.fRt.MC.(myplots.histo.fRt.inputs{i})];
    if ismember(myplots.histo.fRt.inputs{i}, {'C_U235', 'C_U238', 'C_Th232','C_K40'})
        myplots.histo.vsuccess.concData = myplots.histo.vsuccess.concData ./ ([processed.fwd.conc.fRt.MC.fRGinMantle]./100);
        myplots.histo.vsuccess.fugData = myplots.histo.vsuccess.fugData ./ ([processed.fwd.fug.fRt.MC.fRGinMantle]./100);
    end

    % print out appropriate ranges to use
    disp(['PD successful range for ',myplots.histo.fRt.inputs{i},': [',num2str(min([myplots.histo.vsuccess.concData(myplots.histo.fRt.vsuccessConc) myplots.histo.vsuccess.fugData(myplots.histo.fRt.vsuccessFug)])),' - ', ...
                                                                       num2str(max([myplots.histo.vsuccess.concData(myplots.histo.fRt.vsuccessConc) myplots.histo.vsuccess.fugData(myplots.histo.fRt.vsuccessFug)])),']'])
end
disp(['TOTAL PD SUCCESSES: CONC -- ', num2str(sum(myplots.histo.fRt.vsuccessConc)),', FUG -- ', num2str(sum(myplots.histo.fRt.vsuccessFug))])

% ADJUSTING tRminusD0
myplots.adjustedRminusD.conc = NaN(sum(processed.fwd.conc.fRt.successes), 2);
myplots.adjustedRminusD.fug = NaN(sum(processed.fwd.fug.fRt.successes), 2);
myplots.adjustedRminusD.modelViscDep = {'conc', 'fug'};
for vd = 1:length(myplots.adjustedRminusD.modelViscDep)
    for i = 1:size(myplots.adjustedRminusD.(myplots.adjustedRminusD.modelViscDep{vd}), 1)
        tRminusD = processed.fwd.(myplots.adjustedRminusD.modelViscDep{vd}).fRt.MC(i).tRminusD0;
        myplots.adjustedRminusD.(myplots.adjustedRminusD.modelViscDep{vd})(i, 1) = tRminusD;
    end
end

processed = [];
processed.settings.dirPath = [pwd,'\models\fwd\'];
processed.fwd.conc.fRt.TS = [];
processed.fwd.fug.fRt.TS = [];
%processed.settings.baseDirInfo = dir([processed.settings.dirPath,'var\', '*set*VarfR_EssentialTargetedTS.mat']);
%processed.settings.baseDirInfo = dir([processed.settings.dirPath,'var\', '*set*VarfR_VarphiRum_EssentialTargetedTS.mat']);
%processed.settings.baseDirInfo = dir([processed.settings.dirPath,'var\', '*set*VarfR_FullTargetedTS.mat']);
processed.settings.baseDirInfo = dir([processed.settings.dirPath,'var\', '*set*VarfR_VarphiRum_FullTargetedTS.mat']);

% loading those cases
for j = 1:length(processed.settings.baseDirInfo)
    % display that we're loading things
    disp(['--- ', processed.settings.baseDirInfo(j).name])
    if contains(processed.settings.baseDirInfo(j).name, 'FUG')
        processed.temp = load([pwd, '\models\fwd\var\', processed.settings.baseDirInfo(j).name]);
        processed.fwd.fug.fRt.TS = [processed.fwd.fug.fRt.TS, processed.temp.ECs.solved.TS];
    else
        processed.temp = load([pwd, '\models\fwd\var\', processed.settings.baseDirInfo(j).name]);
        processed.fwd.conc.fRt.TS = [processed.fwd.conc.fRt.TS, processed.temp.ECs.solved.TS];
    end
end
for vd = 1:length(myplots.adjustedRminusD.modelViscDep)
    for i = 1:size(myplots.adjustedRminusD.(myplots.adjustedRminusD.modelViscDep{vd}), 1)
        RminusD = processed.fwd.(myplots.adjustedRminusD.modelViscDep{vd}).fRt.TS(i).RminusD;
        RminusDat2Gyr = interp1(processed.fwd.(myplots.adjustedRminusD.modelViscDep{vd}).fRt.TS(i).t, RminusD, 2);
        myplots.RminusDat2Gyr.(myplots.adjustedRminusD.modelViscDep{vd})(i, 2) = RminusDat2Gyr/(3.85e11)*100;
    end
end
myplots.adjustedRminusD.adjustedFig = figure;
hold on; grid on;
fontsize(gca, scale=1.33)
xlabel('R-D at 2 Gyr as Fraction of Present-Day Average [%]'); ylabel('Number of Successes [-]');
histogram(myplots.RminusDat2Gyr.conc(:, 2), 'BinWidth', 0.1, 'FaceColor', [.85 .325 .098], 'FaceAlpha', 0.75, 'EdgeColor', [.85 .325 .098], 'EdgeAlpha', 0.75)
histogram(myplots.RminusDat2Gyr.fug(:, 2), 'BinWidth', 0.1, 'FaceColor', [0 .447 .741], 'FaceAlpha', 0.75, 'EdgeColor', [0 .447 .741], 'EdgeAlpha', 0.75)

%% LOAD IN ALL-CONSTRAINTS TARGETED VARIANT CASES FOR TIME-SERIES ANALYSIS
SectionBreak()
disp('NOW LOADING: TARGETED fR(t) TS')
% processed.settings.baseDirInfo = dir([pwd, '\models\fwd\var\', '*set*VarfR_EssentialTargetedTS.mat']);
% processed.settings.baseDirInfo = dir([pwd, '\models\fwd\var\', '*set*VarfR_VarphiRum_EssentialTargetedTS.mat']);
% processed.settings.baseDirInfo = dir([pwd, '\models\fwd\var\', '*set*VarfR_FullTargetedTS.mat']);
processed.settings.baseDirInfo = dir([pwd, '\models\fwd\var\', '*set*VarfR_VarphiRum_FullTargetedTS.mat']);

processed.fwd.conc.fRt.successes = [];
processed.fwd.fug.fRt.successes = [];
processed.fwd.conc.fRt.TS = [];
processed.fwd.fug.fRt.TS = [];
% loading those cases
for j = 1:length(processed.settings.baseDirInfo)
    % display that we're loading things
    disp(['--- ', processed.settings.baseDirInfo(j).name])
    if contains(processed.settings.baseDirInfo(j).name, 'FUG')
        processed.temp = load([pwd, '\models\fwd\var\', processed.settings.baseDirInfo(j).name]);
        processed.fwd.fug.fRt.successes = [processed.fwd.fug.fRt.successes, processed.temp.ECs.solved.successes];
        processed.fwd.fug.fRt.TS = [processed.fwd.fug.fRt.TS, processed.temp.ECs.solved.TS];
    else
        processed.temp = load([pwd, '\models\fwd\var\', processed.settings.baseDirInfo(j).name]);
        processed.fwd.conc.fRt.successes = [processed.fwd.conc.fRt.successes, processed.temp.ECs.solved.successes];
        processed.fwd.conc.fRt.TS = [processed.fwd.conc.fRt.TS, processed.temp.ECs.solved.TS];
    end
end

myplots.TS.avgFig = figure;
myplots.TS.labels = ('a':'z');
myplots.TS.params = {'T', 'Xum', 'u', 'Qs', 'oceanmass', 'mantleH2Omass', 'eta_um', 'UreyRatio', 'RminusD', 'fH2Oum'};
myplots.TS.timeInterval = linspace(processed.temp.ECs.range.t(1) .* processed.temp.ECs.conversion.sec2yr ./ 1e9, ...
                                   processed.temp.ECs.range.t(end) .* processed.temp.ECs.conversion.sec2yr ./ 1e9, ...
                                   500);
myplots.ageEarth = processed.temp.ECs.pd.ageEarth;

for i = 1:length(myplots.TS.params)
    % interpolating the successful CONC data
    myplots.TS.interpConc = NaN(length(processed.fwd.conc.fRt.TS), length(myplots.TS.timeInterval));
    for j = 1:length(processed.fwd.conc.fRt.successes)
        myplots.TS.interpConc(j, :) = interp1(processed.fwd.conc.fRt.TS(j).t, ...
                                              processed.fwd.conc.fRt.TS(j).(myplots.TS.params{i}), ...
                                              myplots.TS.timeInterval);
    end
    % interpolating the successful FUG data
    myplots.TS.interpFug = NaN(length(processed.fwd.fug.fRt.TS), length(myplots.TS.timeInterval));
    for j = 1:length(processed.fwd.fug.fRt.successes)
        myplots.TS.interpFug(j, :) = interp1(processed.fwd.fug.fRt.TS(j).t, ...
                                             processed.fwd.fug.fRt.TS(j).(myplots.TS.params{i}), ...
                                             myplots.TS.timeInterval);
    end

    % median, 25/75 percentiles, min/max
    myplots.TS.avgConc = median(myplots.TS.interpConc, 'omitnan');
    myplots.TS.prctile25Conc = prctile(myplots.TS.interpConc, 25);
    myplots.TS.prctile75Conc = prctile(myplots.TS.interpConc, 75);
    myplots.TS.minConc = min(myplots.TS.interpConc, [], 'omitnan');
    myplots.TS.maxConc = max(myplots.TS.interpConc, [], 'omitnan');
    myplots.TS.avgFug = median(myplots.TS.interpFug, 'omitnan');
    myplots.TS.prctile25Fug = prctile(myplots.TS.interpFug, 25);
    myplots.TS.prctile75Fug = prctile(myplots.TS.interpFug, 75);
    myplots.TS.minFug = min(myplots.TS.interpFug, [], 'omitnan');
    myplots.TS.maxFug = max(myplots.TS.interpFug, [], 'omitnan');

    % setting up uncertainty bands
    myplots.TS.timeBand = [myplots.TS.timeInterval, myplots.TS.timeInterval(end:-1:1)];
    myplots.TS.concBand = [myplots.TS.prctile75Conc, myplots.TS.prctile25Conc(end:-1:1)];
    myplots.TS.concExtremes = [myplots.TS.maxConc, myplots.TS.minConc(end:-1:1)];
    myplots.TS.fugBand = [myplots.TS.prctile75Fug, myplots.TS.prctile25Fug(end:-1:1)];
    myplots.TS.fugExtremes = [myplots.TS.maxFug, myplots.TS.minFug(end:-1:1)];

    % conditional by parameter
    if strcmp('T', myplots.TS.params{i})
        subplot(3, 4, [1,2]);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel(['T [',char(176),'C]']); ylim([1150, 1850]);
        errorbar(myplots.ageEarth .* processed.temp.ECs.conversion.sec2yr ./ 1e9, errors.T(1), errors.T(2), errors.T(3), 'k-', 'linewidth', 1)
    elseif strcmp('Xum', myplots.TS.params{i})
        subplot(3, 4, [3,4]);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('\chi_{um} [ppm]'); set(gca, 'YScale', 'log'); ylim([0.5 750]); yticks([1 10 100])
        errorbar(myplots.ageEarth .* processed.temp.ECs.conversion.sec2yr ./ 1e9, errors.Xum(1), errors.Xum(2), errors.Xum(3), 'k-', 'linewidth', 1)
        myplots.TS.concBand(myplots.TS.concBand < 0) = 1e-10;
        myplots.TS.fugBand(myplots.TS.fugBand < 0) = 1e-10;
    elseif strcmp('u', myplots.TS.params{i})
        subplot(3, 4, 5);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('u [cm yr^{-1}]'); set(gca, 'YScale', 'log'); ylim([2 15]);
        myplots.TS.concBand(myplots.TS.concBand < eps) = eps;
        myplots.TS.fugBand(myplots.TS.fugBand < eps) = eps;
        errorbar(myplots.ageEarth .* processed.temp.ECs.conversion.sec2yr ./ 1e9, errors.u(1), errors.u(2), errors.u(3), 'k-', 'linewidth', 1)
    elseif strcmp('Qs', myplots.TS.params{i})
        subplot(3, 4, 6);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('Q_s [TW]'); set(gca, 'YScale', 'log'); ylim([25 70]);
        errorbar(myplots.ageEarth .* processed.temp.ECs.conversion.sec2yr ./ 1e9, errors.Qs(1), errors.Qs(2), errors.Qs(3), 'k-', 'linewidth', 1)
        myplots.TS.concBand(myplots.TS.concBand < 0) = 1e-10;
        myplots.TS.fugBand(myplots.TS.fugBand < 0) = 1e-10;
    elseif strcmp('oceanmass', myplots.TS.params{i})
        subplot(3, 4, 7);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('Surface Water [OM]'); ylim([0.25 1.5]);
        errorbar(myplots.ageEarth .* processed.temp.ECs.conversion.sec2yr ./ 1e9, errors.OM(1), errors.OM(2), errors.OM(3), 'k-', 'linewidth', 1)
    elseif strcmp('mantleH2Omass', myplots.TS.params{i})
        subplot(3, 4, 8);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('Mantle Water [OM]'); set(gca, 'YScale', 'log'); ylim([0 2.5]); yticks([0.001 0.01 0.1 1])
        errorbar(myplots.ageEarth .* processed.temp.ECs.conversion.sec2yr ./ 1e9, errors.mantleOM(1), errors.mantleOM(2), errors.mantleOM(3), 'k-', 'linewidth', 1)
    elseif strcmp('eta_um', myplots.TS.params{i})
        subplot(3, 4, 9);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('\eta_{um} [Pa s]'); set(gca, 'YScale', 'log'); ylim([5e19 2e22]); yticks([1e20 1e21 1e22])
        errorbar(myplots.ageEarth .* processed.temp.ECs.conversion.sec2yr ./ 1e9, errors.eta_um(1), errors.eta_um(2), errors.eta_um(3), 'k-', 'linewidth', 1)
    elseif strcmp('UreyRatio', myplots.TS.params{i})
        subplot(3, 4, 10);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('Urey Ratio [-]'); ylim([0 1.75])
        errorbar(myplots.ageEarth .* processed.temp.ECs.conversion.sec2yr ./ 1e9, errors.UreyRatio(1), errors.UreyRatio(2), errors.UreyRatio(3), 'k-', 'linewidth', 1)
    elseif strcmp('RminusD', myplots.TS.params{i})
        subplot(3, 4, 11);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('R - D [kg yr^{-1}]'); ylim([-1e12 2.5e12])
        errorbar(myplots.ageEarth .* processed.temp.ECs.conversion.sec2yr ./ 1e9, errors.RminusD_Korenaga(1), errors.RminusD_Korenaga(2), errors.RminusD_Korenaga(3), 'k-', 'linewidth', 1)
    elseif strcmp('fH2Oum', myplots.TS.params{i})
        subplot(3, 4, 12);
        hold on; grid on;
        fontsize(gca, scale=1.33)
        ylabel('f_{H_2O, um} [%]'); ylim([0 65]);
        errorbar(myplots.ageEarth .* processed.temp.ECs.conversion.sec2yr ./ 1e9, errors.fH2Oum(1), errors.fH2Oum(2), errors.fH2Oum(3), 'k-', 'linewidth', 1)
    end
    if ~strcmp('T', myplots.TS.params{i}) && ~strcmp('Xum', myplots.TS.params{i})
        xticks([0 2 4]);
    end

    % plotting
    xlabel('Time [Gyr]'); xlim([myplots.TS.timeInterval(1) myplots.TS.timeInterval(end)]);
    title(['(',myplots.TS.labels(i),')'])
    if strcmp('T', myplots.TS.params{i})
        scatter(geochem.herzbergAgeData, geochem.herzbergTempData, 40, [46/255, 201/255, 123/255], 'filled')
        errorbar(geochem.condieUniqueAges, geochem.condieTempAvgs, geochem.condieTempMinDiffs, geochem.condieTempMaxDiffs, ...
                 "Marker", 'square', "MarkerSize", 6.5, ...
                 "MarkerEdgeColor", [219/255 186/255 39/255], ...
                 "MarkerFaceColor", [219/255 186/255 39/255], ...
                 'CapSize', 0, "Color", [219/255 186/255 39/255], "LineStyle", "none")
    end
    fill(myplots.TS.timeBand, myplots.TS.concExtremes, [.85 .325 .098], 'FaceAlpha', 0.15, 'EdgeColor', 'none')
    fill(myplots.TS.timeBand, myplots.TS.fugExtremes, [0 .447 .741], 'FaceAlpha', 0.15, 'EdgeColor', 'none')
    fill(myplots.TS.timeBand, myplots.TS.concBand, [.85 .325 .098], 'FaceAlpha', 0.35, 'EdgeColor', 'none')
    fill(myplots.TS.timeBand, myplots.TS.fugBand, [0 .447 .741], 'FaceAlpha', 0.35, 'EdgeColor', 'none')
    plot(myplots.TS.timeInterval, myplots.TS.avgConc, 'r-')
    plot(myplots.TS.timeInterval, myplots.TS.avgFug, 'b-')
end
% show plot
figure(myplots.TS.avgFig)

%% SUCCESS INPUT FACTOR MAPPING
SectionBreak()
disp('NOW LOADING: TARGETED fR(t) MC')

processed.fum.conc.fRt.successes = [];
processed.fum.fug.fRt.successes = [];
processed.phiRum.conc.fRt.successes = [];
processed.phiRum.fug.fRt.successes = [];
processed.fum.conc.fRt.MC = [];
processed.fum.fug.fRt.MC = [];
processed.phiRum.conc.fRt.MC = [];
processed.phiRum.fug.fRt.MC = [];
processed.settings.timeDirections = {'fwd'};
for i = 1:length(processed.settings.timeDirections)
    % setting the file path to the baseline cases
    processed.settings.dirPath = [pwd,'\models\',processed.settings.timeDirections{i},'\'];

    % investigate the additional targeted MC data
    processed.settings.baseDirInfo = dir([processed.settings.dirPath,'var\', '*FullTargetedMC.mat']);

    % loading those cases
    for j = 1:length(processed.settings.baseDirInfo)
        % display that we're loading things
        disp(['--- ', processed.settings.baseDirInfo(j).name])

        % checking for fum, phiRum
        if contains(processed.settings.baseDirInfo(j).name, 'VarphiRum')
            % checking for fug, conc
            if contains(processed.settings.baseDirInfo(j).name, 'FUG')
                processed.temp = load([processed.settings.dirPath, 'var\', processed.settings.baseDirInfo(j).name]);
                processed.phiRum.fug.fRt.successes = [processed.phiRum.fug.fRt.successes, processed.temp.ECs.solved.successes];
                processed.phiRum.fug.fRt.MC = [processed.phiRum.fug.fRt.MC, processed.temp.ECs.solved.MC];
            else
                processed.temp = load([processed.settings.dirPath, 'var\', processed.settings.baseDirInfo(j).name]);
                processed.phiRum.conc.fRt.successes = [processed.phiRum.conc.fRt.successes, processed.temp.ECs.solved.successes];
                processed.phiRum.conc.fRt.MC = [processed.phiRum.conc.fRt.MC, processed.temp.ECs.solved.MC];
            end
        else
            if contains(processed.settings.baseDirInfo(j).name, 'FUG')
                processed.temp = load([processed.settings.dirPath, 'var\', processed.settings.baseDirInfo(j).name]);
                processed.fum.fug.fRt.successes = [processed.fum.fug.fRt.successes, processed.temp.ECs.solved.successes];
                processed.fum.fug.fRt.MC = [processed.fum.fug.fRt.MC, processed.temp.ECs.solved.MC];
            else
                processed.temp = load([processed.settings.dirPath, 'var\', processed.settings.baseDirInfo(j).name]);
                processed.fum.conc.fRt.successes = [processed.fum.conc.fRt.successes, processed.temp.ECs.solved.successes];
                processed.fum.conc.fRt.MC = [processed.fum.conc.fRt.MC, processed.temp.ECs.solved.MC];
            end
        end
    end
end

myplots.fum.fRt.inputs = {'T0', 'totalH2Omass', 'H2Opartition', 'Xppd', 'fRGinMantle', ...
                          'C_U235', 'C_U238', 'C_Th232','C_K40', 'E', 'eta_pd', 'fR0', 'MagniCScale', 'fum'};
myplots.phiRum.fRt.inputs = {'T0', 'totalH2Omass', 'H2Opartition', 'Xppd', 'fRGinMantle', ...
                             'C_U235', 'C_U238', 'C_Th232','C_K40', 'E', 'eta_pd', 'fR0', 'MagniCScale', 'phiRum'};
modelVars = {'fum', 'phiRum'};

% gathering inputs
for i = 1:length({'fum', 'phiRum'})
    modelVar = modelVars{i};
    myplots.(modelVar).fug.inputParamsMat = NaN(sum(processed.(modelVar).fug.fRt.successes), length(myplots.(modelVar).fRt.inputs)+1);
    myplots.(modelVar).fug.inputParamsMat(:, 1) = 1;
    for p = 1:length(myplots.(modelVar).fRt.inputs)
        myplots.(modelVar).fug.inputData = [processed.(modelVar).fug.fRt.MC.(myplots.(modelVar).fRt.inputs{p})];
        if ismember(myplots.(modelVar).fRt.inputs{p}, {'C_U235', 'C_U238', 'C_Th232','C_K40'})
            myplots.(modelVar).fug.inputData = myplots.(modelVar).fug.inputData ./ ([processed.(modelVar).fug.fRt.MC.fRGinMantle]./100);
        end
        myplots.(modelVar).fug.inputParamsMat(:, p+1) = myplots.(modelVar).fug.inputData;
    end
    myplots.(modelVar).conc.inputParamsMat = NaN(sum(processed.(modelVar).conc.fRt.successes), length(myplots.(modelVar).fRt.inputs)+1);
    myplots.(modelVar).conc.inputParamsMat(:, 1) = 0;
    for p = 1:length(myplots.(modelVar).fRt.inputs)
        myplots.(modelVar).conc.inputData = [processed.(modelVar).conc.fRt.MC.(myplots.(modelVar).fRt.inputs{p})];
        if ismember(myplots.(modelVar).fRt.inputs{p}, {'C_U235', 'C_U238', 'C_Th232','C_K40'})
            myplots.(modelVar).conc.inputData = myplots.(modelVar).conc.inputData ./ ([processed.(modelVar).conc.fRt.MC.fRGinMantle]./100);
        end
        myplots.(modelVar).conc.inputParamsMat(:, p+1) = myplots.(modelVar).conc.inputData;
    end
    myplots.(modelVar).inputParamsMat = cat(1, myplots.(modelVar).fug.inputParamsMat, myplots.(modelVar).conc.inputParamsMat);
end

% let's try gridded heatmap!
myplots.heatmap.bins = 10;
myplots.heatmap.n_inputs = 14;
myplots.heatmap.lower = Inf(myplots.heatmap.n_inputs, 1);
myplots.heatmap.upper = -Inf(myplots.heatmap.n_inputs, 1);
myplots.heatmap.fum.bounds = [];
myplots.heatmap.phiRum.bounds = [];
myplots.heatmap.largecounts = -Inf;
for q = 1:2
    modelVar = modelVars{q};
    for i = 1:myplots.heatmap.n_inputs
        myplots.heatmap.conc_inputData = [processed.(modelVar).conc.fRt.MC.(myplots.(modelVar).fRt.inputs{i})];
        if ismember(myplots.(modelVar).fRt.inputs{i}, {'C_U235', 'C_U238', 'C_Th232','C_K40'})
            myplots.heatmap.conc_inputData = myplots.heatmap.conc_inputData ./ ([processed.(modelVar).conc.fRt.MC.fRGinMantle]./100);
        end
        myplots.heatmap.fug_inputData = [processed.(modelVar).fug.fRt.MC.(myplots.(modelVar).fRt.inputs{i})];
        if ismember(myplots.(modelVar).fRt.inputs{i}, {'C_U235', 'C_U238', 'C_Th232','C_K40'})
            myplots.heatmap.fug_inputData = myplots.heatmap.fug_inputData ./ ([processed.(modelVar).fug.fRt.MC.fRGinMantle]./100);
        end
        if strcmp(myplots.(modelVar).fRt.inputs{i}, 'eta_pd')
            myplots.heatmap.conc_inputData = myplots.heatmap.conc_inputData / 1e20;
            myplots.heatmap.fug_inputData = myplots.heatmap.fug_inputData / 1e20;
        end
        myplots.heatmap.inputData = [myplots.heatmap.conc_inputData, myplots.heatmap.fug_inputData];
        if length(myplots.heatmap.inputData) > myplots.heatmap.largecounts
            myplots.heatmap.largecounts = length(myplots.heatmap.inputData);
        end
        if min(myplots.heatmap.inputData) < myplots.heatmap.lower(i)
            myplots.heatmap.lower(i) = min(myplots.heatmap.inputData);
        end
        if max(myplots.heatmap.inputData) > myplots.heatmap.upper(i)
            myplots.heatmap.upper(i) = max(myplots.heatmap.inputData);
        end
        if strcmp(myplots.(modelVar).fRt.inputs{i}, 'fum')
            myplots.heatmap.fum.bounds = [min(myplots.heatmap.inputData), max(myplots.heatmap.inputData)];
        end
        if strcmp(myplots.(modelVar).fRt.inputs{i}, 'phiRum')
            myplots.heatmap.phiRum.bounds = [min(myplots.heatmap.inputData), max(myplots.heatmap.inputData)];
        end
    end
end
cmap = colormap('parula'); close;
myCmap = [0.25 0.25 0.25; cmap];
myplots.heatmap.maxprob = -Inf;
for q = 1:2
    figure;
    modelVar = modelVars{q};
    myplots.heatmap.grid = NaN(myplots.heatmap.bins, myplots.heatmap.n_inputs);
    for i = 1:myplots.heatmap.n_inputs
        myplots.heatmap.conc_inputData = [processed.(modelVar).conc.fRt.MC.(myplots.(modelVar).fRt.inputs{i})];
        if ismember(myplots.(modelVar).fRt.inputs{i}, {'C_U235', 'C_U238', 'C_Th232','C_K40'})
            myplots.heatmap.conc_inputData = myplots.heatmap.conc_inputData ./ ([processed.(modelVar).conc.fRt.MC.fRGinMantle]./100);
        end
        myplots.heatmap.fug_inputData = [processed.(modelVar).fug.fRt.MC.(myplots.(modelVar).fRt.inputs{i})];
        if ismember(myplots.(modelVar).fRt.inputs{i}, {'C_U235', 'C_U238', 'C_Th232','C_K40'})
            myplots.heatmap.fug_inputData = myplots.heatmap.fug_inputData ./ ([processed.(modelVar).fug.fRt.MC.fRGinMantle]./100);
        end
        if strcmp(myplots.(modelVar).fRt.inputs{i}, 'eta_pd')
            myplots.heatmap.conc_inputData = myplots.heatmap.conc_inputData / 1e20;
            myplots.heatmap.fug_inputData = myplots.heatmap.fug_inputData / 1e20;
        end
        myplots.heatmap.inputData = [myplots.heatmap.conc_inputData, myplots.heatmap.fug_inputData];
        if ismember(myplots.(modelVar).fRt.inputs{i}, {'fum', 'phiRum'})
            myplots.heatmap.edges = linspace(myplots.heatmap.(modelVar).bounds(1), myplots.heatmap.(modelVar).bounds(2), myplots.heatmap.bins+1);
        else
            myplots.heatmap.edges = linspace(myplots.heatmap.lower(i, 1), myplots.heatmap.upper(i, 1), myplots.heatmap.bins+1);
        end
        myplots.heatmap.counts = histcounts(myplots.heatmap.inputData, myplots.heatmap.edges);
        myplots.heatmap.grid(:, i) = 100 * myplots.heatmap.counts / length(myplots.heatmap.inputData);
    end
    %disp(sum(myplots.heatmap.grid, 'omitnan'))
    imagesc(myplots.heatmap.grid)
    colormap(myCmap);
    cb = colorbar();
    set(gca, 'Clim', [0.03 100]);
    set(gca, 'ColorScale', 'log')
    ylabel(cb, 'Fraction of Successes [%]', 'Rotation', 270)
    %cb.Label.String = 'Counts [-]';
    set(gca, 'YDir', 'normal');
    % two x-axes
    ax1 = gca;
    ax1.YTick = [];
    ax1.XAxisLocation = 'bottom';
    ax1_pos = ax1.Position;
    ax2 = axes('Position', ax1_pos);
    ax2.Color = 'none';
    ax2.XAxisLocation = 'top';
    ax2.YAxisLocation = 'right';
    ax2.YTick = [];
    ax2.XLim = ax1.XLim;
    ax2.YLim = ax1.YLim;
    % fill the upper x-axis
    xticks(ax2, 1:myplots.heatmap.n_inputs);
    xticklabels(ax2, arrayfun(@(x) sprintf('%.2f', x), myplots.heatmap.upper, 'UniformOutput', false));
    xtickangle(ax2, 0)
    % fill the lower x-axis
    xticks(ax1, 1:myplots.heatmap.n_inputs);
    xticklabels(ax1, arrayfun(@(x) sprintf('%.2f', x), myplots.heatmap.lower, 'UniformOutput', false));
    xtickangle(ax1, 0)
    % post-processing
    linkaxes([ax1, ax2], 'xy');
end

%% SECTIONBREAK
function SectionBreak()
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
end