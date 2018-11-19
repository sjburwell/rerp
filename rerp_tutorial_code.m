clear 
load('./tutorial_data');
example_trials = 1:10;
% Stimulus-locked ERP trials from Go/Nogo task, +/-2 seconds, 0.5 Hz Highpass filter (nTRIALSxSAMPLES)
data;   figure; plot(data(example_trials,:)'); legend([repmat('Trial',length(example_trials),1) num2str(example_trials')]); 
% event1 and event2 consist of integer samples where 1st and 2nd event occur on a given trial
[event1(example_trials)  event2(example_trials)] 
% Range start-end samples for which to compute rERPs around events (1x2 integer array), here (SR=64), -1 to +1 second
rerpwin
% Clusters for which to separately compute rERPs; here, Go = 3, Nogo = 5
cluster(example_trials)
[rerp1, rerp2, regflag, ntrials] = do_rerp(data, event1, event2, 'rerpwin', rerpwin, 'cluster', cluster);

% Plot (top) standard trial-averaged ERPs and (middle) overlap-corrected rERPs. Note the differences between the top and middle
% panels
respdata = data; 
for ii = 1:size(respdata,1), if ~isnan(event2(ii)), respdata(ii,:) = circshift(respdata(ii,:),[0 event1(ii)-event2(ii)]); end; end
h = figure; set(h, 'Position', [10 10 600 900])
yaxis    = [-10 15];
h = subplot(321); plot(times,mean(    data(cluster==3,event1(1)+rerpwin(1):event1(1)+rerpwin(2))),'-g','LineWidth',2); hold on; grid on; title('Stimulus-locked ERPs');
                  plot(times,mean(    data(cluster==5,event1(1)+rerpwin(1):event1(1)+rerpwin(2))),'-r','LineWidth',2); 
                  legend('Go','Nogo','Location','northwest'); xlim([times(1) times(end)]); ylim(yaxis); ylabel('Potential (uV)'); xlabel('Latency (ms)');
h = subplot(322); plot(times,mean(respdata(cluster==3,event1(1)+rerpwin(1):event1(1)+rerpwin(2))),'-g','LineWidth',2); hold on; grid on; title('Response-locked ERPs');
                  legend('Go','Location','northwest');        xlim([times(1) times(end)]); ylim(yaxis)
h = subplot(323); plot(times,rerp1(1,:),'-g','LineWidth',2); hold on; grid on; title('Overlap-corrected rERPs');
                  plot(times,rerp1(2,:),'-r','LineWidth',2); 
                  legend('Go','Nogo','Location','northwest'); xlim([times(1) times(end)]); ylim(yaxis)
h = subplot(324); plot(times,rerp2(1,:),'-g','LineWidth',2); hold on; grid on; title('Overlap-corrected rERPs');
                  legend('Go','Location','northwest');        xlim([times(1) times(end)]); ylim(yaxis)

% Instead of computing rERPs separately for Go and Nogo trials, we can also view Go trials (the more common of the two) as 
% our "reference" level rERP, and then obtain the Nogo deviation rERP from that Go reference rERP. This can be achieved by 
% dropping the 'cluster' arguement and specifying a 'design' matrix. This approach is akin to "difference waves" commonly 
% used in ERP research (e.g. here, the "Nogo-minus-Go difference wave"), but has the added benefit such that these rERPs
% are OVERLAP-CORRECTED. Although the ERPs, rERPs, and factored rERPs appear similar, I encourage the user to use MATLAB's
% "data cursor" (+ icon) to examine the actual values output by the different methods. 
design =  [ones(size(event1)) double(cluster==5)];
design(example_trials,:)
[ferp1, ferp2, regflag, ntrials] = do_rerp(data, event1, event2, 'rerpwin', rerpwin, 'design', design);
% Plot the (bottom) factored rERPs
h = subplot(325); plot(times,ferp1(1,:),'-g','LineWidth',2); hold on; grid on; title('Factored rERPs');
                  plot(times,ferp1(2,:),'-r','LineWidth',2); 
                  legend('Go (reference)','Nogo (deviation)','Location','northwest'); xlim([times(1) times(end)]); ylim(yaxis)
h = subplot(326); plot(times,ferp2(1,:),'-g','LineWidth',2); hold on; grid on; title('Factored rERPs');
                  legend('Go (reference)','Location','northwest');        xlim([times(1) times(end)]); ylim(yaxis)
saveas(h,'Comparison-ERP-rERP-frERP.pdf')

%% Weighting

% Finally, it is possible to computed "weighted" rERPs, that may enable exploration of other features of the task design. 
% For example, let's say a person wants to re-weight the rERPs accounting for the number of trials going into its computation
% e.g., this could be useful for "equalizing" events that are different in terms of frequency of occurrence. So, the goal
% here would to be ensure that the design matrix coefficients corresponding to each trial type (i.e. Go vs. Nogo) add up to
% the same value (e.g., the number of trials for the Go condition):
design = ones(size(event1)); 
design(cluster==5) = sum(cluster==3).*(1/sum(cluster==5));
design(example_trials,:)
[wrerp1, wrerp2, regflag, ntrials] = do_rerp(data, event1, event2, 'rerpwin', rerpwin, 'cluster', cluster, 'design', design);
h = figure; set(h, 'Position', [10 10 600 900])
h = subplot(321); plot(times, rerp1(1,:),'-g','LineWidth',2); hold on; grid on; title('Overlap-corrected rERPs (unweighted)');
                  plot(times, rerp1(2,:),'-r','LineWidth',2); ylabel('Potential (uV)'); xlabel('Latency (ms)');
                  legend('Go','Nogo','Location','northwest'); xlim([times(1) times(end)]); ylim(yaxis)
h = subplot(322); plot(times, rerp2(1,:),'-g','LineWidth',2); hold on; grid on; title('Overlap-corrected rERPs (unweighted)');
                  legend('Go','Location','northwest');        xlim([times(1) times(end)]); ylim(yaxis)
h = subplot(323); plot(times,wrerp1(1,:),'-g','LineWidth',2); hold on; grid on; title('Weighted rERPs (stim-locked)');
                  plot(times,wrerp1(2,:),'-r','LineWidth',2); 
                  legend('Go (weight=1.00)',['Nogo (weight=' num2str(sum(cluster==3).*(1/sum(cluster==5)),'%.2f') ')'],'Location','northwest'); xlim([times(1) times(end)]); ylim(yaxis)
h = subplot(324); plot(times,wrerp2(1,:),'-g','LineWidth',2); hold on; grid on; title('Weighted rERPs (resp-locked)');
                  legend('Go (weight=1.00)','Location','northwest');        xlim([times(1) times(end)]); ylim(yaxis)

% Or, it's just as easy to code this within the factored rERP design matrix: 
design = [ones(size(event1)) double(cluster==5)];
design = length(event1).*(design./repmat(sum(design),length(event1),1)); 
design(example_trials,:) 
[wferp1, wferp2, regflag, ntrials] = do_rerp(data, event1, event2, 'rerpwin', rerpwin, 'design', design);
h = subplot(325); plot(times,wferp1(1,:),'-g','LineWidth',2); hold on; grid on; title('Weighted-factored rERPs');
                  plot(times,wferp1(2,:),'-r','LineWidth',2); 
                  legend('Go','Nogo','Location','northwest'); xlim([times(1) times(end)]); ylim(yaxis)
h = subplot(326); plot(times,wferp2(1,:),'-g','LineWidth',2); hold on; grid on; title('Weighted-factored rERPs');
                  legend('Go','Location','northwest');        xlim([times(1) times(end)]); ylim(yaxis)
saveas(h,'Comparison-Unweighted-vs-Weighted-rERPs.pdf')

