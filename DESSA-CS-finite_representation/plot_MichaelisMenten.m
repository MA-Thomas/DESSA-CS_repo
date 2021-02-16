% PLOT RESULTS - Michaelis-Menten

indices = timePoints_trials{1,1} > 0;
times = timePoints_trials{1,1}(indices);
E = numE_trials{1,1}(indices);
S = numS_trials{1,1}(indices);
ES = numES_trials{1,1}(indices);
P = numP_trials{1,1}(indices);

figure
plot(times,E)
hold on
plot(times,S)
plot(times,ES)
plot(times,P)
xlabel('Time (s)')
ylabel('Molecules')
legend('E','S','ES','P')
ylim([0 1000])
xlim([0 100])
title(['Diffusion spheres can exceed boundaries by ',num2str(exceedBoundaryDist),' um [runtime(sec)=',num2str(runtime),']'])
ylim([0 numReactants])
xlim([0 100])
hold off
% title(['THARD-preCompute,THARD-initializ,THARD-main: ',num2str(8),', ',num2str(8),', ',num2str(8),...
%     ' Code to restrict partners uncommented'])

%%

% load('allAssembliesXYZ__SAVED.mat','cell_saved')
% for i = 1:length(cell_saved)
%    mat = cell_saved{i};
%    X = mat(1,:);
%    Y = mat(2,:);
%    Z = mat(3,:);
%    
%    figure
%    scatter3(X,Y,Z)
%     
% end

%% PLOT runtime scaling if run_spatial_simulation.m was run with multiple trials at increasing #molecules

figure
for i = 1:trials
    loglog( numReactants(i), [initialization_time_trials(i) + mainLoop_time_trials(i)],'*')
    
    hold on
    % 1000 molecules takes roughly 20s
    %loglog( 1000, 37, '*' )
    
end
yticks([0,60,3600,2*3600,3*3600,4*3600,7*3600, 32*3600])
yticklabels({'0','1 minute','1 hour','2 hours','3 hours','4 hours','7 hours', '35 hours'})
ylim([0 35*3600])
xticks([0,10,10^2,10^3,10^4])
xticklabels({'0','10','10^2','10^3','10^4'})
xlim([0 11^4])
title('V = 90 um^3 -- Simulation Time 100s')
ylabel('RUN TIME (s)')
xlabel('# Molecules')