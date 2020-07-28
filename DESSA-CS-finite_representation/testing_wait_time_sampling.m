% I want to investigate the error involved in precomputing integrated
% propensities (finite particle / Green's function model) at discrete r0
% values rather than computing real time integrated propensities at the
% actual r0 values.
% That is, using the precomputed values at the closest r0 available as a
% proxy for the actual value at the actual r0.


Da = 1;
Db = 1;
t_offset_b = 3%.3;
maxDiffusionTime = 10;
R = 0.01;
c = 10;


numSeparations = 1000;
numtrials = 300;
initial_separation_list = linspace(2e-2,16,numSeparations);
wait_time_matrix = zeros(numtrials,numSeparations);
parfor i = 1:numtrials
    
    for d_ind = 1:numSeparations
        
        d = initial_separation_list(d_ind);
        wait_time_matrix(i,d_ind) = computeWaitTime_finiteParticleSize(Da,Db,t_offset_b,d,...
            maxDiffusionTime,R,c);
    end    
    savefunction(i);


end

mean_wait_times = mean(wait_time_matrix);
std_wait_times = std(wait_time_matrix);
errorbar(initial_separation_list,mean_wait_times,std_wait_times,'-s','MarkerSize',1,'MarkerEdgeColor','red','MarkerFaceColor','red');
xlabel('r0');
ylabel('<t-wait>')

save(['workspace__t_offset_',num2str(t_offset_b),'.mat'])
% The blue region corresponds to 2 standard deviations around the mean in
% red.
 %%% -----------------------------------
%%
% % 
% % t_offset_b = 3;
% % maxDiffusionTime = 10;
% % R = 0.01;
% % c = 10;
% % 
% % 
% % numSeparations = 1000;
% % numtrials = 300;
% % initial_separation_list = linspace(2e-2,16,numSeparations);
% % wait_time_matrix2 = zeros(numtrials,numSeparations);
% % for i = 1:numtrials
% %     
% %     for d_ind = 1:numSeparations
% %         
% %         d = initial_separation_list(d_ind);
% %         wait_time_matrix2(i,d_ind) = computeWaitTime_finiteParticleSize(Da,Db,t_offset_b,d,...
% %             maxDiffusionTime,R,c);
% %     end
% %     if mod(i,10)==1; disp(['Trial ',num2str(i)]); end
% % 
% % end
% % 
% % mean_wait_times2 = mean(wait_time_matrix2);
% % std_wait_times2 = std(wait_time_matrix2);
% % figure
% % errorbar(initial_separation_list,mean_wait_times2,std_wait_times2,'-s','MarkerSize',1,'MarkerEdgeColor','red','MarkerFaceColor','red');
% % xlabel('r0');
% % ylabel('<t-wait>')
% % 
% % %%
% % t_offset_b = 3e-4;
% % maxDiffusionTime = 10;
% % R = 0.01;
% % c = 10;
% % 
% % 
% % numSeparations = 1600;
% % numtrials = 300;
% % initial_separation_list = linspace(2e-2,16,numSeparations);
% % wait_time_matrix3 = zeros(1,numSeparations);
% % for i = 1:numtrials
% %     
% %     for d_ind = 1:numSeparations
% %         
% %         d = initial_separation_list(d_ind);
% %         wait_time_matrix3(i,d_ind) = computeWaitTime_finiteParticleSize(Da,Db,t_offset_b,d,...
% %             maxDiffusionTime,R,c);
% %     end
% %     if mod(i,10)==1; disp(['Trial ',num2str(i)]); end
% % 
% % end
% % 
% % mean_wait_times3 = mean(wait_time_matrix3);
% % std_wait_times3 = std(wait_time_matrix3);
% % figure
% % errorbar(initial_separation_list,mean_wait_times3,std_wait_times3,'-s','MarkerSize',1,'MarkerEdgeColor','red','MarkerFaceColor','red');
% % xlabel('r0');
% % ylabel('<t-wait>')