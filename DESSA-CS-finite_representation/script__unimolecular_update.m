% IF REACTION IS UNIMOLECULAR
% This is v2. PriorityQ stores absolute times of next reactions rather than
% wait time. 


unimolecular_count = unimolecular_count+1;

% % What are the products? 
%   (For a dimer, this is trivial. Not so for larger oligomers).
%   For now, choose a molecule from the assembly at random to
%   dissociate.
assemID_to_break = max(a1,a2);
avail_subunitIDs = find(subunitIDsAssemblyIDs == assemID_to_break);
assemblySize = length(avail_subunitIDs);
inds = randperm(assemblySize);
subunitID_p1 = avail_subunitIDs(inds(1)); % a randomly chosen subunitID for product 1

% Get the diffusion constant of the reactant so we can sample how far it
% has diffused before undergoing unimolecular reaction.
j = find(currentAssemblyIDs == assemID_to_break);
D_reactant = Dmonomer; % as in Chew et al.

% % -----------------------------------------------------------------------
% % -----------------------------------------------------------------------

% We need to assign a new location to the new assembly (consisting
% of subunitID_p1), and to the initial assembly.        


% % 1) TRANSLATE THE REACTANT COORDINATES GLOBALY, BASED ON HOW FAR 
% IT DIFFUSED BEFORE UNDERGOING THE UNIMOLECULAR REACTION
% % 2) ADD THE LOCAL POSITIONS CHANGES S.T. PRODUCTS ARE "IN CONTACT"
% (RELATIVE TO THE THE UNIMOLECULAR REACTION LOCATION).


% % KEEP SAMPLING NEW POSITIONS UNTIL THERE'S NO OVERLAP WITH ANY OTHER
% % PARTICLES.
overlap = 1;
while overlap
    % % 1. NOW PERFORM THE GLOBAL TRANSLATION ---------------------------------
    r_global = sample_POS_ONLY_update( t_elapsed_a,D_reactant,n_sigma );
    azimuth = (2*pi)*rand(1,1); %unif sampling on [0,2pi]
    elevation = (pi)*rand(1,1); %unif sampling on [0,pi] 
    % we want sph2cart to produce both positive and negative z values (at
    % random)
    if rand(1,1)>0.5; r_global = -r_global;end 
    [gLOBAL_TRANSLATION_x,gLOBAL_TRANSLATION_y,gLOBAL_TRANSLATION_z] = sph2cart(azimuth,elevation,r_global);

    % For use in applying reflective BC if necessary (unlikely to be necessary)
    xyz_Global_A = allAssembliesXYZ(:,currentAssemblyIDs == assemID_to_break);

    xyz_Global = allAssembliesXYZ(:,currentAssemblyIDs == assemID_to_break) + ...
        [gLOBAL_TRANSLATION_x,gLOBAL_TRANSLATION_y,gLOBAL_TRANSLATION_z]';

    % % CHECK THAT NO PARTICLE OVERLAP OCCURS.
    % % Overlap in the context of DESSA-CS means, for some particle p:
    % 1. distance(product,p) < R_contact 
    %    -- AND --
    % 2. t_offset(p) is "small"

    % Compute list of distances to all other molelcules
    closeDistVector = real( (bsxfun(@plus,dot(allAssembliesXYZ,allAssembliesXYZ,1)',...
        dot(xyz_Global,xyz_Global,1))-...
        2*(allAssembliesXYZ'*xyz_Global) ).^(1/2) )';
    close_ind = find( closeDistVector < sqrt(contact_sphere_radiusSqd) );
    if ~isempty(close_ind)
        smallMeanSquareDistance = contact_sphere_radiusSqd;
        smallTime = smallMeanSquareDistance/(6*Dmonomer); 
        if ~any(t_absolute - assemblies_lastUpdate_absTime(close_ind) < smallTime)
            overlap = 0;
        else
            % Sample another position.
        end
    else
        overlap = 0;
    end


end




% Apply Reflective Boundary Conditions to ensure xyz_Global is in
% simulation box
if any(abs(xyz_Global) > containerLength/2)
    %xyz_Global = testing_reflective_boundary_conditions(xyz_Global,containerLength);
    xyz_Global = testing_reflective_boundary_conditions(xyz_Global,xyz_Global_A,containerLength);
end


% % 2. --------------------------------------------------------------------
% Each product will be placed a distance equal 
% to the contact radius away from the unimolecular reaction location. 
% (i.e. the products are placed essentially in contact with each 
% other after dissociating). eGFRD places products in contact as 
% well (see their Supplement section S1.2.1 "Single domains".)

xyz_Global_original = xyz_Global; 

radialDist_p1 = sqrt(contact_sphere_radiusSqd) + contact_sphere_radiusSqd/100;

% we want sph2cart to produce both positive and negative z values (at
% random)
randnum = rand(1,1);
if randnum>0.5; radialDist_p1 = -radialDist_p1;end 

azimuth = (2*pi)*rand(1,1); %unif sampling on [0,2pi]
elevation = (pi)*rand(1,1); %unif sampling on [0,pi]

[x_prime,y_prime,z_prime] = sph2cart(azimuth,elevation,radialDist_p1);
cols_inf = find(~(currentAssemblyIDs<inf));
xyz_Global_1 = xyz_Global_original + [x_prime,y_prime,z_prime]'; % Position of dissociated molecule 

% Apply Reflective Boundary Conditions to ensure xyz_Global is in
% simulation box
if any(abs(xyz_Global_1) > containerLength/2)
   
    xyz_Global = xyz_Global_1; % testing_periodic_boundary_conditions acts on 'xyz_Global'
    
    xyz_Global = testing_reflective_boundary_conditions(xyz_Global,xyz_Global_A,containerLength);
    xyz_Global_1 = xyz_Global;
end
allAssembliesXYZ(:,cols_inf(1)) = xyz_Global_1;

% -------------------------------------------------------------------------
% For Michaelis-Menten model, all assemblies have D ~ 1
% -------------------------------------------------------------------------

azimuth = (2*pi)*rand(1,1); %unif sampling on [0,2pi]
elevation = (pi)*rand(1,1); %unif sampling on [0,pi]
% we want sph2cart to produce both positive and negative z values (at
% random)
radialDist_p1 = -radialDist_p1;
[x_prime,y_prime,z_prime] = sph2cart(azimuth,elevation,radialDist_p1);
xyz_Global_2 = xyz_Global_original + [x_prime,y_prime,z_prime]'; % Position of dissociated molecule 

% Apply Reflective Boundary Conditions to ensure xyz_Global is in
% simulation box
if any(abs(xyz_Global_2) > containerLength/2)
  
    xyz_Global = xyz_Global_2; % testing_periodic_boundary_conditions acts on 'xyz_Global'
    
    xyz_Global = testing_reflective_boundary_conditions(xyz_Global,xyz_Global_A,containerLength);
    xyz_Global_2 = xyz_Global;
end
allAssembliesXYZ(:,currentAssemblyIDs == assemID_to_break) = xyz_Global_2;

% % -----------------------------------------------------------------------
% % -----------------------------------------------------------------------



% Assign a new assemblyID to the dissociated molecule.
% (The original assembly, now one subunit less, will keep original
% assemblyID.)
% UPDATE: ORIGINAL ASSEMBLY MUST BE GIVEN A NEW assemblyID ALSO SO THAT
% EVENTS LEFT IN THE PQ FOR THIS ASSEMBLY ID ARE DISREGUARDED.
assemID_to_break = currentMaxID + 1; % create new assembly ID for original assembly 
currentAssemblyIDs(j) = assemID_to_break; % record new assembly ID
subunitIDsAssemblyIDs(avail_subunitIDs) = assemID_to_break; % update assemblyID belonged to by subunits in original assembly 

currentMaxID = currentMaxID + 2;
subunitIDsAssemblyIDs(subunitID_p1) = currentMaxID; % update assemblyID belonged to by subunit that dissociated  

% We can't have unimolecular reactions unless there was a prev 
% bimolecular reaction leading to at least one 'inf' entry.
assert( ~all(currentAssemblyIDs < inf) ) 

% Record the newly created assemblyID.
currentAssemblyIDs(cols_inf(1)) = currentMaxID;   



% Update absolute time, the time the most recent reaction occured.
t_absolute = t_reaction;


% Keep track of the absolute time these two products were born, i.e. the
% time their states were last updated.
assemblies_lastUpdate_absTime(currentAssemblyIDs==currentMaxID) = t_absolute;
assemblies_lastUpdate_absTime(currentAssemblyIDs==assemID_to_break) = t_absolute;


% --------------------------------------------------
% Begin Update Michaelis Menten benchmark parameters
% --------------------------------------------------
% Which unimolecular reaction happened? (2) ES -> E + S  or (3) ES -> E + P
if reactionRule == 2
    assert(currentAssemblyTYPES(cols_inf(1)) == -1)
    currentAssemblyTYPES(j) = 1; % 1 is E
    currentAssemblyTYPES(cols_inf(1)) = 2; % 2 is S

elseif reactionRule == 3
    
    assert(currentAssemblyTYPES(cols_inf(1)) == -1)
    currentAssemblyTYPES(j) = 1; % 1 is E
    currentAssemblyTYPES(cols_inf(1)) = 4; % 4 is P

end
% ------------------------------------------------
% End Update Michaelis Menten benchmark parameters
% ------------------------------------------------













% % ----------------------------------------------------
% % Update the Priority Queue with potential new events.
% % ----------------------------------------------------
% 1. Bimolecular event involving the two new products. 
%    If none, temporarily store position only updates info for this channel.

% 2. Bimolecular events involving 1 new product and 1 other assembly
%    For each channel, if none, temporarily store position only updates for
%    the channel.

% Let rows cover the two new products.
% Let cols cover all other assemblies.
rows_two_products = [find(currentAssemblyIDs==currentMaxID),...
    find(currentAssemblyIDs == assemID_to_break)];
numProducts = 2;


    
for row_number = 1:numProducts
    row = rows_two_products(row_number);
    
    % Compute list of distances to all other molelcules
    assemblyDistVector = real( (bsxfun(@plus,dot(allAssembliesXYZ,allAssembliesXYZ,1)',...
        dot(allAssembliesXYZ(:,row),allAssembliesXYZ(:,row),1))-...
        2*(allAssembliesXYZ'*allAssembliesXYZ(:,row)) ).^(1/2) )';
    
   
    % If this reactant is not allowed to react with any others, skip. 
    if ~any(assemblyDistVector < inf);continue;end 

   
    
    % Iterate over upperTriangular part of assemblyDistMatrix only.
    % And, only consider columns with:
    % 1. distances less than <maxRxDistance>
    % 2. distances greater than <maxRxDistance> only if partner has a
    %    nonzero t_offset and d < maxRxDistance+R(t_offset)
    %    i.e. partner has been diffusing long enough to reach distance
    %    of <maxRxDistance>.
    columns = 1:length(assemblyDistVector);
    col_indices1 = assemblyDistVector(columns) <= maxRxDistance;
     
    Ds = Dmonomer*ones(1,length(columns));
    Rs = real( n_sigma*(6*Ds*...
        t_absolute - assemblies_lastUpdate_absTime(columns)).^(1/2) );       
    
    col_indices2 = assemblyDistVector(columns) <= maxRxDistance+Rs;    
    col_indices = col_indices1 | col_indices2;  
    
    % Don't over count. I.e. don't consider the same pairs in successive
    % iterations of the for loop.
    col_indices( rows_two_products(1:row_number) ) = false;
    
    % 3. Michaelis-Menten allows 1 bimol Rx: E + S -> ES, i.e. 1 + 2 -> ...
    if currentAssemblyTYPES(row) == 1
        col_indices_ReactionType = currentAssemblyTYPES(columns) == 2;
        col_indices = col_indices & col_indices_ReactionType;
    elseif currentAssemblyTYPES(row) == 2
        col_indices_ReactionType = currentAssemblyTYPES(columns) == 1;
        col_indices = col_indices & col_indices_ReactionType;
    else
        col_indices = col_indices & 0; % currentAssemblyTYPES(row) is not E or S
    end
    
    if sum(col_indices == 0); continue; end
    % ********** BEGIN NEW CODE ************
    cL = (containerLength/2);
    
    Cols = columns(col_indices); % i.e. [ 3 4 11 13 ]
    distList__currRow = assemblyDistVector(Cols); % i.e. [ 1.2 3.01 5.7 1.3 ]
    
    % Calculate distance to the nearest boundary.
    curr_XYZ_colVec = allAssembliesXYZ(:,row);
    assert( sum( abs(curr_XYZ_colVec) > cL*ones(size(curr_XYZ_colVec)) ) ==0) 
    
    partners_XYZ_colVecs =  allAssembliesXYZ(:,Cols); 
    
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------

    three_smallest_distances_to_reflective_boundaries = sort(cL - abs(curr_XYZ_colVec));
    d_reflectiveBoundary = three_smallest_distances_to_reflective_boundaries(1); % Distance to current molecule's nearest reflective boundary
    
    
    % Particles can always exceed boundary. 
    d_partners_nearestBoundary = min( cL - abs(partners_XYZ_colVecs) ) + ...
        exceedBoundaryDist;
    d_nearestBoundary =  d_reflectiveBoundary + exceedBoundaryDist;
    
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------  
    D = Dmonomer; % Michaelis Menten
    
    % Evaluate the max waiting time for current product to reach its
    % nearest relevant boundary.
    maxAllowed_t_wait_curr = (d_nearestBoundary)^2 / (6*D*n_sigma^2);

    
    % Evaluate the max waiting time for partners to reach their
    % nearest relevant boundary.
    T_LAST_partners = assemblies_lastUpdate_absTime(Cols); % time of last update
    t_offset_partners = t_absolute - T_LAST_partners; % past diffusing time (until now)
    D_partners = D*ones(1,length(Cols));
    maxAllowed_t_wait_partners = (d_partners_nearestBoundary - ....
        6.*D_partners.*n_sigma^2.*t_offset_partners).^2 ./ (6*D*n_sigma^2);

    assert( all(maxAllowed_t_wait_partners > 0) );
    % For each pairing (i.e. product & partner), list the max waiting
    % time.
    maxAllowed_t_wait_list = min(maxAllowed_t_wait_curr,maxAllowed_t_wait_partners);
    inds_partner_sets_maxWait = maxAllowed_t_wait_curr > maxAllowed_t_wait_partners;
    
    % Which pairs of diff spheres intersect by the max allowed waiting time?
    % i.e., R_curr(max_t) + R_partner(max_t) >= d ??
    % First, remember that partners may have been diffusing already, so we
    % want to look at diffusion sphere intersections after *their* ELAPSED
    % diffusion times.
    t_offset_list = zeros(1, length(maxAllowed_t_wait_list));
    t_offset_list(inds_partner_sets_maxWait) = t_offset_partners(inds_partner_sets_maxWait);
    
    % When the partner sets the max waiting time, we need to consider its
    % total elapsed diffusion time in determining if there is diffusion 
    % sphere overlap by the max waiting time.
    maxAllowed_t_elapsed_list = maxAllowed_t_wait_list + t_offset_list; 
    
    % Get all pairs with overlapping diffusion spheres by their max waiting time.
    % Also, to be valid, the wait time for a given pair must be long enough
    % bc we integrate the propensity for a duration at most this wait time.
    % We don't want to integrate the propensity from 0 to 1e-9 for example
    % as that's just a waste of resources.
    valid_pair_indices = n_sigma*sqrt(6.*D.*maxAllowed_t_elapsed_list) + ...
        n_sigma*sqrt(6.*D_partners.*maxAllowed_t_elapsed_list) ...
        >= distList__currRow & maxAllowed_t_wait_list > minPropensityDuration;
    
    % SET HARD GLOBAL UPPER LIMIT ON TIME CURRENT ASSEMBLY CAN DIFFUSE.
    % THIS HARD LIMIT IS INDEPENDENT OF NEAREST BOUNDARIES AND INDEPENDENT
    % OF THE DEFINITION OF NEARBY.
    % THIS HARD GLOBAL UPPER LIMIT IS ALSO APPLIED TO POS-ONLY-UPDATES.
    % THIS HARD GLOBAL UPPER LIMIT ALSO APPLIES TO PROPENSITY INTEGRATION
    % IN computeWaitTime_____v2_alt.m
    valid_pair_indices_HARDLIMIT = maxAllowed_t_wait_list <= T_HARDLIMIT_DIFFUSE;
    valid_pair_indices = valid_pair_indices & valid_pair_indices_HARDLIMIT;
    
    % This will be an input to computeWaitTime______v2_alt.
    % Only keep valid entries. 
    maxAllowed_t_wait_list = maxAllowed_t_wait_list(valid_pair_indices);
    
    % Update <Cols> and <distList__currRow> for valid pairs
    Cols = Cols(valid_pair_indices); % i.e. [ 3 4 ]
    distList__currRow = assemblyDistVector(Cols); % i.e. [ 1.2 3.01 ]
    
    numCols = length(Cols);
    lists1 = zeros(numCols, 7);    
 
%     % Determine which, if any, Integrated Propensity curves will be sent to the
%     % parfor workers below (i.e. only those corresponding to the distances
%     % just computed.
%     if numCols > 0
%     Integrated_propensities_parfor = zeros(length(distList__currRow), size(Integrated_propensities,2));
%     for di = 1:length(distList__currRow)
%         pairwisedist = distList__currRow(di);
%         [~,di_index] = min( abs(pairwisedist - distances) );
%         Integrated_propensities_parfor(di,:) = Integrated_propensities(di_index,:);
%     end
%     end

    % ********** END NEW CODE ************
    
    % INNER LOOP.
    if numCols > 0
    %parfor c_ind = 1:numCols 
    for c_ind = 1:numCols
        
    d = distList__currRow(c_ind);    
    col = Cols(c_ind);
    if d < inf % This entry is a valid reaction pair.
        
        T_LASTs = [assemblies_lastUpdate_absTime(row), assemblies_lastUpdate_absTime(col) ];        

        Ds = [Dmonomer Dmonomer]; % to match Chew et. al
        diffSphereRadii = [ n_sigma*(6*Ds(1)* (t_absolute - T_LASTs(1)) )^(1/2),...
            n_sigma*(6*Ds(2)* (t_absolute - T_LASTs(2)) )^(1/2) ];      
        
        
        % Assign particle labels. A -> row (just updated particle).
        % B -> partner
        DaIndex=1; DbIndex=2;
        Da = Ds(DaIndex);
        Db = Ds(DbIndex);
        T_LASTa = T_LASTs(DaIndex); % time reactant A's state last updated.
        T_LASTb = T_LASTs(DbIndex); % time reactant B's state last updated.
        t_offset_a = t_absolute - T_LASTa; % Elapsed time since A update.
        t_offset_b = t_absolute - T_LASTb; % Elapsed time since B update.

        assert(t_offset_a>=0)
        assert(t_offset_b>=0)
        
        % ------------------------
        % BEGIN Compute wait_time
        % ------------------------
        % % Online wait_time computation.
%         maxDiffusionTime = maxAllowed_t_wait_list(c_ind);
%      
%         wait_time = computeWaitTime_finiteParticleSize(Da,Db,t_offset_b,d,...
%             maxDiffusionTime,contact_sphere_radiusSqd,rate_constant);
        
        % % wait_time computation from pre-computed Integrated Propensities.
        % Find the appropriate CDF curve in the pre-computed data (the CDF
        % curve at the distance closest to our currently considered <d>.
%         IntegratedPropensity_curve = Integrated_propensities_parfor(c_ind,:);
       
%         assert(t_offset_a==0)
%         [wait_time] = computeWaitTime_preComputed_IntegratedPropensity(...
%             IntegratedPropensity_curve, times_list, t_offset_b);

        % Find the appropriate CDF curve in the pre-computed data (the CDF
        % curve at the distance closest to our currently considered <d>.
        [~,d_index] = min( abs(d - distances) );        
        
        % This condition added for efficiency. No need to call the wait
        % time function if we know immediately that there isn't enough
        % integrated propensity for the reaction to occur in the allowed
        % duration.
        Pk = log( 1/rand );     
        if Integrated_propensities(d_index,end) < Pk
            % No reaction in allowed max diffusion time
            wait_time = -1;

        else            
            [wait_time] = computeWaitTime_preComputed_IntegratedPropensity(...
                Integrated_propensities(d_index,:), times_list, t_offset_b, Pk);   
        end
        % ------------------------
        % END Compute wait_time
        % ------------------------
     

        % ADD [WAITTIME+t_absolute, slower diffusing ASSEMBLY, faster diffusing ASSEMBLY, updateType, t_elapsed_a,t_elapsed_b] 
        % TO PRIORITY QUEUE. updateType==1 means reaction&position(s). updateType==2 means                        % 
        % position(s) only. 
        % t_elapsed_a, t_elapsed_b are the durations each particle had been
        % diffusing.
        
        % -----------------------------------------------------------------
        if wait_time > 0 % Add Bimolecular Event to PQ
        % -----------------------------------------------------------------    
            % t_elapsed_x is the total duration the diffusion sphere has been
            % expanding, starting from when x's state was last updated, ending
            % when x undergoes its next reaction. 
            % These variables are needed for correctly sampling the reaction
            % LOCATION. Also needed for position only updates.
            t_elapsed_a = wait_time + t_offset_a; 
            t_elapsed_b = wait_time + t_offset_b;           

            % Add the bimolecular reaction to the PQ.

            % The smaller diffusion (A) sphere will be used as the origin
            % when sampling a reaction location.
            reactant_list = [currentAssemblyIDs(row),currentAssemblyIDs(col)];
            list = [wait_time + t_absolute,reactant_list(DaIndex),reactant_list(DbIndex),1,t_elapsed_a,t_elapsed_b,1];
            lists1(c_ind,:) = list;
                
        end
    end
    % -----------------------------------------------------------------    
    end     % ENDS PARFOR
    % -----------------------------------------------------------------
    % Add Bimolecular Events to PQ.
    [logical_inds, ~]=ismember(lists1,zeros(1,7),'rows');
    if sum(logical_inds)>0; lists1(logical_inds,:) = []; end
    PQueue = [PQueue;lists1];   
    
    end
    
    
    % -----------------------------------------------------------------        
    % 2. Add POSITION ONLY UPDATE
    % -----------------------------------------------------------------
    t_wait = T_HARDLIMIT_DIFFUSE;
    list = [t_wait + t_absolute,currentAssemblyIDs(row),-1,2,t_wait,-1,-1];
    PQueue(end+1,:) = list;

    % ---------------------------------------------------------------------
    % 3. Unimolecular events involving the new product.
    % ---------------------------------------------------------------------

    currentAssemID = currentAssemblyIDs( row );
    currentAssemSize = sum(subunitIDsAssemblyIDs == currentAssemID);
    if currentAssemSize >= 2
        
    % Iterate over possible unimolecular reaction types (Michaelis-Menten).
    for rType = 2:3
        % ADD [REACTION TIME, ASSEMBLY to Break, -1, updateType, t_elapsed_a, -1, reactionTpe] TO PRIORITY QUEUE.
        % The slower diffusing assembly will be used as the origin
        % when sampling a reaction location.
        %wait_time = -k_UNI(rType-1)*log(rand(1,1));
        wait_time = (1 / k_UNI(rType-1))*log(1/rand(1,1));
        T_LAST = assemblies_lastUpdate_absTime(row);
        t_offset = t_absolute - T_LAST;
        assert(t_offset==0)

        % The amount of time A has diffused before it undergoes unimolecular reaction
        t_elapsed_a = wait_time;% + t_offset;

        % UPDATE
        % Unimolecular event does not depend on diff sphere encountering
        % boundary
        list_UNI = [wait_time + t_absolute, currentAssemID, -1, 1, t_elapsed_a -1, rType];
        PQueue(end+1,:) = list_UNI;


    end    
    end
    
[logical_inds_PQueue, ~]=ismember(PQueue,zeros(1,7),'rows');
PQueue(logical_inds_PQueue,:) = [];
% assert( sum(logical_inds_PQueue) == 0 )
end



