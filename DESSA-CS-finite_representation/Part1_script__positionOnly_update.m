% POSITION ONLY UPDATE - Parallel Version

%{
When a Pos-Only-Update (POU) appears in the PQ, keep removing next events
as long as they share the same reaction time. Execute these reactions in
parallel. 

Data representing significant parallelization overhead:
allAssembliesXYZ
assemblyDistMatrix
assemblies_lastUpdate_absTime

Next - after parallel event execution - select new events (normal
parallelized loop with outerloop over just-updated assemblies.

%}
positionOnlyUpdate_count = positionOnlyUpdate_count + numParallel;
assert( size(event_Matrix,1) == numParallel )

new_locations = zeros(3,numParallel);
appropriateCols = zeros(1,numParallel);

% parfor pf = 1:numParallel
for pf = 1:numParallel
    
event_vector = event_Matrix(pf,:);
t_reaction = event_vector(1);
a1 = event_vector(2); % This is the assemblyID a1 (smaller diffusion sphere at t_reaction if bimolecular)
a2 = event_vector(3); % This is the assemblyID a2 (larger diffusion sphere at t_reaction if bimolecular)
updateType = event_vector(4); % 1 means reaction&postion, 2 means only position
t_elapsed_a = event_vector(5); % time elapsed since A started diffusing
t_elapsed_b = event_vector(6); % time elapsed since B started diffusing
reactionRule = event_vector(7); % which (of possibly many) reaction rules to apply within bi/uni update script

row_col_logicals_index_a1 = currentAssemblyIDs==a1;
appropriateCols(pf) = find(row_col_logicals_index_a1);t

Da = Dmonomer; % to match Chew et al.

% % KEEP SAMPLING NEW POSITIONS UNTIL THERE'S NO OVERLAP WITH ANY OTHER
% % PARTICLES.
overlap = 1;
while overlap
    rA = sample_POS_ONLY_update( t_elapsed_a,Da,n_sigma ); 

    % % -----------------------------------------------------------------------
    % % PARTICLE A ------------------------------------------------------------
    % % Convert, using random angular variables, to cartesian displacements
    % we want sph2cart to produce both positive and negative z values (at
    % random)
    if rand(1,1)>0.5; rA = -rA;end 

    [xA,yA,zA] = sph2cart(2*pi*rand,pi*rand,rA);


    XYZ_storage_POU1(xyz_counter,:) = [xA,yA,zA];
    xyz_counter = xyz_counter+1;

    % For use in applying reflective BC if necessary.
    xyz_Global_A = allAssembliesXYZ(:,row_col_logicals_index_a1);

    % Add displacement vector to original position vector.
    xyz_Global = allAssembliesXYZ(:,row_col_logicals_index_a1) + [xA,yA,zA]';  

    

    % Apply Reflective Boundary Conditions to ensure xyz_Global is in
    % simulation box
    if any(abs(xyz_Global) > containerLength/2)

        xyz_Global = testing_reflective_boundary_conditions(xyz_Global,xyz_Global_A,containerLength);
    end

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

% Store new positions so we can update <allAssembliesXYZ> outside the 
% parfor loop because MATLAB won't allow it within parfor. 
new_locations(:,pf) = xyz_Global; 

% % END PARTICLE A --------------------------------------------------------
% % -----------------------------------------------------------------------


end
    



% Update absolute time, new value used to update Tk.
t_absolute = event_Matrix(1,1); 

% Store new positions in the appropriate columns of <allAssembliesXYZ>.
allAssembliesXYZ(:,appropriateCols) = new_locations;


% Keep track of the absolute time these states were last updated.
assemblies_lastUpdate_absTime(appropriateCols) = t_absolute; % updated
