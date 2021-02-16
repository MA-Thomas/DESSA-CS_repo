function [CDFs,Integrated_CDFs] = compute_Integrated_prePropensity_curves(...
    contact_sphere_radiusSqd,distances,var_list,times_list)

% Evaluate CDF = Pr(Q < R^2) at each gridPoint (i.e. at each element
% (dist,var) )

% CDFs: one element for each variance, at each distance.
% Integrated_CDFs: 1 curve for each of the <numDistances> initial particle
% separations. (same size matrix as CDFs)

numDistances = length(distances);
numVariances = length(var_list);

% Each element will be obtained by examining an approximation to the
% infinite power series formulation (up to 200 terms included in the
% approx)
CDFs = zeros(numDistances,numVariances);
% number of consecutive terms that have to be approx equal to consider the
% infinite series converged.
numConsecutive = 10;
%%
nan_reached = 0;
limit_found = 0;

tic
parfor d_ind = 1:numDistances
    Dist = distances(d_ind);
    
    nan_reached = 0;
    limit_found = 0;
    
for v_ind = 1:numVariances



% Var = gridPoints(i,2);    
% Dist = gridPoints(i,1);   
Var = var_list(v_ind);

[x,y,z] = sph2cart(2*pi*rand,pi*rand,Dist);
mu_colVec = [x,y,z]';


% Output: CDF - Probability sample from N(mu_colVec,E) has magnitude less than contact_sphere_radius
% CDF is Pr(Q(X)<R) and thus 0 <= CDF <= 1.

% Inputs (3 dimensional space):
% vector of variances along each dimension: E = [sigma^2_NEW(1),sigma^2_NEW(2),sigma^2_NEW(3)]
% (Cov matrix: E = sigma^2_NEW * IdentityMatrix)

% contact_sphere_radiusSqd: the distance below which two assemblies are
% considered to be "in contact". (distance squared, that is)


%{

Gaussian molecule A: a ~ N(mu_A, sigma^2_A)
Gaussian molecule B: b ~ N(mu_B, sigma^2_B)

Gaussian for difference vector: 
z = a-b ~ N( [mu_A - mu_B], [sigma^2_A + sigma^2_B] )
  = N(mu_NEW, sigma^2_NEW)

We're interested in the distribution of squared distances, i.e. the 
quadratic form:
 Q(z) = z(1)^2 + z(2)^2 + z(3)^2
 
More specifically, we're interested in the CDF(y) of Q(z), Pr{ Q(z) < y }

Also, the variances are time dependent, AND, there may be time
offsets.
e.g.
sigma^2_NEW(t) = sigma^2_A(t) + sigma^2_B(t+offset)
OR
sigma^2_NEW(t) = sigma^2_A(t+offset) + sigma^2_B(t)

%}


%  Theorem 4.2b.1 Mathai-Provost p.95

% Propensity integration cannot start from exactly t==0 and t_offset==0.
assert(Var ~= 0); % This would lead to inf/NaN results.

E = Var*[1 1 1];
p = 3; % 3d space

y = contact_sphere_radiusSqd; % CDF(y)


%{
% Steps:
0. compute b_vector
1. compute d(k) for k = 1:5
2. compute c(k-1)
 
Follow steps on p.96

A = Identity, eye(3)
[P,D] = eig(CovMatrix)

Since the covariance matrix is diagonal, 
P <- 3x3 Identity
D <- Diagonal matrix of eigenVals. All eigenVals equal to sigma^2_NEW
i.e. THE EIGENVALUES ARE JUST THE VARIANCES ALONG EACH DIMENSION AND IN OUR
ISOMETRIC CASE, ALL VARIANCES ARE EQUAL

%}
% Write all the formulas using E (variances).
% In the final formula, replace E with time dependent variance.
% B = chol( E.*eye(3) ); % Matlab v2018
% B = chol( E(1).*eye(3) ); % Matlab v2016 temp
B = ( 1/sqrt(E(1)) ) * eye(3); % Since E = Var*[1 1 1], we're taking the chol
%                        decomp of a multiple of the identity matrix.

% b_colVec = P' * inv(B) * mu   
% But P == Identity
% b_colVec = B \ mu_colVec; %inv(B) * mu_colVec;
b_colVec = B * mu_colVec; % Since B is a scalar multiple of identity, and inv(I)==I

c0 = exp( sym( -.5*(b_colVec'*b_colVec)) ) * sym( (2*E(1))^(-1/2) * ...
    (2*E(2))^(-1/2) * (2*E(3))^(-1/2) );

numTermsInApprox = 200; % Include up to this many terms in approximation. 
d = zeros(1,numTermsInApprox);
c = zeros(1,numTermsInApprox); % Does not include c0. Starts with c1.


sign_terms = ones(1,numTermsInApprox);
sign_terms(1:2:end) = -1*sign_terms(1:2:end);

y_exponents = (p/2)+1:(p/2)+numTermsInApprox;
y_list = y.^y_exponents;


step = 3;
chunks = [1, 11, 22:step:numTermsInApprox+1]; % i.e. [1 21 22 23 ... 151]

Fs = zeros(1,numTermsInApprox);
k=1;
d(k) = sym( (1/2) * ( ( 1-k*(b_colVec(1))^2 )/( (2*E(1))^k ) + ...
    ( 1-k*(b_colVec(2))^2 )/( (2*E(2))^k ) + ...
    ( 1-k*(b_colVec(3))^2 )/( (2*E(3))^k )  ) );
c(k) = sym( d(k)*c0 );



for mP = 1:length(chunks)-1 % i.e. [1 2 3 4 5]

    % If an element of the series Fs is NaN, no need to compute further
    % elements. We can assume the limit of the series is 0.
    if nan_reached || limit_found      
        break
    end
    
    % % Evaluate c and d for the current chunk.
    for k = chunks(mP) : chunks(mP+1)-1 % i.e. 1:20,21:21,22:22,23:23
        % Compute d vector
        d(k) = sym( (1/2) * ( ( 1-k*(b_colVec(1))^2 )/( (2*E(1))^k ) + ...
            ( 1-k*(b_colVec(2))^2 )/( (2*E(2))^k ) + ...
            ( 1-k*(b_colVec(3))^2 )/( (2*E(3))^k )  ) );

        c_inds = 1:k-1; 
        d_inds = k-1:-1:1;
        c(k) = sym( (1/k)* (c0*d(k) + sum( d(d_inds).*c(c_inds) ) ) );

        % Fill in the next chunk of F evaluations     
        Fs(k) = vpa( (c0*sym(y^(3/2 + 0))/gamma(sym(3/2 + 1))) + ...
            sum( sym( sign_terms(1:k).*c(1:k) ).*...
            sym( y_list(1:k) )./ gamma(sym( y_exponents(1:k)+1) ) ) );

        if k > numConsecutive + 1
           % TODO: Test whether last 10 consecutive terms are roughly 
           %       equal. If so, then this is the limit of the series.
           %       I.e. this is CDF(var(i) | d)
           % Also, must be less than/equal to 1 as its a CDF.

           if all( diff( Fs(k-(numConsecutive+1):k) ) < 1e-5 ) && ...
               all(  abs( Fs(k-(numConsecutive+1):k) ) ) <= 1
               
               CDFs(d_ind,v_ind) = Fs(k);
               limit_found = 1;
               
%                figure
%                plot(Fs(1:k))
%                legend([num2str(d_ind),' ',num2str(v_ind)])
%                title('convergence')
               break
           end
        end
        if isnan(Fs(k))

           nan_reached = 1;
           limit_found = 1;
           CDFs(d_ind,v_ind) = 0;
           
%            figure
%            plot(Fs(1:k))
%            legend([num2str(d_ind),' ',num2str(v_ind)])
%            title('no convergence')
           break
        end
        
    end
    
end

% If an element of the series Fs is NaN, no need to compute further
% elements. We can assume the limit of
if nan_reached || limit_found
    nan_reached = 0; % reset for evaluation of next Fs series
    limit_found = 0;
    continue
end 

end
end



%plot(times_list,CDFs(34,:))

toc
tic
Integrated_CDFs = zeros(numDistances,numVariances);
t = times_list;
% for d_ind = 1:numDistances 
% for v_ind = 1:numVariances
%    
%     if v_ind == 1
%        % Integrated Prop at the first var (/time) point is 0.
%     else
%        Integrated_CDFs(d_ind,v_ind) = Integrated_CDFs(d_ind,v_ind-1) + ...
%            trapz( [t(v_ind-1),t(v_ind)], [CDFs(d_ind,v_ind-1),CDFs(d_ind,v_ind)] );
%     end
% end
% end

parfor d_ind = 1:numDistances    

   Integrated_CDFs(d_ind,:) = cumtrapz(t,CDFs(d_ind,:))  

end
toc


end