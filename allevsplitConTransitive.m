function [grps,Qmax,grpscon,Qcon,ctr,maxQ,varargout] = allevsplitConTransitive(A,varargin)
   
% ALLEVSPLITCONTRANSITIVE partition graph eigenvectors of modularity matrix (with consensus)
%   [C,Qmax,Ccon,Qc,N,Q] = ALLEVSPLITCONTRANSITIVE(A) splits the vertices of the graph in adjacency matrix
%   A into multiple groups. 
%   The groups defined in vector C are for the parition with maximum
%   modularity; each element takes an integer value to indicate group membership 
%   of the vertex at that index. Qmax is the corresponding modularity score. 
%   The groups defined in vector Ccon are the result of consensus
%   clustering [see Notes]; Qcon is the modularity score for the consensus
%   partition.
%   N is the number of iterations until consensus was reached. Q is the
%   vector of maxium Q scores for each tested number of groups
%   
%   ...= ALLEVSPLITCONTRANSITIVE(...,DIST,N) sets k-means distance metrics to the strings in 
%   cell array DIST. Set to '' to omit - default is: 'sqEuclidean' (squared Euclidean).
%   Other options include: 'cityblock', 'correlation', 'cosine'. Runs the k-means
%   clustering N times for each specified metric (default is 50); 
%
%   [C,Qmax,Ccon,Qc,N,Q,CLU] = ALLEVSPLITCONTRANSITIVE(...) where CLU is an optional output argument, 
%   returns every single clustering of the adjacency matrix A in the first
%   passs (i.e. before the consensus) - this is useful for further
%   post-processsing.
%
%   Notes: 
%   (0) Adjacency matrix: this can be weighted and directed. When analysing time-series, 
%   this can be the similarity matrix. However, when
%   starting from a similarity  matrix, ensure: 
%       (i) no self-loops - diagonal of A is all zeros; 
%       (ii) it's a similiarity matrix, not a correlation matrix: no
%       negative values
%   Warnings for both of these will be given
%
%   (1) This is a one-step multiple partition method, following up a
%   suggestion in Newman (2006) that all eigenvectors corresponding to positive
%   eigenvalues of the modularity matrix contain information about group
%   structure. The algorithm implemented takes the C such eigenvectors, and
%   uses k-means clustering on those eigenvectors to cluster the nodes into k = C+1 groups. 
%   A value for Q is computed for each k-means clustering (using the defined distance metrics).
%
%   (2) This is repeated for each C in 1:M, where M is number of positive
%   eigenvalues
%
%   (3) Consensus: this attempts to extract a stable set of groups that are robust to repeats
%   of the clustering process. All clusterings with Q>0 across all k-means variants and numbers of groups are
%   pooled. A consensus matrix is computed (Lancichinetti & Fortunato 2012): entry p_ij gives the
%   proportion of clusterings that placed nodes i and j in the same group.
%   The consensus matrix is then run through the one-step multiple parition
%   method (eigenvectors and k-means clustering). A new consensus matrix is
%   created. This is repeated until the distribution of p_ij has become
%   sufficiently bimodal, indicating that the groupings are stable. 
%   
%   (4) Applies suggestion of Reichardt & Bornhaldt (2006) to work with directed networks.
%
%   (5) For the community-detection algorithm, kmeans centres are initialised using the kmeans++ algorithm (Arthur & Vassilvitskii, 2007)
%
%   (6) Detection of bimodal distribution of consensus matrix now sets the
%   k-means initial centres dynamically to the 5th and 95th percentile of
%   the entries on the consensus matrix (not fixed to [0.4 0.9] as per
%   Bruno et al results).
%
%   References: 
%   (1) Newman, M. E. J. (2006) "Finding community structure in
%   networks using the eigenvectors of matrices". Phys Rev E, 74, 036104.
%
%   (2) Reichardt & Bornhaldt (2006) "Statistical mechanics of community detection".
%   Phys Rev E. 74, 016110
%   
%   (3) Lancichinetti, A. & Fortunato, S. (2012) Consensus clustering in complex networks.
%   Scientific Reports, 2, 336
%   
%   (4) Arthur, D. & Vassilvitskii, S. (2007) k-means++: the advantages of careful seeding. 
%   SODA '07: Proceedings of the eighteenth annual ACM-SIAM symposium on Discrete algorithms, Society for Industrial and Applied Mathematics, 1027-1035
%
%   Mark Humphries 19/3/2015

strDist = {'sqEuclidean'}; % 'cityblock','correlation','cosine'};
nreps = 50; % of each distance metric

%% check if the passed matrix is a graph: catch common errors when passing a similarity matrix

% (1) no self-loops allowed
if ~all(diag(A)==0) 
    warning('Results likely unreliable: adjacency matrix has self-loops. Set diagonal to zero if no self-loops are needed.')
end

% (2) no negative links
x = sum(sum(A < 0));
if x > 0
    warning('Results likely unreliable: adjacency matrix has negative values')
end

%% set up options
blnSave = 0; % internal flag for setting saving of data
if nargin >= 2
    if ~isempty(varargin{1}) 
        strDist = varargin{1}; 
        if ~iscell(strDist)
            error('Distance strings should be passed as a cell array')
        end
    end
    
    if ~isempty(varargin{2}) nreps = varargin{2}; end
end

% set up saving of each iteration
if blnSave  
    fname = ['Consensus_Iterations_' datestr(now,30)];
    save(fname,'strDist','nreps');  % save initial data to allow -append to work below
end

%% internal parameters
% total number of clusterings for a given number of groups
Treps = nreps .* numel(strDist);

[nIDs,c] = size(A);

blnConverged = 0;
ctr = 1;
Gstar = ones(nIDs,1); % current clustering
mQ = []; stdQ = [];

Aorig = A;  % save a copy of the original matrix for Q calculations

%% main clustering loop
while ~blnConverged
    [allgrps,Vpos,B,Q] = dosplit(A,nreps,Treps,strDist);
   
    % get max Q from initial set of clusterings
    if ctr == 1
        if isempty(allgrps) || all(Q(:) <= 0)
            % then no groups detected; return empty
            grps = zeros(nIDs,1);
            Qmax = 0; 
            grpscon = zeros(nIDs,1);
            Qcon = 0;
            ctr = 0; 
            maxQ = []; 
            blnConverged = 1;  % no answer
            return
        else
            maxQ = max(Q);
            Qmax = max(maxQ);
            numgrps = find(maxQ == Qmax);
            rpt = find(Q(:,numgrps) == max(Q(:,numgrps)));
            cl = (numgrps-1)*Treps + rpt(1);  % if more than 1, choose first arbitrarily...
            grps = allgrps(:,cl);  
            varargout{1} = allgrps;
        end
    end
    
     
    % compute consensus BETWEEN numbers of groups
    Allowed = (Q(:) > 0);       % only take consensus using Q-positive clusterings...

    A = consensus(allgrps(:,Allowed)); % get consensus matrix 
         

    % threshold: Lanchinetti & Fortunato use a threshold for discarding
    % weakly paired nodes; but this threshold depends strongly on the
    % algorithm being used. So instead we auto-threshold based on
    % convergence to a bimodal distribution for the consensus matrix
    % weights
    
    % auto-threshold
    allWs = triu(A); allWs = allWs(allWs > 0);
    
    % convergence of weights?
    if all(allWs == 1)
        theta = -inf;   % converged to all 1s, so no need to check for bimodality
    else
        try
            % initialised as per Bruno et al paper:    
%             idx = kmeans(allWs,2,'Start',[0.4;0.9]); % use k-means to separate into mode groups
            
            % initialised dynamically according to spread of data
            s =  prctile(allWs,[5,95]);
            idx = kmeans(allWs,2,'Start',s'); % use k-means to separate into mode groups
           
        catch
            keyboard
        end
        
        if blnSave
            eval(['A', num2str(ctr),' = A;']);  % store consensus matrix
            eval(['Groups', num2str(ctr),' = idx;']); 
            eval(['AllGroups', num2str(ctr),' = allgrps;']); 
            eval(['AllQ', num2str(ctr),' = Q;']);       
            save(fname,['A', num2str(ctr)],['Groups', num2str(ctr)],['AllGroups', num2str(ctr)],['AllQ', num2str(ctr)],'ctr','-append');
        end

        m1 = mean(allWs(idx==1)); mx1 = max(allWs(idx==1)); mn1 =  min(allWs(idx==1));
        m2 = mean(allWs(idx==2)); mx2 = max(allWs(idx==2)); mn2 =  min(allWs(idx==2));        
        if m1 > m2  % work out which is lower and which is higher distribution
            theta = mn1;
        elseif m2 > m1 
            theta = mn2;
        end
    end
    
    A_upper = A;
    A_upper(A < theta) = 0;  % remove all of lower mode links

    % find groups deterministically...
    grpscon = [0 0]; grpctr = 1; blnTrans = 1;
    for iN = 1:nIDs
        if ~any(iN == grpscon(:,1))  % if not already sorted
            thisgrp = [iN; find(A_upper(iN,:) > 0)']; % all nodes connected to this one are in the same group
            % if any of this group are already in the storage matrix,
            % then this is not transitive...
            for iT = 1:numel(thisgrp)
                if any(thisgrp(iT) == grpscon(:,1))
                    blnTrans = 0; break
                else
                    blnTrans = 1;
                end
            end
            if blnTrans == 0
                break
            else
                grpscon = [grpscon; thisgrp ones(numel(thisgrp),1)*grpctr]; % save this...
                grpctr = grpctr + 1;
            end
        end
    end

    if blnTrans
        grpscon(1,:) = []; % remove padding zeros
        [Nsrt,ixsrt] = sort(grpscon(:,1));
        grpscon = grpscon(ixsrt,2);
        ngrps = numel(unique(grpscon));

        % Q of converged answer: given original data!!
        P = expectedA(Aorig);
        B = Aorig - P;
        m = sum(sum(Aorig))/2; 

        S = zeros(nIDs,ngrps);
        for loop = 1:ngrps
            S(:,loop) = grpscon == loop;
        end
        % compute modularity
        Qcon = trace(S' * B * S) / (2*m);
        blnConverged = 1;
    else
        % carry on
        ctr = ctr+1;
        if ctr > 50
            warning('Did not converge in 50 iterations')
            try
                % return what we do have...    
                % find groups deterministically...
                grpscon = [0 0]; grpctr = 1;
                for iN = 1:nIDs
                    if ~any(iN == grpscon(:,1))  % if not already sorted
                        thisgrp = [iN; find(A(iN,:) > 0)']; % all nodes connected to this one are in the same group
                        grpscon = [grpscon; thisgrp ones(numel(thisgrp),1)*grpctr]; % save this...
                        grpctr = grpctr + 1;
                    end
                end

                grpscon(1,:) = []; % remove padding zeros
                [Nsrt,ixsrt] = sort(grpscon(:,1));
                grpscon = grpscon(ixsrt,2);
                ngrps = numel(unique(grpscon));

                % Q of converged answer: given original data!!
                P = expectedA(Aorig);
                B = Aorig - P;
                m = sum(sum(Aorig))/2;

                S = zeros(nIDs,ngrps);
                for loop = 1:ngrps
                    S(:,loop) = grpscon == loop;
                end
                % compute modularity
                Qcon = trace(S' * B * S) / (2*m);

                blnConverged = 1;
            catch
                keyboard
            end

        end
    end
end
      


function [allgrps,Vpos,B,Q] = dosplit(A,nreps,Treps,strDist)
    % this function implements the main eigenspectra-decomposition
    % algorithm
    
    [n c] = size(A);
    Abar = (A + A') / 2;
    m = sum(sum(A))/2; % edges
    
    % generate expected graph
    P = expectedA(A);
    Pbar = (P + P') / 2;

    % create modularity matrix
    B = Abar - Pbar;

    % find eigenvalues (D) and eigenvectors (V) of B
    try
        [V,D] = eig(B);         % assumes symmetry
    catch
        keyboard
    end

    eg = diag(D);
    if ~isreal(eg)
        warning('Complex eigenvalues have occurred')
        ns = sum(V > 0);
    end
    
    % eigenvectors corresponding to positive eigenvalues
    egpos = find(eg > 1e-3);    % allow for rounding error
    ngrps = numel(egpos) + 1;   % upper bound is one more group than vectors
    Vpos = V(:,egpos);


    ndivs = numel(2:ngrps);
    Q = zeros(Treps,ndivs);
    allgrps = zeros(n,Treps*ndivs);
    C = cell(Treps,1);

    % run over all possible numbers of groups to max, and record groups with
    % max Q....
    for j = 1:ndivs
        Dctr = 1;
        thisEigVector = Vpos;
        
        for rep = 1:Treps
            ixNow = (j-1)*Treps + rep;
            cpos = kmeansplus(thisEigVector ,j+1); % initialise centers
            % keyboard
            try            
                % j+1 is the number of groups this time...                
                [allgrps(:,ixNow),C{rep},sumD{rep},allD{rep}] = kmeans(thisEigVector,j+1,'Distance',strDist{Dctr},'Start',cpos);
                
                
                % construct S matrix of group membership: each column is a group
                % See: Newman (2006) Eq 31            
                S = zeros(n,j+1);

                for loop = 1:j+1
                    S(:,loop) = (allgrps(:,ixNow) == loop);
                end
                % compute modularity
                Q(rep,j) = trace(S' * B * S) / (2*m);

            catch
                % if kmeans throws a wobbly, set to "no groups"...
                allgrps(:,ixNow) = zeros(n,1);
                Q(rep,j) = -1;
                % warning('kmeans wobbly')
            end
            
            % use next distance metric if reached N repetitions of this one...
            if rep ./ nreps == Dctr Dctr = Dctr + 1; end

        end % end k-means loop
        Qmax(j) = max(Q(:,j));

        % safety-catch: if consecutive attempts at groupings find no groups 
        % after multiple k-means, then extremely unlikely that adding
        % increasing k further (recruiting smaller +ve eigenvalues) will
        % find anything
        if j > 2 
            if Qmax(j) <= 0 & Qmax(j-1) <= 0
                break
            end
        end

    end % end groups loop

function Sc = consensus(Grps)
    % Pass: NxC matrix of C clusterings of N objects; and initial grouping  
    [nIDs nreps] = size(Grps);
    Sc = zeros(nIDs);
    
    % pair-wise similarity matrix
    for nr = 1:nIDs
        for nC = nr:nIDs
            if nr ~= nC
                Gi = Grps(nr,:); Gj = Grps(nC,:);
                Sc(nr,nC) = sum(Gi == Gj) / nreps;
                Sc(nC,nr) = Sc(nr,nC);
            end
        end
    end
     

function m = kmeansplus(X,k)

[n,d] = size(X); 
m = zeros(k,d);
v = inf(n,1); % current minimum distance of data-point to closest centre
m(1,:) = X(ceil(n*rand),:); % put is co-ordinates into 

for i = 2:k
    % find distance from all data-points to closest center
    D = sqrt(sum((X - repmat(m(i-1,:),n,1)).^2,2));  % distance of all data-points from last picked center
    v = min(v,D);  % update minimum distance to closest centre
    p = cumsum(v/sum(v));  % ECDF
    m(i,:) = X(find(rand < p,1),:);  % pick next point proportional to p
end




