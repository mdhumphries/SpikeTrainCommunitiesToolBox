function [Gmax,Gcon,varargout] = consensus_cluster_spike_data_binless(spkdata,Didxs,T,bin,varargin)

% CONSENSUS_CLUSTER_SPIKE_DATA_BINLESS consensus  cluster spike-data using binless-based metrics
% [GMAX,GCON] = CONSENSUS_CLUSTER_SPIKE_DATA_BINLESS(S,I,T,B,OPTS) clusters spike trains into
% groups, based on similarities between pairs of "binless" spike train representations.
% Similarities between all pairs are stored in similarity matrix, which forms
% the basis for the clustering.
%
% Two forms of clusterings are returned: (1) the clustering that maximises
% the benefit function (modularity), in structure GMAX; (2) the consensus clustering created
% from across all clusterings of the data-set, in structure GCON. The
% whole clustering algorithm is run to completion for every time-scale of
% correlation passed in optional argument OPTS.BLpars (see below).
%
% S is a 2-column vector: first column is ID stamps; second column is
% time-stamps. IDs can be either for events (e.g. stimulus presentations
% for single unit data; time-stamps are then relative to stimulus presentation) 
% or neurons (for multi-unit data; time-stamps are then relative to recording period).
%
% I is a column vector of all possible ID stamps - note that they need not all appear in S
% because the ID event may have caused no spikes or been a silent neuron!
% If a neuron is silent, it will not appear in the groupings: its ID stamp
% will be omitted from G.
%
% T is 2-element vector of [start end] for all time-stamps (e.g. if using zero as 
% time-stamp for stimulus presentation: for 2 seconds of data starting at stimulus presentation T = [0 2]; for 2
% seconds starting 0.5 second before stimulus T = [-0.5 1.5]). The function
% will only use spike-times with T(1) <= t <= T(2)
%
% B is the quantisation value for the convolution function - try 1ms (B=0.001)  
%
% ... = CLUSTER_SPIKE_DATA_BINS(...,OPTS) sets all options contained in the OPTS structure - omitting a field will use the default option, indicated
% by the parentheses:
%       OPTS.BLmeth = {('Gaussian') | 'exponential'}
%               Convolve each spike with a Gaussian (Fellous et al, 2004; used in Humphries, 2011); or 
%               forward exponential function (van Rossum, 2001; included for future development of algorithm)           
%       OPTS.BLpars : array of values for the binless metric parameter - each value will be tested,
%               'Gaussian' - standard deviation of Gaussian in seconds (default: 0.01)
%               'exponential' - decay constant of exponential in seconds 
%       OPTS.Dmeth = {'cosine'| {'corrcoef'} | 'corr'} 
%           sets the comparison-between-pairs method to: cosine of angle between vectors (Fellous et al 2004; used in Humphries, 2011); 
%           rectified correlation coefficient (corrcoef or corr) [default]
%       OPTS.modopts ({}) - passes anything in this cell array to the corresponding optional
%           arguments of the clustering functions (ALLEVSPLITCON) - see
%           their help for details. For example, for ALLEVSPLITCON, the first cell will be
%           passed into the first optional argument (a string), 
%           the second cell into the second optional argument (a number) 
%       OPTS.blnS = {(0)|1} - if set to 1, does not do the clustering algorithm, but only computes 
%           the similarity matrices for each value in opts.BLpar: useful for just creating a set of such matrices for further processing.
%           [Default is 0]
%
% Returns: two structures GMAX and GCON each containing fields:
%               grps: the group membership of each node as a two-column vector, first column 
%                        is ID stamp, second column is group membership (integer value) 
%                        - note that there is no guarantee all IDs will be assigned a group. 
%               grpsizes: vector of group sizes
%   `           ngrps: the number of groups
%               Q: modularity score for the grouping
%
% There is one struct for each tested binless time-scale specified.
%
% [...,Sxy,BD,I] = CLUSTER_SPIKE_DATA_BINLESS(...) are optional outputs useful for further 
% processing: 
%       Sxy is a cell array of the similarity matrices; one cell per tested binless time-scale. 
%       BD is a cell array of matrices containing  the convolved spike-train vectors (in spikes per second) - 
%               each column is one spike train; one cell per tested binless time-scale. 
%       I is a cell array of the retained spike-train IDs; one cell per tested binless time-scale. 
% 
% Notes:
% (1) Choice of Gaussian widths: typically, we choose a set of discrete binsizes based on 
%  statistics of spike-train data-set, then follow Kruskal et al (2007)
%  and convert to Gaussian widths by s = w / sqrt(12)
%
% (2) If an ID stamp has no associated spikes then the train with that ID
% is omitted from the analysis: the groupings in G will not contain that ID
% stamp.
%
% (3) This helper function can be easily extended with (a) more choices of
% convolution windows and (b) more choices of pairwise similarity metric
%
% References:
% (1) Humphries, M. D. (2011) Spike-train communities: finding groups of similar
% spike trains J Neurosci 31, 2321-2336
% (2) Humphries, Wood & Gurney (2009) Dopamine-modulated dynamic cell assemblies generated by the GABAergic striatal
% microcircuit, Neural Networks, 22, 1174-1188. 
% (3) van Rossum M. C. (2001) A novel spike distance. Neural Comput, 13, 751-763
% (4) Fellous, J. M.; Tiesinga, P. H.; Thomas, P. J. & Sejnowski, T. J.
%       (2004) Discovering spike patterns in neuronal responses J Neurosci, 24, 2989-3001
% (5) Kruskal, P. B., Stanis, J. J., McNaughton, B. L., Thomas, P. J. (2007). A binless correlation measure reduces the variability
% of memory reactivation estimates. Stat Med 26: 3997–4008.
% 
% Mark Humphries 6/2/2017

[r c] = size(Didxs);
if r==1 Didxs = Didxs'; end   % make column vector

% defaults  
% bin = 0.001;    % 1ms bins
opts.BLmeth = 'Gaussian';   % use Gaussian window around spikes 
opts.BLpars = 0.01;        % std dev is 10 ms
opts.Dmeth = 'corrcoef';     % use correlation coefficient
opts.nlimit = 6;       % minimum number of nodes  
opts.modopts = '';   % set no options for the all eigenvector method...
opts.blnS = 0;       % run clustering algorithm by default

if nargin >= 5
    if isstruct(opts) 
        tempopts = varargin{1}; 
        fnames = fieldnames(tempopts);
        for i = 1:length(fnames)
            opts = setfield(opts,fnames{i},getfield(tempopts,fnames{i}));
        end
    end
end

% storage
Rn = cell(numel(opts.BLpars),1);        % retained neurons...
Gmax = struct('grps',[],'grpsizes',[],'ngrps',[],'Q',[]);
Gcon = struct('grps',[],'grpsizes',[],'ngrps',[],'Q',[]);

spkfcn = cell(numel(opts.BLpars),1);

% bins = T(1)+bin/2:bin:T(2);
bins = T(1):bin:T(2);


%% analyse data
for loop = 1:numel(opts.BLpars)
    
    % set up convolution window if using...
    sig = round(opts.BLpars(loop)/bin);  % 1 SD in counts of bins
    switch opts.BLmeth
        case 'Gaussian'
            x = [-5*sig:1:5*sig]';  % x-axis values of the discretely sampled Gaussian, out to 5xSD
            h = (1/(sqrt(2*pi*sig^2)))*exp(-((x.^2*(1/(2*sig^2))))); % y-axis values of the Gaussian
            shiftbase = floor(length(x)/2);
            % keyboard
        case 'exponential';
            x = [0:1:10*sig]';   % spread to 10 times the time constant
            h = exp(-x/sig);
            shiftbase = length(x);
    end      
    h = h ./sum(h); % make sure kernel has unit area, then can use for rate functions
    
    
    % convolve window with data-spike trains
    [spkfcn{loop},idxs] = convolve_spiketrains(spkdata,h,shiftbase,Didxs,bins,bin,T,opts,x);
    Nidxs = numel(idxs);
    Rn{loop} = idxs;  % list of retained spike-train indices
    

    %% now compute selected distance between those functions
    [Sxy{loop}] = constructS(spkfcn{loop},Nidxs,opts); 

    % keyboard
    % do modularity on that graph
    edges(loop) = sum(sum(Sxy{loop}~=0)); 
    nodes(loop) = Nidxs;
    
    if ~opts.blnS & nodes(loop) >= opts.nlimit & edges(loop) > log(nodes(loop))   
        % only group if: (a) asked to do so; 
        % (b) there are enough nodes and
        % (c) the number of edges ensures a likely fully-connected graph
        
        if numel(opts.modopts) == 1
            [grps,Gmax(loop).Q,grpscon,Gcon(loop).Q,ctr] = allevsplitConTransitive(Sxy{loop},opts.modopts{1}); 
        elseif numel(opts.modopts) == 2
            [grps,Gmax(loop).Q,grpscon,Gcon(loop).Q,ctr]  = allevsplitConTransitive(Sxy{loop},opts.modopts{1},opts.modopts{2}); 
        else
            [grps,Gmax(loop).Q,grpscon,Gcon(loop).Q,ctr]  = allevsplitConTransitive(Sxy{loop}); 
        end
        
        % Qmax answers
        Gmax(loop).ngrps = max(grps);
        % get sizes 
        siz = [];
        for i = 1:Gmax(loop).ngrps
            siz = [siz numel(find(grps == i))];   
        end
        Gmax(loop).grpsizes = siz;
        
        % Qcon answers
        Gcon(loop).ngrps = max(grpscon);
        % get sizes 
        siz = [];
        for i = 1:Gcon(loop).ngrps
            siz = [siz numel(find(grpscon == i))];   
        end
        Gcon(loop).grpsizes = siz;

    else
        grps = zeros(Nidxs,1);  % nothing to group!
        Gmax(loop).ngrps = 0;
        Gmax(loop).grpsizes = [];
        Gmax(loop).Q = 0;
        grpscon = zeros(Nidxs,1);  % nothing to group!
        Gcon(loop).ngrps = 0;
        Gcon(loop).grpsizes = [];
        Gcon(loop).Q = 0;
    end

    % remap from retained index count to ID stamps...
    Gmax(loop).grps = [Rn{loop} grps];  % [ID stamp; Group membership]
    Gcon(loop).grps = [Rn{loop} grpscon];  % [ID stamp; Group membership]
      
end


varargout{1} = Sxy;
varargout{2} = spkfcn;
varargout{3} = Rn;


function [spkfcn,idxs] = convolve_spiketrains(spkdata,h,shiftbase,Didxs,bins,bin,T,opts,x)
    Nidxs = numel(Didxs);
    
    %% go round and compute spike-train binless functions
    spkfcn = zeros(numel(bins),Nidxs);
    nspikes = NaN; % just in case there are no spikes....
        
    for iN = 1:Nidxs
        % for each spike-train
        spkts = spkdata(spkdata(:,1) == Didxs(iN),2); % all spike-times passed to function
        spkts((spkts < T(1)) | (spkts > T(2))) = [];  % set of requested spike times
        nspikes(iN) = numel(spkts);  % count how many spikes
        spkts = spkts - T(1);  % time-stamp with respect to the first bin
        try
        for iS =  1:numel(spkts)
            spk = round(spkts(iS) / bin);
            ixs = spk+x;  % shift to spike-time
            spkfcn(ixs(ixs>0 & ixs<=numel(bins)),iN) = spkfcn(ixs(ixs>0 & ixs<=numel(bins)),iN) + h(ixs>0 & ixs<=numel(bins)); % convolve, truncating outside range
        end
        catch
            keyboard
        end
    end
    
%     for j = 1:Nidxs
%         currix = find(spkdata(:,1) == Didxs(j));
%         nspikes(j) = numel(currix);
%         if nspikes(j) > 0   % only bother doing convolution if there's something to convolve!!
%             [spk,bts] = spike_train_from_times(spkdata(currix,2),bin,T);
%             switch opts.BLmeth
%                 case 'Gaussian'
%                     try
%                         y = conv(h,spk);
%                         [r c] = size(y); if c>1 y = y'; end  % for reasons best known to Matlab, certain convolutions will return this as a row vector rather than a column vector
%                         shifty = y(shiftbase+1:end-shiftbase);   % shift convolved signal to line up with spike-times
%                         if numel(shifty) < numel(bins)
%                             % pad with zeros
%                             diffbins = numel(bins) - numel(shifty);
%                             shifty = [zeros(diffbins,1); shifty]; 
%                         end % can occasionally happen with width pars that are not integer multiples of step-size 
%                         spkfcn(:,j) = shifty;
%                     catch
%                         disp('I ran into a problem convolving the Gaussian')
%                         keyboard
%                     end
%                     
%                 case 'exponential'
%                     try
%                         y = conv(h,spk);
%                         [r c] = size(y); if c>1 y = y'; end  % for reasons best known to Matlab, certain convolutions will return this as a row vector rather than a column vector
%                         % keyboard
%                         spkfcn(:,j) = y(1:numel(bins));   % truncate convolved signal to line up with recording time 
%                     catch
%                         disp('I ran into a problem convolving the exponential')
%                         keyboard
%                     end
%             end
%         end
%         % keyboard
%     end
   
   %  keyboard
    
    idxs = Didxs; 
    try
        if any(nspikes == 0)
            % if not firing, then strip out cells
            spkfcn(:,nspikes==0) = []; 
            idxs(nspikes==0) = [];
        end
    catch
        disp('I ran into a problem removing non-firing cells')
        keyboard
    end
    
    % convert to firing rate (spikes/s)
    spkfcn = spkfcn ./ bin; % sums to 1 if spike in every bin

function [Sxy] = constructS(spkfcn,Nidxs,opts)

    switch opts.Dmeth            
        case {'corr','corrcoef'}
            Sxy = eval([opts.Dmeth '(spkfcn);']);     % compute correlation
            if any(isnan(Sxy))
                keyboard
            end
            Sxy(Sxy < 0) = 0;   % rectify correlation coefficient = similarity....
            % place zeros on the diagonal: no self-connections allowed...
            Sxy(eye(Nidxs)==1) = 0;
         case 'cosine'
            try 
                % how many entries in Sxy to compute? 
                inds = 1:Nidxs^2;
                uniqueinds = triu(reshape(inds,Nidxs,Nidxs)); % indexes of all unique pairs
                uniqueinds = uniqueinds - diag(diag(uniqueinds)); % remove self-indices
                [ip,jp] = ind2sub([Nidxs,Nidxs],uniqueinds(uniqueinds>0)); % get all [i,j] pairs...
                Sxy = zeros(Nidxs);
                for prs = 1:numel(ip)
                    A = spkfcn(:,ip(prs)); B = spkfcn(:,jp(prs));   % current pair
                    Sxy(ip(prs),jp(prs)) = (A'*B)/(norm(A)*norm(B));
                    Sxy(jp(prs),ip(prs)) = Sxy(ip(prs),jp(prs));
                end

            catch
                keyboard
            end

        otherwise
            error('Unknown pairwise similarity metric specified') 
    end
