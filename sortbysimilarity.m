function [G,Sgrp,Sin,Sout] = sortbysimilarity(G,Sxy)

% SORTBYSIMILARITY compute intra- and inter-group similarity and sort group order by similarity
% [NG,Sg,Sin,Sout] = SORTBYSIMILARITY(G,S), given the 2-column group structure vector G ([Cell ID, Group #]), 
% and the similarity matrix S from which that group structure was found (as
% returned by CLUSTER_SPIKE_DATA_X...), returns:
%   NG: the groups indexed by order of mean intra-group similarity 
%   Sg: the mean intra-group similarity of the whole group (defines NG)
%   Sin: the mean intra-group similarity of each spike-train (ensemble average defines Sg) 
%   Sout: the mean inter-group similatiry of each spike-train (useful for
%   detecting weak group membership)
%
% Mark Humphries 18/10/2011

IDs = unique(G(:,1));
Grps = unique(G(:,2));
nIDs = numel(IDs);
Ngrps = numel(Grps);
Sin = zeros(nIDs,1); Sout = zeros(nIDs,1); Sgrp = zeros(Ngrps,1); 

% compute each neurons similarity within their group, and outside of their group...
for j = 1:nIDs
    thisgrp = G(j,2);
    ingrp = find(G(:,2) == thisgrp); 
    outgrp = find(G(:,2) ~= thisgrp); 
    Sin(j) = sum(Sxy(j,ingrp)) ./ (numel(ingrp)-1);   % mean intra-group similarity for that neuron
    Sout(j) = sum(Sxy(j,outgrp)) ./ (numel(outgrp)-1);   % mean inter-group similarity for that neuron
    % Sin(j) = median(Sxy(j,ingrp));   % median intra-group similarity for that neuron
    % Sout(j) = median(Sxy(j,outgrp));   % median inter-group similarity for that neuron
end

% compute each group's similarity measure
for j = 1:Ngrps
    ingrp = find(G(:,2) == j); % find all group members 
    Sgrp(j) = sum(Sin(ingrp))/numel(ingrp); % mean intra-group similarity
    % Sgrp(j) = median(Sin(ingrp)); % median intra-group similarity
end

% sort group order by intra-group similarity measure
[srted,s_idx] = sort(Sgrp);    % low-to-high

% reassign group membership indices: group 1: lowest; group N: highest 
tempG = G; 
for j = 1:Ngrps  G(tempG(:,2) == s_idx(j),2) = j; end  % remap group membership indices
