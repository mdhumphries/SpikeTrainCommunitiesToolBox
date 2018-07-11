function [G,Sgrp,Sin,Sout] = sortbysimilarity(G,Sxy)

% SORTBYSIMILARITY compute intra- and inter-group similarity and sort group order by similarity
% [NG,Sg,Sin,Sout] = SORTBYSIMILARITY(G,S), given the 2-column group structure vector G ([Cell ID, Group #]), 
% and the similarity matrix S from which that group structure was found (as
% returned by CLUSTER_SPIKE_DATA_X...), returns:
%   NG: the groups indexed by order of mean intra-group similarity: 1 - lowest; N - highest 
%   Sg: the mean intra-group similarity of the whole group (defines NG)
%   [note: is returned in order of the original groups]
%   Sin: the mean intra-group similarity of each node (ensemble average defines Sg) 
%   Sout: the mean inter-group similarity of each node (useful for
%   detecting weak group membership)
%
%   30/7/2016: fixed incorrect normalisation in calculation of mean
%   similarities (adjusted by -1 incorrectly)
%   23/9/2016: fixed divide by zero bug when group size = 1   
%   11//7/2018: fixed bug where did not compute similarity for groups of size 2
%
% Mark Humphries

IDs = unique(G(:,1));
Grps = unique(G(:,2));
nIDs = numel(IDs);
Ngrps = numel(Grps);
Sin = zeros(nIDs,1); Sout = zeros(nIDs,1); Sgrp = zeros(Ngrps,1); 

% compute each node's similarity within their group, and outside of their group...
for j = 1:nIDs
    thisgrp = G(j,2);
    ingrp = find(G(:,2) == thisgrp); ingrp(ingrp == IDs(j)) = []; % exclude itself
    outgrp = find(G(:,2) ~= thisgrp); 
    if numel(ingrp) > 0  % otherwise leave these as zero
        Sin(j) = sum(Sxy(j,ingrp)) ./ numel(ingrp);   % mean intra-group similarity for that node
    end
    Sout(j) = sum(Sxy(j,outgrp)) ./ numel(outgrp);   % mean inter-group similarity for that node
    % Sin(j) = median(Sxy(j,ingrp));   % median intra-group similarity for that node
    % Sout(j) = median(Sxy(j,outgrp));   % median inter-group similarity for that node
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
