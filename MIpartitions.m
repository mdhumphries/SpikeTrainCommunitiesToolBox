function MI = MIpartitions(A,B)

% MIPARTITIONS compute normalised mutual information of two graph partitions
% I = MIPARTITIONS(A,B) computes the normalised mutual information of the
% two graph partitions given by the column vectors A and B: each is a
% vector of integers indicating community membership (i.e. entries in A are 1..i, 
% entries in B are 1..k) with i communties in A and k communities in B
%
% Returns: I, the normalised mutual information of these two partitions
% I = 1 indicates identical partitions
% I = 0 indicates total independence of the partitions
%
% Notes:
% (1) For comparison between a "real" partition and a "found" partition,
% it is conventional to assign the real partition to A and the found partition to B.
%
% References:
% Danon, L., Diaz-Guilera, A., Duch, J. & Arenas, A. (2005) Comparing community 
% structure identification J Stat Mech, P09008  
%
% Mark Humphries 2/12/09

nA = max(A); nB = max(B);

% make confusion matrix
N = zeros(nA,nB);

for i = 1:nA
    grpA = find(A == i);
    for j = 1:nB
        grpB = find(B == j);
        common = intersect(grpA,grpB);
        N(i,j) = numel(common);
    end
end

Ni = sum(N,2);   % sum over rows
Nj = sum(N,1);    % sum over columns
NT = sum(Nj);   % sum over whole matrix

% keyboard

%%%% Dongen et al form
NormA = sum(Ni .* log(Ni./NT));
NormB = sum(Nj .* log(Nj./NT));

rMI = 0;
for i = 1:nA
    for j = 1:nB
        if N(i,j)
            % if zero, then adds nothing here....
            rMI = rMI + N(i,j) * log(N(i,j)*NT/(Ni(i)*Nj(j)));
        end
    end
end

MI = -2 * rMI ./ (NormA + NormB);




