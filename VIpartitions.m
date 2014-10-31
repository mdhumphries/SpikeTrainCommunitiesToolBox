function VI = VIpartitions(A,B)

% VIPARTITIONS compute variational information of two graph partitions
% VI = VIPARTITIONS(A,B) computes the variational information of the
% two graph partitions given by the column vectors A and B: each is a
% vector of integers indicating community membership of each node (i.e. entries in A are 1..i, 
% entries in B are 1..k) with i communities in A and k communities in B
%
% Returns: VI, the variational information of these two partitions
% VI = 0 : indicates identical partitions
% VI = ln(n) : indicates total independence of the partitions (where n is the
% total number of nodes in the original network)
%
% Notes:
% (1) For comparison between a "real" partition and a "found" partition,
% it is conventional to assign the real partition to A and the found partition to B.
%
% (2) For comparison between networks, it is useful to normalise VI by
% ln(n)
%
% References:
% (1) Merila, M. (2007) Comparing clusterings--an information based distance 
%       Journal of Multivariate Analysis, 98, 873-895 
% (2) Karrer, B.; Levina, E. & Newman, M. E. J. (2008) Robustness of community
% structure in networks. Phys Rev E, 77, 046119
%
% Mark Humphries 4/12/09

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
NT = sum(sum(N));   % sum over whole matrix

% keyboard


HAB = 0;
HBA = 0;

% keyboard
for x = 1:nA
    nx = numel(find(A == x));
    for y = 1:nB
        if N(x,y)
            % if zero, then adds nothing here....
            ny = numel(find(B == y));           
            pxy = N(x,y) / NT;
            px = nx/NT;
            py = ny/NT;
            
            HAB = HAB + pxy * log(pxy/py);
            HBA = HBA + pxy * log(pxy/px);
        end
    end
end

VI = -HAB - HBA;
