function P = expectedA(A)
% EXPECTEDA computes the expected number of connections between each vertex
%   P = EXPECTEDA(A) computes the matrix P of the expected number of
%   connections between each vertex, given the adjacency matrix A.
%
%   Notes: 
%   (1) the expected number of connection is computed based on the
%   assumption of the "null model": i.e. random connections between nodes
%   of the same degree sequence as A. Thus, for undirected graphs, Pij =
%   (ki*kj) / 2m. For directed graphs, we have to make the additional
%   assumption that the in-degree and out-degree distributions are
%   themselves not correlated - i.e. that we can place an edge between any
%   randomly selected pair of available out-nodes and in-nodes.
%
%   (2) In addition, the underlying "null model" allows multiple and
%   self-edges, and hence has values along the diagonal.
%
%   References: Newman, M. E. J. (2006) "Finding community structure in
%   networks using the eigenvectors of matrices". Phys Rev E, 74, 036104.
%
%   Mark Humphries 24/8/2006

[n c] = size(A);
inks = sum(A);
outks = sum(A');
m = sum(inks);  % number of edges

P = (outks' * inks) ./ m;