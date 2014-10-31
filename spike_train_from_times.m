function [y,bins] = spike_train_from_times(times,bin_size,T)

% SPIKE_TRAIN_FROM_TIMES create spike train array
%
%   SPIKE_TRAIN_FROM_TIMES(A,BINSIZE,T) converts the spike train represented by event times A (in seconds) into
%   a binary spike train array using bins BINSIZE seconds wide, and where T is a two-element array specifying the start
%   and end times in seconds. 
%   
%   [Y,BINS] = SPIKE_TRAIN_FROM_TIMES(...) returns the spike train Y and the bin centers (times) BINS
%
%   Mark Humphries 22/4/2004

time_seconds = T(2) - T(1);
bins = T(1)+bin_size/2:bin_size:T(2)-bin_size/2;
y = hist(times(times<T(2)),bins);

