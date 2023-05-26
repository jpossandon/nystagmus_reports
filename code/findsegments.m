function [segments] = findsegments(nums,minsiz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [segments] = findsegments(nums,min)
%
% finds continous segments of 1 in vector of 0 and 1, including segments at
% start and end of sequence, if min provided, it finds only segment of
% length >= min
%
% 19.07.19 JPO
% Hamburg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
    minsiz = 0;
end
[r,c]               = size(nums);
if r>1 && c>1
    error('findsgements only works with vectors')
elseif c>1
    nums   = nums';
    [r,c]  = size(nums);
end
clusaux             = zeros(r+2,c);    % two more rows to find cluster at the begining and end
clusaux(2:end-1,:)  = nums;
clusaux             = diff(clusaux);   % find start and end of continuous segments of 1s
[~,j]               = ind2sub(size(clusaux),find(clusaux==1));
[~,jj]              = ind2sub(size(clusaux),find(clusaux==-1));
segments            = [find(clusaux==1)-j+1,find(clusaux==-1)-jj+1-1];

if minsiz>1
    keepseg = diff(segments')>=minsiz;
    segments = segments(keepseg,:);
end