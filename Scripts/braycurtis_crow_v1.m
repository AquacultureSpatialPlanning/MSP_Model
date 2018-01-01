function d=braycurtis(x,y)

%takes two vectors, like this from Filter_Seed_v6_BrayCurtis.m
% x=Policy_i_a_filter_unique(1,:)';
% y=Policy_i_a_filter_unique(2,:)';

%%% Calculate the Bray-Cortis distance between samples x and y

%%% Check: x and y must be column vectors of the same size
if length(x) ~= length(y)
    fprintf('Error: vectors x and y must have the same size.\n\n')
    return
end    

if iscolumn(x)==0 | iscolumn(y)==0
    fprintf('Error: x and y must be column vectors.\n\n')
    return
end

%% Computation of the Bray-Curtis distance (dissimilarity index):
% http://homepages.ulb.ac.be/~dgonze/INFO/matlab/braycurtis.m
%%% Author: Didier Gonze
%%% Created: 3/8/2013
%%% Updated: 3/8/2013
% Also see http://www.umass.edu/landeco/teaching/multivariate/readings/McCune.and.Grace.2002.chapter6.pdf
% m=sum(min([x y],[],2));
% s=sum(x)+sum(y);
% d=1-2*m/s;

% http://www.code10.info/index.php%3Foption%3Dcom_content%26view%3Darticle%26id%3D46:articlebray-curtis-dissim%26catid%3D38:cat_coding_algorithms_data-similarity%26Itemid%3D57
numerator=sum(abs(x-y));
denominator=sum(x+y);
d=numerator/denominator;

