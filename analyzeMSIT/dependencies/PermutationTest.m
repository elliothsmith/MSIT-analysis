function [pvalue] = PermutationTest(data,groups,permutations,alpha)
%{
Function to compute the Permutation test for different data grouped into
several groups. This test allow to assess the degree of statistical
significance between different samples or populations. It avoids the
multicomparison problems with may arise.

INPUT:
    data: MATRIX with the data that one wants to analyze with the 
    permutation test. The format of the matrix is the following:
        * rows: each observation to analyze. For instace, the result of
        each control subject or patient. 
        * columns: as many columns as sensors employed to retrieve the
        information or signal in each observation. For example: MEG of 148
        sensors.
    groups: COLUMN VECTOR which sets how the 'data' are grouped for the
    combination that one wants to analyze. Example:
    groups = [1,1,1,0,0,0].'; means that the three first belong to one
    group (control subjects) and the latest three belong to other (patients).
    permutations: the number of permutations that one wants to perfom.
    Recommended value: 5000.
    alpha: the significance threshold. By default, 0.05

OUTPUT:
    pvalue: ROW VECTOR which holds the p-values computed by means of the
    permutations test. (Tip: same length as sensors or channels).

PROJECT:  Research Master in signal theory and bioengineering - University of Valladolid

DATE: 12/11/2012 and 28/03/2015

AUTHOR: Jess Poza Crespo and Jess Monge lvarez.
%}
%% Checking the required ipunt parameters:
control = ~isempty(data);
assert(control,'The user must introduce a matrix with the data (first inpunt).');
control = ~isempty(groups);
assert(control,'The user must introduce a vetor for grouping the data (second inpunt).');
control = ~isempty(permutations);
assert(control,'The user must indicate the number of permutations that they want (third inpunt).');

%% Configuration of the parameters:
% Intializatons: 
pvalue = [];
if (nargin == 3)
    alpha = 0.05;
end

% The number of permutation must allow to apply the significance threshold,
% that is to say, the number of permutations must be greater than a
% recommended minimum:
control = ~((1/permutations) > alpha);
assert(control,'The number of permutations is not valid for the selected significance threshold.');

% The number of permutations must be lesser than the maximum possible,
% N1+N2+...+NN)!/(N1!xN2!x...xNM!):
n = hist(groups); n = n(logical(n)); den = 1;
for i = 1:length(n)
    den = den*factorial(n(i)); 
end 
max_permutations = factorial(sum(n))/den;
control = ~(permutations > max_permutations);
assert(control,'The number of permutations exceeds the maximum possible value.'),

% We check if the number of observations matches up with the number of 
% total elements in the groups:
control = ~(size(data,1) ~= length(groups));
assert(control,'The number of observations differs from grouping elements.'),

% We create the matrix of N permutations, where the first one is what we
% want to analyze:
indexes_permuted = (1:size(data,1))';
num_indexes_permuted = size(indexes_permuted,2); % Initialization to 1
while (num_indexes_permuted < permutations)
    % Permutation of the indixes:
    index = randperm(length(groups))';
    % We check if the permutation is repetead:
    if ~any(all(indexes_permuted == repmat(index,1,size(indexes_permuted,2))))
        % If it is not repeated, we put it in a new column:
        indexes_permuted = [indexes_permuted,index];
        num_indexes_permuted = num_indexes_permuted + 1;
    end 
end 

%% Computation:
% First, we compute the statistic for all the permutations:
maximal_distribution = zeros(1,permutations);
for i = 1:permutations
    % Permutated data:
    data_permuted = data(indexes_permuted(:,i),:);
    % We initialize the statisticfor the analyzed sensors or channels:
    statistic = zeros(1,size(data,2));
    % We go over the sensors or channels:
    for sensor = 1:size(data,2)
        % We compute the statistic which quantifies the magnitud of the
        % effect in the observations, as the statistic 'F':
        [p,table,~] = anova1(data_permuted(:,sensor),groups,'off');
        statistic(1,sensor) = table{2,5};
        % One can also use the difference between averaged values:
        % statistic(1,sensor) = mean(data_permuted(1:n(1),sensor),1) - ...
        % mean(data_permuted(n(1)+1:n(2),sensor),1); 
    end
    % We select the statistic with the maximum value between all the
    % sensors or channles and the maximum values distribution is formed:
    maximal_distribution(1,i) = max(statistic);
    % In the case of differences between averaged values, the more negative
    % the greater effect:
    % maximal_distribution(1,i) = min(statistic);
    
    % We inform the user:
    fprintf('%.0f de %.0f \n',i,permutations);
end

% Secondly, we generate the histogram of the maximum values distribution:
maximal_distribution = sort(maximal_distribution);
% In the case of differences between averaged values, the most negative
% values appears on the right:
% maximal_distribution = sort(maximal_distribution,'descend');

% Finally, for the observation of interest, we compute the p-values
% associated to each sensor or channel:
for sensor = 1:size(data,2)
    [p,table,~] = anova1(data(:,sensor), groups,'off');
    statistic = table{2,5};
    c = min(find(maximal_distribution > statistic));
    % In the case of differences between averaged values, the more negative
    % the greater effect:
    % statistic = mean(data(1:n(1),sensor),1)-mean(data(n(1)+1:n(2),sensor),1);
    % c = min(find(maximal_distribution < statistic));
    if isempty(c)
        c = 1; 
    end
    
    % We compute the p-value as the ratio between the number of elements
    % with greater statistics thant the current one and the total number of
    % permutations:
    pvalue(sensor) = (permutations-c)/permutations;
end 

end % End of the 'PermutationTest.m' function