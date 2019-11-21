%% first testing for proportions of temporal coding neurons vs. rate coding neurons. 
% rate coding vs temporal coding contingency tables
X = [8 6; 65 136-65-14];
mcnemar(X)



%% [20190325] Some of the numbers below are incorrect. 
%The numbers above correctly acount for the non-independence of samples. 


%% [20170907] this script runs tests for proportions in the MSIT paper


%% testing for differences in proportions of neurons classified as task selective
correct = true;

% interference vs. response
X = [14 12];
N = [136 136];

display('results for Ho: same proportion of conflict and response cells')
[h,p,chi2stat,df] = prop_test(X,N,correct)
clearvars -except correct

% interference vs. block
X = [14 11];
N = [136 136];

display('results for Ho: same proportion of conflict and block cells')
[h,p,chi2stat,df] = prop_test(X,N,correct)
clearvars -except correct

% interference vs. response
X = [12 11];
N = [136 136];

display('results for Ho: same proportion of response and block cells')
[h,p,chi2stat,df] = prop_test(X,N,correct)
clearvars -except correct


%% testing for proportions of task-selective neurons in dACC vs. dlPFC

% conflict 
X = [14 15];
N = [136 367];

display('results for Ho: same proportion of conflict neurons in dACC and dlPFC')
[h,p,chi2stat,df] = prop_test(X,N,correct)
clearvars -except correct

% block 
X = [11 24];
N = [136 367];

display('results for Ho: same proportion of block neurons in dACC and dlPFC')
[h,p,chi2stat,df] = prop_test(X,N,correct)
clearvars -except correct

% response 
X = [12 18];
N = [136 367];

display('results for Ho: same proportion of response neurons in dACC and dlPFC')
[h,p,chi2stat,df] = prop_test(X,N,correct)
clearvars -except correct


%% testing for proportions of coherent dACC neurons w/ dACC and dlPFC LFPS 

% beta 
X = [54 43];
N = [136 136];

display('results for Ho: same proportion of dACC neurons cohere with dACC and dlPFC beta')
[h,p,chi2stat,df] = prop_test(X,N,correct)
clearvars -except correct

% theta 
X = [37 49];
N = [136 136];

display('results for Ho: same proportion of dACC neurons cohere with dACC and dlPFC theta')
[h,p,chi2stat,df] = prop_test(X,N,correct)
clearvars -except correct


 


