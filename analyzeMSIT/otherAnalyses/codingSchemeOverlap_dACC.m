%% [20170426] compare interference sleective unit populaiton with beta and theta coherent populations. 

%% load data
dataDir = '~/Dropbox/MachineInvariantData';
load(fullfile(dataDir,'acc_dlpfc_units_results.mat'))
beta = load(fullfile(dataDir,'acc_neuronCoherenceResults_beta.mat'));
theta = load(fullfile(dataDir,'acc_neuronCoherenceResults_theta.mat'));
results = acc_dlpfc_units_results.model{7}.significant_classification_mat;
[acc_dlpfc_units,acc_units,dlpfc_units] = parseMSITGLMresults();


%% indices for rate-selective neurons
conflictNeurons = logical(results(logical(acc_units),1));
responseNeurons = logical(results(logical(acc_units),2));
feedbackNeurons = logical(results(logical(acc_units),3));


%% indices for coherent neurons
idx = 2;
B = beta.coherenceResults(idx).neuronIndices;
T = theta.coherenceResults(idx).neuronIndices; 


%% Now doing venn diagrams
% first inidividual areas
CONarea = sum(conflictNeurons);
Barea = sum(B);
Tarea = sum(T);

% then intersections
BTarea = sum((B+T)==2);
CONBarea = sum((B+conflictNeurons)==2);
CONTarea = sum((T+conflictNeurons)==2);
CONBTarea = sum((B+T+conflictNeurons)==3);

% then vectors 
A = [CONarea Barea Tarea]
I = [CONBarea CONTarea BTarea CONBTarea]


% % plotting diagrams. 
% figure 
% venn(A,I)

%% indices for coherent neurons
for id = 1:2;
B(:,id) = beta.coherenceResults(idx).neuronIndices;
T(:,id) = theta.coherenceResults(idx).neuronIndices; 

%% Now doing venn diagrams
% % first inidividual areas
% Barea = sum(B);
% Tarea = sum(T);

end




sum(Barea,2)
sum(Tarea,2)




