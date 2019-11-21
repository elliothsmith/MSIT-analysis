function S=normlogspec(S,varargin)
% NORMLOGSPEC normalize the log of the spectrogram
%    S=NORMLOGSPEC(S) normalizes the log of the input spectrogram S. S
%    should be in the format time x frequency (x channel or trial), e.g., 
%    as returned from the Chronux MTSPECGRAMC function. Frequencies are 
%    normalized across time so that the range of each frequency is [0,1].
%
%    S=NORMLOGSPEC(S,FLAG), when FLAG=1, removes outliers before 
%    calculating maxima and minima used to normalize S.  This method is 
%    substantially slower but may help to reduce the effect of artifact on 
%    the normalization.  The default behavior is specified when FLAG=0.
%
%    S=NORMLOGSPEC(S,FLAG,DIM) specifies that the spectrogram should be
%    normalized along dimension DIM.  Default is along the first dimension,
%    which is time for Chronux output.  *NOTE: This option may only be used
%    if FLAG=0.

% defaults
FLAG=0; % do not remove outliers before calculating max/min (faster)
DIM=1; % default dimension: time if the input is from chronux

% user inputs
if(nargin>=2)
    FLAG=varargin{1};
end
if(nargin==3)
    DIM=varargin{2};
end

% check size
if(numel(size(S))<2)
    error('normlogspec:baddim','Input should contain at least two dimensions');
end

% check input consistency
if(DIM~=1 && FLAG~=0)
    error('normlogspec:badarg','Outlier removal only supported for normalizing along the first dimension');
end

% log-transformation
S=10*log10(S);

% normalize
if(FLAG==1) % discard outliers

    for ch=1:size(S,3) % loop over channels
        for fr=1:size(S,2)% loop over frequencies

            % subtract minima
            tmp=S(:,fr,ch);
            tmp(outliers(tmp))=[];
            S(:,fr,ch)=S(:,fr,ch)-repmat(min(tmp),[size(S,1) 1]);
            S(S(:,fr,ch)<0,fr,ch)=0; % hard limit the lower end to zero

            % divide by maxima
            tmp=S(:,fr,ch);
            tmp(outliers(tmp))=[];
            S(:,fr,ch)=S(:,fr,ch)./repmat(max(tmp),[size(S,1) 1]);
            S(S(:,fr,ch)>1,fr,ch)=1; % hard limit the upper end to one

        end
    end

else % use absolute maxima, minima

    % subtract minima
    S=S-repmat(min(S,[],DIM),[ones(1,DIM-1) size(S,DIM) ones(1,numel(size(S))-DIM)]);

    % divide by maxima
    S=S./repmat(max(S,[],DIM),[ones(1,DIM-1) size(S,DIM) ones(1,numel(size(S))-DIM)]);

end