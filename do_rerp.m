% do_rerp() - Compute overlap-corrected regression-Event-Related Potentials ("rERPs") using the OLS method.
%             Currently, only capable of handling two event latencies for a given trial (e.g., stimulus and 
%             response; or StimA and StimB). However, several features of these two events (i.e. 'difference 
%             waves', 'weighting') may be flexibly achieved via the optional input "design."
%
% Usage:
%   >> [rerp1, rerp2, rflag, ntrial] = do_rerp( data, event1, event2, varargin)
%
% Inputs (required): 
%   data      - NxP matrix of trial-level EEG data. Rows of N may correspond to N trials of data for one 
%               subject for one channel for one trial type (default), or the optional arguement 'cluster' 
%               (see below) may be used to accomodate multiple subjects, channels, trial types, etc. The
%   event1    - Nx1 column specifying the positive-valued sample in "data" (integer) where the 1st event 
%               occurs, if no "event1" for a given row, should be NaN
%   event2    - Nx1 column specifying the positive-valued sample in "data" (integer) where the 2nd event 
%               occurs, if no "event2" for a given row, should be NaN
%
% Inputs (optional key-value pairs):
%   'cluster' - Nx1 numeric column array specifying ways in which to "cluster" the rERP computation into 
%               C iterations. For example, {'cluster' [1 2 1 2 1 2]'} will compute separate rERPs for two
%               channels (or subjects, or trial types, etc.) having 3 trials of data each. 
%               default = ones(size(data,1),1)
%   'design'  - NxF numeric matrix that may be used to compute F overlap-corrected "factored" rERPs (cf. 
%               overlap-corrected "difference waves") and/or apply weighting to rERPs (see below examples).
%               default = ones(size(data,1),1)
%   'rerpwin' - Range start-end samples for which to compute rERPs around events (1x2 integer array)
%               default = round([-1 1] .* (size(data,2)/4))
%   'verbose' - Print to screen messages from the process (e.g., warnings)
%               default = 1
%
% Outputs:
%   rerp1     - CxRERPLEN array containing rERPs for the first event, where C corresponds to the number of 
%               "clusters" identified in optional inputs (default = 1) and RERPLEN corresponds to the number
%               of samples between rerpwin(1) and rerpwin(2)
%   rerp2     - CxRERPLEN array containing rERPs for the first event, where C corresponds to the number of 
%               "clusters" identified in optional inputs (default = 1) and RERPLEN corresponds to the number
%               of samples between rerpwin(1) and rerpwin(2)
%   rflag     - Cx1 logical array consisting of 1 if some values in rerp1 or rerp2 could not be estimated
%               using OLS regression (e.g., if design was rank-deficient, or if all event2 for a given trial
%               type were missing). For values in rerp1 and rerp2 not estimatable by the algorithm, will be 
%               output as NaNs. 
%   ntrial    - Cx1 numeric array consisting of the number of trials used to compute each rERP within rerp1 
%               and rerp2.
%
% Examples:
%  
%  [rerp1, rerp2] = do_rerp( data, event1, event2, )
%
%
%
% Citation(s):
% 
%
% Scott Burwell, 2018-11-18
function [rerp1, rerp2, rflag, ntrial] = do_rerp(data, event1, event2, varargin); 

%[rerp1, rerp2] = do_rerp(data, event1, event2, 'rerpwin', rerpwin, 'cluster', cluster, 'design', design);
% uses OLS (consistent w/ standard avereage ERPs), currently only handles overlap of two events, but can model multiple factors for each event

if isempty(keyval('verbose',varargin)),
   verbose= 1;                     else, verbose= keyval('verbose',varargin);
end
if isempty(keyval('cluster',varargin)),        
   cluster = ones(size(data,1),1); else, cluster = keyval('cluster',varargin);        
end
if isempty(keyval('rerpwin',varargin)),
   rerpwin= round([-1 1] * (size(data,2)/4)); 
   if verbose>0, disp(['   do_rerp; using default "rerpwin" of length ' num2str(rerpwin)]); end
   else, rerpwin= keyval('rerpwin',varargin);
end
if isempty(keyval('design',varargin)),
   design = ones(size(data,1),1);  else, design = keyval('design',varargin);
end
%if length(unique(design(:,1)))>1, 
%   design = [ones(size(data,1)) design];
%   disp(['   do_rerp; pre-pending column of ones to "design" for intercept']);
%end

droprows = find((event1+rerpwin(1))>size(data,2));
droprows = unique([droprows, find((event2+rerpwin(2))>size(data,2))]);
if ~isempty(droprows),
   if verbose>0, disp(['   do_rerp; WARNING some instances of event+rerpwin exceed the limits of the data, removing ' num2str(length(droprows)) ' row(s) from "data"']); end
   data(droprows,:) = ''; cluster(droprows,:) = ''; design(droprows,:) = ''; event1(droprows) = ''; event2(droprows) = '';
end

%make matrix of dummy values 
rerp1mtx = zeros(size(data,1),size(data,length(size(data)))); 
for ii = 1:length(cluster), if ~isnan(event1(ii)), rerp1mtx(ii,event1(ii)+rerpwin(1):event1(ii)+rerpwin(2)) = 1; end; end
rerp2mtx = zeros(size(data,1),size(data,length(size(data)))); 
for ii = 1:length(cluster), if ~isnan(event2(ii)), rerp2mtx(ii,event2(ii)+rerpwin(1):event2(ii)+rerpwin(2)) = 1; end; end

rerps  = [];
rflag  = [];
ntrial = [];
uclust = unique(cluster);

for ii = 1:length(uclust),
  if verbose>0, disp(['   do_rerp; obtaining rERPs for cluster : ' num2str(uclust(ii))]); end
    tmpclustidx = cluster==uclust(ii);
    tmprerp1mtx = rerp1mtx(tmpclustidx,:);
    tmprerp2mtx = rerp2mtx(tmpclustidx,:);
    tmpclustdata= data(tmpclustidx,:,:);
    if ~isempty(tmpclustdata),
      
      % handle tmpdesign, i.e., look for rank-deficiency-issues
      if size(design,2)>1,

        dropcols= [];        
        tmpdesign = design(tmpclustidx,2:end);

        %find perfectly correlated columns, zero-out...
        tmpcor  = corr(tmpdesign).*double(eye(size(corr(tmpdesign)))==0);
	[jr,jc] = find(tmpcor==1); if ~isempty(jc), dropcols= jc(2:end); end

        %find columns w/ zero variance
        dropcols= [dropcols, find(var(tmpdesign)==0) ];

        if ~isempty(dropcols), tmpdesign(:,dropcols) = 0; end; %make these columns all zeroes
        tmpdesign             = [design(tmpclustidx,1) tmpdesign]; 
      else,
        tmpdesign = design(tmpclustidx,:);
      end

      %compile X and Y variables
      X = [];
      Y = [];
      for kk = 1:size(tmprerp1mtx,1),
        tmpX1 = diag(tmprerp1mtx(kk,:));
        tmpX2 = diag(tmprerp2mtx(kk,:));
        tmpX1(:,sum(abs(tmpX1))==0) = '';  
        tmpX2(:,sum(abs(tmpX2))==0) = '';  
        if isempty(tmpX1), tmpX1 = zeros(size(data,2),length(rerpwin(1):rerpwin(2))); end %SJB added 2018-11-18 to flexibly handle trials where no elements >0
        if isempty(tmpX2), tmpX2 = zeros(size(data,2),length(rerpwin(1):rerpwin(2))); end
        tmpX  = [];
        for ll = 1:size(tmpdesign,2), tmpX = [tmpX [tmpX1, tmpX2].*tmpdesign(kk,ll)]; end
        X  = [X; tmpX];
      end

      for kk = 1:size(tmprerp1mtx,1), 
        if     length(size(tmpclustdata))==2,
              Y  = [Y;  tmpclustdata(kk,:)'];
        elseif length(size(tmpclustdata))==3,
              Y = cat(1,Y,squeeze(tmpclustdata(kk,:,:))');
        end
      end

      % remove extraneous rows of all 0s
      Y = Y(sum(abs(X),2)>0,:,:);                   %don't reorder this and next line.
      X = X(sum(abs(X),2)>0,:  ); Xcol = size(X,2); %

      % test for instances of rank-deficiency and log this information...
      rflagidx = zeros(1,Xcol);
      if sum(rank(X)<size(X))==2,
        curregflag          = 1; %cannot estimate full requested model
        rflagidx(sum(X)==0) = 1; %where in resulting Bs is invalid?
        X(:,rflagidx==1)    =''; %nullify redundant (i.e., all zero columns)
      else,
        curregflag          = 0; %can estimate full requested model
      end

      % solve for B (i.e., get the rERP)
      B = []; Yhat = zeros(1,Xcol);
      if size(Y,2)==1,
        B = inv(X' * X) * X' * Y;
        B = B';
      elseif size(Y,2)>1,
        for kk = 1:size(tmpclustdata,2),
          B = [B inv(X' * X) * X' * Y(:,kk)];
        end
        B = reshape(B',[1,size(B,2),size(B,1)]);
      end
      Yhat(rflagidx==0) = B;
      Yhat(rflagidx==1) = NaN;

      % store output data
      rerps = cat(1,rerps,Yhat);
      rflag = cat(1,rflag,curregflag);
      ntrial= cat(1,ntrial,size(tmpclustdata,1));
      %rimtx= cat(1,rimtx,reshape(repmat(1:size(design,2)*2,size(Yhat,2)./(size(design,2)*2),1),size(Yhat)));
      %rfmtx= cat(1,rfmtx,rflagidx);
  end
end

rerp1 = []; 
rerp2 = [];
if size(design,2)==1,
   rerp1 = rerps(:,1:length(rerpwin(1):rerpwin(2))       );
   rerp2 = rerps(:, (length(rerpwin(1):rerpwin(2))+1):end);
else,
  for ii = 1:size(design,2),
    rerp1 = [rerp1; rerps(1:length(rerpwin(1):rerpwin(2)))]; rerps(1:length(rerpwin(1):rerpwin(2))) = '';
    rerp2 = [rerp2; rerps(1:length(rerpwin(1):rerpwin(2)))]; rerps(1:length(rerpwin(1):rerpwin(2))) = '';
  end
end
rflag = logical(rflag);
