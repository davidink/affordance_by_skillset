classdef LogisticRegression < ProbClassifier & handle
%classdef LogisticRegression
    properties
        inputDims;
        featureDims;
        outputDims;
        W;
        phiFunc;
    end
    
    methods
        %Constructor
        function [obj]=LogisticRegression(nrInputDims,nrOutputDims,featureFunctionName)
            obj.phiFunc=str2func(featureFunctionName);
            
            obj.inputDims=nrInputDims;
            obj.outputDims=nrOutputDims;
            obj.featureDims=length(obj.phiFunc(ones(1,obj.inputDims)));
            
            obj.W=zeros(obj.featureDims,obj.outputDims);
            
        end
        
        %Train the model from weighted input-output data. The first theta
        %is the regularization/(weight prior) for the regression
        function [success]=Train(obj, X,Y,theta)
            Phi=obj.computeFeatures(X);
            lambda=theta(1);
            labels=Y;
            impW=sum(labels,2)*ones(1,size(labels,2)); % normalize
            obj.W=zeros(size(obj.W)); % initialize with 0 for weights
            
            H=zeros(obj.featureDims*obj.outputDims); % Not used?
            
            prevError=10^6;
            loopFlag=1;
            y=exp(Phi*obj.W); % compute prediction on current weights
            y=y./repmat(sum(y,2),1,obj.outputDims); % normalize
            y=y.*impW;
           
            counter=0;
            start_loop = tic;
            while(loopFlag)
                counter=counter+1;
              
                
                grad=Phi'*(y-labels)+lambda*obj.W; % compute gradient of error function
 
                
                Hess=zeros(obj.featureDims*obj.outputDims); % Generate Hessian Matrix, MxM blocks with each block size for DxD
                for j=1:obj.outputDims
                    for k=1:obj.outputDims
                        rIndex=(1+(j-1)*obj.featureDims):(j*obj.featureDims);
                        cIndex=(1+(k-1)*obj.featureDims):(k*obj.featureDims);
                        
                        if j==k
                            R=diag(y(:,k).*(sum(labels,2)-y(:,j))); % Rnn = y_n * (1-y_n)
                        else
                            R=diag(y(:,k).*(0-y(:,j)));;
                        end
                        
                        Hess(rIndex,cIndex)=Phi'*R*Phi; % Hessian for logistic regression
                    end
                end
                Hess=Hess+lambda*eye(size(Hess)); % Add regularization term
                
   
                
                while det(Hess)<1e-20 % select for correct size of regularization coefficient
                    Hess=Hess+lambda*eye(size(Hess));
                    grad=grad+lambda*obj.W;
                    lambda=lambda*2;
                end
                Grad=reshape(grad,obj.outputDims*obj.featureDims,1);
                step=reshape(Hess\Grad,obj.featureDims,obj.outputDims); %compute the gradient descent w = w_old - H^-1*DeltaE
                if max(max(abs(step)))>10
                    step=10*step/max(max(abs(step))); %scale the decay rate if needed
                end
                
                obj.W= obj.W-step; % new weights
                
                y=exp(Phi*obj.W); % new prediction
                y=y./repmat(sum(y,2),1,obj.outputDims);
                y=y.*impW;
                newError=-sum(sum(labels.*log(y+1e-8)))+0.5*lambda*sum(sum(obj.W.^2)); % new error function value
                %newError=-sum(sum(labels.*log(y)))+0.5*lambda*sum(sum(obj.W.^2)); % new error function value
                
                if prevError-newError<=0 && counter>20
                    loopFlag=0;
                    success = true;
                else
                    prevError=newError;
                end
                
                if findstr(lastwarn,'Matrix is singular')
                    disp('Logistic Regression aborting due to singular matrix..');
                    success = false;
                    loopFlag=0;
                end
                if (toc(start_loop) > 5) % here, 2 seconds
                    %error('Took toooo loooong');
                    disp('Logistic Regression aborting due to timeout..');
                    success = false;
                    loopFlag=0;
                end                
            end
        end
        
        %Sample the output given an input
        function [Y]=Sample(obj,X)
            prob=cumsum(obj.computeLikelihood(X,[]),2);
            Y=zeros(size(X,1),1);
            for i=1:size(X,1)
                Y(i)=find(rand<prob(i,:),1);
            end
            
        end
        
        %Compute probability of output given input
        function [Y]=computeExpected(obj,X)
            prob=obj.computeLikelihood(X,[]);
            Y=zeros(size(X,1),1);
            for i=1:size(X,1)
                Y(i)=find(prob(i,:)==max(prob(i,:)),1);
            end
        end
        
        %Compute probability of output given input
        function [prob]=computeLikelihood(obj,X,Y)
            Phi=obj.computeFeatures(X);
            nrSamples=size(X,1);
            
            PROB=exp(Phi*obj.W);
            PROB=PROB./repmat(sum(PROB,2),1,obj.outputDims);
            
            if length(Y)==0
                prob=PROB;
            else
                prob=zeros(nrSamples,1);
                for i=1:nrSamples
                    prob(i)=PROB(i,Y(i));
                end
            end
        end
        
        %Compute feature vector for given input
        function [Phi]=computeFeatures(obj,X)
            Phi=zeros(size(X,1),obj.featureDims);
            for i=1:size(X,1)
                Phi(i,:)=obj.phiFunc(X(i,:));
            end
        end
    end
end

