function [parameter,state] = getStateParam(X,P,options,functionType)

   
    %
    options = check_struct(options,'inF',struct,'inG',struct);
    
    switch functionType
        case 'evolution'    
            inS = options.inF;
        case 'observation'
            inS = options.inG;
    end
  
    

if  isfield(inS,'paramLabel') 
   parameter = priorPrettifyer(inS.paramLabel,P) ;
else
   parameter = [] ; 
end

if  isfield(inS,'stateLabel') 
   state = priorPrettifyer(inS.stateLabel,X) ;
else
   state = [] ; 
end

