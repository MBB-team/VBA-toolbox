function [g,dgdx,dgdP] = g_ttest(x,P,u,in)
    % scenarios
    switch in.hypothesis
        case 'null'
            % prediction
            g = [P(1) ; P(1)]; % H0-prediction
            % derivative
            dgdP = zeros(2);
            dgdP(1,1) = 1;
            dgdP(2,1) = 1;
            
        case 'alternative'
            % prediction
            g = [P(1) ; P(2)]; % H1-prediction
            % derivative
            dgdP = eye(2);
    end
    
    
    % derivative
    dgdx = [];
end