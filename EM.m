function [lambdas] =  EM(words,speakers,Pr_Background,epsilon,niterations)
% [lambdas] =  EM(words,speakers,Pr_Background,epsilon,niterations)
%
% Function to fit model described in:
% T. Iwata and S. Watanabe, “Influence relation estimation based on lexical entrainment in conversation,”
% Speech Communication, vol. 55, no. 2, pp. 329–339, Feb. 2013.
%
% Uses EM for inference as described in Iwata paper.
%  
% INPUTS
%  words: 1 x T matrix, where each entry is a number representing a word.
%  Contains W unique values
%  speakers: 1 x T matrix, where each entry is a number representing the
%  speaker who spoke the corresponding word.
%  Pr_Background: 1 x W matrix, where each entry is the probability of a
%  word appearing (perhaps based on corpus as a whole or on an outside
%  corpus).
%  epsilon: convergence threshold used to detect convergence
%  niterations: maximum number of iterations to perform
%
% RETURNS
% lambdas: M x (M+1) matrix. The first column represents the influence of the background
% on each speaker. Entry (i,j+1) tells you the influence of speaker j on speaker i. Each row sums to 1.

BETA = 10^-8;
ALPHA = 1;

[~,~,words] = unique(words);
[~,~,speakers] = unique(speakers);

T = length(speakers);
M = length(unique(speakers));
W = length(unique(words));
lambdas = ones(M,M+1)./(M+1);
weights = zeros(T,M+1);

Pr_C = zeros(M,T);
for i = 1:M
    for t=1:T
        if t==1
            numerator = BETA;
            denominator = BETA*W;
        else
            w = words(t);
            numerator = sum(words(1:t-1)==w & speakers(1:t-1)==i) + BETA;
            denominator = sum(speakers(1:t-1)==i) + BETA*W;
        end;
        Pr_C(i,t) = numerator./denominator;
    end;
end;

prev_CDLL = -Inf;
for i=1:niterations
    % E-step
    for t=1:T
        s_t = speakers(t);
        w_t = words(t);
        Pr_m_given_t = zeros(M+1,1);
        for m=0:M
            if m==0
                Pr_m_given_t(m+1) = lambdas(s_t,m+1) .* Pr_Background(w_t);
            else
                Pr_m_given_t(m+1) = lambdas(s_t,m+1) .* Pr_C(m,t);
            end;
        end;
        weights(t,:) = (Pr_m_given_t./sum(Pr_m_given_t))';
    end;
    
    % M-step
    for m=0:M
        for n=1:M
            lambdas(n,m+1) = sum(weights(speakers==n,m+1)) + ALPHA;
        end;
    end;
    lambdas = lambdas ./ repmat(sum(lambdas,2),1,M+1); % normalize so that each row adds up to 1
    
    % Calculate complete-data log likelihood:
    first_term=0;
    for t=1:T
        s_t = speakers(t);
        w_t = words(t);
        total = 0;
        for m=0:M
            if m==0
                total = total + lambdas(s_t,m+1) .* Pr_Background(w_t);
            else
                total = total + lambdas(s_t,m+1) .* Pr_C(m,t);
            end;
        end;
        first_term = first_term + log(total);
    end;
    second_term = sum(sum(ALPHA .* log(lambdas),1),2);
    CDLL = first_term + second_term;
    
    if CDLL - prev_CDLL < epsilon
        break
    end;
    prev_CDLL = CDLL;
end;

end

