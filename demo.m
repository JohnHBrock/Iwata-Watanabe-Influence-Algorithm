words = [1 2 3 4 5 6 1 2 3 4 7 8 1 2 3 4 5 6 7 8];
speakers = [1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 3 3];
Pr_Background = [.125 .125 .125 .125 .125 .125 .125 .125];
EPSILON = 10^-9;
MAX_ITERATIONS = 200;

lambdas = EM(words, speakers, Pr_Background, EPSILON, MAX_ITERATIONS);
disp(lambdas);

visualize_influence(lambdas);
