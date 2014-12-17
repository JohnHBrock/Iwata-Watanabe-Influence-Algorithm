function [] = visualize_influence(lambdas)
%visualize_influence Plot the influence network

M = size(lambdas,1);

difference = triu(lambdas(:,2:end) - triu(lambdas(:,2:end)',1));

radius = 10;
theta = linspace(.5*pi,2.5*pi,M+1);
theta = theta(1:length(theta)-1);
rho = ones(1,M)*radius;
[X,Y] = pol2cart(theta,rho);

colors = repmat('rgbcmyk',1,10);

hold on;
axis equal;

THRESHOLD = .05; % arbitrary
line_width_multiplier = 100;
arrow_length_multiplier = 400;
for r=1:M
    for c=r+1:M % don't show background influences
        if c==r % don't show self-transitions
            continue;
        end;
        if abs(difference(r, c)) < THRESHOLD
            continue;
        end;
        line_width = abs(difference(r,c))*line_width_multiplier;
        arrow_length = abs(difference(r,c))*arrow_length_multiplier;
        if difference(r,c) > 0
            start = [X(1,c) Y(1,c)];
            finish = [X(1,r) Y(1,r)];
            arrow(start, finish, 'Width', line_width, 'Length', arrow_length, 'BaseAngle', 30, 'FaceColor', colors(1,c), 'EdgeColor', colors(1,c));
        else
            start = [X(1,r) Y(1,r)];
            finish = [X(1,c) Y(1,c)];
            arrow(start, finish, 'Width', line_width, 'Length', arrow_length, 'BaseAngle', 30, 'FaceColor', colors(1,r), 'EdgeColor', colors(1,r));
        end;
    end;
end;
plot(X,Y,'r*');
labels = cellstr(num2str([1:M]'));  %' # labels correspond to their order
text(X, Y, labels, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'red')

hold off;

end

