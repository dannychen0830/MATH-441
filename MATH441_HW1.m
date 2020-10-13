%% Problem 3

% prey per capita birth rate
alpha = 7.5;

% per capita growth per prey
gamma = 0.003;

% per capita death per predator
delta = 0.01;

% predator per capita death rate per predator
beta = [0.1, 0.2, 0.5];

% prey per capita death rate per prey
phi = [0.5, 1, 2];

counter = 1;
% Model
for i = 1: length(beta)
    for j = 1: length(phi)
        f = @(t, p) [gamma*p(1)*p(2) - beta(i)*p(1)*p(1);
        alpha*p(2) - delta*p(1)*p(2) - phi(j)*p(2)*p(1)];

        x0 = [75, 850];
        [time, sol] = ode45(f, [0,10], x0);

        subplot(3,3, counter)
        yyaxis left
        plot(time, sol(:,1),'linewidth',2)
        ylim([0 inf])
        
        yyaxis right 
        plot(time, sol(:,2),'linewidth',2)
        xlim([0 10])
        ylim([0 inf])
        title("beta =" + beta(i) + ", phi = " + phi(j))
        counter = counter + 1;
    end
end

%% Problem 4

% rate of susceptible to exposed
alpha = 1/56000;

% rate of exposed to infected
beta = 1/15;

% rate of infected to recovery
gamma = 1/21;

% natural death rate
phi = 1/27375;

t = [1, 50, 100, 200];

%model
f = @(t, p) [phi*(p(2)+p(3)+p(4)) - alpha*p(1)*p(3);
            alpha*p(1)*p(3) - beta*p(2) - phi*p(2);
            beta*p(2) - gamma*p(3) - phi*p(3);
            gamma*p(3) - phi*p(4)];
p0 = [39999, 0, 1, 0];

for i = 1:length(t)  
    [time, sol] = ode45(f, [0, 365*t(i)], p0);
    
    subplot(2,4,i)
    plot(time, sol,'linewidth',1.5)
    ylim([0 inf])
    legend('S','E','I','R')
    title(t(i) + " years")
    xlabel('days')
    ylabel('population')
    
    subplot(2,4,i+4)
    plot(time,sol,'linewidth',1.5)
    ylim([0 50])
    xlabel('days')
    ylabel('population')
    legend('S','E','I','R')
end


    

