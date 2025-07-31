clc; clear; close all;

% Parameters
total_time = 1000;
num_channels = 50;
jam_prob = 0.3; % 100% chance a channel is jammed
seed = 42;
trials = 1000;  % increase from 100 to 1000 or more


% DASH Hopping (Formula)
rng(seed);
dash_sequence = mod(seed .* (1:total_time) + floor((1:total_time).^1.1), num_channels) + 1;

% Random Hopping
rng(seed + 100);
random_sequence = randi([1, num_channels], 1, total_time);

% Jamming Simulation: 1 channel jammed every time unit
jammed_channels = 7;

% Success Check
success_dash = (dash_sequence ~= jammed_channels);
success_random = (random_sequence ~= jammed_channels);

% Visualization
figure;
subplot(3,1,1);
plot(dash_sequence, 'g', 'LineWidth', 2); hold on;
plot(random_sequence, 'b--', 'LineWidth', 2);
ylabel('Channel');
title('DASH (green) vs Random (blue) Hopping');
legend('DASH', 'Random');

subplot(3,1,2);
stem(jammed_channels, 'r', 'filled');
ylabel('Jammed Channel');
title('Jammer Activity');

subplot(3,1,3);
bar(1:total_time, [success_dash' success_random'], 1);
ylim([0 1.2]);
legend('DASH Success', 'Random Success');
ylabel('Success (1=OK, 0=Jammed)');
xlabel('Time');
title('Success Comparison');

% Print Success Rates
fprintf('DASH Success Rate: %.2f%%\n', mean(success_dash) * 100);
fprintf('Random Success Rate: %.2f%%\n', mean(success_random) * 100);
