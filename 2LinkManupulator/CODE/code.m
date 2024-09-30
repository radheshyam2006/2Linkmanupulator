clc,clear vars,close all;
% Simulation settings
t_span = 0:0.01:20;
y0 = [0.1, 0.1, 0, 0, 0, 0]; % Initial joint angles and velocities

% Desired joint angles (stop at [theta1_des, theta2_des])
theta1_des = 0;  % Desired angle for joint 1
theta2_des = 0;  % Desired angle for joint 2

% Run simulation
[t, s] = ode45(@(t, state) func(t, state, theta1_des, theta2_des), t_span, y0);

% Plotting
figure;
subplot(2, 1, 1);
plot(t, s(:, 1), 'LineWidth', 2);
title('Joint Angle q1 vs Time');
xlabel('Time (s)');
ylabel('q1 (rad)');

subplot(2, 1, 2);
plot(t, s(:, 2), 'LineWidth', 2);
title('Joint Angle q2 vs Time');
xlabel('Time (s)');
ylabel('q2 (rad)');

function der_S = func(t, state, theta1_des, theta2_des)
    % System parameters
    m1 = 10;
    m2 = 5;
    l1 = 0.2;
    l2 = 0.1;
    g = 9.81;

    % Controller gains
    kp1 = 200;
    kd1 = 5;
    ki1 = 500;
    kp2 = 400;
    kd2 = 150;
    ki2 = 600;

    % Extracting state variables
    q1 = state(1);
    q2 = state(2);
    q1_dot = state(3);
    q2_dot = state(4);
    neg_int_q1 = state(5);
    neg_int_q2 = state(6);

    % Equations of motion
    m11 = (m1 + m2) * (l1^2) + m2 * l2 * (l2 + 2 * l1 * cos(q2));
    m12 = m2 * l2 * (l2 + l1 * cos(q2));
    m22 = m2 * (l2^2);

    M = [m11, m12; m12, m22];

    c11 = -m2 * l1 * l2 * sin(q2) * q2_dot;
    c12 = -m2 * l1 * l2 * sin(q2) * (q1_dot + q2_dot);
    c21 = 0;
    c22 = m2 * l1 * l2 * sin(q2) * q1_dot;

    C = [c11, c12; c21, c22];

    g1 = m1 * l1 * g * cos(q1) + m2 * g * (l2 * cos(q1 + q2) + l1 * cos(q1));
    g2 = m2 * g * l2 * cos(q1 + q2);

    G = [g1; g2];

    % PDI control law (with desired angles)
    q1_error = q1 - theta1_des;
    q2_error = q2 - theta2_des;

    % PDI control law with integral terms
    tau1 = -kp1 * q1_error - kd1 * q1_dot + ki1 * neg_int_q1;
    tau2 = -kp2 * q2_error - kd2 * q2_dot + ki2 * neg_int_q2;

    % Control input vector
    Tau = [tau1; tau2];

    % Solve for q1_ddot and q2_ddot
    q_ddot = M \ (Tau - C * [q1_dot; q2_dot] - G);

    % State derivatives
    der_S = [q1_dot; q2_dot; q_ddot(1); q_ddot(2); -q1_error; -q2_error];
end