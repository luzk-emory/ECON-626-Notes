clear all;clc;
% luzhikun
beta = 0.98; alpha = 0.3; delta = 1;

f_ss = @(k)alpha*k^(alpha-1)+1-delta-1/beta
k_ss = fsolve(f_ss,1)
k_initial = 0.5*k_ss

num_state = 1000;
phi = 2*k_ss/num_state;
k_state = phi:phi:(2*k_ss);

[K_x, K_y] = meshgrid(k_state,k_state);

v = zeros(1,num_state);
epsilon = 10^(-20)

num_iter = 100

xx = [1, 10, 50, 100]
figure(1)
hold on 
for ii = 1:num_iter
	v_primitive = v;
	c =  max( K_x.^alpha + (1 - delta)*K_x - K_y, epsilon);
	v_matrix =  log( c ) + beta*v'*ones(1,num_state);
	[v_improved, k_choice_vector] = max(v_matrix);
	v = v_improved;
	%error_ = max(v_improved - v_primitive)
	if (ii == xx(1))|(ii == xx(2))|(ii == xx(3))|(ii == xx(4))
		plot(k_state, v, 'linewidth', 2)
	end
end

% exact solution
v_exact = v;
a = 1/(1-beta)*(log(1-alpha*beta)+alpha*beta/(1-alpha*beta)*log(alpha*beta));
b = alpha/(1-alpha*beta);
v_exact = a+b*log(k_state);
plot(k_state, v_exact, 'linewidth', 2)

hold off

%title('Value function iteration'); 
xlabel('K_t'); ylabel('V(K_t)');
legend('V_1', 'V_{10}', 'V_{50}', 'V_{100}','V exact')



%{

error_ = max(v_improved - v_primitive)

figure(2)
%value function
x = k_state; y = v; plot(x,y);
title('value function'); xlabel('k_t'); ylabel('v(k)');

k_conv_path = zeros(1,1000);
k_conv_path(1) = k_initial;
t = 0; h = 1;
while h > 0.0001
	k_current = k_conv_path(t+1)
	k_next = interp1(k_state,k_pol,k_current,'linear')
	t = t+1
	k_conv_path(t+1) = k_next;
	h = abs(k_next - k_ss)/k_ss 
end

figure(3)
plot(0:t,k_conv_path(1:t+1))

% (b)

A = 1/(1-beta)*( alpha*beta/(1-alpha*beta)*log(alpha*beta) + log(1-alpha*beta) )
B = alpha/(1-alpha*beta)

pf = @(k) beta*B*k.^alpha/(1+beta*B)
vf = @(k) A + B*log(k)

figure(4)
k_ = phi:10^(-6):2*k_ss;
subplot(1,2,1)
plot(k_state,k_pol,k_,pf(k_))
legend('numerical solution','exact solution')
title('policy function'); xlabel('k_t'); ylabel('k_{t+1}');
subplot(1,2,2)
plot(k_state,v,k_,vf(k_))
legend('numerical solution','exact solution')
title('value function'); xlabel('k_t'); ylabel('v(k)');

c_current = k_state.^alpha - k_pol;
c_next = ( k_pol.^alpha - interp1(k_state,k_pol,k_pol,'linear') );
c_implied_by_euler = c_next./(beta*alpha*k_pol.^(alpha-1));
error_c = c_implied_by_euler./c_current - 1;
figure(5)
plot(k_state,error_c*100)
title('relative error of c'); xlabel('k'); ylabel('percentage');

% How to improve accuracy?
% expand the state space

% Where do these error come from?
% discretization of the state space

%}

