% This script plots basic figures, including q(eta), psi(eta),
% sigma^q(eta), CN(eta), S(eta), MU(eta)
% Platform: MATLAB R2019a

%% Plotting basic figures

figure(1)
subplot(3,2,1);  hold on
plot(Eta,Q);                 
ylabel('q');

subplot(3,2,2);  hold on
plot(Eta,Psi_e);             
ylabel('\psi');

subplot(3,2,3);  hold on
plot(Eta,SSQ - sigma);       
ylabel('\sigma^q');    
ylim([0 .12]);

subplot(3,2,4);  hold on
plot(Eta,CN_e);              
ylabel('expert consumption share');

subplot(3,2,5);  hold on
plot(Eta,S);  
ylabel('\sigma^\eta \eta');

subplot(3,2,6);  hold on
plot(Eta,MU);                
ylabel('\mu^\eta \eta');     
color = 'b';

figure(2)

subplot(2,2,1);  hold on
plot(Eta,Q,color);           
ylabel('q');                
xlabel('\eta');

subplot(2,2,2); hold on
plot(Eta,S,color);          
ylabel('\sigma^\eta \eta'); 
xlabel('\eta');

subplot(2,2,3);  hold on
plot(Eta,SSQ - sigma,color); 
ylabel('\sigma^q');         
xlabel('\eta'); 
ylim([0 .12])

subplot(2,2,4); hold on
plot(Eta,MU,color);          
ylabel('\mu^\eta \eta');    
xlabel('\eta');