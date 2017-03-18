%Net Present Value
figure; 
plot([1:T],discount_factor); 
xlabel('Future year'); 
ylabel('Discounted value'); 
title([num2str(discount_rate*100),'% discount rate']); 
set(gcf,'color','white'); box off

PVpayoff_iy=Payoffiy.*discount_rate_iy; %present value  in each patch in each year
PVpayoff_y=sum(PVpayoff_iy,1); %present value of system in each year
NPVpayoff=sum(PVpayoff_y); %Net present value of system over time horizon

PVyield_iy=Yiy.*discount_rate_iy; %present value in each patch in each year
PVyield_y=sum(PVyield_iy,1); %present value of system in each year
NPVyield=sum(PVyield_y); %Net present value of system over time horizon

figure
plot(1:T,squeeze(sum(sum(Bijy,1),2)))
xlabel('Year')
ylabel('Biomass [kg]')
set(gcf,'color','white'); 

figure
plot(1:T,squeeze(sum(Bijy,1)))
xlabel('Year')
ylabel('Biomass [kg]')
set(gcf,'color','white'); 
% legend([num2str(1:length(Bijy(1,:,1)))'])

figure
plot(1:T,sum(Yiy,1))
xlabel('Year')
ylabel('Yield [kg]')
set(gcf,'color','white'); 

figure
plot(1:T,sum(Payoffiy,1))
xlabel('Year')
ylabel('Profit [$]')
set(gcf,'color','white'); 

figure
plot(1:T,PVpayoff_y)
xlabel('Year')
ylabel('Present value [$]')
set(gcf,'color','white'); 

figure
plot(1:T,PVpayoff_y)
xlabel('Year')
ylabel('Present value [$]')
set(gcf,'color','white'); 

figure
plot(1:T,PVyield_y)
xlabel('Year')
ylabel('Present value [yield]')
set(gcf,'color','white'); 
