function [Vvopt_mc,Vhopt_mc,Dopt_mc] = diameteroptimization(X,Yv,Yh,Z_Th,Z_Peh,Z_Eh,Z_ethah,Z_Tv,Z_Pev,Z_Ev,Z_ethav)
       
Ndiam=length(X);
%MVel=length(max(size(Yv),size(Yh)));
%size(Z_Th);
%size(Z_Tv);
%% Calculate the maximun values for every variable as a function of speed 
Z_Thprima=max(Z_Th);
Z_ethahprima=max(Z_ethah);
Z_Pehprima2=max(Z_Peh);
Z_Ehprima2=max(Z_Eh);

Z_Tvprima=max(Z_Tv);
Z_ethavprima=max(Z_ethav);
Z_Pevprima2=max(Z_Pev);
Z_Evprima2=max(Z_Ev);


%% Calculate region of EPower and Energy that is zero, as it won't be useful to minimization problem (not defined zone)
Z_Peh_zeros=Z_Peh==0;
Z_Eh_zeros=Z_Eh==0;

Z_Pev_zeros=Z_Pev==0;
Z_Ev_zeros=Z_Ev==0;

%% Given the boolean matrices of zeros, get rid of them to calculate minimum
%Multiplying by maximum

Z_Pehprima=min(Z_Peh+Z_Peh_zeros.*Z_Pehprima2);
Z_Ehprima=min(Z_Eh+Z_Eh_zeros.*Z_Ehprima2);

Z_Pevprima=min(Z_Pev+Z_Pev_zeros.*Z_Pevprima2);
Z_Evprima=min(Z_Ev+Z_Ev_zeros.*Z_Evprima2);



%% Normalize matrices with the max(or min) value for a certain Diameter


Z_Thnorm=Z_Th./Z_Thprima;
Z_Pehnorm= Z_Peh./Z_Pehprima;
Z_Ehnorm= Z_Eh./Z_Ehprima;
Z_ethahnorm =Z_ethah./Z_ethahprima ;

Z_Tvnorm=Z_Tv./Z_Tvprima;
Z_Pevnorm=Z_Pev./Z_Pevprima;
Z_Evnorm=Z_Ev./Z_Evprima;
Z_ethavnorm=Z_ethav./Z_ethavprima;

ihest=zeros(Ndiam,1);
ivest=zeros(Ndiam,1);

for j=1:Ndiam
    
    %Calculate for each diameter the min squared solution normalized by its max/min. 
    
    %As parenthesis get closer to 1, the lower the result. That is what we
    %look for (closer to the local and relative extrema)
    
    [~,ihest(j)]=min((Z_Thnorm(:,j)-1).^2+(Z_Pehnorm(:,j)-1).^2+(Z_Ehnorm(:,j)-1).^2+(Z_ethahnorm(:,j)-1).^2);
    
    [~,ivest(j)]=min((Z_Tvnorm(:,j)-1).^2+(Z_Pevnorm(:,j)-1).^2+(Z_Evnorm(:,j)-1).^2+(Z_ethavnorm(:,j)-1).^2);

    
end

%% For the maximum for each colum (each D), from velocity index (ihest,ivest); normalize all column-wise local max/min with respect to the global
%    Z_Threl=Z_Thnorm(ihest,:).*Z_Thprima; %Thrust
%    Z_Thast=Z_Threl/max(Z_Threl);
indexZ_Thast=sub2ind(size(Z_Th),ihest',1:Ndiam); %Thrust
Z_Thast=Z_Th(indexZ_Thast)/max(Z_Th(:));

indexZ_Pehast=sub2ind(size(Z_Peh),ihest',1:Ndiam); %Epower
Z_Pehast=Z_Peh(indexZ_Pehast)/max(Z_Peh(:));

indexZ_Ehast=sub2ind(size(Z_Eh),ihest',1:Ndiam); %Energy
Z_Ehast=Z_Eh(indexZ_Ehast)/max(Z_Eh(:));

indexZ_ethahast=sub2ind(size(Z_ethah),ihest',1:Ndiam); %Eta
Z_ethahast=Z_ethah(indexZ_ethahast)/max(Z_ethah(:));
%%%%%%%%%%%%%%%%%%%%
indexZ_Tvast=sub2ind(size(Z_Tv),ivest',1:Ndiam); %Thrust
Z_Tvast=Z_Tv(indexZ_Tvast)/max(Z_Tv(:));

indexZ_Pevast=sub2ind(size(Z_Pev),ivest',1:Ndiam);  %Epower
Z_Pevast=Z_Pev(indexZ_Pevast)/max(Z_Pev(:));

indexZ_Evast=sub2ind(size(Z_Ev),ivest',1:Ndiam); %Energy
Z_Evast=Z_Ev(indexZ_Evast)/max(Z_Ev(:));

indexZ_ethavast=sub2ind(size(Z_ethav),ivest',1:Ndiam); %Eta
Z_ethavast=Z_ethav(indexZ_ethavast)/max(Z_ethav(:));


   
%% Solve now for D though the minimun squared technique 
[~,jest]=min((Z_Thast-1).^2+(Z_Pehast-1).^2+(Z_Ehast-1).^2+(Z_ethahast-1).^2+(Z_Tvast-1).^2+(Z_Pevast-1).^2+(Z_Evast-1).^2+(Z_ethavast-1).^2);

%% Particularize for the index I got

Dopt_mc=X(jest);
Vvopt_mc=Yv(ivest(jest));
Vhopt_mc=Yh(ihest(jest));

       
            
            
            