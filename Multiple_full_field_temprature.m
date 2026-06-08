% /Calculate the full field temperature of hot and cold sides under multiple heat source/

clear
clc

a=0.5e-3;  % Half length of heat source
b=0.5e-3;  % Half width of heat source
he=0.2e-3;  % Thickness of encapsulation layer
hi=0.2e-3;  % Thickness of interface layer 
hs=1e-3;    %  Thickness of substrate layer 
ht=1.5e-3;  % Thickness of TEG layer
Q=100e-3;     % Heat power
ke=0.17;   % Thermal conductivity of encapsulation layer
ki=0.17; % Thermal conductivity of interface layer
ks=0.17;   % Thermal conductivity of substrate layer
kt=1; % Thermal conductivity of TEG layer
r=0.1*2*a:0.1*2*a:2*2*a; % Heat source array spacing
M=zeros(3,3); % Array number
N=zeros(3,3); % Array number


for i=1:3
    for j=1:3

        M(i,j)= j-2;
        N(i,j)= 2-i; % Array number

    end
end

lim=12e-3; % Calculation area (Treat as infinite boundary)
xy = -lim:1e-4:lim;
num = length(xy);
center_temp_hot = zeros(length(r),num);
center_temp_cold = zeros(length(r),num);

for idx=1:length(r)

    % iterate heat sources
    for k=1:3
        for q=1:3

            I = M(k,q);
            J = N(k,q); 

            % Calculate full-field temperature
            for i=1:num


                f = @(x,y)2.*Q.*sin(a.*x).*sin(b.*y).*cos((xy(i)+(r(idx)+2*a)*I).*x).*cos((r(idx)+2*a)*J.*y)./(pi()^2.*a.*b.*x.*y.*sqrt(x.^2+y.^2).*((exp(2*sqrt(x.^2+y.^2).*hi))-1)) ...
                ./((exp(2.*sqrt(x.^2+y.^2).*ht).*(kt.*coth(sqrt(x.^2+y.^2).*he)+ke)./(kt.*coth(sqrt(x.^2+y.^2).*he)-ke)+1)...
                .*(ki+ks*coth(sqrt(x.^2+y.^2)*hs).*coth(sqrt(x.^2+y.^2)*hi))+...
                (exp(2*sqrt(x.^2+y.^2)*ht).*(kt*coth(sqrt(x.^2+y.^2)*he)+ke)./(kt*coth(sqrt(x.^2+y.^2)*he)-ke)-1)...
                .*(kt*coth(sqrt(x.^2+y.^2)*hi)+kt*ks/ki*coth(sqrt(x.^2+y.^2)*hs)))...
                .*((exp(2*sqrt(x.^2+y.^2)*ht).*(kt*coth(sqrt(x.^2+y.^2)*he)+ke)./(kt*coth(sqrt(x.^2+y.^2)*he)-ke))+1)...
                .*exp(sqrt(x.^2+y.^2)*(hi)); %  hot side
                 temp=integral2(f,0,10e3,0,10e3);
                center_temp_hot(idx,i) = center_temp_hot(idx,i)+temp; 


                f = @(x,y)2.*Q.*sin(a.*x).*sin(b.*y).*cos((xy(i)+(r(idx)+2*a)*I).*x).*cos((r(idx)+2*a)*J.*y)./(pi()^2.*a.*b.*x.*y.*sqrt(x.^2+y.^2).*((exp(2*sqrt(x.^2+y.^2).*hi))-1)) ...
                ./((exp(2.*sqrt(x.^2+y.^2).*ht).*(kt.*coth(sqrt(x.^2+y.^2).*he)+ke)./(kt.*coth(sqrt(x.^2+y.^2).*he)-ke)+1)...
                .*(ki+ks*coth(sqrt(x.^2+y.^2)*hs).*coth(sqrt(x.^2+y.^2)*hi))+...
                (exp(2*sqrt(x.^2+y.^2)*ht).*(kt*coth(sqrt(x.^2+y.^2)*he)+ke)./(kt*coth(sqrt(x.^2+y.^2)*he)-ke)-1)...
                .*(kt*coth(sqrt(x.^2+y.^2)*hi)+kt*ks/ki*coth(sqrt(x.^2+y.^2)*hs)))...
                .*(exp(2.*sqrt(x.^2+y.^2).*ht).*(kt.*coth(sqrt(x.^2+y.^2).*he)+ke)./(kt.*coth(sqrt(x.^2+y.^2).*he)-ke).*exp(sqrt(x.^2+y.^2)*(hi-ht))...
                +exp(sqrt(x.^2+y.^2)*(hi+ht))); % cold side
                temp=integral2(f,0,5e3,0,5e3);
                center_temp_cold(idx,i) = center_temp_cold(idx,i)+temp; 



            end
        end
    end
end



