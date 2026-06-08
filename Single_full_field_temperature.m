% /Calculate the full field temperature of hot and cold sides according to Eq.(10)/

clear
clc

a=0.5e-3;  % Half length of heat source
b=0.5e-3;  % Half width of heat source
he=0.2e-3;  % Thickness of encapsulation layer
hi=0.2e-3;  % Thickness of interface layer 
hs=1e-3;    %  Thickness of substrate layer 
ht=1.5e-3;  % Thickness of TEG layer
Q=300e-3;     % Heat power
ke=0.17;   % Thermal conductivity of encapsulation layer
ki=0.17; % Thermal conductivity of interface layer
ks=0.17;   % Thermal conductivity of substrate layer
kt=1; % Thermal conductivity of TEG layer
lim=12e-3; % Calculation area (Treat as infinite boundary)
xy = -lim:1e-4:lim;
num = length(xy);
temp_ass=zeros(num,num);
[mesh_x,mesh_y]=meshgrid(xy,xy);

xxx=[];
yyy=[];
for i=1:num

    for j=1:num

        xxx = [xxx; xy(i)]; % x 
        yyy = [yyy; xy(j)]; % y

    end

end

% Full-field temperature of hot side
Hot_temp=[];
for i=1:num

    for j=1:num

            f = @(x,y)2.*Q.*sin(a.*x).*sin(b.*y).*cos(xy(j).*y).*cos(xy(i).*x)./(pi()^2.*a.*b.*x.*y.*sqrt(x.^2+y.^2).*((exp(2*sqrt(x.^2+y.^2).*hi))-1)) ...
                ./((exp(2.*sqrt(x.^2+y.^2).*ht).*(kt.*coth(sqrt(x.^2+y.^2).*he)+ke)./(kt.*coth(sqrt(x.^2+y.^2).*he)-ke)+1)...
                .*(ki+ks*coth(sqrt(x.^2+y.^2)*hs).*coth(sqrt(x.^2+y.^2)*hi))+...
                (exp(2*sqrt(x.^2+y.^2)*ht).*(kt*coth(sqrt(x.^2+y.^2)*he)+ke)./(kt*coth(sqrt(x.^2+y.^2)*he)-ke)-1)...
                .*(kt*coth(sqrt(x.^2+y.^2)*hi)+kt*ks/ki*coth(sqrt(x.^2+y.^2)*hs)))...
                .*((exp(2*sqrt(x.^2+y.^2)*ht).*(kt*coth(sqrt(x.^2+y.^2)*he)+ke)./(kt*coth(sqrt(x.^2+y.^2)*he)-ke))+1)...
                .*exp(sqrt(x.^2+y.^2)*(hi));
             temp=integral2(f,0,10e3,0,10e3);
            Hot_temp = [Hot_temp;temp];

    end

end


% Full-field temperature of cold side
Cold_temp=[];
for i=1:num

    for j=1:num

        f = @(x,y)2.*Q.*sin(a.*x).*sin(b.*y).*cos(xy(j).*y).*cos(xy(i).*x)./(pi()^2.*a.*b.*x.*y.*sqrt(x.^2+y.^2).*((exp(2*sqrt(x.^2+y.^2).*hi))-1)) ...
        ./((exp(2.*sqrt(x.^2+y.^2).*ht).*(kt.*coth(sqrt(x.^2+y.^2).*he)+ke)./(kt.*coth(sqrt(x.^2+y.^2).*he)-ke)+1)...
        .*(ki+ks*coth(sqrt(x.^2+y.^2)*hs).*coth(sqrt(x.^2+y.^2)*hi))+...
        (exp(2*sqrt(x.^2+y.^2)*ht).*(kt*coth(sqrt(x.^2+y.^2)*he)+ke)./(kt*coth(sqrt(x.^2+y.^2)*he)-ke)-1)...
        .*(kt*coth(sqrt(x.^2+y.^2)*hi)+kt*ks/ki*coth(sqrt(x.^2+y.^2)*hs)))...
        .*(exp(2.*sqrt(x.^2+y.^2).*ht).*(kt.*coth(sqrt(x.^2+y.^2).*he)+ke)./(kt.*coth(sqrt(x.^2+y.^2).*he)-ke).*exp(sqrt(x.^2+y.^2)*(hi-ht))...
        +exp(sqrt(x.^2+y.^2)*(hi+ht)));
        temp=integral2(f,0,5e3,0,5e3);
        Cold_temp = [Cold_temp;temp];

    end

end

