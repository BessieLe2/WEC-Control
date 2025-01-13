function z_bar=toZbar(x,u)
z_=[x;u];
z_bar=[];
for i = 1:length(z_)
    for j = i:length(z_)
        z_bar=[z_bar;z_(i)*z_(j)];

    end
end