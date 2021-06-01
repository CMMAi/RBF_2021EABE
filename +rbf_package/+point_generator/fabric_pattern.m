function [x,y]=fabric_pattern(n, alpha)
    
    b = round(alpha*sqrt(n));
    phi = (sqrt(5)+1)/2;
    x=zeros(n,1); y=x;
    for k=1:n
        r = radius(k,n,b);
        theta = 2*pi*k/phi^2;
        x(k)=r*sin(theta); y(k)=r*cos(theta);

    end
end

function r = radius(k,n,b)
    if k>n-b
        r = 1;
    else
        r = sqrt(k-1/2)/sqrt(n-(b+1)/2);
    end
end