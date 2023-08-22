function y=value_at(x,xx,yy)
if x<xx(1) || x>xx(end)
    error("error");
end



for i=1:length(xx)
    if x==xx(i)
        y=yy(i);
    elseif x<xx(i)
        y=yy(i-1)+(yy(i)-yy(i-1))/(xx(i)-xx(i-1))*(x-xx(i-1));
        break;
    end
end


end