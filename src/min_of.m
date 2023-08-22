function [x,y]=min_of(xx,yy,x1,y1,x2,y2,x_max)

x=[];
y=[];
i=1;
while xx(i)<x1
    x(end+1)=xx(i);
    y(end+1)=yy(i);
    i=i+1;
end

flag=0; % if Z_SI smaller：1， else -1
while xx(i)<x2
    if yy(i)>value_at(xx(i),[x1,x2],[y1,y2])   
        if flag==-1
            [is_intersect,x_inter,y_inter]=intersect_at(xx(i-1),yy(i-1),xx(i),yy(i),x1,y1,x2,y2);
            x=[x,x_inter];
            y=[y,y_inter];
        end
        flag=1;
    end

    if yy(i)<value_at(xx(i),[x1,x2],[y1,y2])
        if(flag==1)
            [is_intersect,x_inter,y_inter]=intersect_at(xx(i-1),yy(i-1),xx(i),yy(i),x1,y1,x2,y2);
            x=[x,x_inter,xx(i)];
            y=[y,y_inter,yy(i)];
        end
        flag=-1;
    end
    i=i+1;
    if i==length(xx)

    end
end

if value_at(x2,xx,yy)>y2
    if flag==-1
        [is_intersect,x_inter,y_inter]=intersect_at(xx(i-1),yy(i-1),xx(i),yy(i),x1,y1,x2,y2);
            x=[x,x_inter,x2];
            y=[y,y_inter,y2];
    end
    flag=1;
end
if value_at(x2,xx,yy)<y2
    if(flag==1)
        [is_intersect,x_inter,y_inter]=intersect_at(xx(i-1),yy(i-1),xx(i),yy(i),x1,y1,x2,y2);
        x=[x,x_inter,xx(i)];
        y=[y,y_inter,yy(i)];
    end
    flag=-1;
end
i=i+1;

while xx(i)<x_max

end


















