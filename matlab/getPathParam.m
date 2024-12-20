function [Si,So,Sid,Sod,Sidd,Sodd] = getPathParam(flag,inner,outer,nf,svec);
if flag == 1
    [Si,Sid,Sidd] = cubicfit(inner,svec);
    [So,Sod,Sodd] = cubicfit(outer,svec);
elseif flag ==2 
    [Si,Sid,Sidd] = fourierfit(inner,svec,nf);
    [So,Sod,Sodd] = fourierfit(outer,svec,nf);
else 
    svec= svec';
    aC = 3;
    bC = 4;
    xC = aC*cos(2*pi*svec);
    yC = bC*sin(2*pi*svec);

    xCd = -2*pi*aC*sin(2*pi*svec);
    yCd = 2*pi*bC*cos(2*pi*svec);

    xCdd = -2*pi*2*pi*aC*cos(2*pi*svec);
    ycdd = -2*pi*2*pi*bC*sin(2*pi*svec);


    aInner = aC -1 ;
    bInner = bC - 1;

    aOutter = aC +1 ;
    bOutter = bC + 1;

    centerline = [xC,yC];
    innerline = [aInner*cos(2*pi*svec),bInner*sin(2*pi*svec)];
    outterline =[aOutter*cos(2*pi*svec),bOutter*sin(2*pi*svec)];

    Si = innerline;
    Sid = [-2*pi*aInner*sin(2*pi*svec),2*pi*bInner*cos(2*pi*svec)];
    Sidd = [-2*pi*2*pi*aInner*cos(2*pi*svec),-2*pi*2*pi*bInner*sin(2*pi*svec)];
    So = outterline;
    Sod = [-2*pi*aOutter*sin(2*pi*svec),2*pi*bOutter*cos(2*pi*svec)];
    Sodd = [-2*pi*2*pi*aOutter*cos(2*pi*svec),-2*pi*2*pi*bOutter*sin(2*pi*svec)];
end

    function [S,Sd,Sdd] = cubicfit(waypoints,svec)
        for i = 1 : 2
            S_fit = spline(svec, waypoints(:,i));
            Sieval = ppval(S_fit, svec);
            % First derivatives (first-order derivatives)
            Sideval =  ppval(fnder(S_fit, 1), svec);
            % Second derivatives (second-order derivatives)
            Siddeval = ppval(fnder(S_fit, 2), svec);
            S(:,i) = Sieval;
            Sd(:,i) = Sideval;
            Sdd(:,i) = Siddeval;
        end
    end


    function [S,Sd,Sdd] = fourierfit(waypoints,svec,order)
        s =sym('s');

        Xi = waypoints(:,1)';
        Yi = waypoints(:,2)';
        N = length(waypoints);
        L = svec(end);
        ds = (svec(2)-svec(1))/L;
        Ax0 = 2/L*sum(Xi.*ones(size(svec)))*ds;
        fFSx = Ax0/2;

        Ay0 = 2/L*sum(Yi.*ones(size(svec)))*ds;
        fFSy = Ay0/2;
        Fxeq =fFSx;
        Fyeq =fFSy;

        for k=1:order
            Ax(k) = 2/L*sum(Xi.*cos(2*pi*k*svec/L))*ds;
            Bx(k) = 2/L*sum(Xi.*sin(2*pi*k*svec/L))*ds;
            fFSx = fFSx + Ax(k)*cos(k*2*pi*svec/L) + Bx(k)*sin(k*2*pi*svec/L);

            Ay(k) = 2/L*sum(Yi.*cos(2*pi*k*svec/L))*ds;
            By(k) = 2/L*sum(Yi.*sin(2*pi*k*svec/L))*ds;
            fFSy = fFSy + Ay(k)*cos(k*2*pi*svec/L) + By(k)*sin(k*2*pi*svec/L);

            Fxeq = Fxeq +  Ax(k)*cos(k*2*pi*s/L) + Bx(k)*sin(k*2*pi*s/L);
            Fyeq = Fyeq +  Ay(k)*cos(k*2*pi*s/L) + By(k)*sin(k*2*pi*s/L);
        end

        FxeqDot = diff(Fxeq,s);
        FxeqDDot = diff(FxeqDot,s);

        FyeqDot = diff(Fyeq,s);
        FyeqDDot = diff(FyeqDot,s);
        S = [double(subs(Fxeq,s,svec))',double(subs(Fyeq,s,svec))'];
        Sd = [double(subs(FxeqDot,s,svec))',double(subs(FyeqDot,s,svec))'];
        Sdd = [double(subs(FxeqDDot,s,svec))',double(subs(FyeqDDot,s,svec))'];
    end


end
