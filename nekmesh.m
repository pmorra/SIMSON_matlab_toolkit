classdef nekmesh
    methods(Static)
        function [xm,ym] = meshinterp(xp,yp,xd,yd,varargin)
            %[xm,ym] = meshinterp(xp,yp,xd,yd,varargin)
            %varargin > 0 length is based on end of the grid in y dir
            sx = 0.*xp;
            for ii = 1 : size(xp,1)
                sx(ii,:) =[0 cumsum(sqrt(diff(xp(ii,:)).^2+diff(yp(ii,:)).^2))];
                sx(ii,:) = sx(ii,:)/sum(sx(ii,end));
            end
            sy = 0.*xp;
            for ii = 1 : size(xp,2)
                sy(:,ii) =[0;cumsum(sqrt(diff(xp(:,ii)).^2+diff(yp(:,ii)).^2))];
                sy(:,ii) = sy(:,ii)/sum(sy(end,ii));
            end
            %nekmesh.meshplot(sx,sy);
            tmp1 = cumsum(sqrt(diff(xd).^2+diff(yd).^2));
            s1x = [0;tmp1(:)];
            if ~isempty(varargin{1})
                tmp2 = cumsum(sqrt(diff(xp(:,end)).^2+diff(yp(:,end)).^2)) ;
                s1y = [0; tmp2(:)];
            else
                tmp2 = cumsum(sqrt(diff(xp(:,1)).^2+diff(yp(:,1)).^2)) ;
                s1y = [0; tmp2(:)];
            end
            s1x = s1x/s1x(end);
            s1y = s1y/s1y(end);
            
            
            [ssx,ssy]=meshgrid(s1x,s1y);
            
            xm = griddata(sx,sy,xp,ssx,ssy,'linear');
            ym = griddata(sx,sy,yp,ssx,ssy,'linear');
            
        end
        function e = norm(X,idim,varargin)
            % e = norm(X,idim,varargin); varargin = 1(2norm) 'fro' Frobenius
            if nargin == 2 
                varargin{1} = 2;
            end
            switch(idim)
                case  1
                    e = 0.*[1:size(X,2)];
                    for ii = 1 : size(X,2)
                        e(ii) = norm(X(:,ii),varargin{1});
                    end
                case 2
                    e = 0.*[1:size(X,1)];
                    for ii = 1 : size(X,1)
                        e(ii) = norm(X(ii,:),varargin{1});
                    end
                case 3                    
                    e = 0.*[1:size(X,3)];
                    for ii = 1 : size(X,3)
                        e(ii) = norm(squeeze(X(:,:,ii)),varargin{1});
                    end
                case 4
                    nn = size(X);
                    XX = reshape(X,prod(nn)/nn(end),nn(end));
                    clear X;
                    e = 0.*[1:nn(end)];
                    for ii = 1 : nn(end)
                        e(ii) = norm(squeeze(XX(:,ii)),varargin{1});
                    end
            end
        end
        function text(x,y,varargin)
            % text(x,y,tex,format)
            if nargin == 2
                tex = (1:length(x)).';
            elseif nargin ==3
                tex = varargin{1};
            else
                tex = varargin{1};
                fmt = varargin{2};
            end
            if(nargin <= 3)
                for ii = 1 : length(x)
                    text(x(ii),y(ii),num2str(tex(ii)));hold on
                end
            else
                for ii = 1 : length(x)
                    text(x(ii),y(ii),sprintf(fmt,tex(ii)));hold on
                end
            end
            
        end
        function [dd,t,nbc]=load(fname,varargin)
            switch nargin
                case 1
                    col=2;
                otherwise
                    col=varargin{1};
            end
            M=importdata(fname,' ',col);
            dd=M.data;
            if col==2;
                t=str2num(M.textdata{2});
                nbc=str2num(M.textdata{1});
            else
                t=0;
                nbc=0;
            end
        end
        function [ s, curv ] = curvture( x,y )
            %
            %CURVATURE This function computes curvature of an airfoil
            %defined by x and y arrays.
            %
            %[ s, curv] = curvture( x,y )
            %
            %s is the arclength.
            %
            
            s=[0;cumsum(sqrt(diff(x).^2+diff(y).^2))];
            
            [L1,L2]=diffop(s,length(s));
            xt0=L1*x;
            yt0=L1*y;
            xtt0=L2*x;
            ytt0=L2*y;
            curv=abs(xt0.*ytt0-yt0.*xtt0)./(xt0.*xt0+yt0.*yt0).^1.5;
            function [ L1, L2 ] = diffop( y,n )
                %This routine sets up arrays L1 and L2 which contain coefficients
                %     of a standard second order FD scheme for discretization of first
                %     and second derivatives, respectively.
                L1=spalloc(n,n,3*n);%zeros(n,n);
                L2=spalloc(n,n,3*n);%zeros(n,n);
                %L1=zeros(n,n);
                %L2=zeros(n,n);
                %
                % set up 2nd order FD coeff. at lower boundary
                %
                d1 = y(2)-y(1);
                d2 = y(3)-y(1);
                
                L1(1,1) = -(d1+d2)/d1/d2;
                L1(1,2) = -d2/d1/(d1-d2);
                L1(1,3) =  d1/d2/(d1-d2);
                
                L2(1,1) =  2.0/d1/d2;
                L2(1,2) =  2.0/d1/(d1-d2);
                L2(1,3) = -2.0/d2/(d1-d2);
                
                %
                %     set up 2nd order FD coeff. at iner points
                %
                
                for j=2:n-1
                    
                    d1 = y(j-1)-y(j);
                    d2 = y(j+1)-y(j);
                    
                    L1(j,j-1) = -d2/d1/(d1-d2);
                    L1(j,j) = -(d1+d2)/d1/d2;
                    L1(j,j+1) =  d1/d2/(d1-d2);
                    
                    L2(j,j-1) =  2.0/d1/(d1-d2);
                    L2(j,j) =  2.0/d1/d2;
                    L2(j,j+1) = -2.0/d2/(d1-d2);
                    
                end
                %
                %     set up 2nd order FD coeff. at upper boundary
                %
                
                d1 = y(n-1)-y(n);
                d2 = y(n-2)-y(n);
                
                L1(n,n-2) =  d1/d2/(d1-d2);
                L1(n,n-1) = -d2/d1/(d1-d2);
                L1(n,n) = -(d1+d2)/d1/d2;
                
                L2(n,n-2) = -2.0/d2/(d1-d2);
                L2(n,n-1) =  2.0/d1/(d1-d2);
                L2(n,n) =  2.0/d1/d2;
                
            end
        end
        function [xy,xyf,lmult]=wing(varargin) % NOT IN USE
            %[[xy,xyf,lmult]=wing(varargin)
            %varargin{1}>0 multiply eveything by lmult
            %xy and xyf are course and fine representaion of the wing
            %respectively;
            load('MW_135.mat')
            xy(:,1)=MW_135_punkte(:,1)/1000;
            xy(:,2)=MW_135_punkte(:,3)/1000;
            %correct data
            xy(1,2)=0;
            
            
            xxf=fliplr(xy(:,1)');
            yyf=fliplr(xy(:,2)');
            xxf=[ xxf(end) xxf];
            yyf=[ yyf(end) yyf];
            % lower side = 88.89 & upperside 50.37
            if nargin > 1
                is=find(xxf<1.2 ,1,'first');
                ie=find(xxf>0.680 & yyf > 0.0223,1,'first');
            else
                is=1;
                ie=length(xxf);
            end
            x1=xxf(is:ie);
            y1=yyf(is:ie);
            t=[0  cumsum(sqrt(diff(x1).^2+diff(y1).^2))];
            tt=linspace(0,t(end),2e6);
            xyf(:,1)=csaps(t,x1,1,tt);
            xyf(:,2)=csaps(t,y1,1,tt);
            cx=xyf(:,1);
            cy=xyf(:,2);
            xyf(:,3)=[0; cumsum(sqrt(diff(cx).^2+diff(cy).^2))]';
            
            %             sol=load('/scratch/dadfar/wing/solution.dat');
            %             cord=1.350000;
            %             data.xr=sol(:,1)+cord/2;
            %             data.yr=sol(:,2);
            %             data.ur=sol(:,3) ;
            %
            %             ind=find(abs(data.ur)<1e-1);
            %             xyard(:,1)=data.xr(ind);
            %             xyard(:,2)=data.yr(ind);
            if nargin > 0
                lmult =  82.329514741143313  ;
                xy=xy*lmult;
                xyf=xyf*lmult;
            else
                lmult=1;
            end
            if nargin>1
                %varargin{1}=lower and varargin{2} upperside
                ind=find(abs(xyf(:,1)-0)<=1e-8,1,'first');
                indc=0.*xyf(:,1);
                indc(1:ind)=0;
                indc(ind+1:end)=1;
                
                
                
                is=find(xyf(:,1)<=varargin{2}/100*1.35*lmult,1,'first');
                ie=find(xyf(:,1)>=varargin{3}/100*1.35*lmult & indc>0,1,'first');
                %                 is=1;
                %                 ie=length(xyf);
                
                
                ttf(:,1)=xyf(is:ie,1);
                ttf(:,2)=xyf(is:ie,2);
                ttf(:,3)=xyf(is:ie,3);
                clear xyf indc
                xyf=ttf;
                clear ttf xy
                n1=50;n2=30;n3=30;n4=10;
                %                 n1=80;n2=20;n3=20;n4=50;
                ind=find(abs(xyf(:,1)-0)<=1e-8,1,'first');
                indc=0.*xyf(:,1);
                indc(1:ind)=0;
                indc(ind+1:end)=1;
                
                is=find(xyf(:,1)<=3/100*1.35*lmult,1,'first');
                ie=find(xyf(:,1)>=3/100*1.35*lmult & indc>0,1,'first');
                
                xpts1=nekmesh.georefine(xyf(1,1),xyf(is,1),n1,0.94);
                ypts1=interp1(xyf(1:is,1),xyf(1:is,2),xpts1,'cubic');
                lpts1=interp1(xyf(1:is,1),xyf(1:is,3),xpts1,'cubic');
                
                lpts2=nekmesh.findR(xyf(is,3),xyf(ind,3),[0.9 1.2],n2,lpts1(end)-lpts1(end-1));
                xpts2=interp1(xyf(is:ind,3),xyf(is:ind,1),lpts2,'cubic');
                ypts2=interp1(xyf(is:ind,3),xyf(is:ind,2),lpts2,'cubic');
                
                lpts3=nekmesh.findR(xyf(ind,3),xyf(ie,3),[0.9 1.2],n3,lpts2(end)-lpts2(end-1));
                xpts3=interp1(xyf(ind:ie,3),xyf(ind:ie,1),lpts3,'cubic');
                ypts3=interp1(xyf(ind:ie,3),xyf(ind:ie,2),lpts3,'cubic');
                
                lpts4=nekmesh.findR(xyf(ie,3),xyf(end,3),[0.9 1.2],n4,lpts3(end)-lpts3(end-1));
                xpts4=interp1(xyf(ie:end,3),xyf(ie:end,1),lpts4,'cubic');
                ypts4=interp1(xyf(ie:end,3),xyf(ie:end,2),lpts4,'cubic');
                
                
                xy(:,1)=[xpts1 xpts2(2:end) xpts3(2:end) xpts4(2:end)];
                xy(:,2)=[ypts1 ypts2(2:end) ypts3(2:end) ypts4(2:end)];
                
            end
        end
        %****************************************************************
        function writerea2D(EL,ElmConn,FluidBC,pdd)
            Nelm=length(EL);
            Dim=2;
            %-------------------------------------------------------------------------
            %write core of the rea file
            %-----------------------------------------------------------------------
            
            disp(' ================================================================')
            disp('                       WRITING REA FILE                          ')
            disp(' ================================================================')
            fid = fopen('nekmesh.dat','w');
            format short
            
            %Write node/mesh data
            
            disp(' ')
            disp(['write ' num2str(Nelm) ' elments'])
            fprintf(fid,' **MESH DATA** 6 lines are X,Y,X,Y. Columns corners 1-4;5-8\n');
            fprintf(fid,'%10i%10i%10i NEL,NDIM,NELV \n',Nelm,Dim,Nelm);
            for e=1:Nelm;
                fprintf(fid,'      ELEMENT%11i [    1A]    GROUP     0\n',e);
                fprintf(fid,'%11.6f%15.6f%15.6f%15.6f\n',EL(e).nodes(1,:));
                fprintf(fid,'%11.6f%15.6f%15.6f%15.6f\n',EL(e).nodes(2,:));
            end
            %Write curved side data
            Ncurv = 0;
            Curve(1,1) = 0; %no curved sides are considered in this problem
            disp(' ')
            disp(['write ' num2str(Ncurv) ' curved sides'])
            if Curve(1,:) == 0
                n = 0;
            else
                n = length(Curve(:,1));
            end
            fprintf(fid,'  ***** CURVED SIDE DATA ***** \n');
            fprintf(fid,'%10i Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n',n);
            if n ~= 0
                for j = 1:length(Curve(:,1))  % 3D ONLY for now
                    if Nelm < 1000;
                        fprintf(fid,'%3i%3i%14.6f%14.6f%14.6f%14.6f%14.6f m\n',Curve(j,:));
                    else
                        fprintf(fid,'%2i%6i%14.6f%14.6f%14.6f%14.6f%14.6f m\n',Curve(j,:));
                    end
                end
            end
            %Write boundary conditions
            
            fprintf(fid,'  ***** BOUNDARY CONDITIONS ***** \n');
            fprintf(fid,'  ***** FLUID BOUNDARY CONDITIONS ***** \n');
            
            %write periodic boundary conditions
            C = ElmConn;
            for i = 1:length(ElmConn)
                switch FluidBC(i)
                    case 'E'
                        BCstring = 'E';
                    case 'W'
                        BCstring = 'W';
                    case 'O'
                        BCstring = 'O';
                    case 'N'
                        BCstring = 'ON';
                    case 'P'
                        BCstring = 'P';
                    case 'S'
                        BCstring = 'SYM';
                    case 'v'
                        BCstring = 'v';
                    case 'a'
                        BCstring = 'v';
                    case 'j'
                        BCstring = 'j';
                    otherwise
                        disp('Unknown BC')
                        break
                end
                if Nelm < 1000
                    fprintf(fid,[' ' BCstring '%5i%3i%10.5f%14.5f%14.5f%14.5f%14.5f\n'],[ElmConn(i,:) zeros(1,3)]);
                else
                    dummystring = sprintf('%.20f',C(i,3));
                    
                    fprintf(fid,' %-3s%5i%2i  %.8s%14.6f%14.6f%14.6f%14.6f\n',BCstring,...
                        C(i,1),C(i,2),dummystring,C(i,4),zeros(1,3));
                    
                end
            end
            
            fclose(fid);
            
            eval(['cd ' pdd]);
            disp('Generate complete rea-file')
            system('cat /scratch/nicolo/work/201407_Wing/reza/nlfgrid_rough/temp1.rea > nekmesh.rea');
            system('cat nekmesh.dat >> nekmesh.rea');
            system('cat /scratch/nicolo/work/201407_Wing/reza/nlfgrid_rough/temp2.rea >> nekmesh.rea');
            system('rm nekmesh.dat')
            disp(' ')
            disp(' ======================================  FINISHED')
            
        end
        %****************************************************************
        function [x,y] = viewgrid(prm, v)
            % VIEWGRID Visualises grids generated by `gridgen'.
            % Usage:
            %   viewgrid('<gridgen prm file>')       -- display grid only
            %   viewgrid('<gridgen prm file>', 'v')  -- display grid and image on one
            %                                           chart, vertically one under another
            %   viewgrid('<gridgen prm file>', 'h')  -- display grid and image on one
            %                                           chart, horizontally one aside other
            %   viewgrid('<gridgen prm file>', '<any other character than 'v' or 'h'>) --
            %                                           display grid and image on two
            %                                           different charts
            %
            %   For the image display to work, an entry `rectangle' must be
            %   specified in the `gridgen' parameter file.
            
            if (exist('v','var'))
                warning off;
                [geom,data,rect] = rfnames(prm);
                warning on;
                if (~exist('rect','var'))
                    fprintf('  Warning: the image file not specified in \"%s\"\n', prm);
                end
            else
                [geom,data] = rfnames(prm);
            end
            
            f = fopen(geom);
            count = 0;
            if (f >= 0)
                count = 0;
                str = fgetl(f);
                while (str >= 0)
                    if (str(1) ~= '#')
                        count = count + 1;
                        xy = sscanf(str, '%f %f %f');
                        xc(count) = xy(1);
                        yc(count) = xy(2);
                        if (length(xy) > 2)
                            zc(count) = xy(3);
                        else
                            zc(count) = 0;
                        end
                    end
                    str = fgetl(f);
                end
                
                xc_pos = xc(zc > 0);
                xc_neg = xc(zc < 0);
                yc_pos = yc(zc > 0);
                yc_neg = yc(zc < 0);
                
                if (xc(count) ~= xc(1) | yc(count) ~= yc(1))
                    count = count + 1;
                    xc(count) = xc(1);
                    yc(count) = yc(1);
                end
                fclose(f);
            end
            
            f = fopen(data);
            if (f > 0)
                count = 0;
                str = fgetl(f);
                while (str >= 0)
                    if (str(1) ~= '#')
                        count = count + 1;
                        xy = sscanf(str, '%f %f');
                        x(count) = xy(1);
                        y(count) = xy(2);
                    end
                    str = fgetl(f);
                end
                fclose(f);
            end
            
            %figure1 = figure;
            if (exist('rect','var'))
                if (exist('v','var') & v == 'v')
                    %subplot(2,1,2);
                elseif (exist('v','var') & v == 'h')
                    %subplot(1,2,2);
                else
                    %figure2 = figure;
                end
                r = load(rect);
                xr = r(:, 1);
                yr = r(:, 2);
                n = length(xr);
                xr(n + 1) = xr(1);
                yr(n + 1) = yr(1);
                %plot(xr, yr, 'k-', xr, yr, 'k.');
                %     title('Image');
                %     axis equal;
                %     xmin = min(xr);
                %     xmax = max(xr);
                %     ymin = min(yr);
                %     ymax = max(yr);
                %     dx = (xmax - xmin) * 0.05;
                %     dy = (ymax - ymin) * 0.05;
                %     axis([(xmin - dx) (xmax + dx) (ymin - dy) (ymax + dy)]);
                %     grid on;
            end
            
            %figure(figure1);
            if (exist('rect','var'))
                if (exist('v','var') & v == 'v')
                    %subplot(2,1,1);
                elseif (exist('v','var') & v == 'h')
                    %subplot(1,2,1);
                end
            end
            if (exist('xc','var'))
                if (length(xc_neg) > 0)
                    %plot(x, y, '.', xc, yc, 'k-', xc, yc, 'k.', xc_pos, yc_pos, 'go', xc_neg, yc_neg, 'rd');
                    %legend('grid nodes', 'boundary', 'boundary vertices', 'positive turn points', 'negative turn points');
                else
                    %plot(x, y, '.', xc, yc, 'k-', xc, yc, 'k.', xc_pos, yc_pos, 'go');
                    %legend('grid nodes', 'boundary', 'boundary vertices', 'turn points');
                end
            else
                %plot(x, y, '.');
            end
            %title(strcat('Original polygon and generated nodes ("', data, '")'));
            %axis equal;
            %xmin = min(xc);
            %xmax = max(xc);
            %ymin = min(yc);
            %ymax = max(yc);
            %dx = (xmax - xmin) * 0.05;
            %dy = (ymax - ymin) * 0.05;
            %axis([(xmin - dx) (xmax + dx) (ymin - dy) (ymax + dy)]);
            
            return;
            
            function [geom,data,rect] = rfnames(prm);
                f = fopen(prm);
                if (f < 0)
                    return;
                end
                
                seps = char([32,92]);
                str = fgetl(f);
                while (str >= 0)
                    if (str(1) ~= '#')
                        [first,str] = strtok(str, seps);
                        if (strcmpi(first, 'input'))
                            geom = strtok(str, seps);
                        elseif (strcmpi(first, 'output'))
                            data = strtok(str, seps);
                        elseif (strcmpi(first, 'rectangle'))
                            rect = strtok(str, seps);
                        end
                    end
                    str = fgetl(f);
                end
                
                fclose(f);
            end
        end
        %****************************************************************
        function [EL,ElmConn,FluidBC]=readrea(filename,varargin)
            %[EL,ElmConn,FluidBC]=readrea(filename,varargin)
            %varargin{1} folder default ./
            if nargin==1
                dir='./';
            else
                dir=varargin{1};
            end
            pdd=pwd;
            eval(['cd ' dir])
            %             filename=[filename '.rea'];
            fid=fopen(filename,'r+');
            
            while  ~feof(fid)
                tline = fgetl(fid);
                if(strcmp(tline(1:14),' **MESH DATA**') || strcmp(tline(1:14),'**MESH DATA**'))
                    tline = fgetl(fid);
                    dum = sscanf(tline,'%f');
                    if dum(2) == 3 ; dum(2) = 6 ; end
                    nn=1;
                    for i=1:dum(1)
                        tline=fgetl(fid);
                        for j=1:dum(2)
                            tline=fgetl(fid);
                            info=sscanf(tline,'%f');
                            EL(nn).nodes(j,:)=info';
                        end
                        nn=nn+1;
                    end
                end
                if(strcmp(tline(1:13),'  ***** FLUID'))
                    ElmConn =zeros(dum(1)*dum(2)*2,dum(2)*2);
                    FluidBC(1:dum(1)*dum(2)*2) = 'E';
                    if (dum(1) > 1000)
                        disp('element is greater than 1000')
                        for k=1:dum(1)*dum(2)*2
                            tline=fgetl(fid);
                            con=sscanf(tline(6:end),'%f');
                            ElmConn(k,1)=con(1);
                            ElmConn(k,2)=con(2);
                            ElmConn(k,3)=con(3);
                            ElmConn(k,4)=con(4);
                            FluidBC(k)=tline(2);
                        end
                    elseif(dum(1)<1000)
                        disp('element is less than 1000')
                        for k=1:dum(1)*dum(2)*2
                            tline=fgetl(fid);
                            con=sscanf(tline(6:end),'%f');
                            ElmConn(k,1)=con(1);
                            ElmConn(k,2)=con(2);
                            ElmConn(k,3)=con(3);
                            ElmConn(k,4)=con(4);
                            FluidBC(k)=tline(2);
                        end
                    end
                    if(strcmp(tline(1:13),'  ***** DRIVE'))
                        break
                    end
                    
                end
                
                eval(['cd ' pdd])
                
            end
        end
        function plotrea(EL)
            %plotrea(EL) ; just plot the xy mesh
            for ii=1:length(EL)
                for jj=1:4
                    x(jj)=EL(ii).nodes(1,jj);
                    y(jj)=EL(ii).nodes(2,jj);
                end
                plot([x(1) x(2)],[y(1) y(2)],'r');hold on
                plot([x(2) x(3)],[y(2) y(3)],'r');hold on
                plot([x(3) x(4)],[y(3) y(4)],'r');hold on
                plot([x(4) x(1)],[y(4) y(1)],'r');hold on
            end
            axis equal
        end
        function ELplot(EL,im,varargin)
            %ELplot(EL,im,varargin)
            %varargin the list of elemets to be ploted.
            is = im(1);
            ie = im(2);
            if ~isempty(varargin )
                iel = varargin{1};
                is = 1 ;
                ie = length(iel);
            end
            
            n = 1 ;
            for i = is : ie
                j = i;
                if ~isempty(varargin )
                    j = iel(i);
                end
                try
                    x(n) = mean(EL(j).nodes(1,:)) ;
                    y(n) = mean(EL(j).nodes(2,:)) ;
                end
                try
                    N1= sqrt(length(EL(1).GLL(:,1)));
                    x1 = EL(j).GLL(1,1) ;
                    y1 = EL(j).GLL(1,2) ;
                    
                    x2 = EL(j).GLL(N1,1) ;
                    y2 = EL(j).GLL(N1,2) ;
                    
                    x3 = EL(j).GLL(N1 * (N1 -1) + 1,1) ;
                    y3 = EL(j).GLL(N1 * (N1 -1) + 1,2) ;
                    
                    
                    x4 = EL(j).GLL(N1.^2 ,1) ;
                    y4 = EL(j).GLL(N1.^2 ,2) ;
                    
                    
                    xy = nekmesh.intersect([x1 y1],[x4 y4],[x2 y2],[x3 y3]) ;
                    x(n) = xy(1);
                    y(n) = xy(2);
                    %                     plot(x1, y1,'rs')
                    %                     hold on
                    %                     plot(x2, y2,'bs')
                    %                     plot(x3,y3,'gs')
                    %                     plot(x4,y4,'ms')
                    %                     plot(x,y,'rs')
                    %
                    %debug
                end
                hold on
                text(x(n),y(n),num2str(j)) ;
                n = n + 1 ;
            end
        end
        function BCplot_nodes(Eg,FluidBC,is,ie,type,varargin)
            % BCplot_nodes(Eg,FluidBC,is,ie,type,varargin{1})
            %varargin{1} shows you want to plot elements or not
            for iel=is:ie
                M=FluidBC(4*(iel-1)+1:4*iel);
                x(1)=(Eg(iel).nodes(1,1)+Eg(iel).nodes(1,2))/2;
                y(1)=(Eg(iel).nodes(2,1)+Eg(iel).nodes(2,2))/2;
                
                x(2)=(Eg(iel).nodes(1,2)+Eg(iel).nodes(1,3))/2;
                y(2)=(Eg(iel).nodes(2,2)+Eg(iel).nodes(2,3))/2;
                
                x(3)=(Eg(iel).nodes(1,3)+Eg(iel).nodes(1,4))/2;
                y(3)=(Eg(iel).nodes(2,3)+Eg(iel).nodes(2,4))/2;
                
                x(4)=(Eg(iel).nodes(1,4)+Eg(iel).nodes(1,1))/2;
                y(4)=(Eg(iel).nodes(2,4)+Eg(iel).nodes(2,1))/2;
                
                for ii=1:4
                    if strcmp(M(ii),type(1))
                     text(x(ii),y(ii),M(ii));hold on
                    end
                end
                if ~isempty(varargin{1})
                    if ~isempty(strfind(M,type(1)))
                        text(mean(x),mean(y),num2str(iel));hold on
                    end
                end
            end
        end
        %******************************************************************
        function [EL,ElmConn,box,Nelm] = conectivity(box)
            %Define faces
            F(1,:) = [1,2];
            F(2,:) = [2,3];
            F(3,:) = [4,3];
            F(4,:) = [1,4];
            
            
            elnum = 0;
            
            for n = 1:length(box)
                
                box(n).size = size(box(n).x);
                box(n).nelx = box(n).size(2)-1;
                box(n).nely = box(n).size(1)-1;
                box(n).elnumstart = elnum+1;
                box(n).nelz = 0;
                
                box(n).nodes = 1:1:(box(n).nelx+1)*(box(n).nely+1)*(box(n).nelz+1);
                box(n).nodes = reshape(box(n).nodes,box(n).nelx+1,box(n).nely+1,box(n).nelz+1);
                
                
                
                for j = 1:box(n).nely
                    for i = 1:box(n).nelx
                        elnum = elnum+1;
                        
                        EL(elnum).nodenum = [elnum box(n).nodes(i,j) box(n).nodes(i+1,j) box(n).nodes(i+1,j+1) box(n).nodes(i,j+1)];
                        
                        EL(elnum).nodes(1,:) = [box(n).x(j,i) box(n).x(j,i+1) box(n).x(j+1,i+1) box(n).x(j+1,i)];
                        EL(elnum).nodes(2,:) = [box(n).y(j,i) box(n).y(j,i+1) box(n).y(j+1,i+1) box(n).y(j+1,i)];
                        EL(elnum).BC = ['E' 'E' 'E' 'E'];
                        
                    end
                end
                
                
                elnum
            end
            Nelm = elnum;
            ElmConn = zeros(Nelm*4,4);
            for i = 1:Nelm
                ElmConn((i-1)*4+1:i*4,1)=ones(4,1)*i;
                ElmConn((i-1)*4+1:i*4,2)=[1:1:4];
            end
            
        end
        %******************************************************************
        function EL=gll(EL,N)
            Nelm=length(EL);
            n=N+1;
            gll = nekmesh.gllnodes(N);
            gllx=gll';
            glly=gll';
            for i=1:Nelm
                x= EL(i).nodes(1,:)';
                y= EL(i).nodes(2,:)';
                f1=zeros(n^2,1);
                f2=zeros(n^2,1);
                f3=zeros(n^2,1);
                f4=zeros(n^2,1);
                f1(:)=(1-gllx)*(1-glly)';
                f2(:)=(gllx)*(1-glly)';
                f3(:)=(gllx)*(glly)';
                f4(:)=(1-gllx)*(glly)';
                
                G=[f1 f2 f3 f4];
                
                EL(i).GLL=G*[x y];
            end
            
        end
        function [Eg,box,nnx,nny]=easygll(xx,yy,N1,varargin)
            if nargin>3
                box=nekmesh.gridbox(xx,yy,size(xx,2),size(xx,1));
                xy = nekmesh.getgeolead(20,1,6,[]);
                Arc.x=xy(:,1);
                Arc.y=xy(:,2);
                [Arc,box,npx,npy,np2]=nekmesh.arccalc(box,Arc,'spline');
                [E,~,box,~] = nekmesh.conectivity(box);
                Eg=nekmesh.gll(E,N1-1);
                Eg=nekmesh.gllcorrect(Eg,N1-1,Arc,npx,npy,np2);
                [nnx,nny]=nekmesh.meshinfo(Eg);
            else
                box(1).x=xx;
                box(1).y=yy;
                [E,~,box,~] = nekmesh.conectivity(box);
                Eg=nekmesh.gll(E,N1-1);
                [nnx,nny]=nekmesh.meshinfo(Eg);
            end
            for ii=1:length(Eg)
                Eg(ii).Vinit=Eg(ii).GLL;
                Eg(ii).fields='XU';
            end
            
        end
        function easyrea(Eg,varargin)
            %easyrea(Eg or box ,varargin)
            % Eg is a structure or box 
            %varagin{1} = 1 ; load from the BC files 
            %varargin{1} = ['W' 'O' 'O' 'v'];
            if ~isempty(cell2mat(strfind(fieldnames(Eg),'GLL')))
                [xtt,ytt]=nekmesh.meshmap(Eg);
                box(1).x=xtt;
                box(1).y=ytt;
            else
                box=Eg;
            end
            % nekmesh.meshplot(xtt,ytt,'r.',2)
            [EE,ElmConn,box,~] = nekmesh.conectivity(box);
            if nargin == 1 
                try
                    !mkdir info
                end
                FluidBC=nekmesh.Boundary(box,EE);
                save ('./info/FluidBC.mat','FluidBC');
            elseif (nargin>1 && ~iscell(varargin{1}))
                load('./info/FluidBC.mat')
            else
                bb = varargin{1};
                mx = size(box(1).x,2) - 1;
                my = size(box(1).x,1) - 1;
                
                FluidBC=reshape(char((ones(1,4*mx*my)*69)),4,mx*my);
                FluidBC(1,1:mx)=bb{1};
                FluidBC(2,mx:mx:mx*my)=bb{2};
                FluidBC(3,mx*(my-1)+1:mx*my)=bb{3};
                FluidBC(4,1:mx:mx*(my-1)+1)=bb{4};
                FluidBC=FluidBC(:);
                try
                    !mkdir info
                end
                save ('./info/FluidBC.mat','FluidBC');
            end
            nekmesh.writerea2D(EE,ElmConn,FluidBC,pwd);
            
        end
        function easyinit(fb,fi,casename,Eg,time,fii )
            %fb is the filename for extracting hte boundary condition and
            %fi is the filename for extracting the initial condition
            
            [~,~,N1]=nekmesh.meshinfo(Eg);
            
            if(isempty(fi))
                if(isempty(fb))
                    disp('two file names are empty!!!')
                else
                    fi=fb;
                end
            else
                if(isempty(fb))
                    fb=fi;
                end
            end
            
            
            
            ELr=nekmesh.readfile(fb);
            
            if ~isempty(fii)
                EE=nekmesh.readfile(fii);
                for ii=1:length(EE);
                    ELr(ii).GLL=EE(ii).GLL;
                    ELr(ii).fields=EE(ii).fields;
                end
            end
            
            %             for ii=1:length(ELr)
            %                 plot(ELr(ii).GLL(:,1),ELr(ii).GLL(:,2),'r.'); hold on
            %             end
            
            [nnx nny NN1]=nekmesh.meshinfo(ELr);
            [uu,vv,~,xx,yy]=nekmesh.results(ELr,nnx,nny,NN1);
            
            % pcolor(xx,yy,uu);shading interp
            % nekmesh.meshplot(xx,yy,'r.',2)
            % nekmesh.meshplot(box(1).x,box(1).y,'b',2)
            
            data.xr=xx(:);
            data.yr=yy(:);
            data.ur=uu(:);
            data.vr=vv(:);
            load('./info/FluidBC.mat')
            BC=nekmesh.calcBC(N1-1,FluidBC,Eg,data,'v');
            nekmesh.writeBC('./',casename,BC,'.BC');
            %**************************************************************************
            %                        initial condition
            %**************************************************************************
            clear ELr clear nnx nny NN1 uu vv xx yy
            
            ELr=nekmesh.readfile(fi);
            if ~isempty(fii)
                for ii=1:length(EE);
                    ELr(ii).GLL=EE(ii).GLL;
                    ELr(ii).fields=EE(ii).fields;
                end
            end
            
            [nnx nny NN1]=nekmesh.meshinfo(ELr);
            [uu,vv,~,xx,yy]=nekmesh.results(ELr,nnx,nny,NN1);
            
            Fu= TriScatteredInterp(xx(:),yy(:),uu(:));
            Fv= TriScatteredInterp(xx(:),yy(:),vv(:));
            
            for i=1:length(Eg)
                x=Eg(i).GLL(:,1);
                y=Eg(i).GLL(:,2);
                u=Fu(x,y);
                v=Fv(x,y);
                if(norm(double(isnan(u)))~=0);
                    disp(['you have nan in u IC el##' num2str(i)] );
                end
                if(norm(double(isnan(v)))~=0);disp('you have nan in v IC');end
                u(isnan(u))=1;
                v(isnan(v))=0;
                
                Eg(i).Vinit(:,1)=u;
                Eg(i).Vinit(:,2)=v;
                Eg(i).Pinit(:,1)=zeros(size(u));
            end
            
            nekmesh.writeic(Eg,N1-1,'./',casename,'.IC',[],time,8);
            
            
        end
        function [box,EL,elm,std,time,nnx,nny]=easyread(varargin)
            %[box,EL,elm,std]=easyread(varargin)
            %input: (1) input file name or EL
            %input: (2) optional : another file name with coordinate infomation
            %output:(1) box  (2) EL (3) elm:element information
            f=varargin{1};
            if ~isstruct(f)
                [EL,~,time,~,std]=nekmesh.myread(f);
            else
                EL=f;
                clear f;
            end
            if nargin<=1
                [nnx,nny,N1,~,elm]=nekmesh.meshinfo(EL);
            else
                [E1,~,time,~,std]=nekmesh.myread(varargin{2});
                [nnx,nny,N1,~,elm]=nekmesh.meshinfo(E1);
                b=nekmesh.results_wing(E1,nnx,nny,N1,elm);
            end
            box=nekmesh.results_wing(EL,nnx,nny,N1,elm);
            if length(varargin)>1
                for ii=1:length(box)
                    box(ii).x=b(ii).x;
                    box(ii).y=b(ii).y;
                    for jj=1:length(EL)
                        EL(jj).GLL=E1(jj).GLL;
                    end
                end
            end
            
        end
        %********************wing *****************************************
        function [Eg,box]=easygll_wing(box,N1,curv,xyf)
            np2=1;
            [Arc,box,npx,npy,len]=nekmesh.arccalc_wing(box,np2,xyf);
            [E,~,box,~] = nekmesh.conectivity(box);
            Eg=nekmesh.gll(E,N1-1);
            
            if ~isempty(curv)
                Eg=nekmesh.gllcorrect_wing(Eg,N1-1,Arc,npx,npy,np2-1,len);
            end
            %             [nnx,nny]=nekmesh.meshinfo(Eg(1:npx*npy));
            for ii=1:length(Eg)
                Eg(ii).Vinit=Eg(ii).GLL;
                Eg(ii).fields='XU';
                for jj=1:length(box)
                    Eg(ii).box(jj,:)=box(jj).size-1;
                end
            end
            
        end
        function box = easyrea_wing(Eg,lod)
            %easyrea_wing(Eg,lod)
            %lod is flag for FLUIDBC var. it saves it if it empty and
            %loads it if it is not empty
            
            elm=[];
            nn=0;
            for ii=1:length(fieldnames(Eg))
                dname=fieldnames(Eg);
                if isempty(strfind(dname(ii),'box'))
                    nn=nn+1;
                end
            end
            if nn==0
                [~,~,~,~,elm]=nekmesh.meshinfo(Eg);
            else
                for ii=1:size(Eg(1).box,1)
                    elm(ii)=Eg(1).box(ii,1)*Eg(1).box(ii,2);
                end
                elm=cumsum(elm);
                elm=[0 elm];
            end
            for ii=1:length(elm)-1
                [xtt,ytt]=nekmesh.meshmap(Eg(elm(ii)+1:elm(ii+1)));
                box(ii).x=xtt;
                box(ii).y=ytt;
                clear xtt ytt
            end
            % nekmesh.meshplot(box(1).x,box(1).y);
            [EE,ElmConn,box,Nelm] = nekmesh.conectivity(box);
            
            if isempty(lod)
                FluidBC=nekmesh.Boundary(box,EE);
                try
                    !mkdir info
                end
                save ('./info/FluidBC.mat','FluidBC');
            else
                load('./info/FluidBC.mat')
            end
            nekmesh.writerea2D(EE,ElmConn,FluidBC,pwd);
            
            
        end
        %******************************************************************
        function EL=gllcorrect(EL,N,Arc,npx,npy,np)
            %element which should be correct up to np element. it  is the number
            %which fc
            %proceed along i=const direction
            Nelm=length(EL);
            n=N+1;
            gll = nekmesh.gllnodes(N);
            gllx=gll';
            glly=gll';
            % 1 2 3 4 is not like element numbering
            ielem=1;
            for i=1:npy-1
                for j=1:npx
                    if j< np
                        
                        xy1=EL(ielem).GLL(1,:);
                        xy2=EL(ielem).GLL(n,:);
                        xy3=EL(ielem).GLL(n*(n-1)+1,:);
                        xy4=EL(ielem).GLL(n*n,:);
                        %----- length as a funciton of x
                        t1=interp1(Arc.x(i,:),Arc.L(i,:),xy1(1));
                        t2=interp1(Arc.x(i,:),Arc.L(i,:),xy2(1));
                        
                        s1=interp1(Arc.x(i+1,:),Arc.L(i+1,:),xy3(1));
                        s2=interp1(Arc.x(i+1,:),Arc.L(i+1,:),xy4(1));
                        if(sum(isnan([t1 t2 s1 s2]))~=0);disp(['you have nanin ' num2str(ielem)]);end
                        %------ dividing the length and x ,y as a funciotn of L
                        tt=t1+gllx*(t2-t1);
                        g1x=interp1(Arc.L(i,:),Arc.x(i,:),tt);
                        g1y=interp1(Arc.L(i,:),Arc.y(i,:),tt);
                        
                        ss=s1+gllx*(s2-s1);
                        g2x=interp1(Arc.L(i+1,:),Arc.x(i+1,:),ss);
                        g2y=interp1(Arc.L(i+1,:),Arc.y(i+1,:),ss);
                        st=1:n^2;
                        st=reshape(st,n,n);
                        
                        px=zeros(n^2,1);
                        py=zeros(n^2,1);
                        
                        for kk=1:n
                            px(st(kk,:))=g1x(kk)+glly*(g2x(kk)-g1x(kk));
                            py(st(kk,:))=g1y(kk)+glly*(g2y(kk)-g1y(kk));
                        end
                        EL(ielem).GLL(:,1)=px;
                        EL(ielem).GLL(:,2)=py;
                        
                        
                    end
                    ielem=ielem+1;
                end
            end
            
            %------ fix any possible error
            
            for i=1:npy
                st=(i-1)*npx+1:i*npx;
                stu=(i+1)*npx+1:(i+2)*npx;
                for j=1:length(st)
                    try
                        ie=st(j);
                        iep=st(j+1);
                        jep=stu(j);
                        
                        EL(ie).GLL(n:n:n^2,:) = EL(jep).GLL(1:n,:);
                        EL(ie).GLL(n:n:n^2,:) = EL(iep).GLL(1:n:n*(n-1)+1,:);
                        
                    end
                end
            end
            
            %-------- correcting for sym BCs
            st= 1:npx:(npy-1)*npx+1 ;
            for i=1:length(st)
                EL(st(i)).GLL(1:n:n*(n-1)+1,2)=0;
            end
            
            
            
            
        end
        function EL=gllcorrect_first_line(EL,Arc,nx,ny,nps,npe)
            npx=nx;
            npy=ny;
            %element which should be correct
            Nelm=length(EL);
            N=sqrt( length( EL(1).GLL(:,1))) - 1;
            n = N + 1;
            gll = nekmesh.gllnodes(N);
            gllx=gll';
            glly=gll';
            % 1 2 3 4 is not like element numbering
            ielem=1;
            i=1;
            for j=1:npx
                if (j<= npe && j>= nps)
                    xy1=EL(ielem).GLL(1,:);
                    xy2=EL(ielem).GLL(n,:);
                    xy3=EL(ielem).GLL(n*(n-1)+1,:);
                    xy4=EL(ielem).GLL(n*n,:);
                    %------------------------
                    %         3    4
                    %         1    2
                    %----- length as a funciton of x,t1 start of xy1
                    ind = find(Arc.x(i,:) == 0);
                    %---------------- correct the mesh point
                    if xy1(2)< 0
                        t1=interp1(Arc.x(i,1:ind),Arc.L(i,1:ind),xy1(1));
                        t2=interp1(Arc.x(i,1:ind),Arc.L(i,1:ind),xy2(1));
                    else
                        t1=interp1(Arc.x(i,ind:end),Arc.L(i,ind:end),xy1(1));
                        t2=interp1(Arc.x(i,ind:end),Arc.L(i,ind:end),xy2(1));
                        
                    end
                    
                    if(sum(double(isnan([t1 t2])))~=0);
                        disp(['you have nan @ ', num2str(ielem)]);
                    end
                    %------ dividing the length and x ,y as a funciotn of L
                    tt=t1+gllx*(t2-t1);
                    g1x=interp1(Arc.L(i,:),Arc.x(i,:),tt);
                    g1y=interp1(Arc.L(i,:),Arc.y(i,:),tt);
                    
                    
                    %g2x=xy3(1)+gllx*(xy4(1)-xy3(1));
                    %g2y=xy3(2)+gllx*(xy4(2)-xy3(2));
                    g2x = EL(ielem).GLL(n*(n-1)+1:n^2,1);
                    g2y = EL(ielem).GLL(n*(n-1)+1:n^2,2);
                    
                    
                    st=1:n^2;
                    st=reshape(st,n,n);
                    
                    px=zeros(n^2,1);
                    py=zeros(n^2,1);
                    
                    for kk=1:n
                        px(st(kk,:))=g1x(kk)+glly*(g2x(kk)-g1x(kk));
                        py(st(kk,:))=g1y(kk)+glly*(g2y(kk)-g1y(kk));
                    end
                    EL(ielem).GLL(:,1)=px;
                    EL(ielem).GLL(:,2)=py;
                    
                    
                    
                end
                ielem=ielem+1;
            end
            
            %------ fix any possible error
            
            for i=1:npy
                st=(i-1)*npx+1:i*npx;
                stu=(i+1)*npx+1:(i+2)*npx;
                for j=1:length(st)
                    try
                        ie=st(j);
                        iep=st(j+1);
                        jep=stu(j);
                        
                        EL(ie).GLL(n:n:n^2,:) = EL(jep).GLL(1:n,:);
                        EL(ie).GLL(n:n:n^2,:) = EL(iep).GLL(1:n:n*(n-1)+1,:);
                        
                    end
                end
            end
            
            %-------- correcting for sym BCs
            st= 1:npx:(npy-1)*npx ;
            for i=1:length(st)
                EL(st(i)).GLL(1:n:n*(n-1)+1,2)=0;
            end
            
            
            
            
        end
        function EL=gllcorrect_wing(EL,N,Arc,npx,npy,np2,len)
            %element which should be correct up to np element. it  is the number
            %which fc
            %proceed along i=const direction
            Nelm=length(EL);
            n=N+1;
            gll = nekmesh.gllnodes(N);
            gllx=gll';
            glly=gll';
            % 1 2 3 4 is not like element numbering
            ielem=1;
            for i=1:np2+1
                for j=1:npx
                    t1=len(i,j);
                    t2=len(i,j+1);
                    s1=len(i+1,j);
                    s2=len(i+1,j+1);
                    %------ dividing the length and x ,y as a funciotn of L
                    tt=t1+gllx*(t2-t1);
                    
                    g1x=interp1(Arc.L(:,i),Arc.x(:,i),tt);
                    g1y=interp1(Arc.L(:,i),Arc.y(:,i),tt);
                    
                    if(i==np2+1)
                        g1x=EL(ielem-npx).GLL(n*(n-1)+1:n^2,1);
                        g1y=EL(ielem-npx).GLL(n*(n-1)+1:n^2,2);
                    end
                    ss=s1+gllx*(s2-s1);
                    g2x=interp1(Arc.L(:,i+1),Arc.x(:,i+1),ss);
                    g2y=interp1(Arc.L(:,i+1),Arc.y(:,i+1),ss);
                    if(i==np2+1)
                        g2x=EL(ielem+npx).GLL(1:n,1);
                        g2y=EL(ielem+npx).GLL(1:n,2);
                    end
                    
                    
                    st=1:n^2;
                    st=reshape(st,n,n);
                    
                    px=zeros(n^2,1);
                    py=zeros(n^2,1);
                    
                    for kk=1:n
                        px(st(kk,:))=g1x(kk)+glly*(g2x(kk)-g1x(kk));
                        py(st(kk,:))=g1y(kk)+glly*(g2y(kk)-g1y(kk));
                    end
                    EL(ielem).GLL(:,1)=px;
                    EL(ielem).GLL(:,2)=py;
                    
                    
                    
                    ielem=ielem+1;
                    disp(['element # = ' num2str(ielem-1)]);
                end
                
            end
            %------ fix any possible error
            for i=1:npy
                st=(i-1)*npx+1:i*npx;
                stu=(i+1)*npx+1:(i+2)*npx;
                for j=1:length(st)
                    try
                        ie=st(j);
                        iep=st(j+1);
                        jep=stu(j);
                        
                        EL(ie).GLL(n:n:n^2,:) = EL(jep).GLL(1:n,:);
                        EL(ie).GLL(n:n:n^2,:) = EL(iep).GLL(1:n:n*(n-1)+1,:);
                        
                    end
                end
            end
            
        end
        %******************************************************************
        function [x,w,P] = gllnodes(N)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % lglnodes.m
            %
            % Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde
            % matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
            % integration and spectral methods.
            %
            % Reference on LGL nodes and weights:
            %   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
            %   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
            %
            % Written by Greg von Winckel - 04/17/2004
            % Contact: gregvw@chtm.unm.edu
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Truncation+1
            N1 = N+1;
            
            % Use the Gauss-Lobatto-Chebyshev nodes as the first guess
            x = cos(pi*(0:N)/N)';
            
            % glc = x;
            % glc = 1-glc;
            % glc = 0.5*glc;
            
            % The Legendre Vandermonde Matrix
            P = zeros(N1,N1);
            
            % Compute P_(N) using the recursion relation
            % Compute its first and second derivatives and
            % update x using the Newton-Raphson method.
            
            xold = 2;
            
            while max(abs(x-xold))>eps
                
                xold = x;
                
                P(:,1) = 1;
                P(:,2) = x;
                
                for k = 2:N
                    P(:,k+1) = ( (2*k-1)*x.*P(:,k) - (k-1)*P(:,k-1) )/k;
                end
                
                x = xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
                
            end
            
            % Get the GLL points on the interval [0,1]
            x = 1-x;
            x = 0.5*x;
            
            % Integration weights for [0,1]
            w = 1./(N*N1*P(:,N1).^2);
            
            % Both x and w as row vector
            x = x';
            w = w';
            
            % figure,
            % plot(x,0.99.*ones(size(x)),'or',glc,1.01.*ones(size(glc)),'ok','MarkerSize',8);
            % set(gca,'DataAspectRatio',[2.5 1 1]);
            % axis off;
            % cd ~/
            % print('-dpng','-r300','GLL-GLC');
        end
        %******************************************************************
        function FluidBC=Boundary(box,EL)
            %Define faces
            F(1,:) = [1,2];
            F(2,:) = [2,3];
            F(3,:) = [4,3];
            F(4,:) = [1,4];
            
            
            for n = 1:length(box)
                color = ['b','g','m','y','k','r'];
                hbc(n) = figure(110+n);
                set(hbc(n),'UserData',[])
                hold on
                for i = 1:box(n).nely+1
                    plot(box(n).x(i,:),box(n).y(i,:),color(n))
                end
                for i = 1:box(n).nelx+1
                    plot(box(n).x(:,i),box(n).y(:,i),color(n))
                end
                axis equal
                grid on
                set(hbc(n),'WindowButtonDownFcn',{@defineBC,box(n).nelx+1,box(n).nely+1})
                evalresponse = input('Boundary conditions defined?');
                
                box(n).BCpoints = get(hbc(n),'UserData');
                box(n).BCtype = evalresponse;
                
                for i = 1:2:length(box(n).BCpoints)
                    
                    if box(n).BCpoints(i,2) == box(n).BCpoints(i+1,2)
                        signtest = sign(box(n).BCpoints(i+1,1)-box(n).BCpoints(i,1));
                        xind = box(n).BCpoints(i,1):signtest*1:box(n).BCpoints(i+1,1);
                        yind = box(n).BCpoints(i,2);
                        box(n).BC{ceil(i/2)} = squeeze(box(n).nodes(xind,yind));
                    else
                        signtest = sign(box(n).BCpoints(i+1,2)-box(n).BCpoints(i,2));
                        yind = box(n).BCpoints(i,2):signtest*1:box(n).BCpoints(i+1,2);
                        xind = box(n).BCpoints(i,1);
                        box(n).BC{ceil(i/2)} = squeeze(box(n).nodes(xind,yind))';
                    end
                    
                    %plot(box(n).BC{ceil(i/2)}(1,:),box(n).BC{ceil(i/2)}(2,:),'ks')
                end
                
                for elind = box(n).elnumstart:box(n).elnumstart+(size(box(n).x,1)-1)*(size(box(n).x,2)-1)-1
                    for bcind = 1:length(box(n).BCpoints)/2
                        for faceind = 1:size(F,1)
                            indtrue = ismember(EL(elind).nodenum(F(faceind,:)+1),box(n).BC{bcind});
                            testel = EL(elind).nodenum(indtrue);
                            
                            if length(testel)==2
                                EL(elind).BC(faceind) = box(n).BCtype(bcind);
                            end
                        end
                    end
                    FluidBC((elind-1)*size(F,1)+1:elind*size(F,1),1) = EL(elind).BC;
                end
            end
            function defineBC(src,evnt,npx,npy)
                
                if strcmp(get(src,'SelectionType'),'normal')
                    set(src,'pointer','circle')
                    cp = get(gca,'CurrentPoint');
                    xinit = cp(1,1);yinit = cp(1,2);
                    
                    ploth = get(gca,'Children');
                    S = length(ploth)-npy+1;
                    
                    is = 0;
                    for i = S:length(ploth)
                        is = is+1;
                        x(is,:) = get(ploth(i),'XData');
                        y(is,:) = get(ploth(i),'YData');
                    end
                    
                    %find closest point
                    dist = sqrt((x-xinit).^2 + (y-yinit).^2);
                    
                    [val,ind] = min(dist(:));
                    Sd = size(dist);
                    [indj,indi] = ind2sub(Sd,ind);
                    
                    indi1 = indi;
                    indj1 = S+indj-1;
                    x1v = get(ploth(indj1),'XData');
                    y1v = get(ploth(indj1),'YData');
                    x1 = x1v(indi1);
                    y1 = y1v(indi1);
                    
                    
                    S1 = length(ploth)-npy-npx+1;
                    indi2 = S1+npx-indi;
                    indj2 = npy-indj+1;
                    x2v = (get(ploth(indi2),'XData'));
                    y2v = (get(ploth(indi2),'YData'));
                    x2 = x2v(indj2);
                    y2 = y2v(indj2);
                    plot(x1,y1,'g*',x2,y2,'bo')
                    
                    %set(src,'WindowButtonMotionFcn',{@wbmcb,indi1,indj1,indi2,indj2})
                    indj = npy-indj+1;
                    set(src,'WindowButtonUpFcn',{@wbucb,indj,indi})
                end
                
                function wbucb(src,evnt,indj,indi)
                    points = get(src,'Userdata');
                    if isempty(points)
                        points(1,1) = indi;
                        points(1,2) = indj;
                    else
                        S = size(points);
                        points(S(1)+1,1) = indi;
                        points(S(1)+1,2) = indj;
                    end
                    set(src,'Userdata',points)
                    set(src,'pointer','arrow')
                    set(src,'WindowButtonUpFcn','')
                end
            end
            
        end
        %******************************************************************
        function BC2 = convertBC(BC,Ngp)
            %BC2 = convertBC(FluidBC,BC,Ngp)
            Fgll(4,:) = 1:Ngp:(Ngp-1)*Ngp+1;
            Fgll(1,:) = 1:1:Ngp;
            Fgll(2,:) = Ngp:Ngp:Ngp^2;
            Fgll(3,:) = (Ngp-1)*Ngp+1:1:Ngp*Ngp;
            pp = zeros(size(Fgll,1),1);
            for ii =1:length(BC)/Ngp
                    BC2(ii).data = BC((ii-1)*Ngp+1 : ii*Ngp , 3:5); 
                    for kk = 1 : size(Fgll,1) ;
                        pp(kk)  = sum(Fgll(kk,:)' - BC((ii-1)*Ngp+1 : ii*Ngp , 2));
                    end
                    ip = find(abs(pp)< 1e-1, 1 , 'first') ;
                    if ip > 4 || isempty(ip)
                        disp('face number is wrong')
                    end
                    BC2(ii).face = ip ;
                    BC2(ii).ielm = BC((ii-1)*Ngp+1 , 1);
            end
        end
        function BC = convert_revBC(BC2,Ngp)
            Fgll(4,:) = 1:Ngp:(Ngp-1)*Ngp+1;
            Fgll(1,:) = 1:1:Ngp;
            Fgll(2,:) = Ngp:Ngp:Ngp^2;
            Fgll(3,:) = (Ngp-1)*Ngp+1:1:Ngp*Ngp;
            BC = zeros(length(BC2)*Ngp,5);
            for ii = 1 : length(BC2)
                BC((ii-1)*Ngp+1:ii*Ngp,1) = BC2(ii).ielm;
                BC((ii-1)*Ngp+1:ii*Ngp,2)  =  Fgll(BC2(ii).face,:);
                BC((ii-1)*Ngp+1:ii*Ngp,3:5)  =  BC2(ii).data ;
            end
        end
        function BC = bc_read(fname,Ngp)
            fid = fopen(fname,'r');
            item = fscanf(fid,'%20i\n',1);
            tr = fscanf(fid,'%20i\n',1);
            for ii =1 : item
                tmp1s = fgets(fid) ;
                tmp1n = str2num(tmp1s) ;
                BC(ii).ielm = tmp1n(1) ;
                BC(ii).face = tmp1n(2) ;
                for jj = 1 : Ngp
                    tmp3s = fgets(fid);
                    tmp3n = str2num(tmp3s) ;
                    BC(ii).data(jj,1:3) = tmp3n;
                end
            end
        end
        function [BC,iface]=calcBC(FluidBC,EL,data,t,varargin)
            
            %define faces for GLL points in each element
            Ngp = sqrt(length(EL(1).GLL(:,1)));
            %only faces 1,3 and 4 are needed in this case
            
            Fgll(4,:) = 1:Ngp:(Ngp-1)*Ngp+1;
            Fgll(1,:) = 1:1:Ngp;
            Fgll(2,:) = Ngp:Ngp:Ngp^2;
            Fgll(3,:) = (Ngp-1)*Ngp+1:1:Ngp*Ngp;
            
            ELbc = find(FluidBC==t);
            for i = 1:length(ELbc)
                face = mod(ELbc(i),4);
                if face == 0
                    face = 4;
                end
                elnum = (ELbc(i)-face)/4+1;
                coord((i-1)*Ngp+1:i*Ngp,:) = EL(elnum).GLL(Fgll(face,:),1:2);
                eln((i-1)*Ngp+1:i*Ngp) = elnum;
                index((i-1)*Ngp+1:i*Ngp) = Fgll(face,:);
                iface((i-1)*Ngp+1:i*Ngp) = face ;
            end
            if ~isempty(data)
                disp('Write user-defined boundary conditions')
                xr=data.xr;
                yr=data.yr;
                ur=data.ur;
                vr=data.vr;
                %                 if nargin>0
                %                     Fu=TriScatteredInterp(xr,yr,ur);
                %                     Fv=TriScatteredInterp(xr,yr,vr);
                %                     ub=Fu(coord(:,1),coord(:,2));
                %                     vb=Fv(coord(:,1),coord(:,2));
                %                 else
                
                ub = griddata(xr,yr,ur,coord(:,1),coord(:,2),'cubic');
                vb = griddata(xr,yr,vr,coord(:,1),coord(:,2),'cubic');
                %wb = griddata(xr,yr,wr,coord(:,1),coord(:,2),'cubic');
                %                 end
                
                wb = zeros(size(ub));
                
                nn=norm(double(isnan(ub)),1)+norm(double(isnan(vb)),1)+norm(double(isnan(wb)),1);
                if nn~=0
                    disp('you have NaN in your interpolation of Bcs');
                end
                BC = [eln' index' ub vb wb];
            else
                for i=1:length(eln)
                    BC(i,:)= [eln(i) index(i) EL(eln(i)).Vinit(index(i),1) ...
                    EL(eln(i)).Vinit(index(i),2) EL(eln(i)).Vinit(index(i),3)];
                end
            end
         end
        function writeBC(BC,casename,varargin)
            % writeBC(BC,casename,varargin)
            % BC: data
            % casename: name of the file
            % varargin{1} : prefix (do not need to add dot .)
            % varargin{2} : time
            % varargin{3} : path of the file
            switch nargin
                case 2
                    path='./';
                    prefix='.BC';
                    time=0;
                case 3
                    prefix=['.' varargin{1}];
                    path='./';
                    time=0;
                case 4
                    prefix=['.' varargin{1}];
                    path='./';
                    time=varargin{2};
                otherwise
                    prefix=['.' varargin{1}];
                    path=varargin{3};
                    time=varargin{2};
            end
            fid = fopen([path casename prefix],'w');
            numBC = length(BC);
            fprintf(fid,'%20i\n',numBC);
            fprintf(fid,'%20i\n',time);
            switch size(BC,2)
                case 4
                    fprintf(fid,'%7i %1i %20.15f %20.15f\n', BC');
                case 5
                    fprintf(fid,'%7i %1i %20.15f %20.15f %20.15f\n', BC');
                case 9
                    fprintf(fid,'%7i %1i %20.15f %20.15f %20.15f %20.15f %20.15f %20.15f %20.15f\n', BC');
            end
            fclose(fid);
            disp('Bcs is written to the file');
        end
        function writeBC_adam(BC,casename,varargin)
            %writeBC_adam(BC,casename,varargin)
            prefix='.BC';
            time=0;
            switch nargin
                case 3
                    prefix=['.' varargin{1}];
                    time=0;
                case 4
                    prefix=['.' varargin{1}];
                    time=varargin{2};
            end
            fid = fopen([casename prefix],'w');
            numBC = length(BC);
            fprintf(fid,'%20i\n',numBC);
            fprintf(fid,'%20i\n',time);
            for ii = 1 : numBC
                fprintf(fid,'%7i %1i\n', BC(ii).ielm,BC(ii).face);
                fprintf(fid,'%20.15f %20.15f %20.15f\n', BC(ii).data');
            end
            fclose(fid);
            disp('Bcs is written to the file');
         end
        function EL=calcIC(EL,data)
            xr=data.xr;
            yr=data.yr;
            ur=data.ur;
            vr=data.vr;
            name=fieldnames(data);
            if length(name) > 4
                pr=data.pr;
            else
                pr=0.*ur;
            end
            
            N1=sqrt(length(EL(1).GLL(:,1)));
            N=N1-1;
            
            gll2D = zeros(length(EL)*(N+1)^2,2);
            for i = 1:length(EL)
                gll2D((i-1)*(N+1)^2+1:i*(N+1)^2,:) = EL(i).GLL(1:(N+1)^2,1:2);
            end
            
            %interpolate RANS solution onto gll points
            uinit = griddata(xr,yr,ur,gll2D(:,1),gll2D(:,2));
            vinit = griddata(xr,yr,vr,gll2D(:,1),gll2D(:,2));
            pinit = griddata(xr,yr,pr,gll2D(:,1),gll2D(:,2));
            
            % check for the NaN
            ind=find(isnan(uinit(:))==1);
            disp(['norm of the NaN in uinit = ' num2str(norm(ind)) ]);
            uinit(ind)=1.0;ind=find(isnan(uinit(:))==1);
            disp(['norm of the NaN in uinit = ' num2str(norm(ind)) ]);ind=[];
            
            ind=find(isnan(vinit(:))==1);
            disp(['norm of the NaN in uinit = ' num2str(norm(ind)) ]);
            vinit(ind)=0.0; ind=find(isnan(vinit(:))==1);
            disp(['norm of the NaN in vinit = ' num2str(norm(ind)) ]);ind=[];
            
            ind=find(isnan(pinit(:))==1);
            disp(['norm of the NaN in pinit = ' num2str(norm(ind)) ]);
            pinit(ind)=0.0;  ind=find(isnan(pinit(:))==1);
            disp(['norm of the NaN in pinit = ' num2str(norm(ind)) ]);ind=[];
            
            data = zeros(length(EL),(N+1)^2,5);
            for i = 1:length(EL)
                for j = 1:1 %2D fields
                    EL(i).Vinit((j-1)*(N+1)^2+1:j*(N+1)^2,1) = uinit((i-1)*(N+1)^2+1:i*(N+1)^2);
                    EL(i).Vinit((j-1)*(N+1)^2+1:j*(N+1)^2,2) = vinit((i-1)*(N+1)^2+1:i*(N+1)^2);
                    EL(i).Pinit((j-1)*(N+1)^2+1:j*(N+1)^2,1) = pinit((i-1)*(N+1)^2+1:i*(N+1)^2);
                end
                
                %write data-file to be written out
            end
        end
        %*****************************************************************
        function gllwrite(EL,N,path,casename)
            plotgll = zeros(length(EL)*(N+1)^2,2);
            for i = 1:length(EL)
                plotgll((i-1)*(N+1)^2+1:i*(N+1)^2,:) = EL(i).GLL;
            end
            
            
            [fid,message] = fopen([path casename '.grid'],'w','ieee-le');
            count = fwrite(fid,plotgll','double');
            fclose(fid);
            
            disp([num2str(count/2) ' gll points written to gll.dat']);
            
            
        end
        %******************************************************************
        function pts = georefine(pt0,pt1,Npt,pratio)
            
            if pratio == 1.0
                pratio = 1.0+1.0e-10;
            end
            
            pts = zeros(1,Npt);
            pts(1) = pt0;
            
            fc  = (pt1-pt0)*(1.0-pratio)/(1.0-pratio^(Npt-1));
            for p = 2:Npt
                pts(p) = pts(p-1)+fc*pratio^(p-2);
            end
            
        end
        %******************************************************************
        function [x,pratio] = findgeoP(xs,xe,pratioinit,N,ddx)
            
            %define anonymous function
            f = @(pratio)findP(xs,xe,pratio,N,ddx);
            
            %find optimal pratio2
            options = optimset('GradObj','off','Display','off','LargeScale','off');
            [pratio,fval] = fminunc(f,pratioinit,options);
            
            
            x = nekmesh.georefine(xs,xe,N,pratio);
            
            
            function f = findP(x0,x1,pratio,N,dxgiven)
                
                xd = nekmesh.georefine(x0,x1,N,pratio);
                
                dx = xd(2)-xd(1);
                
                f = abs(dxgiven - dx);
            end
            
        end
        %******************************************************************
        function [xdist,ratio,nnf]=findR2(x1,x2,dx1,dx2)
            %[xdist,ratio,nnf,tol]=findR2(x1,x2,dx1,dx2
            %find a distribution between x1 and x2 with dx1 at the beginnig
            %and dx2 at the end.
            nmin=abs(round(min(abs(x2-x1)/dx1,abs(x2-x1)/dx2)));
            if(nmin==0);
                nmin=1;disp('your original nmin is too low');
            end
            nmax=abs(round(max(abs(x2-x1)/dx1,abs(x2-x1)/dx2)));
            
            nn=1;
            for ii=min(nmin,nmax):max(nmin,nmax)
                [xdist,~]=nekmesh.findR(x1,x2,[0.8 1.4],ii,dx1);
                tmp=diff(xdist);
                res(nn)=dx2-tmp(end);
                nf(nn)=ii;
                nn=nn+1;
            end
            [~,ind]=min(abs(res));
            nnf=nf(ind);
            [xdist,ratio]=nekmesh.findR(x1,x2,[0.8 1.4],nnf,dx1);
            tmp=diff(xdist);
        end
        function [xdist,ratio]=findR(x1,x2,r1,N,dx,varargin)
            % [xdist ratio]=findR(x1,x2,r1,N,dx,varargin)
            % varargin{1}
            % find a distribution beween x1 and x2 with N element and dx at
            % the beginning if varargin{1} > empty dx would be calculated
            % at the end of the domain
            L=x2-x1;
            n=N-1;
            
            fb = f(r1(2),dx,n,L);
            fa = f(r1(1),dx,n,L);
            if fa <=0 && fb >=0
                if(isinf(fb) || isinf(fa));disp('you have infinity in your initial ratio');end
                ep=1e-10;
                a=r1(1);
                b=r1(2);
                maxit=1e6;
                itt=1;
                while (b-a > ep*b) && (itt<maxit)
                    
                    x = (a+b)/2;
                    fx = f(x,dx,n,L);
                    if sign(fx) == sign(fa)
                        a = x; fa = fx;
                    else
                        b = x; fb = fx;
                    end
                    itt=itt+1;
                end
                if itt>maxit
                    disp('itt exceed the maximum value')
                end
                
            elseif fa > 0 && fb >=0
                disp('you hit the lower bound, decrease elements # or lower bound')
                x=0.9;
            elseif   fb < 0
                disp('you hit the upper bound')
                x=1.2;
            else
                disp('some thing is wrong')
                return
            end
            xdist = nekmesh.georefine(x1,x2,N,x);
            ratio=x;
            if nargin>5
                xdist=x2-fliplr(xdist);
                xdist(1)=x1;
                ratio=1/ratio;
            end
            
            function sum=f(r,dx,n,L)
                sum=0;
                for i=1:n
                    sum=sum+dx*r^(i-1);
                end
                sum=sum-L;
            end
        end
        %******************************************************************
        function [x,pratio2] = refine2(x0,x1,xr,N1,N2,pratio,pratio2init)
            
            xd1 = nekmesh.georefine(x0,xr,N1,pratio);
            
            %define anonymous function
            f = @(pratio2)findprat(xr,x1,N2,pratio2,xd1(end)-xd1(end-1));
            
            %find optimal pratio2
            options = optimset('GradObj','off','Display','off','LargeScale','off');
            [pratio2,fval] = fminunc(f,pratio2init,options);
            
            xd2 = nekmesh.georefine(xr,x1,N2,pratio2);
            
            x = [xd1 xd2(2:end)];
            
            function f = findprat(x0,x1,N,pratio,deltaxinit)
                
                xd = georefine(x0,x1,N,pratio);
                
                dx = xd(2)-xd(1);
                
                f = abs(deltaxinit - dx);
            end
            function Points  = georefine(pt0,pt1,Npt,pratio)
                Npt = floor(Npt);
                
                S = length(pt0);
                
                if size(pt1) ~= S
                    error('pt0 and pt1 should have same size')
                end
                
                for i = 1:S
                    
                    pts = zeros(1,Npt);
                    pts(1) = pt0(i);
                    
                    fc  = (pt1(i)-pt0(i))*(1.0-pratio)/(1.0-pratio^(Npt-1));
                    for p = 2:Npt
                        pts(p) = pts(p-1)+fc*pratio^(p-2);
                    end
                    
                    Points(i,:) = pts;
                end
            end
            
        end
        %******************************************************************
        function [xs,ys]=buildmesh(x1,y1,x2,y2,x3,y3,x4,y4)
            %------ 3 ----
            %   4       2
            %-- ----1 ----
            xs=zeros(length(x2),length(x1));
            ys=zeros(size(xs));
            
            %---- first side
            xs(1,:)=x1;
            ys(1,:)=y1;
            %----- second side
            xs(:,end)=x2;
            ys(:,end)=y2;
            %---- third side
            xs(end,:)=x3;
            ys(end,:)=y3;
            %-- fourth side
            xs(:,1)=x4;
            ys(:,1)=y4;
            [xs,ys]=nekmesh.BLI(xs,ys);
            
        end
        %******************************************************************
        function xr=chebpoint(xs,xe,n)
            %  xr=chebpoint(xs,xe,n)
            xx=fliplr(1- cos((n:-1:0)/n*pi/2));
            xr=xs+(xe-xs)*xx;
        end
        %******************************************************************
        function [x,y]=LAP(x,y,k)
            for kk=1:k
                for i=2:size(x,1)-1
                    for j=2:size(x,2)-1
                        
                        g11=0.25*( (x(i+1,j)-x(i-1,j)).^2 + (y(i+1,j)-y(i-1,j)).^2 );
                        g22=0.25*( (x(i,j+1)-x(i,j-1)).^2 + (y(i,j+1)-y(i,j-1)).^2 );
                        g12=0.25*( (x(i+1,j)-x(i-1,j))*(x(i,j+1)-x(i,j-1)) + (y(i+1,j)-y(i-1,j))*(y(i,j+1)-y(i,j-1)) );
                        
                        xx=1/(2*(g11+g22));
                        x(i,j)=xx* ( ...
                            g22*(x(i+1,j)+x(i-1,j)) ...
                            -0.5*g12*(x(i+1,j+1)-x(i+1,j-1)+x(i-1,j-1)-x(i-1,j+1)) ...
                            +    g11*(x(i,j+1)+x(i,j-1))...
                            );
                        
                        y(i,j)=xx* ( ...
                            g22*(y(i+1,j)+y(i-1,j)) ...
                            -0.5*g12*(y(i+1,j+1)-y(i+1,j-1)+y(i-1,j-1)-y(i-1,j+1)) ...
                            +    g11*(y(i,j+1)+y(i,j-1))...
                            );
                    end
                end
            end
        end
        %******************************************************************
        function [x,y]=TFI(x,y)
            
            imax=size(x,1);
            jmax=size(x,2);
            % ss=linspace(0,imax,imax);
            % tt=linspace(0,jmax,jmax);
            for i=1:imax
                s=(i-1)/(imax-1);
                for j=1:jmax
                    t=(j-1)/(jmax-1);
                    xs=(1-s)*x(1,j)+s*x(imax,j);
                    xt=(1-t)*x(i,1)+t*x(i,jmax);
                    xst=(1-s)*((1-t)*x(1,1)+t*x(1,jmax))+...
                        s* ((1-t)*x(imax,1)+t*x(imax,jmax));
                    x(i,j)=xs+xt-xst;
                    
                    
                    ys=(1-s)*y(1,j)+s*y(imax,j);
                    yt=(1-t)*y(i,1)+t*y(i,jmax);
                    yst=(1-s)*((1-t)*y(1,1)+t*y(1,jmax))+...
                        s* ((1-t)*y(imax,1)+t*y(imax,jmax));
                    y(i,j)=ys+yt-yst;
                end
            end
        end
        %******************************************************************
        function [xs,ys]=BLI(xs,ys)
            npx=size(xs,2);
            npy=size(xs,1);
            for j=2:npx-1
                for i=2:npy-1
                    A1=[xs(1,j) ys(1,j)];
                    A2=[xs(end,j) ys(end,j)];
                    
                    B1=[xs(i,1) ys(i,1)];
                    B2=[xs(i,end) ys(i,end)];
                    
                    ts = [(A2(:) - A1(:))  -(B2(:)-B1(:))]\(B1(:) - A1(:));
                    xy=A1 + (A2 - A1)*ts(1) ;
                    
                    xs(i,j)=xy(1);
                    ys(i,j)=xy(2);
                end
            end
        end
        %******************************************************************
        function [x,y]=WCS(x,y,k)
            
            for kk=1:k
                
                for i=2:size(x,1)-1
                    for j=2:size(x,2)-1
                        
                        
                        
                        g11=0.25*( (x(i+1,j)-x(i-1,j)).^2 + (y(i+1,j)-y(i-1,j)).^2 );
                        g22=0.25*( (x(i,j+1)-x(i,j-1)).^2 + (y(i,j+1)-y(i,j-1)).^2 );
                        g12=0.25*( (x(i+1,j)-x(i-1,j))*(x(i,j+1)-x(i,j-1)) + (y(i+1,j)-y(i-1,j))*(y(i,j+1)-y(i,j-1)) );
                        
                        xx=1/(2*(g11+g22));
                        
                        
                        x(i,j)=xx* ( ...
                            g22*(x(i+1,j)+x(i-1,j)) ...
                            -0.5*g12*(x(i+1,j+1)-x(i+1,j-1)+x(i-1,j-1)-x(i-1,j+1)) ...
                            +    g11*(x(i,j+1)+x(i,j-1))...
                            );
                        
                        y(i,j)=xx* ( ...
                            g22*(y(i+1,j)+y(i-1,j)) ...
                            -0.5*g12*(y(i+1,j+1)-y(i+1,j-1)+y(i-1,j-1)-y(i-1,j+1)) ...
                            +    g11*(y(i,j+1)+y(i,j-1))...
                            );
                        
                        
                        
                    end
                end
                
                
            end
            
        end
        %******************************************************************
        function [x3,y3]=TTM(x,y,kk)
            PP=zeros([size(x,1)+2 size(x,2)+2]);
            QQ=zeros(size(PP));
            omega=0.5;
            
            
            for k=1:kk
                x1=x;
                y1=y;
                [xs,ys]=shadow(x,y);
                
                for i=2:size(xs,1)-1
                    for j=2:size(xs,2)-1
                        g11=0.25*( (xs(i+1,j)-xs(i-1,j)).^2 + (ys(i+1,j)-ys(i-1,j)).^2 );
                        g22=0.25*( (xs(i,j+1)-xs(i,j-1)).^2 + (ys(i,j+1)-ys(i,j-1)).^2 );
                        g12=0.25*( (xs(i+1,j)-xs(i-1,j))*(xs(i,j+1)-xs(i,j-1)) + ...
                            (ys(i+1,j)-ys(i-1,j))*(ys(i,j+1)-ys(i,j-1))         );
                        % source terms
                        
                        PP(i,j)=-1/(2*g11)*( ( xs(i+1,j)-xs(i-1,j) ) * (xs(i+1,j)-2*xs(i,j)+xs(i-1,j) ) ...
                            +( ys(i+1,j)-ys(i-1,j) ) * (ys(i+1,j)-2*ys(i,j)+ys(i-1,j) ) ) ...
                            -1/(2*g22)*( ( xs(i+1,j)-xs(i-1,j) ) * (xs(i,j+1)-2*xs(i,j)+xs(i,j-1) ) ...
                            +( ys(i+1,j)-ys(i-1,j) ) * (ys(i,j+1)-2*ys(i,j)+ys(i,j-1) ) ) ;
                        
                        QQ(i,j)= -1/(2*g22)*( ( xs(i,j+1)-xs(i,j-1) ) * (xs(i,j+1)-2*xs(i,j)+xs(i,j-1) ) ...
                            +( ys(i,j+1)-ys(i,j-1) ) * (ys(i,j+1)-2*ys(i,j)+ys(i,j-1) ) ) ...
                            -1/(2*g11)*( ( xs(i,j+1)-xs(i,j-1) ) * (xs(i+1,j)-2*xs(i,j)+xs(i-1,j) ) ...
                            +( ys(i,j+1)-ys(i,j-1) ) * (ys(i+1,j)-2*ys(i,j)+ys(i-1,j) ) );
                    end
                end
                P=PP(2:end,2:end);
                Q=QQ(2:end,2:end);
                [P,Q]=nekmesh.TFI(P,Q);
                
                
                
                
                
                for i=2:size(x,1)-1
                    for j=2:size(x,2)-1
                        
                        %------------- mesh metrics
                        
                        g11=0.25*( (x(i+1,j)-x(i-1,j)).^2 + (y(i+1,j)-y(i-1,j)).^2 );
                        g22=0.25*( (x(i,j+1)-x(i,j-1)).^2 + (y(i,j+1)-y(i,j-1)).^2 );
                        g12=0.25*( (x(i+1,j)-x(i-1,j))*(x(i,j+1)-x(i,j-1)) + ...
                            (y(i+1,j)-y(i-1,j))*(y(i,j+1)-y(i,j-1))         );
                        % source terms
                        
                        %----main--------
                        xx=1/(2*(g11+g22));
                        
                        x(i,j)=xx * (   g22 * ((1+P(i,j)/2)*x(i+1,j)+(1-P(i,j)/2)*x(i-1,j)) - ...
                            -g12*(x(i+1,j+1)-x(i+1,j-1)+x(i-1,j-1)-x(i-1,j+1)) ...
                            +    g11*((1+Q(i,j)/2)*x(i,j+1)+(1-Q(i,j)/2)*x(i,j-1)));
                        
                        y(i,j)=xx * (   g22 * ((1+P(i,j)/2)*y(i+1,j)+(1-P(i,j)/2)*y(i-1,j)) ...
                            -g12*(y(i+1,j+1)-y(i+1,j-1)+y(i-1,j-1)-y(i-1,j+1)) ...
                            +    g11*((1+Q(i,j)/2)*y(i,j+1)+(1-Q(i,j)/2)*y(i,j-1)));
                    end
                end
                
                
                x=omega*x+(1-omega)*x1;
                y=omega*y+(1-omega)*y1;
                
            end
            
            x3=x;
            y3=y;
            function [x,y]=shadow(xr,yr)
                x=zeros(size(xr,1)+2,size(xr,2)+2);
                y=zeros(size(x));
                
                x(2:end-1,2:end-1)=xr;
                y(2:end-1,2:end-1)=yr;
                
                %line s=0 without ending points correspont to i=2
                for j=2:size(x,2)-1
                    
                    x0=x(3,j)-x(2,j);
                    y0=y(3,j)-y(2,j);
                    switch j
                        case 2
                            yt=(y(2,j+1)-y(2,j));
                            xt=(x(2,j+1)-x(2,j));
                        case size(x,2)-1
                            yt=(y(2,j)-y(2,j-1));
                            xt=(x(2,j)-x(2,j-1));
                        otherwise
                            yt=(y(2,j+1)-y(2,j-1))/2;
                            xt=(x(2,j+1)-x(2,j-1))/2;
                    end
                    g22=xt^2+yt^2;
                    n=[yt,-xt]/sqrt(g22);
                    xx=[x0,y0];
                    vv=dot(xx,n)*n;
                    %----- debug---------
                    dy=abs(y(3,j)-y(2,j));
                    vv=vv/norm(vv)*dy;
                    
                    x(1,j)=x(2,j)-vv(1);
                    y(1,j)=y(2,j)-vv(2);
                end
                %line s=end without ending points correspont to i=end-1
                for j=2:size(x,2)-1
                    
                    x0=x(end-1,j)-x(end-2,j);
                    y0=y(end-1,j)-y(end-2,j);
                    
                    switch j
                        case 2
                            yt=(y(end-1,j+1)-y(end-1,j));
                            xt=(x(end-1,j+1)-x(end-1,j));
                        case size(x,2)-1
                            yt=(y(end-1,j)-y(end-1,j-1));
                            xt=(x(end-1,j)-x(end-1,j-1));
                            
                        otherwise
                            yt=(y(end-1,j+1)-y(end-1,j-1))/2;
                            xt=(x(end-1,j+1)-x(end-1,j-1))/2;
                    end
                    
                    
                    g22=xt^2+yt^2;
                    
                    n=[yt,-xt]/sqrt(g22);
                    xx=[x0,y0];
                    vv=dot(xx,n)*n;
                    dy=abs(y(end-2,j)-y(end-1,j));
                    
                    vv=vv/norm(vv)*dy;
                    x(end,j)=x(end-1,j)+vv(1);
                    y(end,j)=y(end-1,j)+vv(2);
                end
                %line t=0 without ending points correspont to j=2
                for i=2:size(x,1)-1
                    x0=x(i,3)-x(i,2);
                    y0=y(i,3)-y(i,2);
                    switch 2
                        case 2
                            ys=(y(i+1,2)-y(i,2));
                            xs=(x(i+1,2)-x(i,2));
                        case size(x,1)-1
                            ys=(y(i,2)-y(i-1,2));
                            xs=(x(i,2)-x(i-1,2));
                        otherwise
                            ys=(y(i+1,2)-y(i-1,2))/2;
                            xs=(x(i+1,2)-x(i-1,2))/2;
                    end
                    
                    g11=xs^2+ys^2;
                    
                    n=[-ys,xs]/sqrt(g11);
                    xx=[x0,y0];
                    
                    vv=dot(xx,n)*n;
                    
                    
                    x(i,1)=x(i,2)-vv(1);
                    y(i,1)=y(i,2)-vv(2);
                end
                
                %line t=1 without ending points correspont to j=end-1
                for i=2:size(x,1)-1
                    x0=x(i,end-1)-x(i,end-2);
                    y0=y(i,end-1)-y(i,end-2);
                    switch 2
                        case 2
                            ys=(y(i+1,end-1)-y(i,end-1));
                            xs=(x(i+1,end-1)-x(i,end-1));
                        case size(x,1)-1
                            ys=(y(i,end-1)-y(i-1,end-1));
                            xs=(x(i,end-1)-x(i-1,end-1));
                        otherwise
                            ys=(y(i+1,end-1)-y(i-1,end-1))/2;
                            xs=(x(i+1,end-1)-x(i-1,end-1))/2;
                    end
                    
                    
                    ys=(y(i+1,end-1)-y(i-1,end-1))/2;
                    xs=(x(i+1,end-1)-x(i-1,end-1))/2;
                    g11=xs^2+ys^2;
                    
                    n=[-ys,xs]/sqrt(g11);
                    xx=[x0,y0];
                    vv=dot(xx,n)*n;
                    
                    x(i,end)=x(i,end-1)+vv(1);
                    y(i,end)=y(i,end-1)+vv(2);
                end
                
                %--- four point in the boundary-----
                
                
                %point 1,1
                xdata1=[x(2,1) x(3,1)];
                ydata1=[y(2,1) y(3,1)];
                xdata2=[x(1,2) x(1,3)];
                ydata2=[y(1,2) y(1,3)];
                xy=intersect(xdata1,ydata1,xdata2,ydata2);
                x(1,1)=xy(1);
                y(1,1)=xy(2);
                %point end,1
                xdata1=[x(end-1,1) x(end-2,1)];
                ydata1=[y(end-1,1) y(end-2,1)];
                xdata2=[x(end,2) x(end,3)];
                ydata2=[y(end,2) y(end,3)];
                
                xy=intersect(xdata1,ydata1,xdata2,ydata2);
                x(end,1)=xy(1);
                y(end,1)=xy(2);
                
                %point end,end
                
                xdata1=[x(end-1,end) x(end-2,end)];
                ydata1=[y(end-1,end) y(end-2,end)];
                xdata2=[x(end,end-1) x(end,end-2)];
                ydata2=[y(end,end-1) y(end,end-2)];
                
                xy=intersect(xdata1,ydata1,xdata2,ydata2);
                x(end,end)=xy(1);
                y(end,end)=xy(2);
                
                %point 1,end
                
                xdata1=[x(1,end-1) x(1,end-2)];
                ydata1=[y(1,end-1) y(1,end-2)];
                xdata2=[x(2,end) x(3,end)];
                ydata2=[y(2,end) y(3,end)];
                
                xy=intersect(xdata1,ydata1,xdata2,ydata2);
                x(1,end)=xy(1);
                y(1,end)=xy(2);
                
                
                
                
            end
            function xy=intersect(xdata1,ydata1,xdata2,ydata2)
                A1=[xdata1(1) ydata1(1)];
                A2=[xdata1(2) ydata1(2)];
                
                B1=[xdata2(1) ydata2(1)];
                B2=[xdata2(2) ydata2(2)];
                
                ts = [(A2(:) - A1(:))  -(B2(:)-B1(:))]\(B1(:) - A1(:));
                xy=A1 + (A2 - A1)*ts(1) ;
                
                
                
                
                
                %                 %point 1,1
                %                 if abs(xdata1(1)-xdata1(2))<1e-6
                %                     xy(1)=xdata1(1);
                %                     Q = polyfit(xdata2,ydata2,1);
                %                     xy(2)=Q(1)*xdata1(1)+Q(2);
                %                 elseif abs(xdata2(1)-xdata2(2))<1e-6
                %                     xy(1)=xdata2(1);
                %                     P = polyfit(xdata1,ydata1,1);
                %                     xy(2)=P(1)*xdata2(1)+P(2);
                %                 else
                %
                %                     P = polyfit(xdata1,ydata1,1);
                %                     Q = polyfit(xdata2,ydata2,1);
                %
                %                     A=[-P(1) 1;-Q(1) 1];
                %                     B=[P(2);Q(2)];
                %
                %                     xy=A\B;
                %                 end
                
            end
            
        end
        function xy=intersect(A1,A2,B1,B2)
            ts = [(A2(:) - A1(:))  -(B2(:)-B1(:))]\(B1(:) - A1(:));
            xy=A1 + (A2 - A1)*ts(1) ;
        end
        %******************************************************************
        function box=gridbox(xp,yp,npx,npy)
            fig6=figure(6);
            %  set(fig6,'Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2]);
            hold on
            for i = 1:npy
                plot(xp(i,:),yp(i,:),'r')
            end
            for i = 1:npx
                plot(xp(:,i),yp(:,i),'r')
            end
            axis equal
            grid on
            %set(gca,'fontsize',24,'fontname','times','linewidth',2)
            xlabel('x')
            ylabel('y')
            hold on
            %---------------------- Correct mesh------------------
            
            set(fig6,'WindowButtonDownFcn',{@catchpoint,npx,npy})
            evalresponse = input('All points ok?');
            %             ylim([-12 30])
            
            h1 = get(gca,'Children');
            
            S = length(h1)-npy+1;
            indd = length(h1):-1:S;
            indb = 1;
            boxnx(indb) = npx;
            indy = 0;
            for i = 1:length(indd)
                indy = indy+1;
                xdata = get(h1(indd(i)),'XData');
                ydata = get(h1(indd(i)),'YData');
                
                xp(i,:) = xdata;
                yp(i,:) = ydata;
                
                indnan = isnan(xdata);
                xdata = xdata(~indnan);
                ydata = ydata(~indnan);
                if length(xdata) == boxnx(indb)
                    box(indb).x(indy,:) = xdata;
                    box(indb).y(indy,:) = ydata;
                else
                    indb = indb+1;
                    indy = 1;
                    boxnx(indb) = length(xdata);
                    box(indb).x(indy,:) = xp(i-1,~indnan);
                    box(indb).y(indy,:) = yp(i-1,~indnan);
                    indy = indy+1;
                    box(indb).x(indy,:) = xdata;
                    box(indb).y(indy,:) = ydata;
                end
            end
            function catchpoint(src,evnt,npx,npy)
                
                if strcmp(get(src,'SelectionType'),'normal')
                    set(src,'pointer','circle')
                    cp = get(gca,'CurrentPoint');
                    xinit = cp(1,1);yinit = cp(1,2);
                    
                    ploth = get(gca,'Children');
                    S = length(ploth)-npy+1;
                    
                    is = 0;
                    for i = S:length(ploth)
                        is = is+1;
                        x(is,:) = get(ploth(i),'XData');
                        y(is,:) = get(ploth(i),'YData');
                    end
                    
                    %find closest point
                    dist = sqrt((x-xinit).^2 + (y-yinit).^2);
                    
                    [val,ind] = min(dist(:));
                    Sd = size(dist);
                    [indj,indi] = ind2sub(Sd,ind);
                    
                    indi1 = indi;
                    indj1 = S+indj-1;
                    x1v = get(ploth(indj1),'XData');
                    y1v = get(ploth(indj1),'YData');
                    x1 = x1v(indi1);
                    y1 = y1v(indi1);
                    
                    
                    S1 = length(ploth)-npy-npx+1;
                    indi2 = S1+npx-indi;
                    indj2 = npy-indj+1;
                    x2v = (get(ploth(indi2),'XData'));
                    y2v = (get(ploth(indi2),'YData'));
                    x2 = x2v(indj2);
                    y2 = y2v(indj2);
                    %plot(x1,y1,'g*',x2,y2,'bo')
                    
                    set(src,'WindowButtonMotionFcn',{@wbmcb,indi1,indj1,indi2,indj2})
                    set(src,'WindowButtonUpFcn',@wbucb)
                end
                
                function wbmcb(src,evnt,indi1,indj1,indi2,indj2)
                    
                    cp = get(gca,'CurrentPoint');
                    ploth = get(gca,'Children');
                    
                    %set(ploth(1),'XData',cp(1,1),'YData',cp(1,2))
                    xvec1 = get(ploth(indj1),'XData');
                    yvec1 = get(ploth(indj1),'YData');
                    xvec1(indi1) = cp(1,1);
                    yvec1(indi1) = cp(1,2);
                    set(ploth(indj1),'XData',xvec1,'YData',yvec1)
                    
                    xvec2 = get(ploth(indi2),'XData');
                    yvec2 = get(ploth(indi2),'YData');
                    xvec2(indj2) = cp(1,1);
                    yvec2(indj2) = cp(1,2);
                    set(ploth(indi2),'XData',xvec2,'YData',yvec2)
                    drawnow
                end
                
                
                function wbucb(src,evnt)
                    set(src,'pointer','arrow')
                    set(src,'WindowButtonMotionFcn','')
                    set(src,'WindowButtonUpFcn','')
                end
            end
            
        end
        %*****************************************************************
        %relp number of point on coarse, index  = 1, half of the lead
        function [xy_fine,xy_coarse,Arc] = getgeolead(relp,index,ar,extend)
            
            %relp number of point on lead coarse, index  = 1, half of the lead,
            %ar aspect ratio
            b0 = 1 ;
            a0 = ar * b0;
            % relp = 99999;
            pi2 = pi/2.;
            xoa = 1.0-cos((0:relp).*pi2/relp);  % x/a, cosine-like distribution gives refinement in nose region
            m = 2+xoa.^2;
            arg = (1-xoa).^m;
            arg = 1-arg;
            
            x = -cos((0:relp).*pi2/relp);
            x = a0*x;      % The x coordinate of the element vertices along the ellipses
            %x = x-a0-x(1);
            x = x-x(1);
            y =  b0*sqrt(arg);
            
            x1 = fliplr(x);
            y1 = -fliplr(y);
            x1(end)= [];y1(end)= [];
            x = [x1'; x'] ;
            y = [y1'; y'] ;
            xy_coarse = [x  y] ;
            
            if index==1
                ind=find(xy_coarse(:,2)<=0,1,'last');
                xy_coarse=xy_coarse(ind:end,:);
            end
            
            
            xy_fine = [] ;
            
            %------------------------------------ length of curve-------------
            relp = 49999 ;
            xoa = 1.0-cos((0:relp).*pi2/relp);  % x/a, cosine-like distribution gives refinement in nose region
            m = 2+xoa.^2;
            arg = (1-xoa).^m;
            arg = 1-arg;
            
            x = -cos((0:relp).*pi2/relp);
            x = a0*x;      % The x coordinate of the element vertices along the ellipses
            x = x-x(1);
            y =  b0*sqrt(arg);
            x1 = fliplr(x);
            y1 = -fliplr(y);
            x1(end)= [];y1(end)= [];
            x = [x1'; x'] ;
            y = [y1'; y'] ;
            xy_fine = [x y] ;
            if index==1
                ind=find(xy_fine(:,2)<=0,1,'last');
                xy_fine=xy_fine(ind:end,:);
            end
            if ~isempty(extend)
                tn = extend ;
                line = linspace(xy_fine(end,1),tn,100000)' ;
                line(1) = [] ;
                xy_fin(:,1) = [xy_fine(:,1);line ];
                xy_fin(:,2) = [xy_fine(:,2);ones(size(line)) ];
                cx(1,:) = xy_fin(:,1) ;
                cy(1,:) = xy_fin(:,2) ;
                cl=[0 cumsum(sqrt(diff(cx(1,:)).^2+diff(cy(1,:)).^2))];
                Arc.x = cx ;
                Arc.y = cy ;
                Arc.L= cl;
                
                
            else
                cx(1,:) = xy_fine(:,1) ;
                cy(1,:) = xy_fine(:,2) ;
                cl=[0 cumsum(sqrt(diff(cx(1,:)).^2+diff(cy(1,:)).^2))];
                Arc.x = cx ;
                Arc.y = cy ;
                Arc.L= cl;
            end
            
        end
        %*****************************************************************
        function writeic(EL,casename,prefix,std,varargin)
            %----------------------------------------------------------------------
            %writeic(EL,casename,prefix,std,varargin)
            %write initial conditions from RANS solution
            % writeic(EL,casename,prefix,std,varargin)
            %casenaem and prefix are the case naeem and prefix
            %prefix  the prefix of the name of the output for ex. 'IC'
            %std is 4 for sigle and 8 for double precision
            %varargin{1}: the time of the fields
            %varargin{2}: data: interpolation data
            %----------------------------------------------------------------------
            switch nargin
                case 4
                    t=0;
                    intp=0;
                case 5
                    t=varargin{1};
                    intp=0;
                case 6
                    t=varargin{1};
                    intp=1;
                    data=varargin{2};
            end
            path='./';
            prefix=['.' prefix];
            var(1:3)=0;
            var(1)= cell2mat(strfind(fieldnames(EL),'GLL'));
            var(2)= cell2mat(strfind(fieldnames(EL),'Vinit'));
            var(3)= cell2mat(strfind(fieldnames(EL),'Pinit'));
       
            N=sqrt(length(EL(1).GLL(:,1)))-1;
            
            disp('Write initial conditions')
            
            %interpolate RANS solution onto gll points
            if intp
                vp=cell2mat(strfind(fieldnames(data),'pr'));
                gll=zeros(length(EL),(N+1)^2,2);
                for ii=1:length(EL)
                    gll(ii,:,1)=EL(ii).GLL(:,1);
                    gll(ii,:,2)=EL(ii).GLL(:,2);
                end
                 uinit = griddata(data.xr,data.yr,data.ur,gll(:,:,1),gll(:,:,2));
                 vinit = griddata(data.xr,data.yr,data.vr,gll(:,:,1),gll(:,:,2));
                 if isempty(cell2mat(strfind(fieldnames(data),'pr')))
                     pinit=0.*uinit;
                 else
                    pinit = griddata(data.xr,data.yr,data.pr,gll(:,:,1),gll(:,:,2));
                 end
                 for ii=1:length(EL)
                     EL(ii).Vinit(:,1)=uinit(ii,:);
                     EL(ii).Vinit(:,2)=vinit(ii,:);
                     EL(ii).Pinit(:,1)=pinit(ii,:);
                 end
            end
            N1 = N+1;
            %write nekton restart file
            %Header: 2D field
            if isempty(t)
                tt=sprintf('%20.13E',0);
            else
                tt=sprintf('%20.13E',t);
            end
            if(std==4);precision='*float32';else precision='*float64';end
            hdr = ['#std ' num2str(std) '        1                       ' tt '         0      0      1 XUP                                              '];
            % N:
            Nst = num2str(N1);
            Nstz = num2str(1);
            hdr(9-length(Nst)+1:9)   = Nst;
            hdr(12-length(Nst)+1:12) = Nst;
            hdr(15-length(Nstz)+1:15) = Nstz;
            % nel:
            Nelm=length(EL);
            Nelmst = num2str(Nelm);
            hdr(26-length(Nelmst)+1:26) = Nelmst;
            hdr(37-length(Nelmst)+1:37) = Nelmst;
            %%% tag
            tag = 6.54321;
%             stsw = nekmesh.schreib(data,Nelm,hdr,tag,1,2,'le',[path casename prefix],'xup',precision);
            nekmesh.writeneknew(EL,hdr,tag,2,'le',[path casename prefix],precision);
        end
        %*****************************************************************
        function status = schreib(data,nel,header,tag,cnttail,dim,endian,fname,fields,precision)
            cnttail=cnttail + 1;
            %%% This function writes binary data in the nek5000 new file format
            if precision == 1
                
            end
            % Open file
            [fid,message] = fopen(fname,'w+',['ieee-' endian]);
            
            if fid == -1
                disp(message)
                return
            end
            
            % % New header
            % header(20:39) = '                    ';
            % snel = num2str(nel);
            % header(23:22+length(snel)) = snel;
            % header(34:33+length(snel)) = snel;
            
            fwrite(fid,header,'*char');
            
            fwrite(fid,tag,'*float32');
            
            % Element numbers
            lglel = 1:nel;
            lglel = lglel';
            fwrite(fid,lglel,'*int32');
            
            
            % How many fields?
            if strcmp(fields,'xup')
                if dim == 2
                    nfields = 5;
                elseif dim == 3
                    nfields = 7;
                end
            elseif strcmp(fields,'up')
                if dim == 2
                    nfields = 3;
                elseif dim == 3
                    nfields = 4;
                end
            end
            
            
            % Write data
            for iel=1:nel
                for ifld = 1:dim
                    fwrite(fid,data(iel,:,ifld),precision); %1:N^dim
                end
            end
            for iel=1:nel
                for ifld = dim+1:2*dim
                    fwrite(fid,data(iel,:,ifld),precision);
                end
            end
            for iel=1:nel
                fwrite(fid,data(iel,:,nfields),precision);
            end
            % % Not understood: More data to be written in the fld file.
            % % Pad with zeros.
            % cnttail = nel*cnttail;
            % quatsch = zeros(1,cnttail);
            % fwrite(fid,quatsch,'*float32');
            
            status = fclose(fid);
        end
        function writeneknew(EL,header,tag,dim,endian,fname,precision)
            nel=length(EL);
            % Open file
            [fid,message] = fopen(fname,'w+',['ieee-' endian]);
            if fid == -1;disp(message);return;end
                       
            fwrite(fid,header,'*char');
            fwrite(fid,tag,'*float32');
            
            % Element numbers
            lglel = 1:nel;
            lglel = lglel';
            fwrite(fid,lglel,'*int32');
            
            
            % How many fields?
%             if strcmp(fields,'xup')
%                 if dim == 2
%                     nfields = 5;
%                 elseif dim == 3
%                     nfields = 7;
%                 end
%             elseif strcmp(fields,'up')
%                 if dim == 2
%                     nfields = 3;
%                 elseif dim == 3
%                     nfields = 4;
%                 end
%             end
%             
%             
            % Write data
            for iel=1:nel
                
                for ifld = 1:dim
                    fwrite(fid,EL(iel).GLL(:,ifld),precision);%1:N^dim
                end
            end
            for iel=1:nel
                for ifld = 1:dim
                    fwrite(fid,EL(iel).Vinit(:,ifld),precision);%1:N^dim
                end
            end
            for iel=1:nel
                fwrite(fid,EL(iel).Pinit,precision);%1:N^dim
            end
            status = fclose(fid);
        end
        %*****************************************************************
        function [EL,N1,time,nel]=readfile(filename)
            %fields='xup';precision='*float32 and *float64';
            [b,hdr,tag,N1,nel,Lvar,sts,Dim] = nekmesh.readnek_mpi('le',filename);
            tt=sscanf(hdr,'%*s %f %f %f %f %f %f %f %f %f %f %*s');
            time=tt(7);
            fields=hdr(84:86);
            %Loop over all elements
            if fields(1)=='X'
                for i = 1:nel
                    EL(i).GLL(:,1)   = squeeze(reshape(b{i}(:,:,1),N1^2,1,1));
                    EL(i).GLL(:,2)   = squeeze(reshape(b{i}(:,:,2),N1^2,1,1));
                    if(fields(2)=='U')
                        EL(i).Vinit(:,1) = squeeze(reshape(b{i}(:,:,3),N1^2,1,1));
                        EL(i).Vinit(:,2) = squeeze(reshape(b{i}(:,:,4),N1^2,1,1));
                    end
                    if(fields(3)=='P')
                        EL(i).Pinit(:,1) = squeeze(reshape(b{i}(:,:,5),N1^2,1,1));
                    end
                    EL(i).fields=fields;
                    EL(i).std=hdr(6);
                    if Dim==2
                        EL(i).Vinit(:,3) =zeros(size(EL(i).Vinit(:,1)));
                    end
                    
                    
                end
            else
                for i = 1:nel
                    EL(i).Vinit(:,1) = squeeze(reshape(b{i}(:,:,1),N1^2,1,1));
                    EL(i).Vinit(:,2) = squeeze(reshape(b{i}(:,:,2),N1^2,1,1));
                    if(strfind(fields,'P'))
                        EL(i).Pinit(:,1) = squeeze(reshape(b{i}(:,:,3),N1^2,1,1));
                    end
                    EL(i).fields=fields;
                    EL(i).std=hdr(6);
                    if Dim==2
                        EL(i).Vinit(:,3) =zeros(size(EL(i).Vinit(:,1)));
                    end
                    
                end
            end
        end
        function [EL,N1,time,nel,std]=myread(fname)
            %[EL,N1,time,nel,std]=myread(fname)
            [fid,message] = fopen(fname,'r','ieee-le');
            header= fread(fid,132,'*char')';
            tt=sscanf(header,'%*s %f %f %f %f %f %f %f %f %f %f %*s');time=tt(7);
            tag = fread(fid,1,'*float32' );
            std=header(6);
            if str2num(header(6))==4;precision='*float32';else precision='*float64';end
            if str2num(header(14:15))==1;dim =2;else dim=3; end
            N1 = str2double(header(8:9));
            nel = str2double(header(20:30));
            fields=header(84:86);
            var=zeros(1,3);
            if(~isempty(strfind(fields,'X')));var(1)=1;end
            if(~isempty(strfind(fields,'U')));var(2)=1;end
            if(~isempty(strfind(fields,'P')));var(3)=1;end
            lglel = fread(fid,nel,'*int32');
            
            %             if var(1)==1;[xr,c1]=fread(fid,2*nel*N1^dim,precision);end
            %             if var(2)==1;[vr,c2]=fread(fid,2*nel*N1^dim,precision);end
            %             if var(3)==1;[pr,c3]=fread(fid,1*nel*N1^dim,precision);end
            %             disp(num2str(2*nel*N1^dim-c1))
            %             disp(num2str(2*nel*N1^dim-c2))
            %             disp(num2str(1*nel*N1^dim-c3))
            %
            if var(1) == 1
                EL(1:nel) = struct('GLL',zeros(N1^dim,dim),'Vinit',zeros(N1^dim,dim) ...
                    ,'Pinit',zeros(N1^dim,1));
            end
            if var(1)==1
                for iel=1:nel
                    for idim=1:dim
                        tmp=double(fread(fid,N1^dim,precision));
                        EL(lglel(iel)).GLL(:,idim)=tmp;
                    end
                end
            end
            if var(2)==1
                
                for iel=1:nel
                    for idim=1:dim
                        tmp=double(fread(fid,N1^dim,precision));
                        EL(lglel(iel)).Vinit(:,idim)=tmp;
                    end
                end
            end
            %             ss=0;
            if var(3)==1
                for iel=1:nel
                    tmp=double(fread(fid,N1^dim,precision));
                    EL(lglel(iel)).Pinit=tmp;
                end
            end
            for ii=1:nel;EL(ii).fields=fields;end
        end
        %*****************************************************************
        function [xr,bm,x,y,bmm,EM,time,indt]=readmatstr(casename,is,ie,istep,pr,varargin)
            %[xr,bm,x,y,bmm,EM,time,indt]=readstruc(casename,is,ie,istep,pr,varargin)
            %varargin{1}=cut  based on the x locations in the wing
            %varargin{2}
            pdd=pwd;
            eval(['cd   '  pr])
            switch nargin
                case 5
                    cut = [ 1 1];
                case 6
                    cut = varargin{1};
            end
            nsnap=length([is:istep:ie]);
            [bmm,EM,~,std]=nekmesh.easyread(['bm1' casename '0.f00001']);
            [nx,ny,N1]=nekmesh.meshinfo(EM);

            ix=abs(bmm.x(1,:)-cut(1));
            iy=abs(bmm.x(1,:)-cut(2));
            [~,ind(1)]=min(ix);
            [~,ind(2)]=min(iy);
            indt(1)=min(ind);
            indt(2)=max(ind);
            if isnan(cut(1))
                indt(1)=1;
             end
            if isnan(cut(2))
                indt(2)=length(bmm.x(1,:));
            end
            
            bm1=bmm.u(:,indt(1):indt(2));
            x=bmm.x(:,indt(1):indt(2));
            y=bmm.y(:,indt(1):indt(2));
            bm=repmat(bm1(:),2,1);
            n=size(bm1);n=[n 2 nsnap];
           
            fid=fopen('nekbin.i','w','n');
            fprintf(fid,[casename '\n']);
            fprintf(fid,['out.dat' '\n']);
            fprintf(fid,[num2str(nx) '\n']);
            fprintf(fid,[num2str(ny) '\n']);
            fprintf(fid,[num2str(is) '\n']);
            fprintf(fid,[num2str(ie) '\n']);
            fprintf(fid,[num2str(istep) '\n']);
            fprintf(fid,[num2str(indt(1)) '\n']);
            fprintf(fid,[num2str(indt(2)) '\n']);
            fprintf(fid,[num2str(N1) '\n']);
            fclose(fid);
            if str2double(std)==4
                !read_snap_s
            else
                !read_snap_d
            end
            up=zeros(ny*(N1-1)+1,nx*(N1-1)+1);
            up=up(:,indt(1):indt(2));
            xs=zeros(1,numel(up)*2);
            xr=zeros(numel(xs),nsnap);
            
            [fid,~] = fopen('out.dat','r','ieee-le.l64');
            time=fread(fid,nsnap,'float32');
            if str2double(std)==4
                xr=fread(fid,numel(xr),'float32');
            else
                xr=fread(fid,numel(xr),'float64');
            end
            xr=reshape(xr,length(xr)/nsnap,nsnap);
            fclose(fid);
            eval(['cd ' pdd]);
        end
        function [xr,xvec,bm1,tt]=readnek(casename,is,ie,istep,pr,varargin)
            %[xr,xvec,bm1,tt]=readnek(casename,is,ie,istep,pr,varargin)
            % pr ex. pr = './'
            % xr; snapshots
            %xvec = bm(:,1);
            %bm1 = bm1 
            %varagin{1} = cut
            pdd=pwd;
            eval(['cd   '  pr])
            switch nargin
                case 5
                    cut = [ 1 1];
                case 6
                    cut = varargin{1};
            end
            
            nsnap=length([is:istep:ie]);
            fname_out='out_file';
            fid=fopen('nekbin.i','w','n');
            fprintf(fid,[ '''' casename ''''       '             #p1=casename\n']);
            fprintf(fid,[num2str(is)    '             #p2=is (first snapshot number) \n']);
            fprintf(fid,[num2str(ie)    '             #P3=nsnap (last snapshot number) \n']);
            fprintf(fid,[num2str(istep) '             #p4=istep \n']);
            fprintf(fid,[num2str(cut(1)) '	 #p8=first number of element\n']);
            fprintf(fid,[num2str(cut(2)) '	 #p9=last number of element\n' ]);
            fprintf(fid,'''./''     #p11=inpath\n');
            fprintf(fid,'''./''      #p12=outpath\n');
            fclose(fid);
            % if ~isempty(dir) ;eval(['!mv nekbin.i  ' dir ]);eval(['cd   '  dir]);end
            !readnekfile
            
            fid = fopen(fname_out,'r');
            nel=fread(fid,1,'int');
            dim=fread(fid,1,'int');
            N=fread(fid,1,'int');
            tt=fread(fid,nsnap,'float32');
            ntot=N^dim*dim*nel;
            xr  =reshape(fread(fid,ntot*nsnap,'float64'),ntot,nsnap);
            bm  =reshape(fread(fid,ntot*2,'float64'),ntot,2);
            bm1 =bm(:,2);
            xvec=bm(:,1);
            fclose(fid);
            %!rm nekbin.i
            eval(['cd   '  pdd])
        end
        function [xr,xvec,bm1,tt,nelc,dim,N]=readpod(casename,is,ie,istep,pr,readelem,varargin)
            %[xr,xvec,bm1,tt]=readnek(casename,is,ie,istep,pr,varargin)
            % pr ex. pr = './'
            % xr; snapshots
            %xvec = bm(:,1);
            %bm1 = bm1
            %varagin{1} = cut
            
            pdd=pwd;
            eval(['cd   '  pr])
            switch nargin
                case 6
                    cut = [ 1 1];
                case 7
                    cut = varargin{1};
            end
            nsnap=length([is:istep:ie]);
            fname_out='out_file';
            fid=fopen('nekbin.i','w','n');
            fprintf(fid,[ '''' casename ''''       '             #p1=casename\n']);
            fprintf(fid,[num2str(is)    '             #p2=is (first snapshot number) \n']);
            fprintf(fid,[num2str(ie)    '             #P3=nsnap (last snapshot number) \n']);
            fprintf(fid,[num2str(istep) '             #p4=istep \n']);
            fprintf(fid,[num2str(cut(1)) '	 #p8=first number of element\n']);
            fprintf(fid,[num2str(cut(2)) '	 #p9=last number of element\n' ]);
            fprintf(fid,'''./''     #p11=inpath\n');
            fprintf(fid,'''./''      #p12=outpath\n');
            fprintf(fid,[num2str(readelem) '	 #number of elem\n']);
            fclose(fid);
            % if ~isempty(dir) ;eval(['!mv nekbin.i  ' dir ]);eval(['cd   '  dir]);end
%             !pod_read
           
            fid = fopen('out_file','r');
            nelc=fread(fid,1,'int');
            dim=fread(fid,1,'int');
            N=fread(fid,1,'int');
            tt=fread(fid,nsnap,'float32');
            ntot=N^dim*dim*nelc;
            xr  =reshape(fread(fid,ntot*nsnap,'float64'),ntot,nsnap);
            bm  =reshape(fread(fid,ntot*2,'float64'),ntot,2);
            
            bm1 =bm(:,2);
            xvec=bm(:,1);
            fclose(fid);
            
            
        end
        function [EL,N1,time,nel,std]=myreadpod(fname,varargin)
            %[EL,N1,time,nel,std]=myread(fname,varargin)
            %if varargin = 0 skape the coordinate info
            if isempty(varargin)
                varargin{1} = 1;
            end
            
            [fid,message] = fopen(fname,'r','ieee-le');
            header= fread(fid,132,'*char')';
            tt=sscanf(header,'%*s %f %f %f %f %f %f %f %f %f %f %*s');time=tt(7);
            tag = fread(fid,1,'*float32' );
            std=header(6);
            if str2num(header(6))==4;precision='*float32';else precision='*float64';end
            if str2num(header(14:15))==1;dim =2;else dim=3; end
            N1 = str2double(header(8:9));
            nel = str2double(header(20:30));
            fields=header(84:86);
            var=zeros(1,3);
            if(~isempty(strfind(fields,'X')));var(1)=1;end
            if(~isempty(strfind(fields,'U')));var(2)=1;end
            if(~isempty(strfind(fields,'P')));var(3)=1;end
            lglel = fread(fid,nel,'*int32');
            
            %             if var(1)==1;[xr,c1]=fread(fid,2*nel*N1^dim,precision);end
            %             if var(2)==1;[vr,c2]=fread(fid,2*nel*N1^dim,precision);end
            %             if var(3)==1;[pr,c3]=fread(fid,1*nel*N1^dim,precision);end
            %             disp(num2str(2*nel*N1^dim-c1))
            %             disp(num2str(2*nel*N1^dim-c2))
            %             disp(num2str(1*nel*N1^dim-c3))
            if var(1) == 1
                EL(1:nel) = struct('GLL',zeros(N1^dim,dim),'Vinit',zeros(N1^dim,dim) ...
                    ,'Pinit',zeros(N1^dim,1));
            end
            
            %%%%%%%%%%%%%%%%%%%
            
            if varargin{1} == 0
                tpp = zeros(N1^dim*dim*nel,1);
                %     fseek(fid, N1^dim*dim*nel*std, 'cof');
                tpp = fread(fid,N1^dim*dim*nel,precision);
                clear tpp
            else
                if var(1)==1
                    for iel=1:nel
                        for idim=1:dim
                            tmp=double(fread(fid,N1^dim,precision));
                            EL(lglel(iel)).GLL(:,idim)=tmp;
                        end
                    end
                end
            end
            
            
            if var(2)==1
                for iel=1:nel
                    for idim=1:dim
                        tmp=double(fread(fid,N1^dim,precision));
                        EL(lglel(iel)).Vinit(:,idim)=tmp;
                    end
                end
            end
            %             ss=0;
            var(3) = 0;
            if var(3)==1
                for iel=1:nel
                    tmp=double(fread(fid,N1^dim,precision));
                    EL(lglel(iel)).Pinit=tmp;
                end
            end
            for ii=1:nel;EL(ii).fields=fields;end
            
        end
        %*****************************************************************
        function [data1,header,tag,N,nel,nfields,status,dim] = readnek_mpi(endian,fname)
            
            %  This function reads binary data from the nek5000 new file format
            
            %  open the file
            %             precision='*float64';
            
            [fid,message] = fopen(fname,'r',['ieee-' endian]);
            
            if fid == -1, disp(message), return, end
            
            %  read the header
            
            header= fread(fid,132,'*char')';
            
            std=header(6);
            if str2num(std)==4
                precision='*float32';
            else
                precision='*float64';
            end
            fields=header(84:86);
            nnz=header(14:15);
            if str2num(nnz)==1
                dim =2;
            else
                dim=3;
            end
            N = str2double(header(8:9));
            nel = str2double(header(20:30));
            
            tag = fread(fid,1,'*float32' );
            
            lglel = fread(fid,nel,'*int32');
            
            var=zeros(1,3);
            if(~isempty(strfind(fields,'X')))
                var(1)=1;
            end
            if(~isempty(strfind(fields,'U')))
                var(2)=1;
            end
            if(~isempty(strfind(fields,'P')))
                var(3)=1;
            end
            nfields = var(1)*dim+var(2)*dim+var(3);
            
            data = zeros(nel,N^dim,nfields);
            
            %  read the data
            ss=0;
            for ivar=1:sum(var(1:2))
                for iel=1:nel
                    for icomp=1:dim
                        [data(lglel(iel),:,(ivar-1)*dim+icomp) count]= fread(fid,N^dim,precision);
                        ss=ss+count;
                    end
                end
            end
            if var(3)==1
                for iel=1:nel
                    [data(lglel(iel),:,nfields) count] = fread(fid,N^dim,precision);
                    ss=ss+count;
                end
            end
            
            
            status = fclose(fid);
            
            if dim==2
                for iel=1:nel
                    dt{lglel(iel)} = data(lglel(iel),:,:);
                    data1{lglel(iel)} = reshape(dt{lglel(iel)},N,N,nfields);
                end
            elseif dim==3
                for iel=1:nel
                    dt{lglel(iel)} = data(lglel(iel),:,:);
                    data1{lglel(iel)} = reshape(dt{lglel(iel)},N,N,N,nfields);
                end
            end
            
            
        end
        %******************************************************************
        function [nnx,nny,N1,nelm,elm]=meshinfo(EL)
            %[nnx nny N1 nelm,elm]=meshinfo(EL)
            nelm=length(EL);
            try
                N1=sqrt(size(EL(1).GLL,1));
            catch
                try
                    N1=sqrt(size(EL(1).Vinit,1));
                end
                N=[];
            end
            ielem=0;pp=1;ny=0;nx=0;
            
            while(nx*ny+ielem<nelm)
                
                nx=nelemx(ielem+1,nelm,N1);
                ny=nelemy(ielem+nx,nx,N1,nelm);
                nnx(pp)=nx;nny(pp)=ny;
                pp=pp+1;ielem=ielem+nx*ny;
                nx=0;ny=0;
            end
            elm(1)=0;
            for ii=1:length(nnx)
                elm(ii+1)=nnx(ii)*nny(ii);
            end
            elm=cumsum(elm);
            
            function ny=nelemy(is,nx,N1,nelm)
                %      -------  -------   --------
                %          xm1  xim xip   xp1
                %
                
                ny=0;
                for i=is:nx:nelm
                    xip=EL(i+0).GLL(N1,1);
                    xim=EL(i+0).GLL(1,1);
                    if(i>=nelm)
                        ny=ny+1;
                        break;
                    end
                    xp1=EL(i+1).GLL(1,1);
                    xm1=EL(i-1).GLL(N1,1);
                    if nx==1
                        if( abs(xip-xp1)> eps && abs(xim-xm1) > eps)
                            ny=ny+1;
                        else
                            disp('next box');
                            break;
                        end
                    else
                        if( abs(xip-xp1)> eps && abs(xim-xm1) < eps)
                            ny=ny+1;
                        elseif ( abs(xip-xp1) < eps && abs(xim-xm1) < eps)
                            disp('next box');
                            break;
                        else ( abs(xip-xp1)> eps || abs(xim-xm1) < eps);
                            disp('next box');
                            break;
                        end
                    end
                end
            end
            function nx=nelemx(is,ie,N1)
                nn=1;
                for i=is:ie
                    %                     disp(num2str(i))
                    try
                        x1=EL(i).GLL(N1,1);
                        x2=EL(i+1).GLL(1,1);
                        y1=EL(i).GLL(N1,2);
                        y2=EL(i+1).GLL(1,2);
                        if (abs(x2-x1) > eps && abs(y2-y1) > eps)
                            nx=nn;
                            break
                        end
                        nn=nn+1;
                    catch
                        if nn==ie-is+1
                            nx=nn;
                            break;
                        end
                        x1=EL(i).nodes(1,2);
                        x2=EL(i+1).nodes(1,1);
                        if abs(x2-x1) > eps
                            nx=i;
                            break
                        end
                    end
                end
                
            end
        end
        %*****************************************************************
        function [xp,yp,nelm,N1,nx,ny]=meshmap(EL)
            %[xp yp nelm N1 nx ny]=meshmap(EL)
            [nx ny N1 nelm]=nekmesh.meshinfo(EL);
            for iii=1:length(nx)
                st=1:nelm;
                st=reshape(st,nx,ny)';
                xp=zeros(ny+1,nx+1);
                yp=zeros(ny+1,nx+1);
                ppp=0;
                n=1;
                for i=1:size(st,1)
                    for j=1:size(st,2)
                        
                        try
                            xp(i,j)=EL(n).GLL(1,1);
                            xp(i,j+1)=EL(n).GLL(N1,1);
                            xp(i+1,j+1)=EL(n).GLL(N1^2,1);
                            xp(i+1,j)=EL(n).GLL(N1*(N1-1)+1,1);
                            
                            
                            yp(i,j)=EL(n).GLL(1,2);
                            yp(i,j+1)=EL(n).GLL(N1,2);
                            yp(i+1,j+1)=EL(n).GLL(N1^2,2);
                            yp(i+1,j)=EL(n).GLL(N1*(N1-1)+1,2);
                            ppp=1;
                            
                        catch
                            xp(i,j)=EL(n).nodes(1,1);
                            xp(i,j+1)=EL(n).nodes(1,2);
                            xp(i+1,j+1)=EL(n).nodes(1,3);
                            xp(i+1,j)=EL(n).nodes(1,4);
                            
                            
                            yp(i,j)=EL(n).nodes(2,1);
                            yp(i,j+1)=EL(n).nodes(2,2);
                            yp(i+1,j+1)=EL(n).nodes(2,3);
                            yp(i+1,j)=EL(n).nodes(2,4);
                            
                            ppp=0;
                        end
                        
                        n=n+1;
                    end
                end
                if(ppp==1)
                    disp('wrote according ot EL.GLL')
                else
                    disp('wrote according ot EL.nodes')
                end
            end
        end
        %*****************************************************************
        function[up,vp,pp,xp,yp]=results(EL,nx,ny,N1)
            n=N1-1;
            xp=zeros(ny*n+1,nx*n+1);
            yp=zeros(ny*n+1,nx*n+1);
            up=zeros(size(xp));
            vp=zeros(size(xp));
            pp=zeros(size(xp));
            numb=1;
            
            for ielm=1:ny
                for jelm=1:nx
                    ng=1;
                    for s=1:N1
                        for t=1:N1
                            ii=(ielm-1)*n+s;
                            jj=(jelm-1)*n+t;
                            if(~isempty(strfind(EL(1).fields,'X')))
                                xp(ii,jj)=EL(numb).GLL(ng,1);
                                yp(ii,jj)=EL(numb).GLL(ng,2);
                            end
                            %---------------------
                            up(ii,jj)=EL(numb).Vinit(ng,1);
                            vp(ii,jj)=EL(numb).Vinit(ng,2);
                            %------------------------
                            if(~isempty(strfind(EL(1).fields,'P')))
                                %                                 pp(ii,jj)=EL(numb).Pinit(ng,1);
                                %                                   [numb ng]
                                pp(ii,jj)=EL(numb).Pinit(ng);
                                
                            end
                            ng=ng+1;
                        end
                    end
                    numb=numb+1;
                end
            end
        end
        function box=results_wing(EL,nnx,nny,N1,elm)
            
            for ii=1:length(nnx)
                [up,vp,pp,xp,yp]=nekmesh.results(EL(elm(ii)+1:elm(ii+1)),nnx(ii),nny(ii),N1);
                box(ii).x=xp;
                box(ii).y=yp;
                box(ii).u=up;
                box(ii).v=vp;
                box(ii).p=pp;
            end
        end
        %*****************************************************************
        function bm1=weights(EL,nx,ny,N1)
            n=N1-1;
            bm1=zeros(ny*n+1,nx*n+1);
            numb=1;
            for ielm=1:ny
                for jelm=1:nx
                    ng=1;
                    for s=1:N1
                        for t=1:N1
                            ii=(ielm-1)*n+s;
                            jj=(jelm-1)*n+t;
                            bm1(ii,jj)=bm1(ii,jj)+EL(numb).Vinit(ng,1);
                            ng=ng+1;
                        end
                    end
                    numb=numb+1;
                end
            end
        end
        function EL=rev_results(up,vp,pp,EL)
            %EL=rev_results(up,vp,pp,EL)
            %up,vp,pp should be the size of mesh with GLL
            %this is used when you do not have duplicate point on the mesh
            N1=sqrt(length(EL(1).GLL(:,1)));
            n=N1-1;
            ny=(size(up,1)-1)/n;
            nx=(size(up,2)-1)/n;
            numb=1;
            for ielm=1:ny
                for jelm=1:nx
                    ng=1;
                    for s=1:N1
                        for t=1:N1
                            ii=(ielm-1)*n+s;
                            jj=(jelm-1)*n+t;
                            %---------------------
                            EL(numb).Vinit(ng,1)=up(ii,jj);
                            EL(numb).Vinit(ng,2)=vp(ii,jj);
                            %------------------------
                            EL(numb).Pinit(ng,1)=pp(ii,jj);
                            ng=ng+1;
                        end
                    end
                    numb=numb+1;
                end
            end
        end
        function EL=rev_coord(xp,yp,N1)
            n=N1-1;
            ny=(size(xp,1)-1)/n;
            nx=(size(xp,2)-1)/n;
            numb=1;
            for ielm=1:ny
                for jelm=1:nx
                    ng=1;
                    for s=1:N1
                        for t=1:N1
                            ii=(ielm-1)*n+s;
                            jj=(jelm-1)*n+t;
                            %---------------------
                            EL(numb).GLL(ng,1)=xp(ii,jj);
                            EL(numb).GLL(ng,2)=yp(ii,jj);
                            ng=ng+1;
                        end
                    end
                    numb=numb+1;
                end
            end
        end
        function EL=rev_nekresults(up,vp,pp,N1,EL)
            % EL=rev_nekresults(up,vp,pp,N1,EL)
            % this is used when you have duplicate point on the mesh
            n=N1;
            ny=size(up,1)/n;
            nx=size(up,2)/n;
            numb=1;
            for ielm=1:ny
                for jelm=1:nx
                    ng=1;
                    for s=1:N1
                        for t=1:N1
                            ii=(ielm-1)*n+s;
                            jj=(jelm-1)*n+t;
                            %---------------------
                            EL(numb).Vinit(ng,1)=up(ii,jj);
                            EL(numb).Vinit(ng,2)=vp(ii,jj);
                            EL(numb).Vinit(ng,3)=0;
                            %------------------------
                            EL(numb).Pinit(ng,1)=pp(ii,jj);
                            ng=ng+1;
                        end
                    end
                    numb=numb+1;
                end
            end
        end
        %*****************************************************************
        function [xp,yp]=meshplot(xp,yp,varargin)
            % [xp,yp]=meshplot(xp,yp,varargin)
            switch nargin
                case 1
                    Eg=xp;
                    clear xp
                    disp(fieldnames(Eg));
                    bp=nekmesh.easyread(Eg);
                    xp=bp(1).x;
                    yp=bp(1).y;
                    colormark='r';
                    big=2;
                case 2
                    colormark='r';
                    big=2;
                case 3
                    colormark=varargin{1};
                    big=2;
                case 4
                    colormark=varargin{1};
                    big=varargin{2};
            end
            
            plot(xp,yp,colormark, 'MarkerSize',big)
            hold on
            plot(xp',yp',colormark, 'MarkerSize',big)
%             axis equal
            
        end
        function gllplot(EL,N,mark,msize,fnum)
            % fnum is the figure number and msize is marker size, mark is
            % marker
            plotgll = zeros(length(EL)*(N+1)^2,2);
            for i = 1:length(EL)
                plotgll((i-1)*(N+1)^2+1:i*(N+1)^2,:) = EL(i).GLL;
            end
            figure(fnum)
            hold on
            
            
            plot(plotgll(:,1),plotgll(:,2),mark,'Markersize',msize)
            axis equal
            drawnow
        end
        function BCplot(Eg,FluidBC,is,ie)
            N1=sqrt(length(Eg(1).GLL(:,1)));
            for iel=is:ie
                
                M=FluidBC(4*(iel-1)+1:4*iel);
                x(1)=(Eg(iel).GLL(1,1)          +Eg(iel).GLL(N1,1))/2;
                y(1)=(Eg(iel).GLL(1,2)          +Eg(iel).GLL(N1,2))/2;
                
                x(2)=(Eg(iel).GLL(N1,1)         +Eg(iel).GLL(N1^2,1))/2;
                y(2)=(Eg(iel).GLL(N1,2)         +Eg(iel).GLL(N1^2,2))/2;
                
                x(3)=(Eg(iel).GLL(N1*(N1-1)+1,1)+Eg(iel).GLL(N1^2,1))/2;
                y(3)=(Eg(iel).GLL(N1*(N1-1)+1,2)+Eg(iel).GLL(N1^2,2))/2;
                
                x(4)=(Eg(iel).GLL(1,1)          +Eg(iel).GLL(N1*(N1-1)+1,1))/2;
                y(4)=(Eg(iel).GLL(1,2)          +Eg(iel).GLL(N1*(N1-1)+1,2))/2;
                
                for ii=1:4
                    text(x(ii),y(ii),M(ii));hold on
                end
            end
        end
        %*****************************************************************
        function [ELt,nelm,N1,nx,ny]=meshflip(filename,nc)
            EL=nekmesh.readfile(filename);
            [nelm N1 nx ny]=nekmesh.meshinfo(EL);
            %----- try to build reflected elment number from nx to 1
            stm=reshape([1:nelm],nx,ny);
            stm=stm(1:nc,:);
            stm=reshape(stm,numel(stm),1);
            
            st=reshape(fliplr(reshape(stm,nc,ny)')',nc*ny,1);
            
            for i=1:length(st)
                ELr(i)=EL(st(i));
            end
            
            %----- change the direction of the gll point numbering
            stn=reshape(fliplr(reshape(1:N1^2,N1,N1)')',N1*N1,1);
            
            for i=1:length(ELr)
                ELr(i).GLL(:,1)=ELr(i).GLL(stn,1);
                ELr(i).GLL(:,2)=-ELr(i).GLL(stn,2);
                %--------------------
                ELr(i).Vinit(:,1)= ELr(i).Vinit(stn,1);
                ELr(i).Vinit(:,2)=-ELr(i).Vinit(stn,2);
                
            end
            %-----------------------------------------
            % --------------------------  combine two elemnts
            [nelmr N1 nxr nyr]=nekmesh.meshinfo(ELr);
            st1=[1:nelmr];
            st1=reshape(st1,nxr,nyr);
            
            st2=[1:nelm];
            st2=reshape(st2,nx,ny);
            
            st=[st1;st2];
            ind=[ones(size(st1)) ; -1*ones(size(st2))];
            st=reshape(st,numel(st),1);
            ind=reshape(ind,numel(ind),1);
            %------- first and second part numering from each row and
            %---------second colum shows which element it belongs
            
            n=1;
            for i=1:nelm+nelmr
                if(ind(i)==1)
                    ELt(i)=ELr(st(i));
                else
                    ELt(i)=EL (st(i));
                end
            end
            %-----------debug-------------------------
            %             for i=1:N1^2
            %                 xg(i)=ELt(1).GLL(i,1);
            %                 yg(i)=ELt(1).GLL(i,2);
            %                 plot(xg,yg,'bo')
            %                 hold on
            %                 text(xg(i),yg(i),num2str(i));
            %             end
            %             %-----------to inversigat v vel on symetric line------------
            %             st=1:N1:N1*(N1-1)+1;
            %             nm=1:nx:nx*(ny-1)+1;
            %             for i=1:length(nm)
            %                 EL(nm(i)).Vinit(st,2)
            %             end
            %
            %             str=N1:N1:N1^2;
            %             nmr=nxr:nxr:nxr*nyr;
            %             for i=1:length(nm)
            %                 ELr(nmr(i)).Vinit(str,2)
            %             end
            %
            %
            %
            %             %---------- elment number shoes
            %             for i=1:length(nm)
            %                 st=[1 N1 N1*(N1-1)+1 N1^2];
            %                 xx=ELt(nm(i)).GLL(st,1);
            %                 yy=ELt(nm(i)).GLL(st,2);
            %                 plot(xx,yy,'b.')
            %                 hold on
            %                 axis equal
            %                 xa=mean(xx);
            %                 ya=mean(yy);
            %                 text(xa,ya,num2str(i))
            %             end
            %---------end debug---------------------
            nelm=nelm+nelmr;
            nx=nx+nxr;
        end
        %*****************************************************************
        function [ELr,box,EL]=cutmesh(filename,xs,xe,ye,last,f2)
            % in this function xs and xe is the start and end of cut,
            % filenmae in the original file and f2 is needed is a file with
            % gll point
            EL=nekmesh.readfile(filename);
            ic=0;
            mm=fieldnames(EL(1));
            for jj=1: length(mm)
                if strcmp(mm{jj},'Pinit')
                    ic=1;
                end
            end
            
            if  isempty(strfind(EL(1).fields,'X'))
                Eh=nekmesh.readfile(f2);
                %Eh is not a biger file
                for ii=1:length(Eh)
                    EL(ii).GLL=Eh(ii).GLL;
                    EL(ii).fields=Eh(ii).fields;
                    if ic
                        EL(ii).Pinit=Eh(ii).Pinit;
                    else
                        EL(ii).Pinit=0.*Eh(ii).Vinit(:,1);
                    end
                end
                disp('use another mesh to get spatial info')
            end
            if isempty(last)
                last='last';
            end
            [ELr box]=nekmesh.cut(EL,xs,xe,ye,last);
        end
        function [ELr] = cuti(EL,is,nx,ny)
            %[ELr] = cuti(EL,is)
            mm = 1 ;
            nn = 1 ;
            clear ELr
            for jj=1 : ny
                for ii = 1 : nx
                    if (ii >= is )
                        ELr(nn) = EL(mm) ;
                        nn = nn + 1 ;
                    end
                    mm = mm + 1 ;
                end
            end
        end
        function [ELr] = cutie(EL,id,nx,ny)
            %[ELr] = cuti(EL,is)
            mm = 1 ;
            nn = 1 ;
            clear ELr
            is = id(1);
            ie = id(2);
            for jj=1 : ny
                for ii = 1 : nx
                    if (ii <= ie && ii >= is)
                        ELr(nn) = EL(mm) ;
                        nn = nn + 1 ;
                    end
                    mm = mm + 1 ;
                end
            end
        end
        function [ELr,box]=cut(EL,xs,xe,ye,last)
            
            [nx,ny,N1]=nekmesh.meshinfo(EL);
            x=zeros(nx+1,1);
            for i=1:nx
                x(i)=EL(i).GLL(1,1);
            end
            x(nx+1)=EL(nx).GLL(N1,1);
            y=zeros(ny+1,1);
            nn=1;
            for ii=1:nx:length(EL)
                %                 xx(nn)=EL(ii).GLL(1,1);
                y(nn)=EL(ii).GLL(1,2);
                nn=nn+1;
            end
            y(nn)= EL(nx*(ny-1)+1).GLL(N1*(N1-1)+1,2);
            %             xx(nn)=EL(nx*(ny-1)+1).GLL(N1*(N1-1)+1,1);
            
            if isempty(xs);js=1;else js=find(x<xs,1,'last');end
            je=find(x<xe,1,last);
            is=1;
            if isempty(ye);ie=ny;else ie=find(abs(y)<abs(ye),1,'last');end
            
            n=1;
            m=1;
            for i=1:ny
                for j=1:nx
                    if(j>=js && j<je) && (i>=is && i<=ie)
                        ELr(n)=EL(m);
                        n=n+1;
                    end
                    m=m+1;
                end
            end
            [nnx,nny,N1,~,elm]=nekmesh.meshinfo(ELr);
            box=nekmesh.results_wing(ELr,nnx,nny,N1,elm);
            
        end
        %*****************************************************************
        function lamx=sponge(x,xstart,xend,delta_rise,delta_fall,lambda_max)
            %lamx=sponge(x,xstart,xend,delta_rise,delta_fall,lambda_max)
            for ii=1:length(x)
                xx=x(ii);
                lamx(ii)=sponges(xx,xstart,xend,delta_rise,delta_fall,lambda_max);
            end
            function laamm=sponges(x,xstart,xend,delta_rise,delta_fall,lambda_max)
                arg1 = (x-xstart)/delta_rise;
                arg2 = ((x-xend)/delta_fall) + 1;
                
                S1=smooth(arg1);
                S2=smooth(arg2);
                laamm = lambda_max*(S1-S2);
                function S=smooth(x)
                    if (x <= 0.)
                        S = 0.;
                    elseif (x > 0. && x < 1.)
                        S = 1./(1.+exp(1./(x-1.)+1./x));
                    elseif  (x >= 1)
                        S = 1.;
                    else
                        disp('x is wrrong');
                    end
                    
                end
            end
            
        end
        function b1=fringe_write(EE,xu,xb,casename,prefix)
            %b1=fringe_write(EE,xu,xb,casename,prefix)
            % write the fringe: on the top in u and on bottom in v velocity
            % inputs: EL
            %         xu: xlocation of the fringe at top
            %         xb: xlocation of the fring at bottom
            %         casename: the name of the output
            %         prefix: prefix of the file name
            %output: b1: box file
            [nnx,nny,N1,nelm,elm]=nekmesh.meshinfo(EE);
            b1=nekmesh.results_wing(EE,nnx,nny,N1,elm);
            %             [b1 EE elm]=nekmesh.easyread(filename);
            %             N1=sqrt(length(EE(1).Vinit(:,1)));
            for ii=1:length(EE);
                for jj=1:max(1,length(b1)-1)
                    EE(ii).Vinit(:,jj)=0.*EE(ii).Vinit(:,jj);
                end
            end
            
            for ii=1:max(1,length(b1)-1)
                FR(ii).fx=0.*b1(ii).u;
                FR(ii).fy=0.*b1(ii).v;
            end
            II=zeros(size(b1(1).x(1,:)));
            II(1:find(b1(1).x(1,:)<1e-2,1,'first'))=1;
            %------------------ top --------------------------------
            if ~isempty(xu)
                
                %                 ind=find((b1(1).x(1,:)>xu & b1(1).y(1,:) > b1(1).y(1,ind0) ) ,1,'first');
                ind=find(( b1(1).x(1,:)>xu & II <=0  ) ,1,'first');
                L=[0 cumsum(sqrt(diff(b1(1).x(1,ind:end)).^2+diff(b1(1).y(1,ind:end)).^2))];
                lamx=nekmesh.sponge(L,L(1),2*L(end),0.1*L(end),0.1*L(end),1.0);
                
                nn=1;
                for jj=ind:size(b1(1).x,2)
                    FR(1).fx(:,jj)=lamx(nn);
                    nn=nn+1;
                end
            end
            %----------------- bottom---------------------------------
            if ~isempty(xb)
                
                %                 ind=find((b1(1).x(1,:)>xb & b1(1).y(1,:) <b1(1).y(1,end) ) ,1,'last');
                ind=find(b1(1).x(1,:)<xb & II > 0 ,1,'first');
                clear L
                L=[0 cumsum(sqrt(diff(b1(1).x(1,1:ind)).^2+diff(b1(1).y(1,1:ind)).^2))];
                lam=nekmesh.sponge(L,-L(end),L(end),0.1*L(end),0.1*L(end),1.0);
                for ii=1:max(length(b1)-1,1)
                    nn=1;
                    for jj=1:ind
                        FR(ii).fy(:,jj)=lam(nn);
                        nn=nn+1;
                    end
                end
            end
            %             pcolor(b1(1).x,b1(1).y,FR(1).fy);shading interp;hold on
            %             pcolor(b1(2).x,b1(2).y,FR(2).fx);shading interp
            %             %----------- results in EE-------------------------------
            for ii=1:max(length(b1)-1,1)
                EE(elm(ii)+1:elm(ii+1))=nekmesh.rev_results(FR(ii).fx,FR(ii).fy,0.*FR(ii).fx,EE(elm(ii)+1:elm(ii+1)));
            end
            nekmesh.writeic(EE,casename,prefix,8);
        end
        %*****************************************************************
        function [xt,yt]=correct_mesh(xt,yt,xyf,acc,row)
            %             [~,imin]=min(xyf(:,1));
            % plot(xyf(:,1),xyf(:,2))
            for ii=2:size(xt,2)-1
                %                 if ii>=223
                %                     jjjj=1;
                %                 end
                %                 if yt(1,ii)<xyf(imin,2)
                ind1=find((abs(xyf(:,1) - xt(row,ii))<acc & abs(xyf(:,2) - yt(row,ii))<acc) ,1,'first');
                ind2=find((abs(xyf(:,1) - xt(row,ii))<acc & abs(xyf(:,2) - yt(row,ii))<acc) ,1,'last');
                if(isempty(ind1) || isempty(ind2));disp('ind1 or ind2 is empty');end
                %
                %    else
                %                     ind1=find( (xyf(:,1) < xt(row,ii)) & (abs(xyf(:,1) - xt(row,ii))<acc) & ...
                %                         (abs(xyf(:,2) - yt(row,ii))<acc) ,1,'first');%& (xyf(:,2) < yt(row,ii)),1,'first');
                %                     ind2=find( (xyf(:,1) > xt(row,ii)) & (abs(xyf(:,1) - xt(row,ii))<acc) & ...
                %                         (abs(xyf(:,2) - yt(row,ii))<acc) ,1,'last');
                %                 end
                %      ind1=find( (xyf(:,1) < xt(row,ii)) & (abs(xyf(:,1) - xt(row,ii))<acc) & ...
                %                         (abs(xyf(:,2) - yt(row,ii))<acc) ,1,'first');%& (xyf(:,2) < yt(row,ii)),1,'first');
                %                     ind2=find( (xyf(:,1) > xt(row,ii)) & (abs(xyf(:,1) - xt(row,ii))<acc) & ...
                %                         (abs(xyf(:,2) - yt(row,ii))<acc) ,1,'last');
                
                
                %     plot(xyf(:,1),xyf(:,2));hold on ;
                %     plot(xyf(ind1,1),xyf(ind1,2),'ro',xyf(ind2,1),xyf(ind2,2),'bs')
                %     plot(xt(row,ii),yt(row,ii),'g*')
                %     nekmesh.meshplot(xt,yt,'r',2)
                %     plot(box.x(1,ii),box.y(1,ii),'g^')
                %     plot(xt(row+1,ii),yt(row+1,ii),'bs',xt(row-1,ii),yt(row-1,ii),'bs')
                if (isempty(ind1) || isempty(ind2));disp('ind is empty') ; end
                lenn=[];lenn2=[];
                lenn=0.*[ind1:ind2];
                lenn2=0.*[ind1:ind2];
                nn=1;
                for jj=ind1:ind2
                    lenn(nn)=sqrt( (xt(row+1,ii)-xyf(jj,1)).^2 + ...
                        (yt(row+1,ii)-xyf(jj,2)).^2 );
                    
                    if row>1
                        lenn2(nn)=sqrt( (xt(row-1,ii)-xyf(jj,1)).^2 + ...
                            (yt(row-1,ii)-xyf(jj,2)).^2 );
                    end
                    nn=nn+1;
                end
                if row>1
                    [~,iimin1]=min(lenn(:)+lenn2(:));
                else
                    [~,iimin1]=min(lenn(:));
                end
                iimin=iimin1+ind1-1;
                %     plot(xyf(iimin,1),xyf(iimin,2),'ks')
                xt(row,ii)=xyf(iimin,1);
                yt(row,ii)=xyf(iimin,2);
                disp(['# elem number in arc ' num2str(ii) 'in row ' num2str(row)])
            end
        end
        %*****************************************************************
        function [Arc,box,npx,npy,np2]=arccalc(box,Arc,method)
            xt=box(1).x;
            yt=box(1).y;
            npx=size(xt,2)-1;
            npy=size(xt,1)-1;
            np2=find(xt(1,:)==6,1);
            if (isempty(np2))
                disp('you do not have a node exactly on x=6');
            end
            %---- from x=0 to x(np) use arc data or lars data with a high
            %-----accuracy up to 10000 point----------------------------
            np=6; % default value change
            disp('used np=6 as a default value otherwise change the m file');
            for i=2:npy
                x1=xt(i,1:np);
                y1=yt(i,1:np);
                if(strcmp(method,'spline'))
                    %--------------- using spline------------
                    cx(i,:)=linspace(x1(1),x1(end),10000)';
                    cy(i,:)=interp1(x1,y1,cx(i,:),'cubic');
                    %------------------using lars data--------
                else
                    cx(i,:)=linspace(x1(1),x1(end),10000)';
                    cy(i,:)=interp1(Arc.x(i,:),Arc.y(i,:),cx(i,:),'cubic');
                end
                cy(i,1)=0;
            end
            cx(:,end)=[];
            cy(:,end)=[];
            
            % -------just used sp line and forget lars data up to np2------
            for i=2:npy
                x1=xt(i,np:np2);
                y1=yt(i,np:np2);
                %--------------- using spline------------
                cx2(i,:)=linspace(x1(1),x1(end),40001)';
                cy2(i,:)=interp1(x1,y1,cx2(i,:),'cubic');
            end
            cx=[cx cx2];
            cy=[cy cy2];
            for i=2:npy
                cl(i,:)=[0 cumsum(sqrt(diff(cx(i,:)).^2+diff(cy(i,:)).^2))];
            end
            
            %boundary case
            xy_fine = nekmesh.getgeolead(60,1,6,[]);
            cx(1,:)=xy_fine(:,1);
            cy(1,:)=xy_fine(:,2);
            cl(1,:)=[0 cumsum(sqrt(diff(cx(1,:)).^2+diff(cy(1,:)).^2))];
            
            Arc.x=cx;
            Arc.y=cy;
            Arc.L=cl;
            %------------ plotting
            %             for i = 1:npy
            %                plot(Arc.x(i,:),Arc.y(i,:),'m')
            %                hold on
            %             end
            %----------- we donot correct for the first row but others-----
            [xt,yt,Arc]=leadboxcorect(npy,np2,xt,yt,Arc);%nekmesh.meshplot(xt,yt,'b',2)
            box(1).x=xt;
            box(1).y=yt;
            function [xt,yt,Arc]=leadboxcorect(npy,np2,xt,yt,Arc)
                Arc.y(:,1)=0;
                Arc.x(1,1)=0;
                for i=1:npy
                    for j=2:np2-1
                        if(i==1)
                            %correct the grid point based on smallest
                            %distance between arc point and a line
                            v1=[ xt(i,j) yt(i,j) 0]';
                            v2=[xt(i+1,j) yt(i+1,j) 0]';
                            
                            tt2=linspace(xt(i,j-1),xt(i,j+1),1000);
                            y2=interp1(Arc.x(i,:),Arc.y(i,:),tt2,'cubic');
                            
                            
                            
                            a = v1 - v2;
                            d(1)=100;
                            for ii=2:length(tt2)
                                pt=[tt2(ii) y2(ii) 0]';
                                b = pt - v2;
                                d(ii) = norm(cross(a,b)) / norm(a);
                            end
                            
                            [~,jj]=min(d);
                            xt(i,j)=tt2(jj);
                            yt(i,j)=y2(jj);
                        end
                    end
                end
                xt(1,1)=0;
                yt(:,1)=0;
                
            end
        end
        function [Arc,box,npx,npy,len]=arccalc_wing(box,np2,xyf)
            app=2e6;
            acc_mesh_correct=1e-3;
            acc_point=1e-4;
            xt=box(1).x;
            yt=box(1).y;
            npx=size(xt,2)-1;
            npy=size(yt,1)-1;
            
            %             [~,xyf]=nekmesh.wing();
            is=find(abs(xyf(:,1)-xt(1,1))<acc_point,1,'first');
            ie=find(abs(xyf(:,1)-xt(1,end))<acc_point,1,'last');
            tt=linspace(xyf(is,3),xyf(ie,3),app);
            cx(:,1)=interp1(xyf(is:ie,3),xyf(is:ie,1),tt,'cubic');
            cy(:,1)=interp1(xyf(is:ie,3),xyf(is:ie,2),tt,'cubic');
            cl(:,1)=[0;cumsum(sqrt(diff(cx(:,1)).^2+diff(cy(:,1)).^2))];
            [xt,yt]=nekmesh.correct_mesh(xt,yt,[cx(:,1) cy(:,1)],acc_mesh_correct,1);
            
            len=0.*xt;
            for i=2:np2+1
                x1=xt(i,:);
                y1=yt(i,:);
                t=[0 cumsum(sqrt(diff(x1).^2+diff(y1).^2))];
                len(i,:)=t;
                tt=linspace(0,t(end),app);
                cx(:,i)=csaps(t,x1,1,tt);
                cy(:,i)=csaps(t,y1,1,tt);
                cl(:,i)=[0;cumsum(sqrt(diff(cx(:,i)).^2+diff(cy(:,i)).^2))];
                [xt,yt]=nekmesh.correct_mesh(xt,yt,[cx(:,i) cy(:,i)],acc_mesh_correct,i);
            end
            
            len=0.*xt;
            for i=1:np2+1
                for j=1:size(xt,2)
                    is=find(( abs(cx(:,i)-xt(i,j))<acc_point & abs(cy(:,i)-yt(i,j))<acc_point ),1,'first');
                    if isempty(is);disp('start index cannot be find');end
                    %                     if j==337
                    %                         rr=2
                    %                     end
                    xt(i,j)=cx(is,i);
                    yt(i,j)=cy(is,i);
                    len(i,j)=cl(is,i);
                end
            end
            
            Arc.x=cx;
            Arc.y=cy;
            Arc.L=cl;
            
            box(1).x=xt;
            box(1).y=yt;
            
            
        end
        %******************************************************************
        function [u1,v1]=shape(um,st)
            u=um(1:end/2,:)    ;
            u1 = reshape(u,st);
            v=um(end/2+1:end,:);
            v1 = reshape(v,st);
        end
        function veiwsnap(x,xp,yp,vel,titl,shad)
            % up,xp,yp,'u',title'(null for no),shadeing
            if(strcmp(vel,'u'))
                [uu,~]=nekmesh.shape(x,size(xp));
            else
                [~,uu]=nekmesh.shape(x,size(xp));
            end
            pcolor(xp,yp,uu);
            if isempty(shad)
                shading interp;
            end
            title(titl);colorbar;drawnow
        end
        %******************************************************************
        function [var,tt,hpts]=hpts(filename,ipoint,ivar)
            %[var tt hpts]=hpts(filename,ipoint,ivar)
            %inputs: filename,ipoint,ivar
            %outpots: ivar, time,hpts (optional) data in hpts is organized as
            %hpts(ivar+1,ipoint,:) this gives the information of the ivar
            %exluding the time and ipoint in different time
            
            M=importdata(filename,' ',1);
            data=M.data;
            npoints=find(abs(data(:,1)-data(1,1))<1e-7,1,'last');
            ntime=ceil(length(data(:,1))/npoints);
            nvar=size(data,2);
            hpts=reshape(data',nvar,npoints,ntime);
            for ii=1:length(ipoint)
                var(:,ii)=squeeze(hpts(ivar+1,ipoint(ii),:));
            end
            tt=squeeze(hpts(1,1,:));
            %             A =[
            %
            %      1     2     3     4
            %      1     0    -3    -4
            %      1    -5    -2    -1
            %      2    -8    -9    -7
            %      2    -5    -6    -7
            %      2     0     1     0];
            
        end
        %****************************************************************** 
        function [um,vm,wm,ym,zm]=replicate(zff,yff,ut,vt,wt,np)
            %[um,vm,wm,ym,zm]=replicate(zff,yff,ut,vt,wt,np)
            % replicate will replicate the width=+-np*Lz 
            [zt,yt]=meshgrid(zff,yff);
            Ly=2*max(yff);
            rep=[-np:np];
            for ii=1:length(rep)
                E(ii).y=yt+rep(ii)*Ly;
                E(ii).z=zt;
                E(ii).u=ut;
                E(ii).v=vt;
                E(ii).w=wt;
                
            end
            
            nn=1;
            mm=length(E(ii).y(1,:));
            for ii=1:length(E)
                for jj=1:mm-1
                    ym(nn,:)=E(ii).y(jj,:);
                    zm(nn,:)=E(ii).z(jj,:);
                    um(nn,:)=E(ii).u(jj,:);
                    vm(nn,:)=E(ii).v(jj,:);
                    wm(nn,:)=E(ii).w(jj,:);
                    
                    nn=nn+1;
                end
            end
            ym(nn,:)=E(end).y(end,:);
            zm(nn,:)=E(end).z(end,:);
            um(nn,:)=E(end).u(end,:);
            vm(nn,:)=E(end).v(end,:);
            wm(nn,:)=E(end).w(end,:);
            
            
        end
        %simson stuff
        %******************************************************************
        %                        simson files
        %**************************  read  ********************************
        function [velu,velv,velw,xmesh,ymesh,zmesh,ww,velphys,t,wy,wwy]=readsim(filename)
            %[velu,velv,velw,xmesh,ymesh,zmesh,ww,velphys,t,wy,wwy]=readsim(filename)
            % read a velocity field in Fourier space as output by
            % bla
            % in variables u,v,w.
            %
            % this is the version that uses the vector form of fread
            % this is ten times faster than the old version
            %
            % the velocity components u, v, w are concatenated
            % on the third dimension of the output vel
            %
            % the odd ball is removed
            %
            %%%%%% boundary layer flow
            verbose=0;
            %disp(filename)
            fid=fopen(filename,'r','ieee-be.l64');
            eol=fread(fid,1,'int');
            Re=fread(fid,1,'float64');
            pou=fread(fid,1,'int');
            length=fread(fid,4,'float64');
            eol=fread(fid,2,'int');
            storl=fread(fid,4,'int');
            eol=fread(fid,2,'int');
            flowtype=fread(fid,1,'int');
            dstar=fread(fid,1,'float64');
            eol=fread(fid,1,'int');
            if flowtype==-1
                eol=fread(fid,1,'int');
                rlam=fread(fid,1,'float64');
                eol=fread(fid,1,'int');
            elseif flowtype==-2
                eol=fread(fid,1,'int');
                rlam=fread(fid,1,'float64');
                spanv=fread(fid,1,'float64');
                eol=fread(fid,1,'int');
            elseif flowtype>=4
                eol=fread(fid,1,'int');
                boundparam=fread(fid,4,'float64');
                eol=fread(fid,1,'int');
                rlam=boundparam(3);
            end
            
            NNx=storl(1)/2;
            NNy=storl(2);
            NNz=storl(3);
            Re=Re*dstar;
            Lx=length(1)/dstar;
            Lz=length(2)/dstar;
            Ly=2/dstar;     % in bls, dstar is defined as 2/boxsize
            t=length(3)/dstar;
            
                 
            
            realpos=kron(ones(1,NNx),[1 0]);
            vec=0.*realpos';
            rlu=zeros(NNx,NNz,NNy);ilu=0.*rlu;
            rlv=zeros(NNx,NNz,NNy);ilv=0.*rlu;
            rlw=zeros(NNx,NNz,NNy);ilw=0.*rlu;
            
            for indz=1:NNz
                for indy=1:NNy
                    fread(fid,1,'int');
                    vec=fread(fid,NNx*2,'float64');
                    rlu(:,indz,indy)=vec(~~realpos);
                    ilu(:,indz,indy)=vec(~realpos);
                    fread(fid,1,'int');
                end
            end
            
            vec=0.*vec;
            
            for indz=1:NNz
                for indy=1:NNy
                    fread(fid,1,'int');
                    vec=fread(fid,NNx*2,'float64');
                    rlv(:,indz,indy)=vec(~~realpos);
                    ilv(:,indz,indy)=vec(~realpos);
                    fread(fid,1,'int');
                end
            end
            
            vec=0.*vec;
            if NNz~=1
                for indz=1:NNz
                    for indy=1:NNy
                        fread(fid,1,'int');
                        vec=fread(fid,NNx*2,'float64');
                        rlw(:,indz,indy)=vec(~~realpos);
                        ilw(:,indz,indy)=vec(~realpos);
                        
                        fread(fid,1,'int');
                    end
                end
            else
                rlw=0.*rlv;
                ilw=0.*ilv;
            end
            fclose(fid);
            
            scale=1/dstar;
            padx=0;
            padz=0;
            NxF=2*NNx;
            NzF=NNz;
            xF=Lx/NxF*(-NxF/2:1:NxF/2)';
            zF=Lz/NzF*(-NzF/2:1:NzF/2)';
            yF=scale*(1+cos(pi*(0:1/(NNy-1):1)))';
            xF=-xF(1)+xF;
            xF=xF(1:end-1);
            zF=zF(1:end-1);
            
            u=reshape(complex(rlu,ilu),NNx,NNz,NNy);
            v=reshape(complex(rlv,ilv),NNx,NNz,NNy);
            w=reshape(complex(rlw,ilw),NNx,NNz,NNy);
            
            %%%% set to zero the zerozero mode for u component
            %u(1,1,:)=u(1,1,:)*0;
            
            %%%% concatenate the components on the 3rd dimension
            vel=nekmesh.ccat(3,u,v,w);
            
            %%%% remove odd ball
            if (NNz > 1)
                vel=nekmesh.ccat(2,vel(:,1:NNz/2,:),vel(:,NNz/2+2:end,:));
            end
            %%%% this one performs fftshift
            
            velphys=nekmesh.fou2phys(vel);
            storl=size(velphys);
            % $$$ Nx=storl(1);
            % $$$ Nz=storl(2);
            % $$$ Ny=storl(3)/3;
            
            boxsize=yF(1);
            % ymesh=kron(ones(NxF,1),yF');
            % xmesh=kron(xF,ones(1,NNy));
            % zmesh=kron(zF,ones(1,NNy));
            if NNz==1
                [xmesh,ymesh]=meshgrid(xF(:),yF(end:-1:1));
                zmesh=[];
                velu=squeeze(velphys(:,1,1:NNy));
                velv=squeeze(velphys(:,1,1+NNy:2*NNy));
                velw=[];
                velu=velu(:,end:-1:1)';
                velv=velv(:,end:-1:1)';
            else
                [xmesh,ymesh,zmesh]=meshgrid(xF,yF(end:-1:1),zF);
                velu=squeeze(velphys(:,[NNz/2+1:NNz 1:NNz/2],1:NNy));
                velv=squeeze(velphys(:,[NNz/2+1:NNz 1:NNz/2],1+NNy:2*NNy));
                velw=squeeze(velphys(:,[NNz/2+1:NNz 1:NNz/2],1+2*NNy:3*NNy));
                velu=velu(:,:,end:-1:1);
                velv=velv(:,:,end:-1:1);
                velw=velw(:,:,end:-1:1);
                velu=permute(velu,[3 1 2]);
                velv=permute(velv,[3 1 2]);
                velw=permute(velw,[3 1 2]);
                
            end
            %             velu=squeeze(velphys(:,[NNz/2+1:NNz 1:NNz/2],1:NNy));
            %             velv=squeeze(velphys(:,[NNz/2+1:NNz 1:NNz/2],1+NNy:2*NNy));
            %             velw=squeeze(velphys(:,[NNz/2+1:NNz 1:NNz/2],1+2*NNy:3*NNy));
            clear wy
            wy=zeros(size(yF));
            wy(1)  =  -0.5 * (yF(2)   - yF(1));
            wy(numel(yF))=  -0.5 * (yF(end) - yF(end-1));
            for ii = 2:numel(yF)-1
                wy(ii) =  0.5 * (yF(ii-1) - yF(ii+1));
            end
            ww=0.*xmesh;
            if NNz~=1
                dx=xF(2)-xF(1);
                dz=zF(2)-zF(1);
                for ii=1:NNx*2
                    for jj=1:NNz
                        ww(:,ii,jj)=wy*(xF(end)-xF(1)+dx)/numel(xF)*(zF(end)-zF(1)+dz)/numel(zF);
                        wwy(:,ii,jj)=wy;
                    end
                end
            else
                dx=xF(2)-xF(1);
                for ii=1:NNx*2
                    ww(:,ii)=wy*(xF(end)-xF(1)+dx)/numel(xF);
                end
                
            end
            
        end
        function [vel,xF,yF,zF,Lx,Ly,Lz,t,Re,flowtype,dstar,pou,rlam,boundparam]=readdns(filename)
            %
            % Read a velocity field in Fourier space as defined in Simson
            % in variables u,v,w.
            %
            % This version uses the vector form of fread
            %
            % The velocity components u, v, w are concatenated
            % on the third dimension of the output vel
            %
            % The odd ball is removed
            %
            rlam=0.0;
            spanv=0.0;
            
            %
            % Open file
            %
            fid=fopen(filename,'r','ieee-le.l64');
            eol=fread(fid,1,'int');
            
            %
            % If record is not 44 the file is either corrput or on big endian format
            %
            if eol ~= 44
                fclose(fid);
                disp(' ')
                disp(['Reading ' filename ' on big endian format'])
                fid=fopen(filename,'r','ieee-be.l64');
                eol=fread(fid,1,'int');
            else
                disp(' ')
                disp(['Reading ' filename ' on little endian format'])
            end
            Re=fread(fid,1,'float64');
            pou=fread(fid,1,'int');
            length=fread(fid,4,'float64');
            eol=fread(fid,2,'int');
            storl=fread(fid,4,'int');
            eol=fread(fid,2,'int');
            flowtype=fread(fid,1,'int');
            dstar=fread(fid,1,'float64');
            eol=fread(fid,1,'int');
            if flowtype==-1
                eol=fread(fid,1,'int');
                rlam=fread(fid,1,'float64');
                eol=fread(fid,1,'int');
            elseif flowtype==-2
                eol=fread(fid,1,'int');
                rlam=fread(fid,1,'float64');
                spanv=fread(fid,1,'float64');
                eol=fread(fid,1,'int');
                %elseif flowtype==4 || flowtype==5 || flowtype==6
                %  eol=fread(fid,1,'int')
                %  boundparam=fread(fid,2,'float64');
                %  eol=fread(fid,1,'int')
                %  bstart=boundparam(1)
                %  blength=boundparam(2)
            elseif flowtype>=4
                eol=fread(fid,1,'int');
                boundparam=fread(fid,4,'float64');
                eol=fread(fid,1,'int');
                bstart=boundparam(1);
                blength=boundparam(2);
                rlam=boundparam(3);
                spanv=boundparam(4);
            end
            
            %
            % Define parameters based on retreived data
            %
            NNx=storl(1)/2;
            NNy=storl(2);
            NNz=storl(3);
            Re=Re*dstar;
            Lx=length(1)/dstar;
            Lz=length(2)/dstar;
            Ly=2/dstar;           % In bls, dstar is defined as 2/boxsize
            t=length(3)/dstar;
            
            realpos=kron(ones(1,NNx),[1 0]);
            
            disp(' - Reading u');
            for indz=1:NNz
                for indy=1:NNy
                    fread(fid,1,'int');
                    vec=fread(fid,NNx*2,'float64');
                    rlu(:,indz,indy)=vec(~~realpos);
                    ilu(:,indz,indy)=vec(~realpos);
                    fread(fid,1,'int');
                end
            end
            
            disp(' - Reading v');
            for indz=1:NNz
                for indy=1:NNy
                    fread(fid,1,'int');
                    vec=fread(fid,NNx*2,'float64');
                    rlv(:,indz,indy)=vec(~~realpos);
                    ilv(:,indz,indy)=vec(~realpos);
                    fread(fid,1,'int');
                end
            end
            
            disp(' - Reading w');
            for indz=1:NNz
                for indy=1:NNy
                    fread(fid,1,'int');
                    vec=fread(fid,NNx*2,'float64');
                    rlw(:,indz,indy)=vec(~~realpos);
                    ilw(:,indz,indy)=vec(~realpos);
                    fread(fid,1,'int');
                end
            end
            fclose(fid);
            
            scale=1/dstar;
            padx=0;
            padz=0;
            NxF=2*NNx;
            NzF=NNz;
            xF=Lx/NxF*(-NxF/2:1:NxF/2)';
            zF=Lz/NzF*(-NzF/2:1:NzF/2)';
            yF=scale*(1+cos(pi*(0:1/(NNy-1):1)))';
            xF=-xF(1)+xF;
            
            %
            % Shift velocity field in the streamwise direction in order
            % to move the fringe to the end of the domain for spatial
            % flows
            %
            kxvec=linspace(0,2*pi/Lx*(NNx-1),NNx);
            kzvec=linspace(0,2*pi/Lz*(NNz/2-1),NNz/2);
            kzvec=[kzvec -fliplr(kzvec(2:end))];
            
            xs = Lx/2.;
            zs = 0.;
            
            for i=1:NNx
                argx = -xs*kxvec(i);
                cx(i) = cos(argx);
                sx(i) = sin(argx);
            end
            for k=1:NNz-1
                argz = -zs*kzvec(k);
                for i=1:NNx
                    ca(i)=cx(i)*cos(argz)-sx(i)*sin(argz);
                    sa(i)=cx(i)*sin(argz)+sx(i)*cos(argz);
                end
                for j=1:NNy
                    for i=1:NNx
                        hr=rlu(i,k,j)*ca(i)-ilu(i,k,j)*sa(i);
                        ilu(i,k,j)=ilu(i,k,j)*ca(i)+rlu(i,k,j)*sa(i);
                        rlu(i,k,j)=hr;
                        hr=rlv(i,k,j)*ca(i)-ilv(i,k,j)*sa(i);
                        ilv(i,k,j)=ilv(i,k,j)*ca(i)+rlv(i,k,j)*sa(i);
                        rlv(i,k,j)=hr;
                        hr=rlw(i,k,j)*ca(i)-ilw(i,k,j)*sa(i);
                        ilw(i,k,j)=ilw(i,k,j)*ca(i)+rlw(i,k,j)*sa(i);
                        rlw(i,k,j)=hr;
                    end
                end
            end
            %
            %u=fftshift(u,1);
            %v=fftshift(v,1);
            %w=fftshift(w,1);
            
            u=reshape(complex(rlu,ilu),NNx,NNz,NNy);
            v=reshape(complex(rlv,ilv),NNx,NNz,NNy);
            w=reshape(complex(rlw,ilw),NNx,NNz,NNy);
            
            
            % Set to zero the zerozero mode for u component
            %u(1,1,:)=u(1,1,:)*0.0;
            
            %
            % Concatenate the components on the 3rd dimension
            %
            vel=nekmesh.ccat(3,u,v,w);
            
            %
            % Remove odd ball
            %
            if NNz>1;
                vel=nekmesh.ccat(2,vel(:,1:NNz/2,:),vel(:,NNz/2+2:end,:));
            end
        end
        function [snap,xmesh,ymesh,zmesh,ww,xm,ym,zm,indt]=readblamult(casename,is,ie,istep,cut,varargin)
            if(nargin>5)
                pr=varargin{1};
                pdd=pwd;
                eval(['cd ' pr]);
            end
            
            
            param=nekmesh.readparam([casename sprintf('%04i',is) '.u']);
            
            if(~isnan(cut(1)));indxs=find(param.xf(:)>=cut(1),1,'first');else indxs=1;end
            if(~isnan(cut(2)));indxe=find(param.xf(:)>=cut(2),1,'first');else indxe=param.nx;end
            if(~isnan(cut(3)));indys=find(param.yf(:)>=cut(3),1,'first');else indys=1;end
            if(~isnan(cut(4)));indye=find(param.yf(:)>=cut(4),1,'first');else indye=param.ny;end
            if mod(indye-indys+1,2)==0;indye=indye+1;end
            if(~isnan(cut(5)));indzs=find(param.zf(:)>=cut(5),1,'first');else indzs=1;end
            if(~isnan(cut(6)));indze=find(param.zf(:)>=cut(6),1,'first');else indze=param.nz;end
            cut=[indxs indxe indys indye indzs indze];
            nsnap=0;
            for ii=is:istep:ie
                nsnap=nsnap+1;
            end
            
            
            fid=fopen('dmd.i','w','n');
            fprintf(fid,[casename                  '     #p1 =casename\n']);
            fprintf(fid,[num2str(is)               '     #p2 =is (first snapshot number) \n']);
            fprintf(fid,[num2str(nsnap)            '     #P3 =nsnap (last snapshot number) \n']);
            fprintf(fid,[num2str(istep)            '     #p4 =istep \n']);
            fprintf(fid,['.'                       '     #p5 =path\n']);
            fprintf(fid,[num2str(cut(1))           '	 #p6 =first number of element in x\n']);
            fprintf(fid,[num2str(cut(2))           '	 #p7 =last number of element  in x\n' ]);
            fprintf(fid,[num2str(cut(3))           '	 #p8 =first number of element in y\n']);
            fprintf(fid,[num2str(cut(4))           '	 #p9 =last number of element  in y\n' ]);
            fprintf(fid,[num2str(cut(5))           '	 #p10=first number of element in z\n']);
            fprintf(fid,[num2str(cut(6))           '	 #p11=last number of element  in z\n' ]);
            fclose(fid);
            
            indt=[indxs indxe indys indye indzs indze];
            
            
            system('readsimfile');
            
            
            fid=fopen('sensor.dat','r');
            fseek(fid,4,0);
            tmp=fread(fid,7,'int');
            nxc=tmp(1); nyc=tmp(2); nzc=tmp(3);nx =tmp(4); ny =tmp(5); nz =tmp(6);nsnaps=tmp(7);
            fseek(fid,8,0);
            dd=fread(fid,nxc*nyc*nzc*3*nsnaps,'float64');
            fclose(fid);
            snap=reshape(dd,nxc,nyc,nzc,3,nsnaps);
            
            xf=param.xf(indxs:indxe);
            yf=param.yf(indys:indye);
            zf=param.zf(indzs:indze);
            
            [xmesh,ymesh,zmesh]=meshgrid(xf,yf,zf);
            snap=permute(snap,[2 1 3 4 5]);
            clear wy
            wy=zeros(size(yf));
            wy(1)  =  0.5 * (yf(2)   - yf(1));
            wy(length(yf)) =  0.5 * (yf(end) - yf(end-1));
            for ii = 2:length(yf)-1
                wy(ii) =  -0.5 * (yf(ii-1) - yf(ii+1));
            end
            ww=0.*xmesh;
            for ii=1:nxc
                for jj=1:nzc
                    ww(:,ii,jj)=wy*(xf(end)-xf(1))/length(xf)*(zf(end)-zf(1))/length(zf);
                end
            end
            [xm,ym,zm]=meshgrid(param.xf(1:end-1),param.yf,param.zf(1:end-1));
            if(nargin>5)
                eval(['cd ' pdd]);
            end
            
        end
        function [snap,xmesh,ymesh,zmesh,ww,x,y,z,indt]=readsimmultmat(is,ie,istep,cut,varargin)
              %[snap,xmesh,ymesh,zmesh,ww,x,y,z,indt]=readsimmultmat(is,ie,istep,cut,varargin)
              %varargin{1} path otherwise '.'
              switch nargin
                  case 4
                      pp='.';
                  otherwise
                      pp=varargin{1};
              end
              
            t=is:istep:ie;
            [~,~,~,x,y,z,wp]= nekmesh.readsim([pp '/t' sprintf('%04d',is) '.u']);
            xp=squeeze(x(1,:,1));
            yp=squeeze(y(:,1,1));
            zp=squeeze(z(1,1,:));
            
            if ~isnan(cut(1)); indt(1)=find(xp(:)>=cut(1),1,'first'); else indt(1)=1;end
            if ~isnan(cut(2)); indt(2)=find(xp(:)>=cut(2),1,'first'); else indt(2)=length(xp);end
            if ~isnan(cut(3)); indt(3)=find(yp(:)>=cut(3),1,'first'); else indt(3)=1;end
            if ~isnan(cut(4)); indt(4)=find(yp(:)>=cut(4),1,'first'); else indt(4)=length(yp);end
            if ~isnan(cut(5)); indt(5)=find(zp(:)>=cut(5),1,'first'); else indt(5)=1;end
            if ~isnan(cut(6)); indt(6)=find(zp(:)>=cut(6),1,'first'); else indt(6)=length(zp);end
            
            xmesh=x(indt(3):indt(4),indt(1):indt(2),indt(5):indt(6));
            ymesh=y(indt(3):indt(4),indt(1):indt(2),indt(5):indt(6));
            zmesh=z(indt(3):indt(4),indt(1):indt(2),indt(5):indt(6));
            
            n=1;
            snap=zeros(size(xmesh,1),size(xmesh,2),size(xmesh,3),3,length(t));
            for ii=is:istep:ie
                [u,v,w]= ...
                    nekmesh.readsim([pp '/t' sprintf('%04d',ii) '.u']);
                disp(['file number ' num2str(ii) ' is read']);
                snap(:,:,:,1,n)=u(indt(3):indt(4),indt(1):indt(2),indt(5):indt(6));
                snap(:,:,:,2,n)=v(indt(3):indt(4),indt(1):indt(2),indt(5):indt(6));
                snap(:,:,:,3,n)=w(indt(3):indt(4),indt(1):indt(2),indt(5):indt(6));
                n=n+1;
                u=[];
                v=[];
                w=[];
            end
            clear wy
            
            xp=xp(indt(1):indt(2));
            yp=yp(indt(3):indt(4));
            zp=zp(indt(5):indt(6));
            
            
            wy=zeros(size(yp));
            wy(1)  = 0.5*(yp(2)-yp(1));
            wy(length(yp)) =  0.5 * (yp(end) - yp(end-1));
            for ii = 2:length(yp)-1
                wy(ii) =  0.5 * (yp(ii+1) - yp(ii-1));
            end
            ww=0.*xmesh;
            for ii=1:length(xp)
                for jj=1:length(zp)
                    ww(:,ii,jj)=wy*(xp(end)-xp(1))/length(xp)*(zp(end)-zp(1))/length(zp);
                end
            end
            if sum(isnan(cut))==6
                ww=wp;
            end
        end
        %******************************************************************
        function [phys,NNx,NNy,NNz]=fou2phys(fou,varargin)
            % ***********************************************************************
            %
            % $HeadURL: https://www2.mech.kth.se/svn/simson/trunk/matlab/fou2phys.m $
            % $LastChangedDate: 2006-11-16 21:05:30 +0100 (Thu, 16 Nov 2006) $
            % $LastChangedBy: mattias@MECH.KTH.SE $
            % $LastChangedRevision: 336 $
            %
            % ***********************************************************************
            % Fourier transform
            % fou Two-dimensional plane in Fourier space
            % ex  Number of extra zeros in x-direction (0 normally) to
            %     improve interpolation to physical space
            % ez  Number of extra zeros in z-direction (0 normally)
            %     improve interpolation to physical space
            
            s=size(fou);
            Nx=2*s(1);
            Nz=s(2)+1;
            
            if length(s)>=3; Ny=s(3); else Ny=1;end;
            if length(s)>=4; Ns=s(4); else Ns=1;end;
            if length(s)>=5; Nr=s(5); else Nr=1;end;
            
            
            if size(varargin) == 0
                ex=0;
                ez=0;
                ncomp=3;
            elseif size(varargin) == 1
                ex=varargin{1};
                ez=0;
                ncomp=3;
            elseif size(varargin) == 2
                ex=varargin{1};
                ez=varargin{2};
                ncomp=3;
            elseif size(varargin) == 3
                ex=varargin{1};
                ez=varargin{2};
                ncomp = varargin{3};
            end
            NNx=Nx+ex;
            NNy=Ny/ncomp;
            NNz=Nz+ez;
            
            phys=zeros(Nx+ex,Nz+ez,Ny,Ns,Nr);
            
            for inds=1:Ns
                for indy=1:Ny
                    for indr=1:Nr
                        % The z odd ball
                        tp=fou(:,:,indy,inds,indr);
                        tp=[tp(:,1:Nz/2) , zeros(Nx/2,1+ez) , tp(:,Nz/2+1:Nz-1)];
                        
                        % The conjugate and x odd ball
                        X=[ tp ;
                            zeros(1+ex,Nz+ez) ;
                            conj(flipud(tp(2:Nx/2,1))) , conj(flipud(fliplr(tp(2:Nx/2,2:Nz+ez))))];
                        
                        % To physical space
                        x=ifft2(X)*(Nx+ex)*(Nz+ez);
                        x=fftshift(x);
                        
                        phys(:,:,indy,inds,indr)=real(x);
                        
                    end
                end
            end
        end
        function [fou]=phys2fou(phys)
            % Fourier transform
            
            s=size(phys);
            Nx=s(1);
            Nz=s(2);
            if length(s)==3; Ny=s(3); else Ny=1;end;
            
            fou=zeros(Nx/2,Nz-1,Ny);
            
            for indy=1:Ny
                % To Fourier space
                x=phys(:,:,indy);
                
                %x=fftshift(x);
                X=fft2(x)/(Nx*Nz);
                
                % Remove conjugate part and odd balls
                fou(:,:,indy)=[X(1:Nx/2,1:Nz/2) ,  X(1:Nx/2,Nz/2+2:Nz)];
            end
        end
        %******************************************************************
        function res=ccat(dim,varargin)
            res=[];
            for ind=1:nargin-1
                res=cat(dim,res,varargin{ind});
            end
        end
        %***********************  stat ************************************
        function [ur,xp,yp,Data,up,uu]=readstat(filename,varargin)
            %[Data,up,xp,yp]=readstat(filename,scalar)
            %filenmae is the name of the file
            %varargin{1} number of scalars; for the simple case please
            %specify zero or nothing
            %Data; contains all the statistics
            %up; is just the streamwise rms
            %xp,yp; mesh data
            switch nargin
                case 1
                    scalar=0;
                otherwise
                   scalar=varargin{1}; 
            end
            fid=fopen(filename,'r','ieee-be.l64');
            
            eol=fread(fid,1,'int');
            Re=fread(fid,1,'float64');
            pou=fread(fid,1,'int');
            length=fread(fid,4+scalar*2,'float64');
            eol=fread(fid,2,'int');
            AA=fread(fid,1,'char');
            dd=fread(fid,4,'float64');
            eol=fread(fid,2,'int');
            storl=fread(fid,4,'int');
            eol=fread(fid,2,'int');
            flowtype=fread(fid,1,'int');
            dstar=fread(fid,1,'float64');
            eol=fread(fid,1,'int');
            if flowtype==-1
                eol=fread(fid,1,'int');
                rlam=fread(fid,1,'float64');
                eol=fread(fid,1,'int');
            elseif flowtype==-2
                eol=fread(fid,1,'int');
                rlam=fread(fid,1,'float64');
                spanv=fread(fid,1,'float64');
                eol=fread(fid,1,'int');
            elseif flowtype>=4
                eol=fread(fid,1,'int');
                boundparam=fread(fid,4,'float64');
                eol=fread(fid,1,'int');
                bstart=boundparam(1);
                blength=boundparam(2);
                rlam=boundparam(3);
                spanv=boundparam(4);
            end
            
            NNx=storl(1)/2;
            NNy=storl(2);
            NNz=storl(3);
            Re=Re*dstar;
            Lx=length(1)/dstar;
            Lz=length(2)/dstar;
            Ly=2/dstar;           % In bls, dstar is defined as 2/boxsize
            t=length(3)/dstar;
            
            eol=fread(fid,1,'int');
            sumw=fread(fid,1,'float64');
            para=fread(fid,3,'int');
            nxys=para(1);
            nxysth=para(2);
            scalar=para(3);
            eol=fread(fid,1,'int');
            
            A=zeros(NNx*2*NNy,nxys+nxysth);
            for xysi=1:nxys+nxysth*scalar
                eol=fread(fid,1,'int');
                A(:,xysi)=fread(fid,2*NNx*NNy,'float64');
                eol=fread(fid,1,'int');
            end
            
            xys=reshape(A,2*NNx,NNy,nxys+nxysth);
            fclose(fid)
            
            for ii=1:nxys
                nn=1;
                for j=1:NNy
                    for k=1:2*NNx
                        totxys(k,j,ii)=A(nn,ii);
                        nn=nn+1;
                    end
                end
            end
            
            yF=(1+cos([0:(NNy-1)]*pi/(NNy-1)))*Ly/2;
            xF=(0:2*NNx-1)/(2*NNx-1)*Lx;
            
            Data.y=[yF(:)];
            % Data.x=[xF(NNx+1:2*NNx) xF(1:NNx)];
            Data.x=xF;
            Data.u=[xys(NNx+1:2*NNx,:,1); xys(1:NNx,:,1)];
            Data.v=[xys(NNx+1:2*NNx,:,2); xys(1:NNx,:,2)];
            Data.w=[xys(NNx+1:2*NNx,:,3); xys(1:NNx,:,3)];
            Data.uu=[xys(NNx+1:2*NNx,:,4); xys(1:NNx,:,4)];
            Data.vv=[xys(NNx+1:2*NNx,:,5); xys(1:NNx,:,5)];
            Data.ww=[xys(NNx+1:2*NNx,:,6); xys(1:NNx,:,6)];
            Data.uv=[xys(NNx+1:2*NNx,:,13); xys(1:NNx,:,13)];
            Data.uw=[xys(NNx+1:2*NNx,:,14); xys(1:NNx,:,14)];
            Data.vw=[xys(NNx+1:2*NNx,:,15); xys(1:NNx,:,15)];
            
            if scalar>0
                Data.T =[xys(NNx+1:2*NNx,1:(NNy+1)/2,nxys+1);
                    xys(1:NNx,1:(NNy+1)/2,nxys+1)];
                Data.uT=[xys(NNx+1:2*NNx,1:(NNy+1)/2,nxys+3);
                    xys(1:NNx,1:(NNy+1)/2,nxys+3)];
                Data.vT=[xys(NNx+1:2*NNx,1:(NNy+1)/2,nxys+4);
                    xys(1:NNx,1:(NNy+1)/2,nxys+4)];
                Data.wT=[xys(NNx+1:2*NNx,1:(NNy+1)/2,nxys+5);
                    xys(1:NNx,1:(NNy+1)/2,nxys+5)];
                Data.TT=[xys(NNx+1:2*NNx,1:(NNy+1)/2,nxys+2);
                    xys(1:NNx,1:(NNy+1)/2,nxys+2)];
            end
            
            Data.urms=(Data.uu-Data.u.^2);
            Data.vrms=(Data.vv-Data.v.^2);
            Data.wrms=(Data.ww-Data.w.^2);
            ind=(Data.urms(:)<=0 & Data.urms(:) >= -1e-13);
            Data.urms(ind)=0;
            Data.urms=sqrt(Data.urms);
            
            ind=(Data.vrms(:)<=0 & Data.vrms(:) >= -1e-13);
            Data.vrms(ind)=0;
            Data.vrms=sqrt(Data.vrms);    
            
            
            Data.kin=(Data.wrms.^2+Data.vrms.^2+Data.wrms.^2)./2;
            ur=permute(Data.urms(:,end:-1:1),[2,1]);
            up=permute(Data.u(:,end:-1:1),[2,1]);
            uu=permute(Data.uu(:,end:-1:1),[2,1]);
            
            [xp,yp]=meshgrid(Data.x,Data.y(end:-1:1));
            
            
            % %Little shift for the sake of nice representations
            % xF(NNx+3:NNx*2+1)=xF(NNx+2:2*NNx);
            % xys(NNx+3:NNx*2+1,:,:)=xys(NNx+2:2*NNx,:,:);
            % xF(NNx+2)=NaN;xys(NNx+2,:,:)=NaN;
            
        end
        function [ur,vr,wr,up,vp,wp,uu]=readstatu(filename,nx,ny,nz)
            fid=fopen(filename,'r','ieee-be.l64');
            A=zeros(nx*ny*nz,6);
            for ii=1:6
                eol=fread(fid,1,'int');
                xyzs(:,ii)=fread(fid,nx*ny*nz,'float64');
                eol=fread(fid,1,'int');
            end
            xyzs=reshape(xyzs,nx,nz,ny,6);
            fclose(fid);
            up=xyzs([end/2+1:end 1:end/2],:,end:-1:1,1);
            vp=xyzs([end/2+1:end 1:end/2],:,end:-1:1,2);
            wp=xyzs([end/2+1:end 1:end/2],:,end:-1:1,3);
            uu=xyzs([end/2+1:end 1:end/2],:,end:-1:1,4);
            vv=xyzs([end/2+1:end 1:end/2],:,end:-1:1,5);
            ww=xyzs([end/2+1:end 1:end/2],:,end:-1:1,6);
            up=squeeze(permute(up,[3 1 2]));
            vp=squeeze(permute(vp,[3 1 2]));
            wp=squeeze(permute(wp,[3 1 2]));
            uu=squeeze(permute(uu,[3 1 2]));
            vv=squeeze(permute(vv,[3 1 2]));
            ww=squeeze(permute(ww,[3 1 2]));
            ur=(uu-up.^2);
            ind=(ur(:)<=0 & ur(:) >= -1e-13);ur(ind)=0;
            vr=(vv-vp.^2);
            ind=(vr(:)<=0 & vr(:) >= -1e-13);vr(ind)=0;
            wr=(ww-wp.^2);
            ind=(wr(:)<=0 & wr(:) >= -1e-13);wr(ind)=0;
            
            ur=sqrt(ur);
            vr=sqrt(vr);
            wr=sqrt(wr);
            
            
            
        end
        %******************************************************************
        function pcolor3(xmesh,ymesh,zmesh,snap,varargin)
            %pcolor3(xmesh,ymesh,zmesh,snap,ii,varargin)
            %input: xmesh,ymesh,zmesh
            %velocities: velocity
            %varagin{1} : ii - snapshot number
            %varargin(2): magitude in procent
            %varargin(3):aspectration(0 or 1)
            %varargin{4}: constatnt color bar
            mt=size(snap);
            switch numel(mt)
                case 3
                    velu=snap;
                    clear snap
                case 4
                    ii=varargin {1};
                    velu=squeeze(snap(:,:,:,ii));
                    clear snap
                case 5
                    ii=varargin {1};
                    velu=squeeze(snap(:,:,:,1,ii));
                    clear snap
            end
            
            mm=max(max(max(max(velu))))/100;
            switch nargin
                case 4
                    mag=50*mm;
                    dasp=0;
                    sc=0;
                 case 5
                     mag=50*mm;
                     dasp=0;
                     sc=0;
                 case 6
                     mag=varargin{2}*mm;
                     dasp=0;
                     sc=0;
                 case 7
                    mag=varargin{2}*mm;
                    dasp=varargin{3};
                    sc=0;
                case 8
                    mag=varargin{2}*mm;
                    dasp=varargin{3};
                    sc=varargin{4};
            end
             
             p = patch(isosurface(xmesh,ymesh,zmesh,velu,mag));
             isonormals(xmesh,ymesh,zmesh,velu,p)
             set(p,'FaceColor','red','EdgeColor','none');
             if dasp
                 daspect([1,1,1]);
             end
             view(3);
             axis tight
             camlight
             lighting gouraud
             hold on 
             p = patch(isosurface(xmesh,ymesh,zmesh,velu,-mag));
             isonormals(xmesh,ymesh,zmesh,velu,p)
             set(p,'FaceColor','blue','EdgeColor','none');
             if dasp
                 daspect([1,1,1]);
             end
             view(3);
             axis tight
             camlight
             lighting gouraud
             xlabel('x');ylabel('y');zlabel('z');
             title(['iso-surface = ' num2str(mag)]);
             drawnow
             hold off
             if sc;clf;end
             
             
        end
        function pcolor2(xmesh,ymesh,zmesh,snap,plane,val,varargin)
            % pcolor2(xmesh,ymesh,zmesh,snap,plane,val,varargin)
            %xmesh,ymesh,zmesh: grid points
            %snap : all the velocity fields usually u(:,:,:,uvw,nsnap) it
            %can also be the single velocity component
            %plane : it can be 'xy' xz' 'yz'
            %val:the intersection location of the plane
            %varargin{1}: the snapshot number ; default = 1
            %varargin{2}: constant color axis with this value
            %varargin{3}: axis equal or not
            eq=0;
            m=[];
            switch nargin
                case 7
                    ii=varargin{1};
                case 8
                    ii=varargin{1};
                    m=varargin{2};
                case 9
                    ii=varargin{1};
                    m=varargin{2};
                    eq=varargin{3};
                otherwise
                    ii=1;disp('plot is shown for the first snapshot')
            end
            mt=size(snap);
            switch numel(mt)
                case 3
                    velu=snap;
                    clear snap
                case 4
                    velu=squeeze(snap(:,:,:,ii));
                    clear snap
                case 5
                    velu=squeeze(snap(:,:,:,1,ii));
                    clear snap
            end
            if strcmp(plane,'xy')
                ind=find(squeeze(zmesh(1,1,:))<=val,1,'last');
                pcolor(squeeze(xmesh(:,:,ind)),squeeze(ymesh(:,:,ind)),squeeze(velu(:,:,ind)))
                shading interp;colorbar
                xlabel('x');ylabel('y');
                if ~isempty(m);caxis(m);end
                if eq;axis equal ;end
                drawnow
            elseif strcmp(plane,'xz')
                magy=abs(squeeze(ymesh(:,1,1))-val);
                [~,ind]=min(magy);
                pcolor(squeeze(xmesh(ind,:,:)),squeeze(zmesh(ind,:,:)),squeeze(velu(ind,:,:)))
                shading interp;colorbar
                xlabel('x','fontsize',12,'FontName','Times','FontAngle','italic');
                ylabel('z','fontsize',12,'FontName','Times','FontAngle','italic');
                set(findobj('type','axes'),'fontsize',12,'FontName','Times','FontAngle','italic')
                if ~isempty(m);caxis([-m m]);end
                if eq;axis equal ;end
                drawnow
                disp(['the plot is at the location ' num2str(ymesh(ind,1,1))])
            elseif strcmp(plane,'yz')
                ind=find(squeeze(xmesh(1,:,1))<=val,1,'last');
                pcolor(squeeze(ymesh(:,ind,:)),squeeze(zmesh(:,ind,:)),squeeze(velu(:,ind,:)))
                shading interp;colorbar
                xlabel('y');ylabel('z');drawnow
                if ~isempty(m);caxis(m);end
                if eq;axis equal ;end
                drawnow
            else
                disp('the plane is not found')
            end
            
        end
        function contour2(xmesh,ymesh,zmesh,snap,plane,val,varargin)
            % pcolor2(xmesh,ymesh,zmesh,snap,plane,val,varargin)
            %xmesh,ymesh,zmesh: grid points
            %snap : all the velocity fields usually u(:,:,:,uvw,nsnap) it
            %can also be the single velocity component
            %plane : it can be 'xy' xz' 'yz'
            %val:the intersection location of the plane
            %varargin{1}: the snapshot number ; default = 1
            %varargin{2}: constant color axis with this value
            %varargin{3}: axis equal or not
            eq=0;
            m=[];
            switch nargin
                case 7
                    ii=varargin{1};
                case 8
                    ii=varargin{1};
                    m=varargin{2};
                case 9
                    ii=varargin{1};
                    m=varargin{2};
                    eq=varargin{3};
                otherwise
                    ii=1;disp('plot is shown for the first snapshot')
            end
            mt=size(snap);
            switch numel(mt)
                case 3
                    velu=snap;
                    clear snap
                case 4
                    velu=squeeze(snap(:,:,:,ii));
                    clear snap
                case 5
                    velu=squeeze(snap(:,:,:,1,ii));
                    clear snap
            end
            if strcmp(plane,'xy')
                ind=find(squeeze(zmesh(1,1,:))<=val,1,'last');
                pcolor(squeeze(xmesh(:,:,ind)),squeeze(ymesh(:,:,ind)),squeeze(velu(:,:,ind)))
                shading interp;colorbar
                xlabel('x');ylabel('y');
                if ~isempty(m);caxis(m);end
                if eq;axis equal ;end
                drawnow
            elseif strcmp(plane,'xz')
                magy=abs(squeeze(ymesh(:,1,1))-val);
                [~,ind]=min(magy);
                contour(squeeze(xmesh(ind,:,:)),squeeze(zmesh(ind,:,:)),...
                    squeeze(velu(ind,:,:)) )
                
                xlabel('x','fontsize',12,'FontName','Times','FontAngle','italic');
                ylabel('z','fontsize',12,'FontName','Times','FontAngle','italic');
                set(findobj('type','axes'),'fontsize',12,'FontName','Times','FontAngle','italic')
                if ~isempty(m);caxis([-m m]);end
                if eq;axis equal ;end
                drawnow
                disp(['the plot is at the location ' num2str(ymesh(ind,1,1))])
            elseif strcmp(plane,'yz')
                ind=find(squeeze(xmesh(1,:,1))<=val,1,'last');
                pcolor(squeeze(ymesh(:,ind,:)),squeeze(zmesh(:,ind,:)),squeeze(velu(:,ind,:)))
                shading interp;colorbar
                xlabel('y');ylabel('z');drawnow
                if ~isempty(m);caxis(m);end
                if eq;axis equal ;end
                drawnow
            end
            
        end
        function surf2(xmesh,ymesh,zmesh,snap,plane,val,varargin)
            if nargin>6
                ii=varargin{1};
            else
                ii=1;disp('plot is shown for the first snapshot')
            end
            if strcmp(plane,'xy')
                ind=find(squeeze(zmesh(1,1,:))<=val,1,'last');
                surf(squeeze(xmesh(:,:,ind)),squeeze(ymesh(:,:,ind)),squeeze(snap(:,:,ind,1,ii)))
                shading interp;colorbar
                xlabel('x');ylabel('y');drawnow
            elseif strcmp(plane,'xz')
                ind=find(squeeze(ymesh(:,1,1))<=val,1,'last');
                surf(squeeze(xmesh(ind,:,:)),squeeze(zmesh(ind,:,:)),squeeze(snap(ind,:,:,1,ii)))
                shading interp;colorbar
                xlabel('x');ylabel('z');drawnow
            elseif strcmp(plane,'yz')
                ind=find(squeeze(xmesh(1,:,1))<=val,1,'last');
                surf(squeeze(ymesh(:,ind,:)),squeeze(zmesh(:,ind,:)),squeeze(snap(:,ind,:,1,ii)))
                shading interp;colorbar
                xlabel('y');ylabel('z');drawnow
            end
            
        end
        %*****************************************************************
        function writeblai(is,ie,st)
            %writeblai(is,ie,st)
            % inputs:
            %is : the first time you want to outpost
            %ie : the last time you want to outpost
            %st : the stride between the times
             fid=fopen('tmp3.dat','w');
             itt=is:st:ie;
             ns=length(itt);
             
             
             system('head -n47 ~/codes/simson/reza/bla.i > tmp1.dat');
             fprintf(fid,[num2str(ns) '           msave:          number of saved interm. 3-D field\n']);
             for ii=1:ns
                 fprintf(fid,[num2str(itt(ii)) '\n']);
                 fprintf(fid,['t' sprintf('%04i',itt(ii)) '.u\n']);
             end
             system('tail -n4 ~/codes/simson/reza/bla.i > tmp2.dat');
             fclose(fid);
             
             system('cat tmp1.dat tmp3.dat tmp2.dat > bla.i');
             system('rm tmp*.dat')
         end
        %************************************************************
        function param=readparam(filename)
            rlam=0.0;
            spanv=0.0;
            
            %
            % Open file
            %
            fid=fopen(filename,'r','ieee-le.l64');
            eol=fread(fid,1,'int');
            
            %
            % If record is not 44 the file is either corrput or on big endian format
            %
            if eol ~= 44
                fclose(fid);
                disp(' ')
                disp(['Reading ' filename ' on big endian format'])
                fid=fopen(filename,'r','ieee-be.l64');
                eol=fread(fid,1,'int');
            else
                disp(' ')
                disp(['Reading ' filename ' on little endian format'])
            end
            Re=fread(fid,1,'float64');
            pou=fread(fid,1,'int');
            length=fread(fid,4,'float64');
            eol=fread(fid,2,'int');
            storl=fread(fid,4,'int');
            eol=fread(fid,2,'int');
            flowtype=fread(fid,1,'int');
            dstar=fread(fid,1,'float64');
            eol=fread(fid,1,'int');
            if flowtype==-1
                eol=fread(fid,1,'int');
                rlam=fread(fid,1,'float64');
                eol=fread(fid,1,'int');
            elseif flowtype==-2
                eol=fread(fid,1,'int');
                rlam=fread(fid,1,'float64');
                spanv=fread(fid,1,'float64');
                eol=fread(fid,1,'int');
                %elseif flowtype==4 || flowtype==5 || flowtype==6
                %  eol=fread(fid,1,'int')
                %  boundparam=fread(fid,2,'float64');
                %  eol=fread(fid,1,'int')
                %  bstart=boundparam(1)
                %  blength=boundparam(2)
            elseif flowtype>=4
                eol=fread(fid,1,'int');
                boundparam=fread(fid,4,'float64');
                eol=fread(fid,1,'int');
                bstart=boundparam(1);
                blength=boundparam(2);
                rlam=boundparam(3);
                spanv=boundparam(4);
            end
            param.nx=storl(1);
            param.ny=storl(2);
            param.nz=storl(3);
            param.re=Re*dstar;
            param.x=length(1)/dstar;
            param.z=length(2)/dstar;
            param.y=2/dstar;           % In bls, dstar is defined as 2/boxsize
            param.t=length(3)/dstar;
            param.xf=param.x/param.nx*(-param.nx/2:1:param.nx/2)';
            param.zf=param.z/param.nz*(-param.nz/2:1:param.nz/2)';
            param.yf=param.y/2*(1+cos(pi*(0:1/(param.ny-1):1)))';
            param.yf=param.yf(end:-1:1);
            param.xf=-param.xf(1)+param.xf;
        end
        %************************************************************
        function writesim(up,vp,wp,pr,fname)
            %writesim(up,vp,wp,pr,fname)
            %velocities;path,output name
            pdd=pwd;
            eval(['cd ' pr])
            
            [~,nx,~]=size(up);
            up=permute(up,[2 3 1]);up=up(:,:,end:-1:1);up=up([nx/2+1:nx 1:nx/2],:,:);
            vp=permute(vp,[2 3 1]);vp=vp(:,:,end:-1:1);vp=vp([nx/2+1:nx 1:nx/2],:,:);
            wp=permute(wp,[2 3 1]);wp=wp(:,:,end:-1:1);wp=wp([nx/2+1:nx 1:nx/2],:,:);
            %------------------------- write file-----------
            % cd /home/guest/dadfar/codes/simson/bls
            fid=fopen('input.txt','w');
            utot=cat(4,up,vp,wp);
            utot=permute(utot,[4 1 2 3]);
            fprintf(fid,'%20.16f \n',utot);
            fclose(fid);
            system('./getsimson bls');
            system('./bls');
            system(['mv b1.u ' pdd '/' fname]);
            eval(['cd  ' pdd])
        end
        %******************************************************************
        function [modes,xm,ym,zm,xmesh,ymesh,zmesh,S,U,V,t]=podsim(is,ie,istep,cut,casename,nm,t)
            
            [snap,xmesh,ymesh,zmesh,ww,xm,ym,zm,indt]=nekmesh.readbla(casename,is,ie,istep,cut,pp);
            
            n=size(snap);
            snap2=reshape(snap,prod(n(1:end-1)),n(end));
            ww2=repmat(ww(:),3,1);
            for ii=1:n(end)
                snap2(:,ii)=snap2(:,ii).*sqrt(ww2(:));
            end
            [U,S,V]=svd(snap2,'econ');
            for ii=1:n(end)
                U(:,ii)=U(:,ii)./sqrt(ww2(:));
            end
            U=reshape(U,size(snap));
            
            % figure;nekmesh.pcolor2(xmesh,ymesh,zmesh,U,'xz',1,4);
            % nekmesh.pcolor3(xmesh,ymesh,zmesh,U,2,0.01)
            
            semilogy(diag(S),'ro')
            s1=diag(S);
            % disp(sum(s1(1:10).^2)/sum(s1.^2)*100);
            
            
            
            modes=zeros([size(xm) 3 nm]);
            for ii=1:nm
                modes(indt(3):indt(4),indt(1):indt(2),indt(5):indt(6),:,ii)=U(:,:,:,:,ii);
            end
            
        end
        %******************************************************************
        function up=symsim(up,varargin)
            % symetric and anti symetrize an snapshots
            % up=symsim(up,varargin)
            % up is the snapshots default symmetry
            %varargin{1} =1 if symmetry is needed and =2 if antisymmetric
            
            switch nargin
                case 1
                    sym=1;
                otherwise
                    sym=varargin{1};
            end
            if sym
                up(:,:,end:-1:end/2+2)=+up(:,:,2:end/2);
            else
                up(:,:,end:-1:end/2+2)=-up(:,:,2:end/2);
            end
            
        end
        %******************************************************************
        %                     free stream turbulence
        %******************************************************************
        function [x,y,z,u,v,w,t] = incbox_getfort(fpath,fname)
            %[x,y,z,u,v,w,t] = getfort(fpath,fname)
            counter = 1;
            
            load([fpath '/' fname])
            nfield=3;
%             xc=25; yc=8*pi;ra=0.5;
%             theta=[0:0.01:2*pi]; xx=xc+ra*cos(theta); yy=yc+ra*sin(theta);
            
            t=fort(1);
            xl=fort(2);
            yl=fort(3);
            zl=fort(4);
            npx=fort(5);
            npy=fort(6);
            npz=fort(7);
            uvw=reshape(fort(8:7+npx*npy*npz*nfield),npx,npy,npz,nfield);
            vor=reshape(fort(8+npx*npy*npz*nfield:end),npx,npy,npz,nfield);
            
            
            nppx=npx-2;
            x=0:xl*pi/(nppx):xl*pi-xl*pi/(nppx);
            y=0:yl*pi/npy:yl*pi-yl*pi/npy;
            z=0:zl*pi/npz:zl*pi-zl*pi/npz;
            
            
            u = squeeze(uvw(end/2,:,:,1));
            v = squeeze(uvw(end/2,:,:,2));
            w = squeeze(uvw(end/2,:,:,3));
%             urms=1/3*sqrt(rms(u(:))^2+rms(v(:))^2+rms(w(:))^2);
        end
        function [urms,tke,relambda,lambda,t]=incbox_getstat(fpath,varargin)
            if nargin == 1
                varargin{1} = 15 ;
            end
                
            dd = pwd ;
            eval(['cd  ' fpath])
            load fort.7
            t = fort(:,1) ;
            tke = fort(:,2);
            relambda = fort(:,3);
            lambda = fort(:,4) ;
            urms = sqrt(2/3*mean(tke(varargin{1}:end)));
            eval(['cd ' dd])
        end
        function [xx,yy,zz,uu,vv,ww,t] = getnekbox(fpath,fname,x_station)
            vel =load([fpath '/' fname]);
            t = vel.t
            ind = find(vel.x(1,1,:)>=x_station,1,'first');
            ind
            xx  = squeeze(vel.x(:,:,ind));
            yy  = squeeze(vel.y(:,:,ind))+abs(min(min(squeeze(vel.y(:,:,ind)))));
            zz  = squeeze(vel.z(:,:,ind))+abs(min(min(squeeze(vel.z(:,:,ind)))));
            uu  = squeeze(vel.u(:,:,ind)-1);
            vv  = squeeze(vel.v(:,:,ind));
            ww  = squeeze(vel.w(:,:,ind));
            pp  = squeeze(vel.p(:,:,ind));
            
        end
        function [rms,x_station,indx] = get_rms(fpath,fname,istart,iend,tu,loadbox)
            counter = 0;
            
            if loadbox
                disp('reading the freestream pert fields ...')
                for i = istart:iend
                    counter = counter + 1;
                    vel =load([fpath '/' fname '_' num2str(i)]);
                    
                    indy = find(vel.y(1,:,1)>=0,1,'first');
                    indz = find(vel.z(:,1,1)>=0,1,'first');
                    
                    t(counter) = vel.t
                    signal(counter,:) =  (abs(vel.v(indz,indy,:)));
                    xx  = squeeze(vel.x(1,1,:));
                    %  yy  = squeeze(vel.y(:,:,ind))+abs(min(min(squeeze(vel.y(:,:,ind)))));
                    %  zz  = squeeze(vel.z(:,:,ind))+abs(min(min(squeeze(vel.z(:,:,ind)))));
                    %  uu  = squeeze(vel.u(:,:,ind));
                    %  vv  = squeeze(vel.v(:,:,ind));
                    %  ww  = squeeze(vel.w(:,:,ind));
                end
                T = t(end)-t(1);
                rms = sqrt((1/T)*trapz(t,signal.^2,1));
                indx = find(rms<=tu,1,'first');
                
                display('To be extracted x_station of the box')
                x_station = squeeze(vel.x(1,1,indx))
                figure()
                loglog(xx,rms)
                save(['rms_with_' num2str(iend-istart) 'fields'],'rms','t');
                save('rms.mat','rms','x_station','indx','T','signal','t','xx')
            else
                
                load(['rms'],'rms','T');
                size(T)
                T = T(iend)-T(istart);
                rms = sqrt((1/T)*trapz(t,signal.^2,1));
                indx = find(rms<=tu,1,'first');
                indx = find(rms>=tu,1,'first');
                
                display('To be extracted x_station of the box')
                x_station = squeeze(vel.x(1,1,indx))
            end
            loglog(xx,rms)
        end
        function VTK_binary(pressure,u,v,w,x,y,z,filename,varargin)
            % save_velocityAndPressureVTK_binary(pressure,u,v,w,x,y,z,filename)
            nr_of_elements=numel(x);
            fid = fopen(filename, 'w');
            
            %ASCII file header
            fprintf(fid, '# vtk DataFile Version 3.0\n');
            fprintf(fid, 'VTK from Matlab\n');
            fprintf(fid, 'BINARY\n\n');
            fprintf(fid, 'DATASET STRUCTURED_GRID\n');
%             if ~isempty(varargin{1})
%                 fprintf(fid, ['FIELD FieldData ' num2str(1) '\n']);
%                 fprintf(fid, 'TIME 1 1 double\n');
%                 fwrite(fid, varargin{1},'float','b');
%                 fprintf(fid, '\n');
%             end
            fprintf(fid, ['DIMENSIONS ' num2str(size(x,1)) ' ' num2str(size(x,2)) ' ' num2str(size(x,3)) '\n']);
            fprintf(fid, ['POINTS ' num2str(nr_of_elements) ' float\n']);
            fclose(fid);
            
            %append binary x,y,z data
            fid = fopen(filename, 'a');
            fwrite(fid, [reshape(x,1,nr_of_elements);  reshape(y,1,nr_of_elements); reshape(z,1,nr_of_elements)],'float','b');
            
            %append another ASCII sub header
            fprintf(fid, ['\nPOINT_DATA ' num2str(nr_of_elements) '\n']);
            fprintf(fid, 'VECTORS velocity_vectors float\n');
            
            %append binary u,v,w data
            fwrite(fid, [reshape(u,1,nr_of_elements);  reshape(v,1,nr_of_elements); reshape(w,1,nr_of_elements)],'float','b');
            
            %append some scalar data x_velocity
            fprintf(fid, '\nSCALARS x_velocity float\n'); %ASCII header
            fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
            fwrite (fid, reshape(u,1,nr_of_elements),'float','b'); %binary data
            
            %append some scalar data y_velocity
            fprintf(fid, '\nSCALARS y_velocity float\n'); %ASCII header
            fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
            fwrite (fid, reshape(v,1,nr_of_elements),'float','b'); %binary data
            
            %append some scalar data  z_velocity
            fprintf(fid, '\nSCALARS z_velocity float\n'); %ASCII header
            fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
            fwrite (fid, reshape(w,1,nr_of_elements),'float','b'); %binary data
            
            %append some scalar data pressure
            fprintf(fid, '\nSCALARS pressure float\n'); %ASCII header
            fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
            fwrite (fid, reshape(pressure,1,nr_of_elements),'float','b'); %binary data
            
            fclose(fid);


        end
        function VTK_binary_2d(pressure,u,v,x,y,filename,varargin)
            % save_velocityAndPressureVTK_binary(pressure,u,v,w,x,y,z,filename)
            nr_of_elements=numel(x);
            fid = fopen(filename, 'w');
            
            %ASCII file header
            fprintf(fid, '# vtk DataFile Version 3.0\n');
            fprintf(fid, 'VTK from Matlab\n');
            fprintf(fid, 'BINARY\n\n');
            fprintf(fid, 'DATASET STRUCTURED_GRID\n');
%             if ~isempty(varargin{1})
%                 fprintf(fid, ['FIELD FieldData ' num2str(1) '\n']);
%                 fprintf(fid, 'TIME 1 1 double\n');
%                 fwrite(fid, varargin{1},'float','b');
%                 fprintf(fid, '\n');
%             end
            fprintf(fid, ['DIMENSIONS ' num2str(size(x,1)) ' ' num2str(size(x,2)) '\n']);
            fprintf(fid, ['POINTS ' num2str(nr_of_elements) ' float\n']);
            fclose(fid);
            
            %append binary x,y,z data
            fid = fopen(filename, 'a');
            fwrite(fid, [reshape(x,1,nr_of_elements);  reshape(y,1,nr_of_elements)],'float','b');
            
            %append another ASCII sub header
            fprintf(fid, ['\nPOINT_DATA ' num2str(nr_of_elements) '\n']);
            fprintf(fid, 'VECTORS velocity_vectors float\n');
            
            %append binary u,v data
            fwrite(fid, [reshape(u,1,nr_of_elements);  reshape(v,1,nr_of_elements)],'float','b');
            
            %append some scalar data x_velocity
            fprintf(fid, '\nSCALARS x_velocity float\n'); %ASCII header
            fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
            fwrite (fid, reshape(u,1,nr_of_elements),'float','b'); %binary data
            
            %append some scalar data y_velocity
            fprintf(fid, '\nSCALARS y_velocity float\n'); %ASCII header
            fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
            fwrite (fid, reshape(v,1,nr_of_elements),'float','b'); %binary data
            
            
            %append some scalar data pressure
            fprintf(fid, '\nSCALARS pressure float\n'); %ASCII header
            fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
            fwrite (fid, reshape(pressure,1,nr_of_elements),'float','b'); %binary data
            
            fclose(fid);


        end
        %******************************************************************
        %                          temp
        %******************************************************************
        function phi = conv(u,base,varargin)
            if isempty(varargin)
                dx = 1 ; dy = 1 ;
            elseif isempty(varargin{2})
                dx = varargin{1}; dy = dx;
            else  
                dx =  varargin{1}; dy =  varargin{2};
            end            
            [Fx,Fy] = gradient(u,dx,dy);
            phi = (base.u.*Fx) + (base.v.*Fy) ;
        end
    end
end



