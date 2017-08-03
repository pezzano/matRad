n=8;
% radius=4;
% sigma_ray=7;
% sigma_sub=5.5;
% method = 'square';
methodi = ['circle';'square'];
% startingPoint(1) = sqrt(sigma_ray^2-sigma_sub^2);
% startingPoint(2) = sigma_sub^2*startingPoint(1)^2/sigma_ray^2;
startingPoint(1) = 16;
startingPoint(2) = 1;
sigma_rayvec = linspace(3,10,8);
%sigma_rayvec = 6;
for met = 1:1
    if met==1
        method='circle';
    else
        method='square';
    end
    for s=1:8
        sigma_ray = sigma_rayvec(s);
        radiusvec = linspace(sigma_ray/5, sigma_ray/4, 3);
        for k=1:3
            radius = radiusvec(k);
            sigma_subvec = linspace(sigma_ray/3.5,sigma_ray/2,10);
            for m=1:10
                sigma_sub = sigma_subvec(m);
                %for test = 1:10;
                tic;
                
                if ~exist('method','var')
                    method = 'square';
                end
                
                if ~strcmp(method,'square') && ~strcmp(method,'circle')
                    error('method not supported');
                end
                
                % setting positions of sub-beams
                if strcmp(method,'square')
                    numOfSub = (2*n +1)^2;
                    points = linspace(-radius*(sqrt(numOfSub)-1)/2,radius*(sqrt(numOfSub)-1)/2,sqrt(numOfSub));
                    posx = points'*ones(1,sqrt(numOfSub));
                    posy = posx';
                else
                    numOfSub = (2^n -1)*6 +1;
                    ang = zeros(1,1);
                    posx = zeros(1,1);
                    posy = zeros(1,1);
                    radiusShell = zeros(1,1);
                    for i=1:n
                        SubsInShell = (2^i -1)*6 +1 - ((2^(i-1) -1)*6 +1 );
                        % this takes the sub-beams index in one shell
                        ang = cat(2, ang, pi .* linspace(0,2-2/SubsInShell, SubsInShell));
                        radiusShell = cat(2, radiusShell, i.*radius.*ones(1, SubsInShell));
                    end
                    posx = cat(2, posx, posx(1) + radiusShell(2:end).*cos(ang(2:end)));
                    posy = cat(2, posy, posy(1) + radiusShell(2:end).*sin(ang(2:end)));
                end
                
                x0 = -3*sigma_ray:sigma_ray/70:3*sigma_ray;
                y0 = x0;
                
                gaussian2 = @(x, y, mux, muy ,sig) (2*pi*sig^2)^(-1) .* exp(-(x-mux).^2/(2*(sig^2)))' * exp(-(y-muy).^2/(2*(sig^2)));
                f1 =@(xi,yi) gaussian2(xi,yi,0,0,sigma_ray);
                f2 = @(X,xi,yi) -f1(xi,yi);
                for i=1:numOfSub
                    f2 = @(X,xi,yi) f2(X,xi,yi) + X(1) .* ...
                        gaussian2(posx(i),posy(i),0,0,X(2)).*gaussian2(xi,yi,posx(i),posy(i),sigma_sub);
                end
                
                
                f3 = @(X) 0;
                x1 = -3*sigma_ray:sigma_ray/4*3:3*sigma_ray;
                [xf,yf]=meshgrid(x1,x1);
                numOfPoints = numel(xf);
                for i=1:numOfPoints
                    f3 = @(X) f3(X) + (f2(X,xf(i),yf(i)))^2;
                    %             if mod(i,50)==0
                    %                 disp(i)
                    %             end
                end
                
                [X1,fx] = fminsearch(f3, startingPoint);
                X_rad(k,m,s,:) = X1;
                todisp = [sigma_ray radius sigma_sub X1];
                disp(todisp)
                %[X1,fx] = fminsearch(f3, [2,1],optimset('MaxIter',1000))
                
                %timex = toc;
                %tempoi(test) = toc;
                timei_rad(k,m,s) = toc;
                %
                %         subGauss = zeros([length(x0) length(x0) numOfSub]);
                f1test = gaussian2(x0,y0,0,0,sigma_ray);
                f2test = @(X) 0;
                for i=1:numOfSub
                    f2test = @(X) f2test(X) + X(1) .* ...
                        gaussian2(posx(i),posy(i),0,0,X(2)).*gaussian2(x0,y0,posx(i),posy(i),sigma_sub);
                    subGauss(:,:,i) = X1(1) .* ...
                        gaussian2(posx(i),posy(i),0,0,X1(2)).*gaussian2(x0,x0,posx(i),posy(i),sigma_sub);
                end
                
                %maxi_test(test) = max(max(abs((f2test(X1)-f1test)./max(max(f1test)).*100)));
                maxi_rad(k,m,s) = max(max(abs((f2test(X1)-f1test)./max(max(f1test)).*100)));
                pointi_rad(k,m,s) = radius;
            end
        end
    end
    
    % figure
    % subplot(2,2,1)
    % surf(f1test,'LineWidth',.1,'EdgeAlpha',.3)
    % title('single gauss')
    %
    % subplot(2,2,2)
    % surf(f2test(X1),'LineWidth',.1,'EdgeAlpha',.3)
    % title('superpos gaussian')
    %
    % subplot(2,2,3)
    % surf(f2test(X1)-f1test,'LineWidth',.1,'EdgeAlpha',.3)
    % title('abs diff')
    %
    %
    % subplot(2,2,4)
    % surf(((f2test(X1)-f1test)./max(max(f1test)).*100),'LineWidth',.1,'EdgeAlpha',.3)
    % title('rel diff')
    
    %title(strcat('max percentage error =  ',num2str(max(max((f2test(X1)-f1test)./f1test.*100))),...
    %   '%      elapsed time =  ',num2str(timex),'s    for  ',num2str(numOfPoints),'  points'))
    
    
    
    if method == 'square'
        save('pezzWeightsData_square3')
    else
        save('pezzWeightsData_circle8')
    end
    
end
