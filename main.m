function curve=main(Np,F,mu_rate,N_para)
% ind is assumed to be normalized
%%%%%%%%%%%%%%%%%%%%%%%   initialization  %%%%%%%%%%%%%%%%%%%%%%%
disp('the work is initializing....')
if nargin==0
    Np=30;F=0.8;mu_rate=0.5; N_para=225;
end
if nargin==1
    Np=30;F=0.8;mu_rate=0.5;
end
close all
clc;tic;figure(1);
%%%%%%%%%%%%%%%%%%%%%%%   environment parameters  %%%%%%%%%%%%%%%%%%%%%%%
disp('Setting environment parameters......')
% load('d_eku.mat')
gene=0; in_num=0;  ind=zeros(1,N_para);
individual=cell(1,Np);
fitness=1e9*ones(Np,1); least_fit=1e9;

%%%%%%%%%%%%%%%%%%%%%%%   recording results  %%%%%%%%%%%%%%%%%%%%%%%
curve=[];
disp('the individuals are initializing...')
disp('........')

for ii=1:Np
    tmp=rand(size(ind));
    individual(ii)=mat2cell(tmp,1,N_para);   
    fitness(ii)=toyfun(tmp);
end
toc

disp('the initial work has been done!')
disp('........')

while least_fit>0.05*N_para
gene=gene+1;
if mod(gene,1000)==0
    save 'curve.mat' curve;
%     plot(app.UIAxes,curve(:,1),curve(:,2),'-r');   
end
    for ii=1:Np
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% crossover & mutation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        in_num=in_num+1;
        tmp=cell2mat(individual(ii));
        a=myrandperm(Np,3,ii);
        tmp1=cell2mat(individual(a(1))); tmp2=cell2mat(individual(a(2))); tmp3=cell2mat(individual(a(3)));
        for mm=1:N_para
            if rand()<mu_rate
                tmp(mm)=tmp1(mm)+rand()*F*(tmp2(mm)-tmp3(mm));
                tmp(mm)=mod(tmp(mm),1);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmp_fit=toyfun(tmp);
        if tmp_fit<fitness(ii)
            fitness(ii)=tmp_fit;
            individual(ii)=mat2cell(tmp,1,N_para);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(groot,'CurrentFigure',1);
    stem(fitness)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% check the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [least_tmp,~]=min(fitness);
    if least_tmp~=least_fit 
        least_fit=least_tmp;
        curve=[curve;[gene,least_fit]];
        toc
        disp(['generations:',num2str(gene),'  The minimization: ',num2str(least_fit)])
        hold off
        pause(0.001)
    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function x=myrandperm(Np,n,ii)
    a=randperm(Np);
    x=a([1:n]);
    while sum(x==ii)
        a=randperm(Np);
        x=a([1:n]);
    end  
end