
%    Topology optimization with simulated annealing

% This is a code for replicating results of a paper entitlred:
% ''A New Non-Gradient-Based Topology Optimization Algorithm with Black-White Density and Manufacturability Constraints''
% Written by H. R. N. and T. G. G for educational purposes 
% parameters can be used from the paper to replicate the results 
% Just 0 or 1 element (void or solid) without intermediate density
%  Just one elemet change the density with an opposite element in each iteration
% The optimizatin problem is cantilever or MBB-beam
% All functions are embedid in this single file, including FE and generator
% of new solution

function TO(nelx, nely, N, vol) % Define the umber of elements in x or y direction, number of iterations and volume fraction
penal=1; 
  initime=cputime;
  Tmax=100; % starting temperature
  Tmin=0.001; % end temperature
  nov=nelx*nely; %number of elements, it will be the number of variable
  alphaa=0.85; % cooling factor
  te=Tmax;
  min_den=0.001; % Using a very small value for density to avoid singularity in FEM
  %generetae a current random solution
   xi(1:(nely*vol), 1:(nelx))=1;
   xi((nely*vol)+1:nely,1:nelx)=min_den;
% xi=randi([0,1])*ones(nely,nelx);
  colormap(gray); imagesc(-xi); axis equal; axis tight; axis off;pause(1e-5);
  aux=FE(nelx, nely, xi, penal);
  Ji=computeCompliance(nelx, nely, xi, penal, aux); % calculation of objective function
  fun=[];
  temp=[];
  acepted=[];
  rejected=[];
  mm= round((log(Tmin)-log(Tmax))/log(alphaa))+1;
%   hist=ones(nely,nelx,N,mm);
  cont=1;
  while te>Tmin
    a1=0; %acepted solution count for each temperature
    a2=0; %rejected solution for each temperature
    for ii=1:N
[xj,x1,x2,y1,y2,r1,r2]=randomselection(xi,nely,nelx,min_den); % randomly change two opposite elements

      aux2=FE(nelx, nely, xj, penal);
      Jn2=computeCompliance(nelx, nely, xj, penal,aux2);
      J1=dist(xj,nelx,nely); % To see if the neighboring elements are the same solid or void as objective 2nd
      J2=J1/3; % To normalize the continuity criteria
      Jn=Jn2+J2;
      dE=Jn-Ji;
      p=exp(-dE/te);
      if dE<=0
        xi=xj; %the current solution receiver the newsolution
        Ji=Jn;
        a1=a1+1;
      elseif rand()<=p
        xi=xj; %the current solution receiver the newsolution
        Ji=Jn;
        a1=a1+1;
      else
          a2=a2+1;

      end    
      
    end
    acepted(cont)=a1;
    rejected(cont)=a2;

    fun(cont)=Ji;
    temp(cont)=te;
    te=te*alphaa
    cont=cont+1;
    %figure(1);
    colormap(gray); imagesc(-xi); axis equal; axis tight; axis off;pause(1e-5);
    end
  %###### display the best value and plot Function vs Temperature
  
  volfracf=sum(sum(xi))/(nelx*nely);
  fprintf('the best value %f \n', Ji);
  fprintf('The volume frac of the best solution %f \n', volfracf)
  
  endtime=cputime;
  runtime=endtime-initime;
  fprintf('the Running time %f \n', runtime);
  save thiago8.mat
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function to analysis thedeslocation of nods
%
function [U]=FE(nelx,nely,x,penal)
  [KE] = lk; 
  K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
  F = sparse(2*(nely+1)*(nelx+1),1); %1to2 for other loads
  U = zeros(2*(nely+1)*(nelx+1),1);  %here too
  for elx = 1:nelx
    for ely = 1:nely
      n1 = (nely+1)*(elx-1)+ely; 
      n2 = (nely+1)* elx   +ely;
      edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
      K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
    end
  end
  % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
  pp=2*(nelx+1)*(nely+1);
 % F(2,1) = -1;
 %fixeddofs   = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
  % Short Cantilever problem (commetn two previous lines)
  F(2*(nelx+1)*(nely+1),1)= -1;
  % F(2*(nelx)*(nely+1)+2,2)= 1; %second load
  fixeddofs = [1:2*(nely+1)];
  alldofs     = [1:2*(nely+1)*(nelx+1)];
  freedofs    = setdiff(alldofs,fixeddofs);
  % SOLVING
  U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
  U(fixeddofs,:)= 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% function to create the stiffines matrix 
function [KE]=lk
  E = 1.; 
  nu = 0.3;
  k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
     -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
  KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                    k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                    k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                    k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                    k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                    k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                    k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                    k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to compute the compliance
%
function [J]=computeCompliance(nelx,nely,x,penal, U) %elemsIn,voidEps,
  %Compute the compliance of the entire mesh
  [KE] = lk;J = 0;
  for elx = 1:nelx
      for ely = 1:nely
          n1 = (nely+1)*(elx-1)+ely; 
          n2 = (nely+1)* elx   +ely;

          Ue = U([2*n1-1;2*n1;2*n2-1;2*n2;2*n2+1;2*n2+2;2*n1+1;2*n1+2]);
          J = J + x(ely, elx)^penal*Ue'*KE*Ue; % transposta deslocamento* matriz rigidez * deslocamento
      end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xj,x1,x2,y1,y2,r1,r2]=randomselection(xi,nely,nelx,min_den)
xj=xi;
y1=randi(nely);
x1=randi(nelx);
r1=xi(y1,x1);
if r1 == 1
    r1=min_den;
else
    r1=1;
end
xj(y1,x1)=r1;
y2=randi(nely);
x2=randi(nelx);
r2=xi(y2,x2);

while r2 ~= r1;
    y2=randi(nely);
    x2=randi(nelx);
    r2=xi(y2,x2);
end
if r2 == 1
    r2= min_den;
else
    r2=1;
end
xj(y2,x2)=r2;
end

function [J1]=dist(xj,nelx,nely)
J1=0;
for k=1:nely
    for l=1:nelx
        num(k,l)=0;
        dis(k,l)=0;
        r=xj(k,l);
        for i1=max(1,k-1):min(k+1,nely)
            for j1=max(1,l-1):min(l+1,nelx)
                num(k,l)=num(k,l)+1;
                if xj(i1,j1)==r
                    dis(k,l)=dis(k,l)+1;
                end
            end
        end
        d(k,l)=(num(k,l)-dis(k,l))/num(k,l);   
        J1=J1+d(k,l);
    end
end
end