function [res]=terzaghi_solution(x,time,L,c,p_o,varargin)

% SOLUTION OF THE 1D DIFFUSION EQUATION OF TERZAGUI'S CONSOLIDATION PROBLEM

% OUTPUT
% res: pressure

% INPUT
% x: spatial coordinate of the domain
% time: time
% L: size of the domain (thickness of soil layer)
% c: Consolidation coefficient
% p_o: initial pressure through the whole domain
% varargin: number of terms to be considered in the infinite sumation

    if nargin < 3
        n=200;
    else
        n=varargin{1};
    end
    res=0.*(x'*time);
   
    for k=1:2:n
        res_1= sin(pi * k * ( x/L ) /2.) ;
        res_2 = exp(-(c*time/(4*L^2))*((pi*k)^2.));
        res=res + p_o*(4./(pi*k))*(res_1'*res_2);
    end
       
end
