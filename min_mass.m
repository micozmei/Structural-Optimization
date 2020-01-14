clear all; close all; clc

lb=ones(1,10).*0.1;
ub=[];
A=[];
b=[];
Aeq=[];
beq=[];

x0=ones(1,10).*5;

nonlcon=@g;
fun=@mass;
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[solution,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)

function [c,ceq] = g(x)
c = stress(x);
ceq = [];
end

function [g_stress] = stress(x)
% ------- INPUT -------
% [x1 y1 z1;
% [x2 y2 z2]; etc...

Nodes_coordinates=360.*[2 1 0;
                   2 0 0;
                   1 1 0;
                   1 0 0;
                   0 1 0;
                   0 0 0];

% Element card [Node1, Node2, Property] - #rows is element number 
Elements_connectivity=[5 3 1;
                       3 1 2;
                       6 4 3;
                       4 2 4;
                       3 4 5;
                       1 2 6;
                       5 4 7;
                       6 3 8;
                       3 2 9;
                       4 1 10];

% Material Card [E, rho, v, Sigmault, etc...] to be defined
Material=[10e6, 0.1];


% Force Vector [Force,Node,Component]
F_input=[-100e3,2,2;
         -100e3,4,2];

% Boundary Conditions [Node, X, Y, Z] - 1 = constraint - 0 = non constraint
BC=[1 0 0 1;
    2 0 0 1;
    3 0 0 1;
    4 0 0 1;
    5 1 1 1;
    6 1 1 1];

% Non Structural Masses [Node, value]
NS_mass=[1 0];

% Mass Formulation 'consistent' or 'lumped'
Mass_formulation='lumped';

% Optimization Constraint Request
% Node Displacement Constraint [Node, Direction, Min/Max, Value]
% Direction: X-axis = 1 | Y-axis = 2 | Z-axis = 3
% Min/Max: Min = 1 | Max = 2
Node_displacement_constraint=[];
                          
% Element Stress Constraint [Element, Tension/Compression, Value] 
% Tension/Compression: Tension = 1 | Compression = 2
Element_stress_constraint=[1 1 25000;
                           1 2 -25000;
                           2 1 25000;
                           2 2 -25000;
                           3 1 25000;
                           3 2 -25000;
                           4 1 25000;
                           4 2 -25000;
                           5 1 25000;
                           5 2 -25000;    
                           6 1 25000;
                           6 2 -25000; 
                           7 1 25000;
                           7 2 -25000; 
                           8 1 25000;
                           8 2 -25000;                            
                           9 1 75000;
                           9 2 -75000;                            
                           10 1 25000;
                           10 2 -25000];                            
                                            
% Frequency Constraint [Frequency Index, Min/Max, value]
% Min/Max: Min = 1 | Max = 2                       
Frequency_constraint=[];
                           
% -------- OUTPUT REQUEST --------
% S, Stress
% U, Displacement
% ['S',Element] can be 'All'
% ['U',Node] can be 'All'

OutStress=['All'];
OutDisplacement=['All'];

% Property Card [Area, Mat]
Property=[x(1),1;
          x(2),1;
          x(3),1;
          x(4), 1;
          x(5), 1;
          x(6), 1;
          x(7), 1;
          x(8), 1;
          x(9), 1; 
          x(10), 1];

[g_stress,Total_mass]=fea_solver(Nodes_coordinates,Elements_connectivity,...
                            Material,Property,F_input,BC,NS_mass,...
                            Mass_formulation,OutStress,OutDisplacement,...
                            Node_displacement_constraint,...
                            Element_stress_constraint,...
                            Frequency_constraint);
end

function [Total_mass] = mass(x)
% ------- INPUT -------
% [x1 y1 z1;
% [x2 y2 z2]; etc...

Nodes_coordinates=360.*[2 1 0;
                   2 0 0;
                   1 1 0;
                   1 0 0;
                   0 1 0;
                   0 0 0];

% Element card [Node1, Node2, Property] - #rows is element number 
Elements_connectivity=[5 3 1;
                       3 1 2;
                       6 4 3;
                       4 2 4;
                       3 4 5;
                       1 2 6;
                       5 4 7;
                       6 3 8;
                       3 2 9;
                       4 1 10];

% Material Card [E, rho, v, Sigmault, etc...] to be defined
Material=[10e6, 0.1];

% Force Vector [Force,Node,Component]
F_input=[-100e3,2,2;
         -100e3,4,2];

% Boundary Conditions [Node, X, Y, Z] - 1 = constraint - 0 = non constraint
BC=[1 0 0 1;
    2 0 0 1;
    3 0 0 1;
    4 0 0 1;
    5 1 1 1;
    6 1 1 1];

% Non Structural Masses [Node, value]
NS_mass=[1 0];

% Mass Formulation 'consistent' or 'lumped'
Mass_formulation='lumped';

% Optimization Constraint Request
% Node Displacement Constraint [Node, Direction, Min/Max, Value]
% Direction: X-axis = 1 | Y-axis = 2 | Z-axis = 3
% Min/Max: Min = 1 | Max = 2
Node_displacement_constraint=[];
                          
% Element Stress Constraint [Element, Tension/Compression, Value] 
% Tension/Compression: Tension = 1 | Compression = 2
Element_stress_constraint=[1 1 25000;
                           1 2 -25000;
                           2 1 25000;
                           2 2 -25000;
                           3 1 25000;
                           3 2 -25000;
                           4 1 25000;
                           4 2 -25000;
                           5 1 25000;
                           5 2 -25000;    
                           6 1 25000;
                           6 2 -25000; 
                           7 1 25000;
                           7 2 -25000; 
                           8 1 25000;
                           8 2 -25000;                            
                           9 1 75000;
                           9 2 -75000;                            
                           10 1 25000;
                           10 2 -25000];                            
                                           
% Frequency Constraint [Frequency Index, Min/Max, value]
% Min/Max: Min = 1 | Max = 2                       
Frequency_constraint=[];
                                  
% -------- OUTPUT REQUEST --------
% S, Stress
% U, Displacement
% ['S',Element] can be 'All'
% ['U',Node] can be 'All'

OutStress=['All'];
OutDisplacement=['All'];

% Property Card [Area, Mat]

Property=[x(1),1;
          x(2),1;
          x(3),1;
          x(4), 1;
          x(5), 1;
          x(6), 1;
          x(7), 1;
          x(8), 1;
          x(9), 1; 
          x(10), 1];

[g_stress,Total_mass]=fea_solver(Nodes_coordinates,Elements_connectivity,...
                            Material,Property,F_input,BC,NS_mass,...
                            Mass_formulation,OutStress,OutDisplacement,...
                            Node_displacement_constraint,...
                            Element_stress_constraint,...
                            Frequency_constraint);
end