function [g_stress,Total_mass,dmdA]=fea_solver(Nodes_coordinates,Elements_connectivity,...
                            Material,Property,F_input,BC,NS_mass,...
                            Mass_formulation,OutStress,OutDisplacement,...
                            Node_displacement_constraint,...
                            Element_stress_constraint,...
                            Frequency_constraint)
% -------- SOLVER --------
% Calculate total number of nodes
Number_nodes=size(Nodes_coordinates,1);
% Calculate total number of Degrees of Freedom
DOF_total=Number_nodes*size(Nodes_coordinates,2);
% Calculate total number of Elements
Number_elements=size(Elements_connectivity,1);
% Calculate total number of Design Variables
% Number_DV=length(Design_variables);

% Check for errors in the input variables
if size(BC,1) < Number_nodes
    error('--- Not all nodes are properly constrained ---')
end
if size(BC,2) < 4
    error('--- Not all constraint directions are properly defined ---')
end

% Presize Variables to save computational time
counteri=0;
counterj=0;
Total_mass=0;
F=zeros(DOF_total,1);
K_global=zeros(DOF_total);
M_global=K_global;
NS_mass_global=M_global;
Cp=cell(1,Number_elements);
dMdA=cell(1,Number_elements);
dKdA=cell(1,Number_elements);
dMdA_global=cell(1,Number_elements);
dKdA_global=cell(1,Number_elements);
Connectivity=cell(1,Number_elements);
dKdA_global_reduced=cell(1,Number_elements);
dMdA_global_reduced=cell(1,Number_elements);
dmdA=zeros(1,Number_elements);

% figure(1)
% Loop over the element
for i=1:Number_elements
   % Find Elements Length
   N1=Nodes_coordinates(Elements_connectivity(i,1),:);
   N2=Nodes_coordinates(Elements_connectivity(i,2),:);
   % Plot Undeformed shape
   % plot3([N1(1);N2(1)],[N1(2);N2(2)],[N1(3);N2(3)],'ko-'); hold on
   % xlabel('x-dir'); ylabel('y-dir'); zlabel('z-dir')
   
% Find Vector and its norm
   V=[N2(1)-N1(1),N2(2)-N1(2),N2(3)-N1(3)];
   L=norm(V);
   % Find Direction Cosines
   Cx=V(1)/L;
   Cy=V(2)/L;
   Cz=V(3)/L;
   % Find Cp matrix needed to find stresses later on
   Cp{i}=(Material(Property(Elements_connectivity(i,3),2))/L).*...
         [-Cx -Cy -Cz Cx Cy Cz];
   % Find Euler Buckling Constant for later on
   Euler_buckling(i)=4*L^2/(pi*Property(Elements_connectivity(i,3),1)*...
                     Material(Property(Elements_connectivity(i,3),2)));
   % Find Element Mass and count the total mass
   Mass=Material(Property(Elements_connectivity(i,3),2),2)*...
        Property(Elements_connectivity(i,3),1)*L;
   dmdA(i)=Mass/Property(Elements_connectivity(i,3),1);
   Total_mass=Total_mass+Mass;
     
   % Connectivity of each element and each node
   Connectivity{i}=[(Elements_connectivity(i,1)*3)-2,...
                    (Elements_connectivity(i,1)*3)-1,...
                     Elements_connectivity(i,1)*3,...
                    (Elements_connectivity(i,2)*3)-2,...
                    (Elements_connectivity(i,2)*3)-1,...
                     Elements_connectivity(i,2)*3];
   
   % Element Stiffness Matrix in Global Coordinates
   K_element=(Property(Elements_connectivity(i,3),1)*...
           Material(Property(Elements_connectivity(i,3),2))/L).*...
           [Cx Cy Cz 0 0 0;
            0 0 0 Cx Cy Cz]'*...
           [1 -1;
            -1 1]*...
           [Cx Cy Cz 0 0 0;
            0 0 0 Cx Cy Cz];
        
    if strcmp(Mass_formulation,'lumped')==1
        % Element Lumped Mass Matrix in Global Coordinates
        M_element=diag(Mass/2.*ones(1,size(K_element,1)));
    end
    if strcmp(Mass_formulation,'consistent')==1
        % Element Consistent Mass Matrix in Global Coordinates
        M_element=(Mass/6).*[Cx Cy Cz 0 0 0;
                       0 0 0 Cx Cy Cz]'*...
                       [2 1;
                        1 2]*...
                       [Cx Cy Cz 0 0 0;
                        0 0 0 Cx Cy Cz];
    end

   % IF SENSITIVITIES ARE WANTED
   % dMdA for each element
   dMdA{i}=M_element./(Property(Elements_connectivity(i,3),1));
   % dKdA for each element
   dKdA{i}=K_element./(Property(Elements_connectivity(i,3),1));

   % Algorithm to Fill the Global Stiffness Matrix Element by Element
   dMdA_global{i}=zeros(DOF_total,DOF_total);
   dKdA_global{i}=zeros(DOF_total,DOF_total);

   for ii=Connectivity{i}
       counteri=counteri+1;
       for jj=Connectivity{i}
            counterj=counterj+1;
            K_global(ii,jj)=K_global(ii,jj)+K_element(counteri,counterj);
            M_global(ii,jj)=M_global(ii,jj)+M_element(counteri,counterj);
            dMdA_global{i}(ii,jj)=dMdA_global{i}(ii,jj)+dMdA{i}(counteri,counterj);
            dKdA_global{i}(ii,jj)=dKdA_global{i}(ii,jj)+dKdA{i}(counteri,counterj);
       end
       counterj=0;
   end
   counteri=0;
end

% Non Structural Masses
for i=1:size(NS_mass,1)
    NS_mass_global(NS_mass(i,1)*3,NS_mass(i,1)*3)=NS_mass(i,2);
    NS_mass_global(NS_mass(i,1)*3-1,NS_mass(i,1)*3-1)=NS_mass(i,2);
    NS_mass_global(NS_mass(i,1)*3-2,NS_mass(i,1)*3-2)=NS_mass(i,2);
end
M_global=M_global+NS_mass_global;

% Loop to populate F
for i=1:size(F_input,1)
   F((F_input(i,2)*3-3)+F_input(i,3))=F_input(i,1);
end

% Initialize Bookeeping Matrix for Rows and Columns Elimination
BC_to_eliminate=zeros(DOF_total,1);
counter=0;
for i=1:size(BC,1)
   for j=2:size(BC,2)
   counter=counter+1;
   BC_to_eliminate(counter)=BC(i,j);
   end
end
% Determine Index of Rows and Column to keep or eliminate
rowes_to_eliminate=find(BC_to_eliminate)';
rowes_to_keep=find(~BC_to_eliminate)';

% Reduce the Stiffness, Mass, and Force Matrix
for i = 1:DOF_total
    K_global_reduced = K_global;
    K_global_reduced(rowes_to_eliminate, :) = [];
    K_global_reduced(:,rowes_to_eliminate) = [];
    M_global_reduced = M_global;
    M_global_reduced(rowes_to_eliminate, :) = [];
    M_global_reduced(:,rowes_to_eliminate) = [];
    F_reduced=F;
    F_reduced(rowes_to_eliminate) = [];
    
    for ii=1:Number_elements
        dKdA_global_reduced{ii} = dKdA_global{ii};
        dKdA_global_reduced{ii}(rowes_to_eliminate, :) = [];
        dKdA_global_reduced{ii}(:,rowes_to_eliminate) = [];    
        dMdA_global_reduced{ii} = dMdA_global{ii};
        dMdA_global_reduced{ii}(rowes_to_eliminate, :) = [];
        dMdA_global_reduced{ii}(:,rowes_to_eliminate) = [];  
    end
end

% Calculate Eigenvalues and Eigenvector
[V,D]=eig(K_global_reduced,M_global_reduced,'vector');
w=sqrt(D);
w_hz=w./(2*pi);

% Frequency_name=cell(1,length(w_hz));
% for i=1:length(w_hz)
%     Frequency_name(i)={['Frequency ',num2str(i)]};
% end
 
% T=table(Frequency_name',w_hz);
% T.Properties.VariableNames={'Frequency','Hz'} 
    
Displacement_reduced=K_global_reduced\F_reduced;
Displacement(rowes_to_eliminate)=0;
Displacement(rowes_to_keep)=Displacement_reduced;

% Stress Calculation and Buckling Constraint Generation
Stress=zeros(1,Number_elements);
for i=1:Number_elements
    Stress(i)=Cp{i}*Displacement(Connectivity{i})';  
    g_buckling(i)=-1-Stress(i)*Euler_buckling(i);
end

% SENSITIVITY MATRICES
dudA_reduced=zeros(length(Displacement_reduced),Number_elements);
dldA_reduced=zeros(length(Displacement_reduced),Number_elements);
dudA=zeros(length(Displacement),Number_elements);
dldA=zeros(length(Displacement),Number_elements);
dSdA=zeros(1,Number_elements);
for i=1:Number_elements
   F_pseudo=-dKdA_global_reduced{i}*Displacement_reduced;
   dudA_reduced(:,i)=K_global_reduced\F_pseudo; 
   dudA(rowes_to_eliminate,i)=0;
   dudA(rowes_to_keep,i)=dudA_reduced(:,i);
   dSdA(i)=Cp{i}*dudA(Connectivity{i},i);
   for j=1:length(Displacement_reduced)
       dldA_reduced(j,i)=(V(:,j)'*(dKdA_global_reduced{i}-D(j).*...
                          dMdA_global_reduced{i})*V(:,j))/...
                          (V(:,j)'*M_global_reduced*V(:,j));
   end
   dldA(rowes_to_eliminate,i)=0;
   dldA(rowes_to_keep,i)=dldA_reduced(:,i);
end

% DISPLACEMENT CONSTRAINTS
g_disp=zeros(1,size(Node_displacement_constraint,1));
dgdispdA=zeros(length(g_disp),Number_elements);
for i=1:length(g_disp)
    Node=Node_displacement_constraint(i,1); %Node Number
    % Displacement X Y Z of the Node above
    Displacements_node=Displacement(3*Node-2:3*Node);
    % LOWER BOUND
    if Node_displacement_constraint(i,3) == 1
        if Node_displacement_constraint(i,4) >= 0
        g_disp(i)=-Displacements_node(2)/Node_displacement_constraint(i,4)+1;
        dgdispdA(i,:)=(-1/Node_displacement_constraint(i,4)).*...
            dudA(3*Node-3+Node_displacement_constraint(i,2),:);
        elseif Node_displacement_constraint(i,4) <= 0
        g_disp(i)=Displacements_node(2)/Node_displacement_constraint(i,4)-1;
        dgdispdA(i,:)=(1/Node_displacement_constraint(i,4)).*...
            dudA(3*Node-3+Node_displacement_constraint(i,2),:);
        end
    end

    % UPPER BOUND
    if Node_displacement_constraint(i,3) == 2
        if Node_displacement_constraint(i,4) >= 0
        g_disp(i)=Displacements_node(2)/Node_displacement_constraint(i,4)-1;
        dgdispdA(i,:)=(1/Node_displacement_constraint(i,4)).*...
            dudA(3*Node-3+Node_displacement_constraint(i,2),:);
        elseif Node_displacement_constraint(i,4) <= 0
        g_disp(i)=-Displacements_node(2)/Node_displacement_constraint(i,4)+1;
        dgdispdA(i,:)=(-1/Node_displacement_constraint(i,4)).*...
            dudA(3*Node-3+Node_displacement_constraint(i,2),:);
        end   
    end 
end   

% STRESS CONSTRAINTS
g_stress=zeros(1,size(Element_stress_constraint,1));
dgstressdA=zeros(1,length(g_stress));
for i=1:length(g_stress)
    Element=Element_stress_constraint(i,1); %Element Number
    if Element_stress_constraint(i,2) == 1
        g_stress(i)=Stress(Element)/Element_stress_constraint(i,3)-1;
    elseif Element_stress_constraint(i,2) == 2
        g_stress(i)=Stress(Element)/Element_stress_constraint(i,3)-1;
    end
    dgstressdA(i)=(1/Element_stress_constraint(i,3)).*dSdA(Element);
end

% FREQUENCY CONSTRAINTS
g_frequency=zeros(1,size(Frequency_constraint,1));
dgwdA=zeros(size(Frequency_constraint,1),Number_elements);
for i=1:size(Frequency_constraint,1)
    Frequency=Frequency_constraint(i,1); %Frequency Index
    if Frequency_constraint(i,2) == 1 
        g_frequency(i)=-w_hz(Frequency)/Frequency_constraint(i,3)+1;
        dgwdA(i,:)=(-1/(2*w_hz(Frequency))).*dldA_reduced(Frequency,:);
    elseif Frequency_constraint(i,2) == 2 
        g_frequency(i)=w_hz(Frequency)/Frequency_constraint(i,3)-1;
        dgwdA(i,:)=(1/(2*w_hz(Frequency))).*dldA_reduced(Frequency,:);
    end
end       

% % For Loop for Plotting new displacement
% magnification=15;
% % Reshape Displacement
% counter=0;
% Displacement_reshaped=zeros(size(Nodes_coordinates));
% for i=1:3:length(Displacement)
%    counter=counter+1;
%    Displacement_reshaped(counter,1:3)=Displacement(i:i+2);
% end
% 
% New_nodes_coordinates=magnification.*...
%                       Displacement_reshaped+Nodes_coordinates;
% for i=1:Number_elements
%    N1=New_nodes_coordinates(Elements_connectivity(i,1),:);
%    N2=New_nodes_coordinates(Elements_connectivity(i,2),:);
%    plot3([N1(1);N2(1)],[N1(2);N2(2)],[N1(3);N2(3)],'ro-'); hold on
%    xlabel('x-dir'); ylabel('y-dir'); zlabel('z-dir')
% end
% 
% view(0,90)
% grid on; grid minor
% title(['Undeformed and Deformed Plot | Magnification = ',...
%         num2str(magnification),'X'])
%     
% % Loop for plotting Eigen Modes
% Displacement_Eigen=zeros(size(Displacement));
% Modes_to_plot=4; % Please use only even number for now!
% magnification=300;
% figure
% for ii=1:Modes_to_plot
%     Displacement_Eigen(rowes_to_eliminate)=0;
%     Displacement_Eigen(rowes_to_keep)=V(:,ii);
%     counter=0;
%     Displacement_Eigen_reshaped=zeros(size(Nodes_coordinates));
%     for i=1:3:length(Displacement_Eigen)
%        counter=counter+1;
%        Displacement_Eigen_reshaped(counter,1:3)=Displacement_Eigen(i:i+2);
%     end
%     
%     New_nodes_coordinates=magnification.*...
%                       Displacement_Eigen_reshaped+Nodes_coordinates;
%                   
%     for i=1:Number_elements
%        % DEFORMED PLOT
%        N1=New_nodes_coordinates(Elements_connectivity(i,1),:);
%        N2=New_nodes_coordinates(Elements_connectivity(i,2),:);
%        subplot(Modes_to_plot/2,Modes_to_plot/2,ii)
%        title(['Eigen Mode ' num2str(ii)])
%        plot3([N1(1);N2(1)],[N1(2);N2(2)],[N1(3);N2(3)],'ro-'); hold on
%        xlabel('x-dir'); ylabel('y-dir'); zlabel('z-dir')
%        view(0,90)
%        grid on; grid minor
%        % UNDEFORMED PLOT
%        N1=Nodes_coordinates(Elements_connectivity(i,1),:);
%        N2=Nodes_coordinates(Elements_connectivity(i,2),:);
%        % Plot Undeformed shape
%        plot3([N1(1);N2(1)],[N1(2);N2(2)],[N1(3);N2(3)],'ko-');
%        axis([0 800 -100 500]) % It has to be a manual update
%     end
% end
% % Uncomment for bigger size window
% % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.75 0.75]);
% suptitle(['Eigen Modes Plot | Magnification = ' num2str(magnification) 'X'])
% 
% counter=0;
% Displacement_name=cell(1,DOF_total);
% for i=1:size(Nodes_coordinates,2):DOF_total
%     counter=counter+1;
%     Displacement_name(i)={['u',num2str(counter)]};
%     Displacement_name(i+1)={['v',num2str(counter)]};
%     Displacement_name(i+2)={['w',num2str(counter)]};
% end
% 
% if strcmp(OutDisplacement,'All') == 1
%     T=table(Displacement_name',Displacement');
%     T.Properties.VariableNames={'Node','Value'}    
% end
% 
% Element_name=cell(1,Number_elements);
% for i=1:Number_elements
%     Element_name(i)={['Element ',num2str(i)]};
% end
% 
% if strcmp(OutStress,'All') == 1
%     T=table(Element_name',Stress');
%     T.Properties.VariableNames={'Element','Value'}    
% end
end