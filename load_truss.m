% Mechanical Engineering 332 - Computational Mechanics
% Stanford University
% Computing Assignment I - Loads on a Truss Bridge
% July 12, 2018
% Professor Dr. Peter Pinsky, TA Qing Yin
% Erik Lamoureux
% Code adapted from truss template provided. 

clear all;
close all;

%% Coordinate vectors, IEN and ID matrices
X = [0 10 20 30 40 50 60 70 10 20 30 40 50 60;
     0  0  0  0  0  0  0  0 20 20 20 20 20 20];
 
IEN = [1  8 1 2 3 4 5 6 7  9 10 11 12 13 2  3  4  5  6  7 3 4  5  3  6  4  5  6;
       9 14 2 3 4 5 6 7 8 10 11 12 13 14 9 10 11 12 13 14 9 9 10 12 11 13 14 14];
    
ID = [0 1 3 5 7  9 11 0 13 15 17 19 21 23;
      0 2 4 6 8 10 12 0 14 16 18 20 22 24];
    
% number of elements
ne = size(IEN,2);
% number of localnodes
na = size(IEN,1);
% number of space dimensions
ni = size(X,1);

% location matrix
L = zeros(na*ni,ne);

% For loop to iterate through Connectivity Array (IEN) to be applied to 
% the Equation Number Array (ID) to get the Destination Array (L) based on
% the forumla provided in the computing assingment. 
for i = 1:2
    for a = 1:2
        for e = 1:ne
            L(2*(a-1)+i,e) = ID(i,IEN(a,e));
        end
    end
end

%% Compute the stiffness matrix and load vector
% number of equations
neqns = nnz(ID);

% Stiffness matrix
K = zeros(neqns,neqns);
for e = 1:ne
    
    % The truss vector, truss length and inclination angle
    dX = X(:,IEN(2,e)) - X(:,IEN(1,e));
    l = norm(dX,2);
    cos_t = dX(1)/l;
    sin_t = dX(2)/l;
    
    % Determine T_e for the current element based on cos_t and sin_t 
    % found above based on the geometry. 
    T_e = [cos_t sin_t   0      0   ;
          0      0     cos_t  sin_t];
      
    % Determine the element stiffness matrix k_e 
    k_e = transpose(T_e)*((100/l).*[1 -1; -1 1])*T_e;
    
    % Iterate through k_e and assign values into the global stiffness
    % matrix K based on locations found in L (given they are not equal to
    % zero.
    for i = 1:length(k_e)
        if L(i,e) ~= 0
            for j = 1:length(k_e)
                if L(j,e) ~= 0
                    K(L(i,e),L(j,e)) = K(L(i,e),L(j,e)) + k_e(i,j);
                        % Specifically find the internal forces for element
                        % 23 and 27, as per the assignment. 
                        if e == 23
                            T_23 = T_e;
                            l_23 = l;
                        end
                        if e == 27
                            T_27 = T_e;
                            l_27 = l;
                        end
                end          
            end
        end
    end
end

% Load vector
F = zeros(neqns,1);
F(8,1) = -6;

% Solve the linear equation, equivalent to d = inv(K)*F;
d = K \ F;


%% Post-processing: plot the undeformed and deformed truss

% Plot the Undeformed truss
hold on
x_coord = X(1,:);
y_coord = X(2,:);
scatter(x_coord, y_coord)
axis([-10 80 -10 30])
for i = 1:ne
    A = IEN(1,i);
    B = IEN(2,i);
    line([x_coord(A) x_coord(B)],[y_coord(A) y_coord(B)])
end

% Plot the deformed truss
new_X = X(:,[2:7 9:end]);
x2_coord = new_X(1,:);
y2_coord = new_X(2,:);
mag = 1;
for i= 1:12
    temp1(i) = new_X(1,i) + mag*d(2*i-1);
    temp2(i) = new_X(2,i) + mag*d(2*i);
end

x2_coord(1) = 0;
y2_coord(1) = 0;
x2_coord(8) = 70;
y2_coord(8) = 0;
x2_coord(2:7) = temp1(1:6);
y2_coord(2:7) = temp2(1:6);
x2_coord(9:14) = temp1(7:end);
y2_coord(9:14) = temp2(7:end);
scatter(x2_coord, y2_coord, 'MarkerFaceColor', 'red')
for i = 1:ne
    A = IEN(1,i);
    B = IEN(2,i);
    line([x2_coord(A) x2_coord(B)],[y2_coord(A) y2_coord(B)], 'Color','red'); 
end

% Compute and find the internal force in element 23 and element 27    
Delta23 = T_23*[d(7); d(8); d(15); d(16)];
k_23 = (100/l_23).*[1 -1; -1 1];
N_23 = k_23*Delta23;
fprintf('The internal force in element 23 is %2.2f \n',N_23(2))

Delta27 = T_23*[d(7); d(8); d(23); d(24)];
k_27 = (100/l_27).*[1 -1; -1 1];
N_27 = k_27*Delta27;
fprintf('The internal force in element 27 is %2.2f \n',N_27(2))

