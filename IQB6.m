% Group Members
% Sahil Yadav - 2016184
% Sneha Sinha - 2016098
% Ritu Kumari - 2016078

No_Runs=5;                                % Total runs
TimeStep=5000;                            % Monte carlo time-step
L=50;                                     % Lattice size
Total_Molecules=100;                      % Total no. of molecules
nx=0;
ny=0;
mx=0;
my=0;
lattice=zeros(L, L);                      % 50x50 lattice
x_coordinate=zeros(Total_Molecules, 1);
y_coordinate=zeros(Total_Molecules, 1);
x0=zeros(Total_Molecules, 1);             % initial x-coordinates of molecules
y0=zeros(Total_Molecules, 1);             % initial y-coordinates of molecules
px=zeros(No_Runs*Total_Molecules, 1);     % P(x)
py=zeros(No_Runs*Total_Molecules, 1);     % P(y)
r=zeros(No_Runs*Total_Molecules, 1);      % variable for calculating <r>
rsquare=zeros(TimeStep, 1);
r2=0;                                     % variable for calculating mean_squared_displacement <r2>
rx=0;                                     % variable for calculating <x>
ry=0;                                     % variable for calculating <y>
j=0;
for i=1:No_Runs                           % Loop over Number of runs
    mol=1;
    while mol<=Total_Molecules            % Loop over Total molecules
        nx=randi([1, L], 1, 1);           % Generate a random x-coordinate between 1 and 50
        ny=randi([1, L], 1, 1);           % Generate a random y-coordinate between 1 and 50
        
        % 1 specifies if there is a molecule; 0 if the lattice site is not already occupied
        if (lattice(nx, ny)~=1)           
            x_coordinate(mol)=nx;
            y_coordinate(mol)=ny;
            x0(mol)=nx;
            y0(mol)=ny;
            lattice(nx, ny)=1;         
            mol=mol+1;
        end
    end
    for ctr=1:TimeStep                      % Loop over Monte Carlo Time-Steps
        for trial=1:Total_Molecules  
            random_molecule=randi([1, Total_Molecules], 1, 1);    % get a random molecule
            mx=x_coordinate(random_molecule);                     % get the x-coordinate of random molecule
            my=y_coordinate(random_molecule);                     % get the y-coordinate of random molecule
            move=randi([1, 4], 1, 1);                             % generate a random number between 1 and 4
            switch(move)
                case 1                                     % move the molcule down if no other molecule is already there.
                    if (lattice(mx, two(my))~=1)
                        x_coordinate(random_molecule)=mx;
                        y_coordinate(random_molecule)=two(my);
                        lattice(mx, my)=0;
                        lattice(mx, two(my))=1;
                    end
                case 2                                      % move the molcule up if no other molecule is already there.
                    if (lattice(mx, one(my, L))~=1)
                        x_coordinate(random_molecule)=mx;
                        y_coordinate(random_molecule)=one(my, L);
                        lattice(mx, my)=0;
                        lattice(mx, one(my, L))=1;
                    end
                case 3                                          % move the molcule to left if no other molecule is already there.
                    if (lattice(two(mx), my)~=1)
                        x_coordinate(random_molecule)=two(mx);
                        y_coordinate(random_molecule)=my;
                        lattice(mx, my)=0;
                        lattice(two(mx), my)=1;
                    end
                case 4                                           % move the molcule to right if no other molecule is already there.
                    if (lattice(one(mx, L), my)~=1)
                        x_coordinate(random_molecule)=one(mx, L);
                        y_coordinate(random_molecule)=my;
                        lattice(mx, my)=0;
                        lattice(one(mx, L), my)=1;
                    end
            end
        end
        for mol=1:Total_Molecules
        r2 = (x_coordinate(mol)-x0(mol))*(x_coordinate(mol)-x0(mol)) + (y_coordinate(mol)-y0(mol))*(y_coordinate(mol)-y0(mol)) + r2;
        rx = (x_coordinate(mol)-x0(mol)) + rx;
        ry = (y_coordinate(mol)-y0(mol)) + ry;
        
        rsquare(ctr) = (x_coordinate(mol)-x0(mol))*(x_coordinate(mol)-x0(mol)) + (y_coordinate(mol)-y0(mol))*(y_coordinate(mol)-y0(mol)) + rsquare(ctr);
        end
    end
end
N=No_Runs*Total_Molecules;
fprintf('Mean Square Displacement = %f\n', r2/N);
mean_dis=sqrt((rx/N).^2 + (ry/N).^2);              % Mean Displacement
fprintf('Mean Displacement = %f\n', mean_dis);
fprintf('<x> = %f , <y> = %f\n', rx/N, ry/N);
rsquare = rsquare/N;
disp(rsquare);
%histfit(rsquare);
plot(rsquare);
xlabel('Monte Carlo Time Steps (= 5000)');
ylabel('Mean Square Displacement (r2)');