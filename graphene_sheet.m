
clear all
close all
clc
C=[];

% Define lattice constants of Graphene 
latt_c = 1.42; 

% Define Mass 
C_mass = 12.011;
%Define length for one repeatation 
ux = latt_c*3.0;
uy = latt_c*sqrt(3);

% Define coordinate of for one repeatation:

u = [0.0*latt_c 0.5*sqrt(3)*latt_c
    0.5*latt_c 0.0*latt_c
    1.5*latt_c 0.0*latt_c
    2.0*latt_c 0.5*sqrt(3)*latt_c];

n = round(input('Angstrom length in X direction = ')/ux);

m = round(input('Angstrom length in Y direction = ')/uy);

s_z = 0 %z_coordinate

num = 0.0;

for j = 1:m;
    for i = 1:n;
        for c = 1:4;
            num = num+1;
            C(num,1) = u(c,1)+(i-1)*ux;
            C(num,2) = u(c,2)+(j-1)*uy;
            C(num,3) = s_z;            
        end
    end
end

n_atom = length(C); 

atom_dat = [(1:n_atom)' C]

x_lo = min(atom_dat(:,2));
x_hi = max(atom_dat(:,2));
y_lo = min(atom_dat(:,3));
y_hi = max(atom_dat(:,3));
z_lo = min(atom_dat(:,4));
z_hi = max(atom_dat(:,4));

% Write the xyz file :
filename3 = ['Gr_H2O.xyz'];
fid3 = fopen(filename3,'wt');
fprintf(fid3,'%4i\n',n_atom);
fprintf(fid3,'%s\n','Atoms');

for i = 1:n_atom;
    
    fprintf(fid3,'C %12.8f %12.8f %12.8f\n',C(i,1:3));
    
end

fclose(fid3);
filename1 = sprintf('graphene.lammps');
fid1 = fopen(filename1,'wt');

fprintf(fid1,'%s\n','LAMMPS readable file ');
fprintf(fid1,'%s\n','');
fprintf(fid1,'\t\t%4i %s\n',n_atom,'atoms');
fprintf(fid1,'\t\t%s\n','1 atom types');

fprintf(fid1,'\t%12.8f %12.8f \t%s\n',x_lo,x_hi+latt_c,'xlo xhi');
fprintf(fid1,'\t%12.8f %12.8f \t%s\n',y_lo,y_hi+1.22976,'ylo yhi');
fprintf(fid1,'\t%12.8f %12.8f \t%s\n',z_lo-50,z_hi+50,'zlo zhi');
fprintf(fid1,'%s\n','');
fprintf(fid1,'%s\n','Masses');
fprintf(fid1,'%s\n','');
 C_mass = 12;
% B_mass = 10.811 ;
% N_mass = 14.0067;
fprintf(fid1,'\t\t%s %12.8f\n','1',C_mass);
% fprintf(fid1,'\t\t%s %12.8f\n','2',B_mass);
% fprintf(fid1,'\t\t%s %12.8f\n','3',N_mass);

%fprintf(fid1,'\t\t%s %12.8f\n','2',Cu_mass);

fprintf(fid1,'%s\n','');
fprintf(fid1,'%s\n','Atoms');
fprintf(fid1,'%s\n','');

for i = 1:n_atom;
    
        fprintf(fid1,'%i 1 %12.8f %12.8f %12.8f\n',atom_dat(i,1:4));

    % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   
    
end



fclose(fid1);





