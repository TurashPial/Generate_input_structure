%%code to generate a simulations system with two types of monomer. 
clear all
tic
%% initialization 
N=49;
f_charged=1/2;  %ratio of charged vs non-charged monomer
f=1-f_charged;
fid = fopen('system.data','r');
box1_x = textscan(fid,'%f %f','Headerlines',14);
box1_y = textscan(fid,'%f %f','Headerlines',1);
box1_z = textscan(fid,'%f %f','Headerlines',1);
fclose(fid);
box_x = cell2mat(box1_x);
box_y = cell2mat(box1_y);
box_z = cell2mat(box1_z);
xlo = box_x(1,1);xhi = box_x(1,2);
ylo = box_y(1,1);yhi = box_y(1,2);
zlo = box_z(1,1);zhi = box_z(1,2);
fid = fopen('system.data','r'); % read various data entry
A = textscan(fid,'%f %f %f %f %f %f %f','Headerlines',36);
B = textscan(fid,'%f %f %f %f','Headerlines',2);
C = textscan(fid,'%f %f %f %f %f','Headerlines',2);
D = textscan(fid,'%f %f %f %f %f %f','Headerlines',2);
E = textscan(fid,'%f %f %f %f %f %f','Headerlines',2);
fclose(fid);
atoms=cell2mat(A);
bonds=cell2mat(B);
angles=cell2mat(C);
dihedrals=cell2mat(D);
impropers=cell2mat(E);
%%

for i=1:length(atoms(:,1)) %change atom idendity
    if atoms(i,3)==7
        atoms(i,3)=12;
    end
    if atoms(i,3)==8
        atoms(i,3)=13;
    end
    if atoms(i,3)==9
        atoms(i,3)=14;
    end
end

s = RandStream('mlfg6331_64'); 
for mol=1:36  %change monomer type
    mol_sel_init=find(atoms(:,2)==mol & atoms(:,3)==4);
    ran_sites=randsample(s,linspace(1,(N-1)/2,(N-1)/2),f*(N-1)/2);
    mol_sel=mol_sel_init(ran_sites);
    for ionize_iter=1:length(mol_sel)
        atm=mol_sel(ionize_iter);
        bsearch=find(bonds(:,3)==atm | bonds(:,4)==atm);
        for b_iter=1:length(bsearch)
            bnd=bsearch(b_iter);
            if atoms(bonds(bnd,3),3)==3
            bdatm=bonds(bnd,3);
            bkcar=bonds(bnd,3);
            atoms(bdatm,3)=7;
            atoms(bdatm,4)=-0.0600;
            elseif (atoms(bonds(bnd,3),3)==5 & rem(atoms(bonds(bnd,3),1),2)==1)
            bdatm=bonds(bnd,3);
            dboxy=bonds(bnd,3);
            atoms(bdatm,3)=9;
            atoms(bdatm,4)=-0.4400;
            elseif (atoms(bonds(bnd,3),3)==5 & rem(atoms(bonds(bnd,3),1),2)==0)
            bdatm=bonds(bnd,3);
            sboxy=bonds(bnd,3);
            sbhyd=length(atoms(:,1))+1;
            atoms(bdatm,3)=10;
            atoms(bdatm,4)=-0.5300;
            atoms=vertcat(atoms,[length(atoms(:,1))+1 mol 11 0.4500...
                atoms(atm,5) atoms(atm,6)-1.25 atoms(atm,7)+1.25]);
            bonds=vertcat(bonds,[length(bonds(:,1))+1 9 atoms(bonds(bnd,3),1) length(atoms(:,1))]);
            elseif atoms(bonds(bnd,3),3)==4
            bdatm=bonds(bnd,3);
            atoms(bdatm,3)=8;
            atoms(bdatm,4)=0.5200;
            end
            
            if atoms(bonds(bnd,4),3)==3
            bdatm=bonds(bnd,4);
            bkcar=bonds(bnd,4);
            atoms(bdatm,3)=7;
            atoms(bdatm,4)=-0.0600;
            elseif (atoms(bonds(bnd,4),3)==5 & rem(atoms(bonds(bnd,4),1),2)==1)
            bdatm=bonds(bnd,4);
            dboxy=bonds(bnd,4);
            atoms(bdatm,3)=9;
            atoms(bdatm,4)=-0.4400;
            elseif (atoms(bonds(bnd,4),3)==5 & rem(atoms(bonds(bnd,4),1),2)==0)
            bdatm=bonds(bnd,4);
            sboxy=bonds(bnd,4);
            sbhyd=length(atoms(:,1))+1;
            atoms(bdatm,3)=10;
            atoms(bdatm,4)=-0.5300;
            atoms=vertcat(atoms,[length(atoms(:,1))+1 mol 11 0.4500...
                atoms(atm,5) atoms(atm,6)-1.25 atoms(atm,7)+1.25]);
            bonds=vertcat(bonds,[length(bonds(:,1))+1 9 atoms(bonds(bnd,4),1) length(atoms(:,1))]);
            elseif atoms(bonds(bnd,4),3)==4
            bdatm=bonds(bnd,4);
            atoms(bdatm,3)=8;
            atoms(bdatm,4)=0.5200;
            end 
        end
        angles=vertcat(angles, [length(angles(:,1))+1 16 atm sboxy sbhyd]);
        dihedrals=vertcat(dihedrals, [length(dihedrals(:,1))+1 15 dboxy atm sboxy sbhyd]);
        dihedrals=vertcat(dihedrals, [length(dihedrals(:,1))+1 16 bkcar atm sboxy sbhyd]);
    end
end

for i=1:length(bonds(:,1)) %modify bonds
    if ((atoms(bonds(i,3),3) == 8 & atoms(bonds(i,4),3) == 10) | (atoms(bonds(i,3),3) == 10 & atoms(bonds(i,4),3) == 8))
        bonds(i,2)=10;
    end
    if ((atoms(bonds(i,3),3) == 8 & atoms(bonds(i,4),3) == 9) | (atoms(bonds(i,3),3) == 9 & atoms(bonds(i,4),3) == 8))
        bonds(i,2)=11;
    end
    if ((atoms(bonds(i,3),3) == 8 & atoms(bonds(i,4),3) == 7) | (atoms(bonds(i,3),3) == 7 & atoms(bonds(i,4),3) == 8))
        bonds(i,2)=12;
    end
    if ((atoms(bonds(i,3),3) == 1 & atoms(bonds(i,4),3) == 7) | (atoms(bonds(i,3),3) == 7 & atoms(bonds(i,4),3) == 1))
        bonds(i,2)=13;
    end
    if ((atoms(bonds(i,3),3) == 6 & atoms(bonds(i,4),3) == 7) | (atoms(bonds(i,3),3) == 7 & atoms(bonds(i,4),3) == 6))
        bonds(i,2)=14;
    end
    if ((atoms(bonds(i,3),3) == 2 & atoms(bonds(i,4),3) == 7) | (atoms(bonds(i,3),3) == 7 & atoms(bonds(i,4),3) == 2))
        bonds(i,2)=15;
    end
end

for i=1:length(angles(:,1)) %modify angles
  if ((atoms(angles(i,3),3) == 9 & atoms(angles(i,4),3) == 8 & atoms(angles(i,5),3) == 10) | ...
          (atoms(angles(i,3),3) == 10 & atoms(angles(i,4),3) == 8 & atoms(angles(i,5),3) == 9))
        angles(i,2)=17;
  end
  if ((atoms(angles(i,3),3) == 9 & atoms(angles(i,4),3) == 8 & atoms(angles(i,5),3) == 7) | ...
          (atoms(angles(i,3),3) == 7 & atoms(angles(i,4),3) == 8 & atoms(angles(i,5),3) == 9))
        angles(i,2)=18;
  end
  if ((atoms(angles(i,3),3) == 7 & atoms(angles(i,4),3) == 8 & atoms(angles(i,5),3) == 10) | ...
          (atoms(angles(i,3),3) == 10 & atoms(angles(i,4),3) == 8 & atoms(angles(i,5),3) == 7))
        angles(i,2)=19;
  end
  if ((atoms(angles(i,3),3) == 8 & atoms(angles(i,4),3) == 7 & atoms(angles(i,5),3) == 2) | ...
          (atoms(angles(i,3),3) == 2 & atoms(angles(i,4),3) == 7 & atoms(angles(i,5),3) == 8))
        angles(i,2)=20;
  end
  if ((atoms(angles(i,3),3) == 8 & atoms(angles(i,4),3) == 7 & atoms(angles(i,5),3) == 1) | ...
          (atoms(angles(i,3),3) == 1 & atoms(angles(i,4),3) == 7 & atoms(angles(i,5),3) == 8))
        angles(i,2)=21;
  end
  if ((atoms(angles(i,3),3) == 8 & atoms(angles(i,4),3) == 7 & atoms(angles(i,5),3) == 6) | ...
          (atoms(angles(i,3),3) == 6 & atoms(angles(i,4),3) == 7 & atoms(angles(i,5),3) == 8))
        angles(i,2)=22;
  end
  if ((atoms(angles(i,3),3) == 7 & atoms(angles(i,4),3) == 1 & atoms(angles(i,5),3) == 3) | ...
          (atoms(angles(i,3),3) == 3 & atoms(angles(i,4),3) == 1 & atoms(angles(i,5),3) == 7))
        angles(i,2)=23;
  end
  if ((atoms(angles(i,3),3) == 7 & atoms(angles(i,4),3) == 1 & atoms(angles(i,5),3) == 7) | ...
          (atoms(angles(i,3),3) == 7 & atoms(angles(i,4),3) == 1 & atoms(angles(i,5),3) == 7))
        angles(i,2)=24;
  end
  if ((atoms(angles(i,3),3) == 1 & atoms(angles(i,4),3) == 7 & atoms(angles(i,5),3) == 1) | ...
          (atoms(angles(i,3),3) == 1 & atoms(angles(i,4),3) == 7 & atoms(angles(i,5),3) == 1))
        angles(i,2)=25;
  end
  if ((atoms(angles(i,3),3) == 1 & atoms(angles(i,4),3) == 7 & atoms(angles(i,5),3) == 6) | ...
          (atoms(angles(i,3),3) == 6 & atoms(angles(i,4),3) == 7 & atoms(angles(i,5),3) == 1))
        angles(i,2)=26;
  end
  if ((atoms(angles(i,3),3) == 7 & atoms(angles(i,4),3) == 1 & atoms(angles(i,5),3) == 2) | ...
          (atoms(angles(i,3),3) == 2 & atoms(angles(i,4),3) == 1 & atoms(angles(i,5),3) == 7))
        angles(i,2)=27;
  end
  if ((atoms(angles(i,3),3) == 7 & atoms(angles(i,4),3) == 6 & atoms(angles(i,5),3) == 2) | ...
          (atoms(angles(i,3),3) == 2 & atoms(angles(i,4),3) == 6 & atoms(angles(i,5),3) == 7))
        angles(i,2)=28;
  end
  if ((atoms(angles(i,3),3) == 2 & atoms(angles(i,4),3) == 7 & atoms(angles(i,5),3) == 1) | ...
          (atoms(angles(i,3),3) == 1 & atoms(angles(i,4),3) == 7 & atoms(angles(i,5),3) == 2))
        angles(i,2)=29;
  end
  if ((atoms(angles(i,3),3) == 2 & atoms(angles(i,4),3) == 7 & atoms(angles(i,5),3) == 6) | ...
          (atoms(angles(i,3),3) == 6 & atoms(angles(i,4),3) == 7 & atoms(angles(i,5),3) == 2))
        angles(i,2)=30;
  end
end

for i=1:length(dihedrals(:,1)) %modify dihedral
  if ((atoms(dihedrals(i,3),3) == 10 & atoms(dihedrals(i,4),3) == 8 & atoms(dihedrals(i,5),3) == 7 ...
            & atoms(dihedrals(i,6),3) == 2) | (atoms(dihedrals(i,3),3) == 2 & atoms(dihedrals(i,4),3) == 7 ...
            & atoms(dihedrals(i,5),3) == 8 & atoms(dihedrals(i,6),3) == 10))
        dihedrals(i,2)=17;
  end
  if ((atoms(dihedrals(i,3),3) == 10 & atoms(dihedrals(i,4),3) == 8 & atoms(dihedrals(i,5),3) == 7 ...
            & atoms(dihedrals(i,6),3) == 1) | (atoms(dihedrals(i,3),3) == 1 & atoms(dihedrals(i,4),3) == 7 ...
            & atoms(dihedrals(i,5),3) == 8 & atoms(dihedrals(i,6),3) == 10))
        dihedrals(i,2)=18;
  end
  if ((atoms(dihedrals(i,3),3) == 10 & atoms(dihedrals(i,4),3) == 8 & atoms(dihedrals(i,5),3) == 7 ...
            & atoms(dihedrals(i,6),3) == 6) | (atoms(dihedrals(i,3),3) == 6 & atoms(dihedrals(i,4),3) == 7 ...
            & atoms(dihedrals(i,5),3) == 8 & atoms(dihedrals(i,6),3) == 10))
        dihedrals(i,2)=19;
  end
  if ((atoms(dihedrals(i,3),3) == 9 & atoms(dihedrals(i,4),3) == 8 & atoms(dihedrals(i,5),3) == 7 ...
            & atoms(dihedrals(i,6),3) == 2) | (atoms(dihedrals(i,3),3) == 2 & atoms(dihedrals(i,4),3) == 7 ...
            & atoms(dihedrals(i,5),3) == 8 & atoms(dihedrals(i,6),3) == 9))
        dihedrals(i,2)=20;
  end
  if ((atoms(dihedrals(i,3),3) == 9 & atoms(dihedrals(i,4),3) == 8 & atoms(dihedrals(i,5),3) == 7 ...
            & atoms(dihedrals(i,6),3) == 1) | (atoms(dihedrals(i,3),3) == 1 & atoms(dihedrals(i,4),3) == 7 ...
            & atoms(dihedrals(i,5),3) == 8 & atoms(dihedrals(i,6),3) == 9))
        dihedrals(i,2)=21;
  end
  if ((atoms(dihedrals(i,3),3) == 9 & atoms(dihedrals(i,4),3) == 8 & atoms(dihedrals(i,5),3) == 7 ...
            & atoms(dihedrals(i,6),3) == 6) | (atoms(dihedrals(i,3),3) == 6 & atoms(dihedrals(i,4),3) == 7 ...
            & atoms(dihedrals(i,5),3) == 8 & atoms(dihedrals(i,6),3) == 9))
        dihedrals(i,2)=22;
  end
  if ((atoms(dihedrals(i,3),3) == 8 & atoms(dihedrals(i,4),3) == 7 & atoms(dihedrals(i,5),3) == 1 ...
            & atoms(dihedrals(i,6),3) == 2) | (atoms(dihedrals(i,3),3) == 2 & atoms(dihedrals(i,4),3) == 1 ...
            & atoms(dihedrals(i,5),3) == 7 & atoms(dihedrals(i,6),3) == 8))
        dihedrals(i,2)=23;
  end
  if ((atoms(dihedrals(i,3),3) == 8 & atoms(dihedrals(i,4),3) == 7 & atoms(dihedrals(i,5),3) == 1 ...
            & atoms(dihedrals(i,6),3) == 3) | (atoms(dihedrals(i,3),3) == 3 & atoms(dihedrals(i,4),3) == 1 ...
            & atoms(dihedrals(i,5),3) == 7 & atoms(dihedrals(i,6),3) == 8))
        dihedrals(i,2)=24;
  end
  if ((atoms(dihedrals(i,3),3) == 8 & atoms(dihedrals(i,4),3) == 7 & atoms(dihedrals(i,5),3) == 1 ...
            & atoms(dihedrals(i,6),3) == 7) | (atoms(dihedrals(i,3),3) == 7 & atoms(dihedrals(i,4),3) == 1 ...
            & atoms(dihedrals(i,5),3) == 7 & atoms(dihedrals(i,6),3) == 8))
        dihedrals(i,2)=25;
  end
  if ((atoms(dihedrals(i,3),3) == 8 & atoms(dihedrals(i,4),3) == 7 & atoms(dihedrals(i,5),3) == 6 ...
            & atoms(dihedrals(i,6),3) == 2) | (atoms(dihedrals(i,3),3) == 2 & atoms(dihedrals(i,4),3) == 6 ...
            & atoms(dihedrals(i,5),3) == 7 & atoms(dihedrals(i,6),3) == 8))
        dihedrals(i,2)=26;
  end
  if ((atoms(dihedrals(i,3),3) == 7 & atoms(dihedrals(i,4),3) == 1 & atoms(dihedrals(i,5),3) == 3 ...
            & atoms(dihedrals(i,6),3) == 2) | (atoms(dihedrals(i,3),3) == 2 & atoms(dihedrals(i,4),3) == 3 ...
            & atoms(dihedrals(i,5),3) == 1 & atoms(dihedrals(i,6),3) == 7))
        dihedrals(i,2)=27;
  end
  if ((atoms(dihedrals(i,3),3) == 7 & atoms(dihedrals(i,4),3) == 1 & atoms(dihedrals(i,5),3) == 3 ...
            & atoms(dihedrals(i,6),3) == 4) | (atoms(dihedrals(i,3),3) == 4 & atoms(dihedrals(i,4),3) == 3 ...
            & atoms(dihedrals(i,5),3) == 1 & atoms(dihedrals(i,6),3) == 7))
        dihedrals(i,2)=28;
  end
  if ((atoms(dihedrals(i,3),3) == 7 & atoms(dihedrals(i,4),3) == 1 & atoms(dihedrals(i,5),3) == 3 ...
            & atoms(dihedrals(i,6),3) == 1) | (atoms(dihedrals(i,3),3) == 1 & atoms(dihedrals(i,4),3) == 3 ...
            & atoms(dihedrals(i,5),3) == 1 & atoms(dihedrals(i,6),3) == 7))
        dihedrals(i,2)=29;
  end
  if ((atoms(dihedrals(i,3),3) == 7 & atoms(dihedrals(i,4),3) == 1 & atoms(dihedrals(i,5),3) == 7 ...
            & atoms(dihedrals(i,6),3) == 2) | (atoms(dihedrals(i,3),3) == 2 & atoms(dihedrals(i,4),3) == 7 ...
            & atoms(dihedrals(i,5),3) == 1 & atoms(dihedrals(i,6),3) == 7))
        dihedrals(i,2)=30;
  end
  if ((atoms(dihedrals(i,3),3) == 7 & atoms(dihedrals(i,4),3) == 1 & atoms(dihedrals(i,5),3) == 7 ...
            & atoms(dihedrals(i,6),3) == 1) | (atoms(dihedrals(i,3),3) == 1 & atoms(dihedrals(i,4),3) == 7 ...
            & atoms(dihedrals(i,5),3) == 1 & atoms(dihedrals(i,6),3) == 7))
        dihedrals(i,2)=31;
  end
  if ((atoms(dihedrals(i,3),3) == 7 & atoms(dihedrals(i,4),3) == 1 & atoms(dihedrals(i,5),3) == 3 ...
            & atoms(dihedrals(i,6),3) == 6) | (atoms(dihedrals(i,3),3) == 6 & atoms(dihedrals(i,4),3) == 3 ...
            & atoms(dihedrals(i,5),3) == 1 & atoms(dihedrals(i,6),3) == 7))
        dihedrals(i,2)=32;
  end
  if ((atoms(dihedrals(i,3),3) == 7 & atoms(dihedrals(i,4),3) == 1 & atoms(dihedrals(i,5),3) == 7 ...
            & atoms(dihedrals(i,6),3) == 6) | (atoms(dihedrals(i,3),3) == 6 & atoms(dihedrals(i,4),3) == 7 ...
            & atoms(dihedrals(i,5),3) == 1 & atoms(dihedrals(i,6),3) == 7))
        dihedrals(i,2)=33;
  end
  if ((atoms(dihedrals(i,3),3) == 2 & atoms(dihedrals(i,4),3) == 7 & atoms(dihedrals(i,5),3) == 1 ...
            & atoms(dihedrals(i,6),3) == 2) | (atoms(dihedrals(i,3),3) == 2 & atoms(dihedrals(i,4),3) == 1 ...
            & atoms(dihedrals(i,5),3) == 7 & atoms(dihedrals(i,6),3) == 2))
        dihedrals(i,2)=34;
  end
  if ((atoms(dihedrals(i,3),3) == 2 & atoms(dihedrals(i,4),3) == 7 & atoms(dihedrals(i,5),3) == 1 ...
            & atoms(dihedrals(i,6),3) == 3) | (atoms(dihedrals(i,3),3) == 3 & atoms(dihedrals(i,4),3) == 1 ...
            & atoms(dihedrals(i,5),3) == 7 & atoms(dihedrals(i,6),3) == 2))
        dihedrals(i,2)=35;
  end
  if ((atoms(dihedrals(i,3),3) == 2 & atoms(dihedrals(i,4),3) == 7 & atoms(dihedrals(i,5),3) == 6 ...
            & atoms(dihedrals(i,6),3) == 2) | (atoms(dihedrals(i,3),3) == 2 & atoms(dihedrals(i,4),3) == 6 ...
            & atoms(dihedrals(i,5),3) == 7 & atoms(dihedrals(i,6),3) == 2))
        dihedrals(i,2)=36;
  end
  if ((atoms(dihedrals(i,3),3) == 1 & atoms(dihedrals(i,4),3) == 7 & atoms(dihedrals(i,5),3) == 1 ...
            & atoms(dihedrals(i,6),3) == 2) | (atoms(dihedrals(i,3),3) == 2 & atoms(dihedrals(i,4),3) == 1 ...
            & atoms(dihedrals(i,5),3) == 7 & atoms(dihedrals(i,6),3) == 1))
        dihedrals(i,2)=37;
  end
  if ((atoms(dihedrals(i,3),3) == 1 & atoms(dihedrals(i,4),3) == 7 & atoms(dihedrals(i,5),3) == 1 ...
            & atoms(dihedrals(i,6),3) == 3) | (atoms(dihedrals(i,3),3) == 3 & atoms(dihedrals(i,4),3) == 1 ...
            & atoms(dihedrals(i,5),3) == 7 & atoms(dihedrals(i,6),3) == 1))
        dihedrals(i,2)=38;
  end
  if ((atoms(dihedrals(i,3),3) == 6 & atoms(dihedrals(i,4),3) == 7 & atoms(dihedrals(i,5),3) == 1 ...
            & atoms(dihedrals(i,6),3) == 3) | (atoms(dihedrals(i,3),3) == 3 & atoms(dihedrals(i,4),3) == 1 ...
            & atoms(dihedrals(i,5),3) == 7 & atoms(dihedrals(i,6),3) == 6))
        dihedrals(i,2)=39;
  end
  if ((atoms(dihedrals(i,3),3) == 6 & atoms(dihedrals(i,4),3) == 7 & atoms(dihedrals(i,5),3) == 1 ...
            & atoms(dihedrals(i,6),3) == 2) | (atoms(dihedrals(i,3),3) == 2 & atoms(dihedrals(i,4),3) == 1 ...
            & atoms(dihedrals(i,5),3) == 7 & atoms(dihedrals(i,6),3) == 6))
        dihedrals(i,2)=40;
  end
  if ((atoms(dihedrals(i,3),3) == 1 & atoms(dihedrals(i,4),3) == 7 & atoms(dihedrals(i,5),3) == 6 ...
            & atoms(dihedrals(i,6),3) == 2) | (atoms(dihedrals(i,3),3) == 2 & atoms(dihedrals(i,4),3) == 6 ...
            & atoms(dihedrals(i,5),3) == 7 & atoms(dihedrals(i,6),3) == 1))
        dihedrals(i,2)=41;
  end
end
for i=1:length(impropers(:,1)) %modify improper
    if ((atoms(impropers(i,3),3) == 8 & atoms(impropers(i,4),3) == 7 & atoms(impropers(i,5),3) == 9 ...
            & atoms(impropers(i,6),3) == 10) | (atoms(impropers(i,3),3) == 8 & atoms(impropers(i,4),3) == 7 ...
            & atoms(impropers(i,5),3) == 10 & atoms(impropers(i,6),3) == 9) | (atoms(impropers(i,3),3) == 8 ...
            & atoms(impropers(i,4),3) == 9 & atoms(impropers(i,5),3) == 7 & atoms(impropers(i,6),3) == 10) | ...
            (atoms(impropers(i,3),3) == 8 & atoms(impropers(i,4),3) == 9 & atoms(impropers(i,5),3) == 10 ...
            & atoms(impropers(i,6),3) == 7) | (atoms(impropers(i,3),3) == 8 & atoms(impropers(i,4),3) == 10 ...
            & atoms(impropers(i,5),3) == 9 & atoms(impropers(i,6),3) == 7) | (atoms(impropers(i,3),3) == 8 ...
            & atoms(impropers(i,4),3) == 10 & atoms(impropers(i,5),3) == 7 & atoms(impropers(i,6),3) == 9))
        impropers(i,2)=2;
    end
end 

%%
masses=[12.011 1.008 12.011 12.011 15.999 12.011 12.011 12.011 15.999 15.999 1.008 15.999 1.008 22.990 15.00794];
masses=vertcat(linspace(1,15,15),masses);
fileID=fopen('system_ionized.data','w'); %write new data file
fprintf(fileID,'LAMMPS Description\n');
fprintf(fileID,'\n%d atoms\n',numel(atoms(:,1)));
fprintf(fileID,'%d bonds\n',numel(bonds(:,1)));
fprintf(fileID,'%d angles\n',numel(angles(:,1)));
fprintf(fileID,'%d dihedrals\n',numel(dihedrals(:,1)));
fprintf(fileID,'%d impropers\n',numel(impropers(:,1)));
fprintf(fileID,'\n%d atom types\n',numel(unique(atoms(:,3)))+1);
fprintf(fileID,'%d bond types\n',numel(unique(bonds(:,2))));
fprintf(fileID,'%d angle types\n',numel(unique(angles(:,2))));
fprintf(fileID,'%d dihedral types\n',numel(unique(dihedrals(:,2))));
fprintf(fileID,'%d improper types\n',numel(unique(impropers(:,2))));
fprintf(fileID,'\n%0.4f %0.4f xlo xhi\n',xlo,xhi);
fprintf(fileID,'%0.4f %0.4f ylo yhi\n',ylo,yhi);
fprintf(fileID,'%0.4f %0.4f zlo zhi\n',-82,400);
fprintf(fileID,'\nMasses\n');
fprintf(fileID,'\n%d %0.5f',masses);
fprintf(fileID,'\n\nAtoms\n');
fprintf(fileID,'\n%d %d %d %0.4f %0.4f %0.4f %0.4f',atoms');
fprintf(fileID,'\n\nBonds\n');
fprintf(fileID,'\n%d %d %d %d',bonds');
fprintf(fileID,'\n\nAngles\n');
fprintf(fileID,'\n%d %d %d %d %d',angles');
fprintf(fileID,'\n\nDihedrals\n');
fprintf(fileID,'\n%d %d %d %d %d %d',dihedrals');
fprintf(fileID,'\n\nImpropers\n');
fprintf(fileID,'\n%d %d %d %d %d %d',impropers');
fprintf(fileID,'\n');

toc