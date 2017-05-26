function run_miso
%miso(datafile,maxeval, surrogate, n_start, init_design, sampling, own_design)

xlow=-10*ones(1,8); %variable lower bounds
xup=10*ones(1,8);  %variable upper bounds
dim = 8; %problem dimension
integ=(1:4); %ind

np_1 = 5;
owndesign_1 = repmat(xlow, np_1, 1) +repmat(xup-xlow, np_1, 1).*rand(np_1, dim);
owndesign_1(:,integ) = round(owndesign_1(:,integ));

np_2 = 9;
owndesign_2 = repmat(xlow, np_2, 1) +repmat(xup-xlow, np_2, 1).*rand(np_2, dim);
owndesign_2(:,integ) = round(owndesign_2(:,integ));

maxeval=[50, 200];
rbft = {'rbf_t','rbf_l', 'rbf_c'};
nstart = [5,10];
sdesign = {'slhd', 'lhs', 'own'};
sampler = {'cp', 'tv', 'cptv', 'rs', 'ms', 'cptvl'};
%test 1

for ii = 2: 3 %maxeval
    if ii == 3
        nev =[];
    else 
        nev =maxeval(ii);
    end
    for jj = 1:4 %surface types
        if jj == 4
            surface = [];
        else
            surface = rbft{jj};
        end
        for kk = 1:3 %number start points
            if kk == 3
                ns = [];
            else
                ns = nstart(kk);
            end
            for ll = 1:4 %initial design
                if ll == 4
                    des = [];
                else
                    des = sdesign{ll};
                end
                
                for mm = 1:7 %sampling strategy
                    if mm == 7
                        sam = [];
                    else
                        sam = sampler{mm};
                    end
                    ii
                    jj
                    kk
                    ll
                    mm
                    if strcmp(des, 'own')
                        od1 = owndesign_1;
                        miso('datainput_convex', nev, surface, ns, des,sam,od1);
                        od2 = owndesign_2;
                        miso('datainput_convex', nev, surface, ns, des,sam,od2);
                    else
                        miso('datainput_convex', nev, surface, ns, des,sam);       
                    end
                    close all
                end %sampling choice
            end %initial design choice
        end %number start point choice
    end %surface choice
end %maxeval choice