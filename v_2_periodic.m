%code written by Rachith Aiyappa at IISc in September 2018 to help with
%Jitesh's PhD project
%email rachithaiyappa96@gmail.com for any queries

%this is to model collective behaviour of N fish without any explcit
%attraction/repulsion functions.
%Two body clocks for pairwise(not triad) and spontaneous interactions 
%An absolute time clock for the purpose of measurements at equal intervals
%of time. %However, v_1 and v_2 contains iformation at unqueal intervals of
%time as extra information pertaining to when exactly the interaction
%occurs has been included.

%Periodic boundary conditions have been included


clear all
close all

v = 5;
a_pair = 3;
b_pair = 3;
a_spont = 0.1;
b_spont = 0.1;
time_int = 0.12;
N = 20;
L = 20;
samp_time = 0:time_int:12;
fish_time = zeros(2*N,1);
phi = [];
phi(1:1:N,1) = pi/2;
x(1:1:N,1) = [1:1:N]';
y(1:1:N,1) = 0;

%generating second column of fish time
fish_time(1:1:N,1) = generate_pair(rand(N,1),a_pair,b_pair);
fish_time(N+1:1:2*N,1) = generate_spont(rand(N,1),a_spont,b_spont);

%time_temp increases in equal intervals and also has time where change
%occurs
time_temp = zeros(2*N,1);
samp_col = 2;
col = 2;
while samp_col < size(samp_time,2)
    [min_fish_time,idx_new] = min(fish_time(:,col-1));
    %min(fish_time(:,col))
    %samp_time(samp_col);
    phi(:,col-1);
    
    %No interaction occurs. Just continues swimming
    while min_fish_time > samp_time(samp_col) && samp_col < size(samp_time,2)
        samp_col;
        min_t(col) = min_fish_time;
        delta_t = samp_time(samp_col)-time_temp(1,col-1);
        x(:,col) = x(:,col-1) + (v*delta_t*cos(phi(:,col-1)));
        y(:,col) = y(:,col-1) + (v*delta_t*sin(phi(:,col-1)));

%         for i=1:1:N
%             if abs(x(i,col))>L
%                 x(i,col) = x(i,col) - L;
%             end
%             if abs(y(i,col))>L
%                 y(i,col) = y(i,col) - L;
%             end
%         end
        %periodic boundary conditions
        for i = 1:1:N
            x(i,col) = mod(x(i,col),L);
            y(i,col) = mod(y(i,col),L);
        end
        
        %updating angle
        phi(:,col) = phi(:,col-1);
        
        %updating times
        time_temp(:,col) = samp_time(samp_col);
        fish_time(:,col) = fish_time(:,col-1);
        
        samp_col = samp_col + 1;
        col = col+1;
    end
    
    %Some Interaction occurs
    while min_fish_time < samp_time(samp_col)
        col;
        delta_t = min_fish_time - time_temp(col-1);
        
        
        %Spontaneous interaction
        if idx_new > N
            min_t(col) = min_fish_time;
            idx = idx_new-N;
            
            delta_t = min_fish_time - time_temp(1,col-1);
            x(:,col) = x(:,col-1) + (v*delta_t*cos(phi(:,col-1)));
            y(:,col) = y(:,col-1) + (v*delta_t*sin(phi(:,col-1)));
            
%             for i=1:1:N
%                 if abs(x(i,col))>L
%                     x(i,col) = x(i,col) - L;
%                 end
%                 if abs(y(i,col))>L
%                     y(i,col) = y(i,col) - L;
%                 end
%             end
            %periodic boundary conditions
            for i = 1:1:N
                x(i,col) = mod(x(i,col),L);
                y(i,col) = mod(y(i,col),L);
            end
            
           %updating by a gaussain angle centered at current angle and SD=pi/4 
            phi(idx,col) = phi(idx,col-1) + ((pi/4)*randn);
            while (abs((pi/4)*randn)) > pi
                phi(idx,col) = phi(idx,col-1) + ((pi/4)*randn);
            end
            %other angles remain the same
            for i = 1:1:N
                if i ~= idx
                    phi(i,col) = phi(i,col-1);
                end
            end
    
            %additional information of fnn. Not required.
            d(1:1:N,col) = sqrt((x(:,col)-x(idx,col)).^2 + (y(:,col)-y(idx,col)).^2);
            [~,idx_d] = sort(d(:,col));
            min_d(col) = d(idx_d(2),col);
            
            %updating times
            time_temp(:,col) = min_fish_time;
            fish_time(idx+N,col) = fish_time(idx+N,col-1) + generate_spont(rand,a_spont,b_spont); 
            for i = 1:1:2*N
                if i ~= idx+N
                    fish_time(i,col) = fish_time(i,col-1);
                end
            end
            col = col + 1;
            [min_fish_time,idx_new] = min(fish_time(:,col-1));
        
        %Pairwise interaction
        else 
            min_t(col) = min_fish_time;
            idx = idx_new; 
            
            delta_t = min_fish_time - time_temp(1,col-1);
            %bringing all fish to the focal fish time
            x(:,col) = x(:,col-1) + (v*delta_t*cos(phi(:,col-1)));
            y(:,col) = y(:,col-1) + (v*delta_t*sin(phi(:,col-1)));
 
%             for i=1:1:N
%                 if abs(x(i,col))>L
%                     x(i,col) = x(i,col) - L;
%                 end
%                 if abs(y(i,col))>L
%                     y(i,col) = y(i,col) - L;
%                 end
%             end
            %periodic boundary conditions
            for i = 1:1:N
                x(i,col) = mod(x(i,col),L);
                y(i,col) = mod(y(i,col),L);
            end
            %identifying the fnn to make the focal fish align
            d(1:1:N,col) = sqrt((x(:,col)-x(idx,col)).^2 + (y(:,col)-y(idx,col)).^2);
            [~,idx_d] = sort(d(:,col));
            min_d(col) = d(idx_d(2),col);
            
            %updating phi's
            for i = 1:1:N
                if i ~= idx
                    phi(i,col) = phi(i,col-1);
                end
            end
            phi(idx,col) = phi(idx_d(2),col);
            
            time_temp(:,col) = min_fish_time;
            fish_time(idx,col) = fish_time(idx,col-1) + generate_pair(rand,a_pair,b_pair);
            for i = 1:1:2*N
                if i ~= idx
                    fish_time(i,col) = fish_time(i,col-1);
                end
            end
            col = col + 1;
            min(fish_time(:,col-1))
            [min_fish_time,idx_new] = min(fish_time(:,col-1));
        end
        
        
    end
end
for i=1:1:N
    plot(x(i,:),y(i,:))
    hold on
end
