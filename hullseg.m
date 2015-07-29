DEBUG = false;

%hullseg.m
%By:    Sepehr Sadighpour
%==========================================================================
% Generate stochastic data with variable degree of self-similarity, as defined by the hurst exponent.
% 
%value.  In each string, draw the convex hull then find and save:
%       1. the number of points lying on the hull
%       2. The number of colinear (ie lying on the same line) pairs
%       3. The number of ANTIPARALLEL (we actually want anti_parallel, such that 
%          the set of these and the colinear pairs are disjoint) pairs
%       4. The number of segments involved in the MAXIMUM perimeter segment
%               collection, (ie for each segment, find all other colinear
%               and ANTIPARALLEL segments, calculate the percentage of the hull
%               that these make up, and keep the set with the largest %)
%       5. The percentage of the hull surface composed of the MAXIMUM
%               perimeter set.

fprintf('Beginning Hullseg...\n')

clear;clc;close all;

HURST_INCREMENT = 0.01; % increment for the Hurst exponent, or degree of self-similarity in data. %
POPULATION_SIZE = 1000; % number of streams of length SERIES_LENGTH for each Hurst exponent %
COLINEARITY_THRESHOLD_DEGREE = 0.5; % hull segments are defined as colinear if the angle between them is smaller than this number in degrees %
ANTIPARALLEL_THRESHOLD_DEGREE = 180-COLINEARITY_THRESHOLD_DEGREE;

% Create array of possible hurst values %

MIN_H = HURST_INCREMENT;
MAX_H = 1-HURST_INCREMENT;
H = MIN_H:HURST_INCREMENT:MAX_H;

NUM_H_VALUES = length(H);

fprintf('hullseg: Range of Hurst values defined as %f to %f\n',MIN_H,MAX_H)
fprintf('hullseg: Number of increments is %d\n',NUM_H_VALUES)

% Begin generating simulations of variable lengths %
for SERIES_LENGTH = [ 100, 500, 1000, 5000, 10000]
    fprintf('\n\nhullseg: Beginning SIMULATION for Series length %d\n\n',SERIES_LENGTH)
    
    %Pre-allocating variables and arrays
    fprintf('hullseg: Pre-allocating arrays...\n')
    SIMULATION = zeros(POPULATION_SIZE, SERIES_LENGTH, NUM_H_VALUES);
    HULL_INDICES = cell(POPULATION_SIZE,NUM_H_VALUES);
    NUM_HULL_POINTS = zeros(POPULATION_SIZE,NUM_H_VALUES);
    HULL_COORDS = cell(POPULATION_SIZE,NUM_H_VALUES);
    HULL_LENGTHS = zeros(POPULATION_SIZE,NUM_H_VALUES);
    HULL_VECTOR_LENGTHS = cell(POPULATION_SIZE,NUM_H_VALUES);
    HULL_VECTORS = cell(POPULATION_SIZE,NUM_H_VALUES);
    HULL_VECTOR_NORMS = cell(POPULATION_SIZE,NUM_H_VALUES);
    DOT_PRODUCTS = cell(POPULATION_SIZE,NUM_H_VALUES);
    DEGREES = cell(POPULATION_SIZE,NUM_H_VALUES);
    ANTIPARALLELS = cell(POPULATION_SIZE,NUM_H_VALUES);
    COUNT_ANTIPARALLELS = zeros(POPULATION_SIZE,NUM_H_VALUES);
    INDEX_UNIQUE_ANTIPARALLEL = cell(POPULATION_SIZE,NUM_H_VALUES);
    INDEX_UNIQUE_COLINEAR = cell(POPULATION_SIZE,NUM_H_VALUES);
    ANTIPARALLEL_LENGTHS = cell(POPULATION_SIZE,NUM_H_VALUES);
    SUM_ANTIPARALLEL_LENGTHS = zeros(POPULATION_SIZE,NUM_H_VALUES);
    COLINEARS = cell(POPULATION_SIZE,NUM_H_VALUES);
    COUNT_COLINEARS = zeros(POPULATION_SIZE,NUM_H_VALUES);
    COLINEAR_LENGTHS = cell(POPULATION_SIZE,NUM_H_VALUES);
    SUM_COLINEAR_LENGTHS = zeros(POPULATION_SIZE,NUM_H_VALUES);
    MAXIMUM_SUBSET_PERCENT = zeros(POPULATION_SIZE,NUM_H_VALUES);
    MAXIMUM_SUBSET_COUNT = zeros(POPULATION_SIZE,NUM_H_VALUES);
    %END PREALLOCATION%
    
    fprintf('hullseg: Done preallocating.\nhullseg: Beginning loop over Hurst increments...\n')
    
    for ih = 1:NUM_H_VALUES; 
        
        fprintf('hullseg: ih = %d\n',ih)
        
        for n = 1:POPULATION_SIZE; %FOR EACH SAMPLE%
            if(DEBUG)
                fprintf('\n\nhullseg: sample = %d of %d\n',n,POPULATION_SIZE)
            end
            %Run fractional brownian motion simulation
            SIMULATION (n,:,ih) =  wfbm(H(ih),SERIES_LENGTH);
            
            %Scale the x-coordinate sensibly
            X_INCREMENT = SERIES_LENGTH.^H(ih)/SERIES_LENGTH;
            XCOORD_SCALED = (0:SERIES_LENGTH-1).*X_INCREMENT;
            
            %Get the coordinates of the hull points and plot them
            HULL_INDICES (n,ih) = {convhull(XCOORD_SCALED,SIMULATION(n,:,ih))}; %returns the index of the elements in SIMULATION falling on the hull%
            NUM_HULL_POINTS (n,ih) = length(HULL_INDICES{n,ih})-1; %first hull point is double counted, so substract 1%
            HULL_COORDS (n,ih) = {[XCOORD_SCALED(HULL_INDICES{n,ih}) ; SIMULATION(n,HULL_INDICES{n,ih},ih)]'};
            
            if(DEBUG)
                fprintf('hullseg: hull coordinates = \n');
                disp(HULL_COORDS{n,ih})
                plot(HULL_COORDS{n,ih}(:,1),HULL_COORDS{n,ih}(:,2),'-ro', XCOORD_SCALED,SIMULATION(n,:,ih), '-b') %SHOW US THE HULL POINTS
                title(['Hull surface for fractional brownian motion-hurst parameter ',num2str(H(ih)),' Series Length ',num2str(SERIES_LENGTH)])
                pause
            end
            
            %Calculate length of Hull vectors
            HULL_VECTORS (n,ih) = {diff(HULL_COORDS{n,ih})}; %CALCULATE ROW DIFFERENCES BETWEEN SUCCESSIVE HULL POINTS
            HULL_VECTOR_LENGTHS(n,ih) = {sqrt(sum((diff(HULL_COORDS{n,ih},1,1).^2), 2))};%LENGTH OF HULL VECTORS
            HULL_LENGTHS (n,ih) = sum(HULL_VECTOR_LENGTHS{n,ih}); %TOTAL LENGTH OF CONVEX HULL
            
            if(DEBUG)
                fprintf('hullseg: hull vector length = %1.3f\n', HULL_VECTOR_LENGTHS{n,ih});% SHOW US HULL VECTOR LENGTHS
                pause
            end
            
            %Calculate angles between Hull vectors
            for i = 1:length(HULL_VECTORS {n,ih});
                HULL_VECTOR_NORMS {n,ih}(i,:) = HULL_VECTORS{n,ih}(i,:) ./ HULL_VECTOR_LENGTHS{n,ih}(i);
            end
            DOT_PRODUCTS(n,ih) = {HULL_VECTOR_NORMS{n,ih} * HULL_VECTOR_NORMS{n,ih}'};
            DEGREES(n,ih) = {real(acos(DOT_PRODUCTS{n,ih}) * 180/pi)}; %ANGLE BETWEEN HULL VECTORS IN DEGREES%
            
            if(DEBUG)
                %check max and min degree, setting diagonal elements to 90
                %degrees
                fprintf('\nhullseg: maximum degree = %f and minimum degree = %f\n',max(max(DEGREES{n,ih})),min(min(DEGREES{n,ih}+diag(90*ones(NUM_HULL_POINTS(n,ih),1),0))))
            end
            
            %{
            %COUNT COLINEAR%
            COLINEARS(n,ih)={DEGREES{n,ih}<COLINEARITY_THRESHOLD_DEGREE}; %1 FOR COLINEAR, 0 FOR OTHER%
            COLINEARS{n,ih} = COLINEARS{n,ih} - diag(diag(COLINEARS{n,ih})); %SET DIAGONAL ELEMENTS TO 0%
            [rowC,colC] = find(COLINEARS{n,ih});
            INDEX_UNIQUE_COLINEAR (n,ih) = {unique(rowC)};
            COUNT_COLINEARS (n,ih) = length(INDEX_UNIQUE_COLINEAR{n,ih}); %NUMBER OF COLINEAR PAIRS OF VECTORS
            COLINEAR_LENGTHS(n,ih) = {HULL_VECTOR_LENGTHS{n,ih}(INDEX_UNIQUE_COLINEAR{n,ih})};
            SUM_COLINEAR_LENGTHS(n,ih) = sum(COLINEAR_LENGTHS{n,ih});
            
            if(DEBUG)
                fprintf('\n\nhullseg: Found %d pairs of colinear hull vectors\n', COUNT_COLINEARS(n,ih))
                fprintf('hullseg: Giving lengths...\n')
                fprintf('hullseg:')
                for i = 1:length(COLINEAR_LENGTHS{n,ih})
                    fprintf('\t%2.3f', COLINEAR_LENGTHS{n,ih}(i))
                end
                fprintf('\n')
                fprintf('hullseg: Sum of COLINEAR lengths = %2.3f\n',SUM_COLINEAR_LENGTHS(n,ih))
            end
            
            
            
            %COUNT ANTIPARALLELS%
            ANTIPARALLELS(n,ih)={DEGREES{n,ih}>ANTIPARALLEL_THRESHOLD_DEGREE}; %1 FOR ANTI-PARALLLEL, 0 FOR OTHER%
            ANTIPARALLELS{n,ih} = ANTIPARALLELS{n,ih} - diag(diag(ANTIPARALLELS{n,ih})); %SET DIAGONAL ELEMENTS TO 0%
            [rowP,colP] = find(ANTIPARALLELS{n,ih});
            INDEX_UNIQUE_ANTIPARALLEL (n,ih)= {unique(rowP)};
            COUNT_ANTIPARALLELS(n,ih) = length(INDEX_UNIQUE_ANTIPARALLEL{n,ih}); %NUMBER OF ANTI-ANTIPARALLEL PAIRS OF VECTORS
            ANTIPARALLEL_LENGTHS(n,ih) = {HULL_VECTOR_LENGTHS{n,ih}(INDEX_UNIQUE_ANTIPARALLEL{n,ih})};
            SUM_ANTIPARALLEL_LENGTHS(n,ih) = sum(ANTIPARALLEL_LENGTHS{n,ih});
            
            if(DEBUG)
                fprintf('\n\nhullseg: Found %d pairs of ANTI-PARALLEL hull vectors\n', COUNT_ANTIPARALLELS(n,ih))
                fprintf('hullseg: Giving lengths...\n')
                fprintf('hullseg:')
                for i = 1:length(ANTIPARALLEL_LENGTHS{n,ih})
                    fprintf('\t%f', ANTIPARALLEL_LENGTHS{n,ih}(i))
                end
                fprintf('\n')
                fprintf('hullseg: Sum of ANTI-PARALLEL lengths = %f\n',SUM_ANTIPARALLEL_LENGTHS(n,ih))
                fprintf('\n')
            end
            %}
            
            %Find the MAXIMUM subset of Colinear and ANTIPARALLEL vectors.
            for hull_vec = 1:NUM_HULL_POINTS(n,ih);
                COLINEAR_OR_ANTIPARALLEL = DEGREES{n,ih}(:,hull_vec)<COLINEARITY_THRESHOLD_DEGREE | DEGREES{n,ih}(:,hull_vec)>ANTIPARALLEL_THRESHOLD_DEGREE; %For hullvector i, find set of antiparallel and colinear vectors (includes itself by default)
                HULL_SURFACE_PERCENT = (COLINEAR_OR_ANTIPARALLEL'* HULL_VECTOR_LENGTHS{n,ih})/HULL_LENGTHS(n,ih); %Find the percent of the hull surface occupied by this set
                if (sum(COLINEAR_OR_ANTIPARALLEL) >=2 && HULL_SURFACE_PERCENT>MAXIMUM_SUBSET_PERCENT(n,ih))% If the set contains a pair and more of the hullsurface than any previous set, update the values
                    MAXIMUM_SUBSET_COUNT(n,ih) = sum(COLINEAR_OR_ANTIPARALLEL);
                    MAXIMUM_SUBSET_PERCENT(n,ih)= HULL_SURFACE_PERCENT;
                end
            end
            
            
            if(DEBUG)
                fprintf('hullseg: Found a maximum length set of %d vectors which are all antiparallel or colinear\n', MAXIMUM_SUBSET_COUNT(n,ih));
                fprintf('hullseg: Subset occupies %f x 100 percent of the hull surface\n',MAXIMUM_SUBSET_PERCENT(n,ih));
                pause
            end
            
        end %POPULATION SIZE
    end %HURST Values
    
    %COUNT_TOTAL = COUNT_ANTIPARALLELS + COUNT_COLINEARS;
    %PERCENT_ANTIPARALLEL = SUM_ANTIPARALLEL_LENGTHS ./ HULL_LENGTHS;
    %PERCENT_COLINEAR = SUM_COLINEAR_LENGTHS ./ HULL_LENGTHS;
    
    %OUTPUT TO CSV %
    
    %dlmwrite(['L' num2str(SERIES_LENGTH),'_COUNT_ANTIPARALLELS.csv'],COUNT_ANTIPARALLELS,'newline','pc')
    %dlmwrite(['L' num2str(SERIES_LENGTH),'_COUNT_COLINEARS.csv'],COUNT_COLINEARS,'newline','pc')
    %dlmwrite(['L' num2str(SERIES_LENGTH),'_COUNT_TOTAL.csv'],COUNT_TOTAL,'newline','pc')
   
    %dlmwrite(['L' num2str(SERIES_LENGTH),'_PERCENT_ANTIPARALLELS.csv'],PERCENT_ANTIPARALLEL,'newline','pc')
    %dlmwrite(['L' num2str(SERIES_LENGTH),'_PERCENT_COLINEARS.csv'],PERCENT_COLINEAR,'newline','pc')
    dlmwrite(['L' num2str(SERIES_LENGTH),'_MAXIMUM_SUBSET_COUNT.csv'],MAXIMUM_SUBSET_COUNT,'newline','pc')
    dlmwrite(['L' num2str(SERIES_LENGTH),'_MAXIMUM_SUBSET_PERCENT.csv'],MAXIMUM_SUBSET_PERCENT,'newline','pc')
    dlmwrite(['L' num2str(SERIES_LENGTH),'_NUM_HULL_POINTS.csv'],NUM_HULL_POINTS,'newline','pc')
    
end %SERIES_LENGTH


% OUTPUT A SAMPLE FBM SERIES ACROSS ALL H %
%{
A=SIMULATION(1,:,:);
A=reshape(A,SERIES_LENGTH,NUM_H_VALUES);
dlmwrite('L10000_Hinc25_FBM_SIM.csv',A,'newline','pc')
%}
