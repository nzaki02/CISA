% This function takes the input protein sequenec, and the SA parameters and
% outputs the prediction sensitivity (recall) and precision, the optimal
% threshold values, and the predicted domain-linker class for each chunk.
% Author: Maad Shatnawi

function [Final_Sensitivity Final_Precision Final_Threshold Final_compositional_index_output] = simAnneal_WeightedProduct (proteinSequence, trueLabels, w, gama, beta)

L= size (proteinSequence,2);
segmentSize =5;                                                               
No_segments = floor(L/segmentSize);     % divide the sequence into chunks.
threshold = ones(1,L);
[compositional_index, compositional_index_output, EVAL] = averaged_compositional_Index_WeightedProduct (proteinSequence, trueLabels, w, threshold, gama, beta);
thresholdTransition = (max(compositional_index)-min(compositional_index))/10;
threshold_new = mean(compositional_index)* ones(1,L);
    
for n = 1:No_segments               
    fprintf(1,'segments = %3.0f\n', n);    
    segment_lower_bound = (n-1) * segmentSize + 1;
    segment_upper_bound = segment_lower_bound + segmentSize-1;
    %Segmented_Threshold = zeros(1,No_segments);        % the initial threshold values can be zeros as well.
    Segmented_Threshold = threshold(1, segment_lower_bound:segment_upper_bound);
    [compositional_index, compositional_index_output, EVAL] = averaged_compositional_Index_WeightedProduct (proteinSequence, trueLabels, w, threshold_new, gama, beta);
    Sensitivity = EVAL(2);
    Precision = EVAL(4);
    T0 = 0.01; 
    T = T0;
    alpha = 0.9;
    iteration_count = 0;
    iteration_count2 = 0;
    deltaSensitivity = 0;
    deltaPrecision = 0;

    while (iteration_count2 <= 30 ) & (not (Precision<=0)) & (not (Precision> 0))
        random_segment= ceil(No_segments * rand);      % Select this line for random segment selection or select the next line for left to right segments
     %      random_segment= n;
        iteration_count2 = iteration_count2 +1;        
        segment_lower_bound = (random_segment-1) * segmentSize + 1;
        segment_upper_bound = segment_lower_bound + segmentSize-1;
        rand1= rand-0.2;
        if rand1>0
            Segmented_Threshold_new  = Segmented_Threshold + thresholdTransition;
        else
            Segmented_Threshold_new  = Segmented_Threshold - thresholdTransition;
        end
        threshold_new(1, segment_lower_bound : segment_upper_bound)= Segmented_Threshold_new;
        [compositional_index_new, compositional_index_output_new, EVAL_new] = averaged_compositional_Index_WeightedProduct (proteinSequence, trueLabels, w, threshold_new, gama, beta);
        Sensitivity_new = EVAL_new(2);
        Precision_new = EVAL_new(4);
        deltaSensitivity = Sensitivity_new - Sensitivity;
        deltaPrecision = Precision_new - Precision;
        Sensitivity = Sensitivity_new;
        Precision = Precision_new;
        EVAL = EVAL_new;
        threshold = threshold_new;
        compositional_index_output = compositional_index_output_new;
        Segmented_Threshold  = Segmented_Threshold_new;
        
    end
    
    while ((iteration_count <= 20 ))
        T = alpha * T;
        iteration_count = iteration_count+1;
        random_segment= ceil(No_segments * rand);   % Select this line for random segment selection or select the next line for left to right segments
     %      random_segment= n;
        segment_lower_bound = (random_segment-1) * segmentSize + 1;
        segment_upper_bound = segment_lower_bound + segmentSize-1;
        rand2= rand-0.5;
        if rand2>0
            Segmented_Threshold_new  = Segmented_Threshold + thresholdTransition; 
        else
            Segmented_Threshold_new  = Segmented_Threshold - thresholdTransition;
        end
        threshold_new(1, segment_lower_bound : segment_upper_bound)= Segmented_Threshold_new;
        [compositional_index_new, compositional_index_output_new, EVAL_new] = averaged_compositional_Index_WeightedProduct (proteinSequence, trueLabels, w, threshold_new, gama, beta);
        Sensitivity_new = EVAL_new(2);
        Precision_new = EVAL_new(4);
        deltaSensitivity = Sensitivity_new - Sensitivity;
        deltaPrecision = Precision_new - Precision;



        if (((deltaSensitivity > 0) & (deltaPrecision >= 0))| ((deltaSensitivity >= 0) & (deltaPrecision > 0)) & ((Segmented_Threshold_new <= max(compositional_index)+ thresholdTransition) & (Segmented_Threshold_new >= min(compositional_index)-thresholdTransition)))
            Sensitivity = Sensitivity_new;
            Precision = Precision_new;
            EVAL = EVAL_new;
            threshold = threshold_new;
            compositional_index_output = compositional_index_output_new;
            Segmented_Threshold  = Segmented_Threshold_new;
       elseif ((Segmented_Threshold_new <= max(compositional_index)+ thresholdTransition) & (Segmented_Threshold_new >= min(compositional_index)-thresholdTransition))
            r = rand;
            if ((r < 1 * (exp(-deltaSensitivity/T)+ exp(-deltaPrecision/T))) &((deltaSensitivity >= 0) & (deltaPrecision >= 0)))
                Sensitivity = Sensitivity_new;
                Precision = Precision_new;
                EVAL = EVAL_new;
                threshold = threshold_new;
                Segmented_Threshold  = Segmented_Threshold_new;
                compositional_index_output = compositional_index_output_new;
            else
                threshold_new = threshold;
                Sensitivity_new = Sensitivity;
                Precision_new = Precision;
                EVAL_new = EVAL;
                Segmented_Threshold_new  = Segmented_Threshold;
            end

        end
    end

    Final_Sensitivity(n) = Sensitivity;
    Final_Precision(n) = Precision;
    Final_Threshold(n,:) = threshold;
    Final_compositional_index_output(n,:) = compositional_index_output;
    
   end

Final_Sensitivity = Sensitivity;
Final_Precision = Precision;
Final_Threshold = threshold;
Final_compositional_index_output = compositional_index_output;

x = 1: size (proteinSequence,2);
plot (x, Final_Threshold,'r', x, trueLabels,'g', x, compositional_index,'b');
