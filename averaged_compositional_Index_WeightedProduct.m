% This function calculates the AVERAGED COMPOSITIPNAL INDEX (mjw) values
% for a given protein sequence over a window of size w. It predicts the
% class of each Amino Acid component of this sequence. If the
% averaged Linker Index value of an AA is less than the threshold, it will
% be assigned as 1, and it will be assigned 0 otherwise.
% It also returns the true labels of the the given protien sequence.
% and then the function calculates the prediction accuracy by comparing the
% predicted labels with the true labels.
% 
% DATE 20-09-2012
% Author: Maad Shatnawi


function [compositional_index, compositional_index_output, EVAL] = averaged_compositional_Index_WeightedProduct (proteinSequence, trueLabels, w, threshold, gama, beta)

L= size (proteinSequence,2);    % L is the sequence length
epsilon= 1;
Freq_Data = aacount (char( (proteinSequence)));
F3 = struct2cell(Freq_Data);
F13= cell2mat(F3);
AA = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
Comp_Matrix =  F13 * 100 / sum(F13); 
Comp_Cells = mat2cell(Comp_Matrix,[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ],[1]);
CompIndex = cell2struct(Comp_Cells, AA, 1);

%Z is the linker index values for each amino acid 

Z.A= 0;
Z.R= -0.8690;
Z.N= 0.4818;
Z.D= 0.0800;
Z.C= 0.9280;
Z.Q= -0.1892;
Z.E= 0.2744;
Z.G= 1.5445;
Z.H= 1.0561;
Z.I= 0.3677;
Z.L= 0.8473;
Z.K= 0.1823;
Z.M= 3.5264;
Z.F= 0.6190;
Z.P= 0.0870;
Z.S= -0.0476;
Z.T= -0.4595;
Z.W= 2.6391;
Z.Y= 0.6931;
Z.V= 0.7167;

% to calculate the the AVERAGED COMPOSITIPNAL INDEX values

for j = 1 : L
    
    if (j >= 1) && (j <= ((w - 1)/2))
        windowSize = j+((w - 1)/2);
        WindowSequence = proteinSequence (1 : j+((w-1)/2));
    
    elseif (j > ((w - 1)/2)) && (j <= (L - ((w - 1)/2)))  
        windowSize = w;
        WindowSequence = proteinSequence (j-((w-1)/2) : j+((w-1)/2));
       
    else
        windowSize = L-j+((w-1)/2)+1;
        WindowSequence = proteinSequence (j-((w-1)/2) : L);
    
    end
    
    T2 = aacount (WindowSequence);    
    m(j) =100*(T2.A * (Z.A ^gama) * ((CompIndex.A + epsilon) ^beta) + T2.R * (Z.R ^gama) * ((CompIndex.R + epsilon) ^beta) + T2.N * (Z.N ^gama) * ((CompIndex.N + epsilon) ^beta) + T2.D * (Z.D ^gama) * ((CompIndex.D + epsilon) ^beta) + T2.C * (Z.C ^gama) * ((CompIndex.C + epsilon) ^beta) + T2.Q * (Z.Q ^gama) * ((CompIndex.Q + epsilon) ^beta) + T2.E * (Z.E ^gama) * ((CompIndex.E + epsilon) ^beta) + T2.G * (Z.G ^gama) * ((CompIndex.G + epsilon) ^beta) + T2.H * (Z.H ^gama) * ((CompIndex.H + epsilon) ^beta) + T2.I * (Z.I ^gama) * ((CompIndex.I + epsilon) ^beta) + T2.L * (Z.L ^gama) * ((CompIndex.L + epsilon) ^beta) + T2.K * (Z.K ^gama) * ((CompIndex.K + epsilon) ^beta) + T2.M * (Z.M ^gama) * ((CompIndex.M + epsilon) ^beta) + T2.F * (Z.F ^gama) * ((CompIndex.F + epsilon) ^beta) + T2.P * (Z.P ^gama) * ((CompIndex.P + epsilon) ^beta) + T2.S * (Z.S ^gama) * ((CompIndex.S + epsilon) ^beta) + T2.T * (Z.T ^gama) * ((CompIndex.T + epsilon) ^beta) + T2.W * (Z.W ^gama) * ((CompIndex.W + epsilon) ^beta) + T2.Y * (Z.Y ^gama) * ((CompIndex.Y + epsilon) ^beta) + T2.V * (Z.V ^gama) * ((CompIndex.V + epsilon) ^beta)) / windowSize;
    if m(j) < threshold(1,j)
        compositional_index_output(j) = 1;
    else 
        compositional_index_output(j) = 0;
    end
        
end

compositional_index = m;

EVAL = Evaluate2(trueLabels, compositional_index_output);
