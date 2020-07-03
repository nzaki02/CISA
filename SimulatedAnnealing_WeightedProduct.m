% This is the main program that is used to apply Simulated Annealing to the
% protein sequence dataset.
% Author: Maad Shatnawi
% Please select the dataset that you want to work with by commenting one of
% the following two lines. Be sure to import the corresponding class labels csv file
% before running this program.
[Protein_Header, Protein_Seq] = fastaread('DomCut.pep');
% [Protein_Header, Protein_Seq] = fastaread('DSAll.txt');
classLabels = data;
w=25;  % this is the averaging window size
gama = 1;  % this is the exponent of linker index (LI)
beta = -1;  % this is the exponent of the amino acid composition (AAC)

for i = 2: 2 %size(Protein_Seq,2)
    Seq = char(Protein_Seq(:,i));
    Labels= classLabels(1: size(Seq,2), i)';
    fprintf(1,'i= %10.5f\n', i);
    fprintf(1,'SequenceSize= %10.5f\n', size(Seq,2));
    
    [Final_Sensitivity Final_Precision Final_Threshold Final_compositional_index_output] = simAnneal_WeightedProduct (Seq, Labels, w, gama, beta);
    fprintf(1,'Final Recall= %10.5f, Final Precision= %10.5f\n', Final_Sensitivity, Final_Precision);
    Sensitivity (i) = Final_Sensitivity;
    Precision (i) = Final_Precision;
end

SensitivityT=Sensitivity';
PrecisionT = Precision';
Average_Sensitivity = sum(SensitivityT) / size(SensitivityT,1);
Average_Precision = sum(PrecisionT) / size(PrecisionT,1);