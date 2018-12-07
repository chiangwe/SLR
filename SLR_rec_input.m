function [  ] = SLR_rec_input( SLR_Method, InputPath, ModelPath, TrainInputPath, KnowledgePath, OutputPath, Ita, Seed)

% SLR_rec_input( ['SLR_sli_rec'], ['./Dataset/Prescription/A_test/A_test_fold_1.mat'], ['./Dataset/Example/model_example_fold_1.mat'], ['./Dataset/Train/train_fold_1.mat'],['./Dataset/KnowledgePool/A_test/Pool_fold_1.mat'], ['./Dataset/Example/Example_Rec.mat'], 5, 123)

% SLR_rec_input( ['SLR_sl_rec'], ['./Dataset/Prescription/A_test/A_test_fold_1.mat'], ['./Dataset/Example/model_example_fold_1.mat'], ['./Dataset/Train/train_fold_1.mat'],['./Dataset/KnowledgePool/A_test/Pool_fold_1.mat'], ['./Dataset/Example/Example_Rec.mat'], 5, 123)

% SLR_rec_input( ['SLR_si_rec'], ['./Dataset/Prescription/A_test/A_test_fold_1.mat'], ['./Dataset/Example/model_example_fold_1.mat'], ['./Dataset/Train/train_fold_1.mat'],['./Dataset/KnowledgePool/A_test/Pool_fold_1.mat'], ['./Dataset/Example/Example_Rec.mat'], 5, 123)

% SLR_rec_input( ['SLR_s_rec'], ['./Dataset/Prescription/A_test/A_test_fold_1.mat'], ['./Dataset/Example/model_example_fold_1.mat'], ['./Dataset/Train/train_fold_1.mat'],['./Dataset/KnowledgePool/A_test/Pool_fold_1.mat'], ['./Dataset/Example/Example_Rec.mat'], 5, 123)

maxNumCompThreads(1);


	switch SLR_Method
	   case 'SLR_sli_rec'
		run('SLR_sli_rec.m');
	   case 'SLR_sl_rec'
	        run('SLR_sl_rec.m');
	   case 'SLR_si_rec'
		run('SLR_si_rec.m');
	   case 'SLR_s_rec'
		run('SLR_s_rec.m');
	end

end


