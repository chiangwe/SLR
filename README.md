# SLR : Drug Recommendation toward Safe Polypharmacy
The source code is developed under MATLAB version 2016b.

Contact author: Wen-Hao Chiang, chiangwe@iupui.edu

The implementation is to solve the optimization problem, "joint learning of SLIM and Logistic regression model."

The git repositories contains the following:

1. Dataset: Raw Data extracted from FDA Adverse Event Reporting System (FAERS) and processed datasets for training and testing

2. Package: All the packages that are used to solve the optimization problem. 

3. SLR\_input.m: The function to invoke the model training for SLR method and pass the parameters.

4. SLR.m: The main algorithm.

5. SLR\_rec\_input.m: The function to invoke the recommendation for SLR mdoel and pass the parameters.

6. SLR\_s\_rec.m: The recommendation that only use SLIM component.

7. SLR\_si\_rec.m: The recommendation that use SLIM component and indication component.

8. SLR\_sl\_rec.m: The recommendation that use SLIM component and LogR component. 

9. SLR\_sli\_rec.m: the recommendation that use SLIM, LogR, and indication component. 

The details of the input data format is discussed in the next section.

All the function are provided with a **sample code** that can run directly in MATLAB prompt. 


## Model Training

MATLAB function SLR\_input( InputPath, OutputPath, OMEGA, ALPHA, LAMBDA, BETA, GAMMA) is to invoke the 
training of the model.

-  InputPath: The path of the input file. 

  The input file is a format of MAT-file, which contains a n X m matrix, where n is the number of drug combination and 
  m is the number of drugs. 

- OutputPath: The path of the output file.
 
  The input file is a format of MAT-file, which contains two n X m matrix, where n is the number of drug combination and
  m is the number of drugs, for SLIM components. Also, a m x 1 vector for LogR component.

- OMEGA, ALPHA, LAMBDA, BETA, GAMMA: model parameters. 

  They control, trade-off between SLIM and LogR, L2 for SLIM, L1 for SLIM, L2 for LogR, and L1 for LogR respectively.

- Code example:

  `SLR_input( ['./Dataset/Train/train_fold_1.mat'], ['./Dataset/Example/model_example_fold_1.mat'], 0.01, 10, 0.00000001, 0.000001, 0.0001)`

## Recommendation and Evaluation 

SLR\_rec\_input( SLR\_Method, InputPath, ModelPath, TrainInputPath, KnowledgePath, OutputPath, Ita, Seed) is to invoke the recommendation and the evaluation of the model.

- SLR\_Method: The string to specify the recommendation method.
  - SLR\_s\_rec: Method SLR\_s\_rec.m
  - SLR\_si\_rec: Method SLR\_si\_rec.m
  - SLR\_sl\_rec: Method SLR\_sl\_rec.m
  - SLR\_sli\_rec: Method SLR\_sli\_rec.m

- InputPath: The path of the input file.

  The input file is a format of MAT-file, which contains a n X m matrix, where n is the number of testing prescription and
  m is the number of drugs.

- ModelPath: The path of the trained model. 

  This is the "OutputPath" for training model function SLR\_input.

- TrainInputPath: The path of the training data. 

  This is used for remove the low frequency recommendation is SLIM component, which is the same as "InputPath" in training model function SLR\_input.

- KnowledgePath: The path of knowledge pool

  The knowledge pool file is a format of MAT-file, which contains a n X m matrix, where n is the number of drug combinations in the pool and
  m is the number of drugs.

- OutputPath: The path of the output file.

   The output file are in format of MAT-file, which contains two matrix. One is (n\*N) X m matrix, where n is the number of prescription in the pool, N is the number of recommendation and m is the number of drugs. So, this one would be a matrix of new prescriptions with recommended drug. 

   The second one contains a n X N matrix, which are a matrix for recommended drug. 

Code example:

  `SLR_rec_input( ['SLR_sli_rec'], ['./Dataset/Prescription/A_test/A_test_fold_1.mat'], ['./Dataset/Example/model_example_fold_1.mat'], ['./Dataset/Train/train_fold_1.mat'],['./Dataset/KnowledgePool/A_test/Pool_fold_1.mat'], ['./Dataset/Example/Example_Rec.mat'], 5, 123)`

  `SLR_rec_input( ['SLR_sl_rec'], ['./Dataset/Prescription/A_test/A_test_fold_1.mat'], ['./Dataset/Example/model_example_fold_1.mat'], ['./Dataset/Train/train_fold_1.mat'],['./Dataset/KnowledgePool/A_test/Pool_fold_1.mat'], ['./Dataset/Example/Example_Rec.mat'], 5, 123)`

  `SLR_rec_input( ['SLR_si_rec'], ['./Dataset/Prescription/A_test/A_test_fold_1.mat'], ['./Dataset/Example/model_example_fold_1.mat'], ['./Dataset/Train/train_fold_1.mat'],['./Dataset/KnowledgePool/A_test/Pool_fold_1.mat'], ['./Dataset/Example/Example_Rec.mat'], 5, 123)`

  `SLR_rec_input( ['SLR_s_rec'], ['./Dataset/Prescription/A_test/A_test_fold_1.mat'], ['./Dataset/Example/model_example_fold_1.mat'], ['./Dataset/Train/train_fold_1.mat'],['./Dataset/KnowledgePool/A_test/Pool_fold_1.mat'], ['./Dataset/Example/Example_Rec.mat'], 5, 123)`


