% 
%clear all;  clc;
addpath(genpath('./'))
%Seed = 123;
%Ita = 5;
%InputPath = ['./Dataset/Prescription/A_test/A_test_fold_1.mat'];
%TrainInputPath = ['./Dataset/Train/train_fold_1.mat'];
%ModelPath = ['./Dataset/Example/model_example_fold_1.mat'];
%KnowledgePath = ['./Dataset/KnowledgePool/A_test/Pool_fold_1.mat'];
%OutputPath = ['./Dataset/Example/Example_Rec.mat'];

RecTopN = 5;
Thre = Ita;

% get Indication Matrix
Str_out = ['./Dataset/SIDER.mat']
load(Str_out, 'INDI');
INDI_Count = INDI * INDI.';
INDI_Count = INDI_Count - diag( diag( INDI_Count ) );

% Load input
load(InputPath)
load(TrainInputPath)
load(ModelPath)
load(KnowledgePath)
%load(RemoveTrainPath)

[TestRectTopNMatrix, TestRectTopNMatrixValue, TestIndex, RecScrSortIndexTopN ] = RecMatrixSide(TestInput, INDI_Count, Seed, RecTopN);
INDIRec = {RecScrSortIndexTopN};

TRAIN_all = train_data;
inst_num = size(train_data,1);

A_p = train_data(              1 : inst_num/2 , :);
A_n = train_data( inst_num/2 + 1 : inst_num   , :);
Count_Ap = A_p.' * A_p; Count_Ap = Count_Ap - diag(diag(Count_Ap));
Count_An = A_n.' * A_n; Count_An = Count_An - diag(diag(Count_An));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get All possible for Alltest LogR recommandation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drugNum = size(TestInput,2);
Test = []; TestNum = []; TestIdx = [];

for test_case = 1 :size(TestInput,1)
	
	% repeat each combination by drug numbers
	% plug one drug one by one
	IdxNonZeros = find( TestInput( test_case ,:)~=0 );
	drugNumLeft = drugNum-numel(IdxNonZeros);
	[row, column, values] = find( repmat( TestInput(test_case,:) , drugNumLeft,1) );
	
	indexAll = (1:drugNumLeft).';
	indexColumn = setdiff((1:drugNum),IdxNonZeros).';
	
	row    = [row;    indexAll ];
	column = [column; indexColumn ];
	values = [values; ones(drugNumLeft,1)];
	
	% unique same elements
	sparseElement = unique( [row column values] ,'rows' );
	
	% making sparse matrix
	TEST_all_can = sparse(sparseElement(:,1),sparseElement(:,2),sparseElement(:,3));
	
	Test     = [Test    ; {TEST_all_can}];
	TestNum  = [TestNum ; size(TEST_all_can,1) ];
	TestIdx = [TestIdx; {indexColumn}];
end

%
TotalTest = {Test};% TotalTestNum = {TestNum}; 
TotalTestIdx = {TestIdx};

All_test = TotalTest{1};
TotalTestIdx = TotalTestIdx{1};

Indirec = cell2mat(  INDIRec  );
Indirec =  mat2cell( Indirec,   repmat(1,size(Indirec,1),1) , size(Indirec,2)  );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Recommendation posive
[TestRectTopNMatrix, TestRectTopNMatrixValue, TestIndex, RecScrSortIndexTopN] = RecMatrix(TestInput, W_rec_p, Seed, RecTopN);
% remove non
RecScrSortIndexTopN = mat2cell(  RecScrSortIndexTopN,   repmat(1,size(RecScrSortIndexTopN,1),1) , RecTopN  );
TestInput_NonIdx = cellfun( @(vec) find(vec~=0) , mat2cell(  TestInput,   repmat(1,size(TestInput,1),1) , size(TestInput,2)  )   ,'UniformOutput', false);
temp = cellfun( @( Test, Rec) Rec( find( sum( Count_Ap(Test,Rec),1) < Thre ) ) , TestInput_NonIdx, RecScrSortIndexTopN ,'UniformOutput', false);
RecNonZero = cellfun( @( vectOrg, vectRev) setdiff(vectOrg, vectRev ) , RecScrSortIndexTopN, temp ,'UniformOutput', false);
RecNonZeroSize = cellfun( @( vectOrg) size(vectOrg,2 ) , RecNonZero ,'UniformOutput', false);

%%%% get

AllPrd = cellfun( @(Mat) 1./( 1 + exp( -(Mat * x_weight + c_intercept) ) ),  All_test,'UniformOutput', false);
[AllPrd_value, AllPrd_Idx] = cellfun( @(Vec) sort( complex( Vec, rand(size(Vec)) ) ,'descend','ComparisonMethod','real'),  AllPrd,'UniformOutput', false);
LogRec = cellfun( @(Vec,Idx) Vec(Idx), TotalTestIdx, AllPrd_Idx,'UniformOutput', false);
LogRec = cellfun( @(Vec,Vec2) setdiff(Vec,Vec2,'stable'), LogRec, RecNonZero,'UniformOutput', false);
LogRec = cellfun( @(Vec,Num) Vec(1:RecTopN-Num).', LogRec, RecNonZeroSize,'UniformOutput', false);

MethodLogR = cellfun( @(Vec,Vec2) [Vec Vec2] , RecNonZero, LogRec,'UniformOutput', false);
MethodLogR = cell2mat(MethodLogR);

[TestRectTopNMatrix] = RecMatrixMake(TestInput, MethodLogR, RecTopN);

AllRectMatrixPrdTopN_Wp = {TestRectTopNMatrix};


 % recommendation negative
[TestRectTopNMatrix, TestRectTopNMatrixValue, TestIndex, RecScrSortIndexTopN] = RecMatrix(TestInput, W_rec_n,seed, RecTopN);
% remove non
RecScrSortIndexTopN = mat2cell(  RecScrSortIndexTopN,   repmat(1,size(RecScrSortIndexTopN,1),1) , RecTopN  );
TestInput_NonIdx = cellfun( @(vec) find(vec~=0) , mat2cell(  TestInput,   repmat(1,size(TestInput,1),1) , size(TestInput,2)  )   ,'UniformOutput', false);
temp = cellfun( @( Test, Rec) Rec( find( sum( Count_An(Test,Rec),1) < Thre) ) , TestInput_NonIdx, RecScrSortIndexTopN ,'UniformOutput', false);
RecNonZero = cellfun( @( vectOrg, vectRev) setdiff(vectOrg, vectRev ) , RecScrSortIndexTopN, temp ,'UniformOutput', false);
RecNonZeroSize = cellfun( @( vectOrg) size(vectOrg,2 ) , RecNonZero ,'UniformOutput', false);


%%%% get

AllPrd = cellfun( @(Mat) 1./( 1 + exp( -(Mat * x_weight + c_intercept) ) ),  All_test,'UniformOutput', false);
[AllPrd_value, AllPrd_Idx] = cellfun( @(Vec) sort( complex( Vec, rand(size(Vec)) ) ,'ascend','ComparisonMethod','real'),  AllPrd,'UniformOutput', false);

LogRec = cellfun( @(Vec,Idx) Vec(Idx), TotalTestIdx, AllPrd_Idx,'UniformOutput', false);
LogRec = cellfun( @(Vec,Vec2) setdiff(Vec,Vec2,'stable'), LogRec, RecNonZero,'UniformOutput', false);
LogRec = cellfun( @(Vec,Num) Vec(1:RecTopN-Num).', LogRec, RecNonZeroSize,'UniformOutput', false);

MethodLogR = cellfun( @(Vec,Vec2) [Vec Vec2] , RecNonZero, LogRec,'UniformOutput', false);
MethodLogR = cell2mat(MethodLogR);

[TestRectTopNMatrix] = RecMatrixMake(TestInput, MethodLogR, RecTopN);

AllRectMatrixPrdTopN_Wn = {TestRectTopNMatrix};


M_num = size(Pos_pool,1);
N_num = size(Neg_pool,1);

%%%%%% Loading recommendation drug combination

RecMatrixRecTopNWp = cell2mat(AllRectMatrixPrdTopN_Wp);
RecMatrixRecTopNWn = cell2mat(AllRectMatrixPrdTopN_Wn);

[~, TargetHit_M, ~ ] = HitSearch( RecMatrixRecTopNWp,TestInput, [Pos_pool] , RecTopN );
[~, TargetHit_N, ~ ] = HitSearch( RecMatrixRecTopNWp,TestInput, [Neg_pool] , RecTopN );
HitCountWp =[ TargetHit_M TargetHit_N];
%HitCountWp = [HitCountWp ; [TargetHit_M TargetHit_N]];
%
[~, TargetHit_M, ~ ] = HitSearch( RecMatrixRecTopNWn,TestInput, [Pos_pool] , RecTopN );
[~, TargetHit_N, ~ ] = HitSearch( RecMatrixRecTopNWn,TestInput, [Neg_pool] , RecTopN );
HitCountWn =[ TargetHit_M TargetHit_N];
%HitCountWn = [HitCountWn ; [TargetHit_M TargetHit_N]];

save(OutputPath,'AllRectMatrixPrdTopN_Wp','AllRectMatrixPrdTopN_Wn','HitCountWp','HitCountWn');

%%%%%%%%%%%%%%%%%%%%%%%%%% Calcualte Max
%
%PrdTopN = RecTopN;
%TestingDataNum = size(TestInput,1);
%
%[HitMatrix, Pssible_pos_ground_Num, PoolHit] = PossibleHitSearch(TestInput, remove_train_M, RecTopN);
%
%index_larger_topN = (Pssible_pos_ground_Num > RecTopN);
%Pssible_pos_ground_Num(index_larger_topN) = RecTopN;
%Pssible_pos_ground_Num = sum(Pssible_pos_ground_Num);
%Possible_hit_num_pos_ground = Pssible_pos_ground_Num;
%
%%find_index = unique( find( sum(find_index,2)~=0 ) );
%pos_ground_Num = sum(sum(HitMatrix,2),1);
%
%pos_ground = pos_ground_Num ;
%
%% load N dataset in which trainning data is removed
%[HitMatrix, Pssible_neg_ground_Num, PoolHit] = PossibleHitSearch(TestInput,remove_train_N, RecTopN);
%
%index_larger_topN = (Pssible_neg_ground_Num > RecTopN);
%Pssible_neg_ground_Num(index_larger_topN) = RecTopN;
%Pssible_neg_ground_Num = sum(Pssible_neg_ground_Num);
%Possible_hit_num_neg_ground = Pssible_neg_ground_Num;
%
%neg_ground_Num = sum(sum(HitMatrix,2),1);
%
%neg_ground = neg_ground_Num ;
%
%% recall
%RecallMaxPoz = Possible_hit_num_pos_ground./pos_ground; RecallMinPoz = zeros(size(pos_ground,1),1);
%RecallMaxNeg = Possible_hit_num_neg_ground./neg_ground; RecallMinNeg = zeros(size(neg_ground,1),1);
%
%% precision
%PrecisionMaxPoz = Possible_hit_num_pos_ground./(TestingDataNum*PrdTopN); PrecisionMinPoz = zeros(size(pos_ground,1),1);
%PrecisionMaxNeg = Possible_hit_num_neg_ground./(TestingDataNum*PrdTopN); PrecisionMinNeg = zeros(size(neg_ground,1),1);


PrdTopN = RecTopN;
TestDataNumPerfold = size( TestInput , 1) * PrdTopN;

TruePos = full(( (HitCountWp(:,1) == 1) ));
TrueNeg = full(( (HitCountWn(:,2) == 1) ));%PredictPos = (predict_total > threshold_out_total);

TruePosPerFold = mat2cell( TruePos , TestDataNumPerfold, 1);
TrueNegPerFold = mat2cell( TrueNeg , TestDataNumPerfold, 1);

TruePosNumPerFold = cellfun(@sum,TruePosPerFold);
TrueNegNumPerFold = cellfun(@sum,TrueNegPerFold);
[PredictPosNumPerFold,~] = cellfun(@size,TruePosPerFold);

%%%% calculate true positve and negtive

PosGroundNum = pos_ground;
NegGroundNum = neg_ground;

RecP = (TruePosNumPerFold ./ PosGroundNum)./(RecallMaxPoz);
RecN = (TrueNegNumPerFold ./ NegGroundNum)./(RecallMaxNeg);

F1 = (2 .* RecP .* RecN) ./ (RecP + RecN);

string_metrics = [ 'RecallTp: '  num2str( RecP , '%2.4f') '\t' 'RecallTn: ' num2str( RecN , '%2.4f')];
string_metrics = [ string_metrics '\t' 'F1: ' num2str( F1 , '%2.4f') '\n'];
disp(string_metrics)

