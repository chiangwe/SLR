function [TestRectTopNMatrix,TestRectTopNMatrixValue,TestIndex,RecScrSortIndexTopN] = RecMatrix(Test,Wcoef,seed,topN)

%index
TestIndex = (1:size(Test,1));
TestIndex = repmat(TestIndex,1,topN);
TestIndex = TestIndex(:);

%get Rec score
RecScr_all = Test* Wcoef;

%Keep Rec score that is not in the Combination
RecScr = (Test == 0) .* RecScr_all;

% set random seed
rng( seed );

% set image part for random select drugs that has scr zero
RandomImag = reshape( randperm(size(Test,1)*size(Test,2)) , size(Test,1) , size(Test,2) ) ;

% push drugs in the combination in the last order
RuleOut = ( -ones(size(RecScr)) ) .* (Test == 1);
RecScrForRank = complex( full(RecScr + RuleOut ), RandomImag) ;

% sort the score
[RecScrSortValue,RecScrSortIndex] = sort(RecScrForRank,2,'descend','ComparisonMethod','real');
RecScrSortValue = real(RecScrSortValue);

%%%% Rec Top N %%%%

RecScrSortIndexTopN = RecScrSortIndex(:,1:topN);
RecScrSortIndexTopNForCat = reshape(RecScrSortIndexTopN,numel(RecScrSortIndexTopN),1);

RecScrSortValueTopN = RecScrSortValue(:,1:topN);
RecScrSortValueTopNForCat = reshape(RecScrSortValueTopN,numel(RecScrSortValueTopN),1);

%%%% Repeat Tesing Matrix %%%%

TestRep = repmat(Test,topN,1);

TestRecScr = (Test == 1).*RecScr_all;
TestRecScrRep = repmat(TestRecScr,topN,1);

RowIndex = (1:size(Test,1)*topN).';

%%%% get sub index and value from nnz for TestRep
I = find(TestRep~=0);
[X, Y] = ind2sub(size(TestRep), I);
V = TestRecScrRep(find(TestRep~=0));

RowIndexSparse      = [X; RowIndex];
ColumnIndexSparse   = [Y; RecScrSortIndexTopNForCat];
ValueSparse         = [V; RecScrSortValueTopNForCat];
BinSparse           = ones(size(RowIndexSparse,1),1);

TestRectTopNMatrix      = sparse(RowIndexSparse,ColumnIndexSparse,BinSparse  , size(TestRep,1),1210);
TestRectTopNMatrixValue = sparse(RowIndexSparse,ColumnIndexSparse,ValueSparse, size(TestRep,1),1210);
%TestRectTopNMatrix(1,:)
%TestRectTopNMatrix(20,:)
%%%% reorder the output 
ReOrderIndex = reshape( (1:size(Test,1)*topN), size(Test,1),topN);
ReOrderIndex = ReOrderIndex.';
ReOrderIndex = ReOrderIndex(:);

%%%% reoreder before output
TestRectTopNMatrix = TestRectTopNMatrix(ReOrderIndex,:);
TestRectTopNMatrixValue = TestRectTopNMatrixValue(ReOrderIndex,:);
TestIndex = TestIndex(ReOrderIndex);
end

