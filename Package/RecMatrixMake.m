function [TestRectTopNMatrix] = RecMatrixMake(Test, RecIdxMat, topN)

%index
TestIndex = (1:size(Test,1));
TestIndex = repmat(TestIndex,1,topN);
TestIndex = TestIndex(:);


%%%% Rec Top N %%%%

RecScrSortIndexTopN = RecIdxMat;
RecScrSortIndexTopNForCat = reshape(RecScrSortIndexTopN,numel(RecScrSortIndexTopN),1);


%%%% Repeat Tesing Matrix %%%%

TestRep = repmat(Test,topN,1);

RowIndex = (1:size(Test,1)*topN).';

%%%% get sub index and value from nnz for TestRep
I = find(TestRep~=0);
[X, Y] = ind2sub(size(TestRep), I);


RowIndexSparse      = [X; RowIndex];
ColumnIndexSparse   = [Y; RecScrSortIndexTopNForCat];
BinSparse           = ones(size(RowIndexSparse,1),1);

TestRectTopNMatrix      = sparse(RowIndexSparse,ColumnIndexSparse,BinSparse  , size(TestRep,1),1210);

ReOrderIndex = reshape( (1:size(Test,1)*topN), size(Test,1),topN);
ReOrderIndex = ReOrderIndex.';
ReOrderIndex = ReOrderIndex(:);

%%%% reoreder before output
TestRectTopNMatrix = TestRectTopNMatrix(ReOrderIndex,:);

end

