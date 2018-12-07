function [HitMatrix, TargetHit, PoolHit ] = HitSearch( RecMatrix,Target,PoolAll, TopN )
 
 TargetRep = repmat(Target,TopN,1);
 
 %%%% reorder the output
 ReOrderIndex = reshape( (1:size(Target,1)*TopN), size(Target,1),TopN);
 ReOrderIndex = ReOrderIndex.';
 ReOrderIndex = ReOrderIndex(:);

 TargetRep = TargetRep(ReOrderIndex,:);
 
  % get order
 TargetOrder = sum(TargetRep,2);
 
 
 step_size = 5;
 step = floor(size(PoolAll,1)/step_size );
 interval = [ (0:step_size-1)*step+1 ; (1:step_size-1)*step size(PoolAll,1)];

 HitMatrix =sparse(size(TargetOrder,1),size(PoolAll,1));
 for i = 1 : step_size
     Pool = PoolAll( interval(1,i):interval(2,i),:);
     PoolOrder = sum(Pool,2);
     
     OverlapBeforeRecNum = ((TargetRep * Pool.') == TargetOrder);
     OverlapAfterRecNum =  ((RecMatrix * Pool.') == PoolOrder.');
 
     HitMatrixSub = (  (PoolOrder.' == (TargetOrder +1 ) )  &  ( OverlapBeforeRecNum ) & ( OverlapAfterRecNum )   );
     HitMatrix(:,  interval(1,i):interval(2,i) ) = HitMatrixSub;
 end
 

 % hit index
 TargetHit = sum(HitMatrix,2);
 PoolHit   = sum(HitMatrix,1);
 

end

