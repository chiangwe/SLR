function [HitMatrix, TargetHit, PoolHit ] = PossibleHitSearch(Target,Pool,PlusNum)
 
 % get order
 TargetOrder = sum(Target,2);
 PoolOrder = sum(Pool,2);
 
%OverlapNum = Target * Pool.';
 OverlapNum =((Target * Pool.') == TargetOrder ) ;
 


 % check hit
 % HitMatrix = (  (PoolOrder.' == (TargetOrder + PlusNum) )  &  ( OverlapNum == TargetOrder)   );
 HitMatrix = (  (PoolOrder.' == (TargetOrder + PlusNum) )  &  ( OverlapNum ) );
 
 % hit index
 TargetHit = sum(HitMatrix,2);
 PoolHit   = sum(HitMatrix,1);
 
 %TargetHitIndex = find( sum(HitMatrix,2)~=0);
 %PoolHitIndex = find( sum(HitMatrix,1)~=0);
end

