{-# LANGUAGE TupleSections #-}

module PriorityQueue (
    PriorityQueue,
    empty,
    size,
    fromList,
    insert,
    pop,
    next,
    toAscList,
    union,
) where

import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map

data PriorityQueue p a = PQ (Map p [a])
-- invariant: none of the lists are empty.

empty :: Ord p => PriorityQueue p a
empty = PQ Map.empty

-- This could be made O(1) by keeping track of the size, but that's probably not worth the effort.
size :: Ord p => PriorityQueue p a -> Int
size (PQ m) = sum $ fmap length m

fromList :: Ord p => [(p,a)] -> PriorityQueue p a
fromList = PQ . Map.fromListWith (++) . map (\(p,x) -> (p,[x]))

insert :: Ord p => p -> a -> PriorityQueue p a -> PriorityQueue p a
insert p x (PQ m) = PQ $ Map.insertWith (++) p [x] m

pop :: Ord p => PriorityQueue p a -> Maybe (a, PriorityQueue p a)
pop (PQ m) = reInsert <$> Map.minViewWithKey m
    where reInsert ((p,[x]),m') = (x, PQ m')
          reInsert ((p,x:xs),m') = (x, PQ $ Map.insert p xs m')

next :: Ord p => PriorityQueue p a -> Maybe p
next (PQ m) = fst <$> Map.lookupMin m

toAscList :: Ord p => PriorityQueue p a -> [(p,a)]
toAscList (PQ m) = Map.foldrWithKey (\p xs xs' -> map (p,) xs ++ xs') [] m

union :: Ord p => PriorityQueue p a -> PriorityQueue p a -> PriorityQueue p a
union (PQ m) (PQ m') = PQ $ Map.unionWith (++) m m'
