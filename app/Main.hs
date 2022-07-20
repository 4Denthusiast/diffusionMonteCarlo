{-# LANGUAGE TupleSections #-}

module Main where

import Lib

import System.Random

main :: IO ()
main = (continue . resetIteration) =<< prepare 6000 =<< initSystem 0.5 1000 (hydrogen 4)
    where prepare 0 s = pure s
          prepare n s = prepare (n-1) =<< (putStrLn (show n) >> stepAndTrace s)
          continue s = continue =<< (stepAndTrace s)
