{-# LANGUAGE TupleSections #-}

module Main where

import Lib

import System.Random

main :: IO ()
main = (continue . resetIteration) =<< prepare 2000 =<< initSystem 0.01 400 (hydride 3)
    where prepare 0 s = pure s
          prepare n s = prepare (n-1) =<< stepAndTrace s
          continue s = continue =<< (stepAndTrace s)
