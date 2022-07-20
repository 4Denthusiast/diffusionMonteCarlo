{-# LANGUAGE TupleSections #-}

module Main where

import Lib

import System.Environment
import System.Random
import Text.Read

main :: IO ()
main = do
        arguments <- getArgs
        case parseArguments arguments defaultExecutionParameters of
            Left err -> putStrLn err
            Right (ExecutionParameters{dimension=d, atomicNumber=z, charge=q}) -> (continue . resetIteration) =<< prepare 6000 =<< initSystem 0.5 1000 (createAtom d z q)
    where prepare 0 s = pure s
          prepare n s = prepare (n-1) =<< (putStrLn (show n) >> stepAndTrace s)
          continue s = continue =<< (stepAndTrace s)

data ExecutionParameters = ExecutionParameters {
    dimension :: Int,
    atomicNumber :: Int,
    charge :: Int
}

defaultExecutionParameters = ExecutionParameters {
    dimension = 3,
    atomicNumber = 1,
    charge = 0
}

parseArguments :: [String] -> ExecutionParameters -> Either String ExecutionParameters
parseArguments [] xp = Right xp
parseArguments (arg:args) xp = case arg of
        "-d" -> requireInteger "-d" args (\n -> xp{dimension = n})
        "-z" -> requireInteger "-z" args (\n -> xp{atomicNumber = n})
        "-q" -> requireInteger "-q" args (\n -> xp{charge = n})
        x -> Left ("Unrecognised argument: " ++ show x)
    where requireInteger a [] f = Left ("\""++a++"\" option must be followed by an integer.")
          requireInteger a (ns:args') f = case readMaybe ns of
            Nothing -> Left ("\""++a++"\" option must be followed by an integer.")
            Just n -> parseArguments args' (f n)
