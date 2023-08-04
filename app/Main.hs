{-# LANGUAGE TupleSections #-}

module Main where

import Particle
import Configuration
import Lib

import System.Environment
import System.Random
import Text.Read

main :: IO ()
main = do
        arguments <- getArgs
        case parseArguments arguments defaultExecutionParameters of
            Left err -> putStrLn err
            Right (ExecutionParameters{dimension=d, atomicNumber=z, charge=q, requiredError=e}) -> (continue e . resetIteration) =<< prepare 2000 =<< initSystem 0.05 200 e (createAtom d z q)
    where prepare 0 s = pure s
          prepare n s = prepare (n-1) =<< (putStrLn (show n) >> snd <$> stepAndTrace s)
          continue e s = do
            (e',s') <- stepAndTrace s
            if e' < e then return () else continue e s'

dihydrogen3 = Conf [([-0.6,0,0],Nucleus 1),([0.6,0,0],Nucleus 1),([1,0,0],Electron),([-1,0,0],Electron)]

data ExecutionParameters = ExecutionParameters {
    dimension :: Int,
    atomicNumber :: Int,
    charge :: Int,
    requiredError :: Double
}

defaultExecutionParameters = ExecutionParameters {
    dimension = 3,
    atomicNumber = 1,
    charge = 0,
    requiredError = 0.01
}

parseArguments :: [String] -> ExecutionParameters -> Either String ExecutionParameters
parseArguments [] xp = Right xp
parseArguments (arg:args) xp = case arg of
        "-d" -> requireNumber True "-d" args (\n -> xp{dimension = n})
        "-z" -> requireNumber True "-z" args (\n -> xp{atomicNumber = n})
        "-q" -> requireNumber True "-q" args (\n -> xp{charge = n})
        "-r" -> requireNumber False "-r" args (\x -> xp{requiredError = x})
        x -> Left ("Unrecognised argument: " ++ show x)
    where requireNumber int a [] f = Left ("\""++a++"\" option must be followed by " ++ if int then "an integer." else "a number.")
          requireNumber int a (ns:args') f = case readMaybe ns of
            Nothing -> Left ("\""++a++"\" option must be followed by " ++ if int then "an integer." else "a number.")
            Just n -> parseArguments args' (f n)
