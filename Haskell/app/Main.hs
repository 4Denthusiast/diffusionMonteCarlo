{-# LANGUAGE TupleSections #-}

module Main where

import Particle
import Configuration
import Lib
import Graphics
import Walk (normalVar)

import Control.Monad.Random
import Data.Char
import Data.List
import Debug.Trace
import System.Environment
import System.Random
import Text.Read

main :: IO ()
main = do
        arguments <- getArgs
        case parseArguments arguments defaultExecutionParameters of
            Left err -> putStrLn err
            Right (ExecutionParameters{dimension=d, atomSeps=as, requiredError=e, timeStep=dt, useGraphics=ug, ansatzName=a, walkerCount=w, prepSteps = p, measurementNames = ms, measurementFade = mf}) -> do
              initialSystem <- initSystem dt w e (createAtoms d as) a ms mf
              stabilised <- prepare p initialSystem
              mWindow <- if ug then Just <$> createWindow else return Nothing
              continue p e mWindow (resetIteration stabilised)
    where prepare 0 s = pure s
          prepare n s = prepare (n-1) =<< (putStrLn (show n) >> snd <$> stepAndTrace s)
          continue n e mWindow s = do
            (e',s') <- stepAndTrace s
            mapM (updateGraphics $ snd s') mWindow
            if e' < e && n < 0 then return () else continue (n-1) e mWindow s'

data AtomSep = Atom Int Int Double Int | SepF Double | SepR Double | SepN Double Double

data ExecutionParameters = ExecutionParameters {
    dimension :: Int,
    atomSeps :: [AtomSep],
    requiredError :: Double,
    timeStep :: Double,
    useGraphics :: Bool,
    ansatzName :: String,
    walkerCount :: Int,
    prepSteps :: Int,
    measurementNames :: String,
    measurementFade :: Double
}

defaultExecutionParameters = ExecutionParameters {
    dimension = 3,
    atomSeps = [],
    requiredError = 0.01,
    timeStep = 0.05,
    useGraphics = False,
    ansatzName = "",
    walkerCount = 200,
    prepSteps = 1000,
    measurementNames = "",
    measurementFade = 0
}

parseArguments :: [String] -> ExecutionParameters -> Either String ExecutionParameters
parseArguments [] xp = Right xp
parseArguments (arg:args) xp = case arg of
        "-d" -> requireNumber True "-d" args (\n -> xp{dimension = n})
        "-r" -> requireNumber False "-r" args (\x -> xp{requiredError = x})
        "-t" -> requireNumber False "-t" args (\x -> xp{timeStep = x})
        "-g" -> parseArguments args xp{useGraphics = True}
        "-w" -> requireNumber True "-w" args (\n -> xp{walkerCount = n})
        "-p" -> requireNumber True "-p" args (\n -> xp{prepSteps = n})
        "-a" -> requireString "-a" args (\s -> xp{ansatzName = s})
        "-m" -> requireString "-m" args (\s -> xp{measurementNames = s})
        "-f" -> requireNumber False "-f" args (\n -> xp{measurementFade = n})
        x -> case parseAsAtomSep x of
            Nothing -> Left ("Unrecognised argument: " ++ show x)
            Just a -> parseArguments args xp{atomSeps = a:atomSeps xp}
    where requireNumber int a [] f = Left ("\""++a++"\" option must be followed by " ++ if int then "an integer." else "a number.")
          requireNumber int a (ns:args') f = case readMaybe ns of
            Nothing -> Left ("\""++a++"\" option must be followed by " ++ if int then "an integer." else "a number.")
            Just n -> parseArguments args' (f n)
          requireString a [] f = Left ("\""++a++"\" option must be followed by a string.")
          requireString a (s:args') f = parseArguments args' (f s)

parseAsAtomSep :: String -> Maybe AtomSep
parseAsAtomSep s = case reads s of
    [(d,"")] -> Just (SepF d)
    [(d,"*")] -> Just (SepR d)
    [(d,'*':s')] -> case reads s' of
        [(d',"")] -> Just (SepN d d')
        _ -> Nothing
    _ -> do
        let (zs,(as,(qs,rest))) = ((span (flip elem "+-") <$>) . span isDigit) <$> span isAlpha s
        m <- case rest of
            "" -> Just (1/0)
            ('m':ms) -> case reads ms of
                [(n,"")] -> Just n
                _ -> Nothing
        z <- case zs of
            "H" -> Just 1
            "He" -> Just 2
            "Li" -> Just 3
            "Be" -> Just 4
            _ -> Nothing
        let a = case as of
                "" -> if z == 1 then 1 else z*2
                _ -> read as
        let q = length (filter (=='+') qs) - length (filter (=='-') qs)
        return $ Atom z a m q

createAtoms :: Int -> [AtomSep] -> Rand StdGen Configuration
createAtoms d [] = createAtoms d [Atom 1 1 (1/0) 0]
createAtoms d as = (Conf . reverse . sortOn snd) <$> createAtoms' (replicate d 0) (if d>=4 then 1 else 0) as
    where createAtoms' x s (SepF sep : Atom z a m q : as) = createAtoms'' (zipWith (+) x (sep:repeat 0)) s z a m q as
          createAtoms' x s (SepR sep : Atom z a m q : as) = do
              x' <- randomOffset sep x
              createAtoms'' x' s z a m q as
          createAtoms' x s (SepN sep sep' : Atom z a m q : as) = do
              x' <- randomOffset sep' (head x + sep : tail x)
              createAtoms'' x' s z a m q as
          createAtoms' x s (Atom z a m q : as) = createAtoms'' (zipWith (+) x (1:repeat 0)) s z a m q as
          createAtoms' x s _ = pure []
          createAtoms'' x s z a m q as = do
              es <- replicateM (z - q) ((,Electron) <$> randomOffset 4 x)
              rest <- createAtoms' x s as
              pure $ (x,Nucleus (fromIntegral z) (fromIntegral a*s) m) : rest ++ es
          randomVector l = replicateM d $ (l*) <$> normalVar
          randomOffset l x = zipWith (+) x <$> randomVector l
