{-# LANGUAGE TupleSections #-}
module Graphics (
    createWindow,
    updateGraphics
) where

import Control.Concurrent
import Control.Concurrent.MVar
import Data.List
import Data.Word
import GHC.Arr
import Graphics.UI.Gtk
import Graphics.UI.Gtk.Gdk.GC (GC, gcNew, gcSetValues, newGCValues, foreground)

import Particle
import Configuration
import Ansatz
import Walk

data GraphicsState = GraphicsState {
  gWindow :: Window,
  gCanvas :: DrawingArea,
  gData :: MVar (Array (Int,Int) [Double], Array Int [Double])
}

createWindow :: IO GraphicsState
createWindow = do
    initGUI
    window <- windowNew
    windowSetTitle window "DMC simulation"
    windowSetDefaultSize window 2000 2000
    canvas <- drawingAreaNew
    containerAdd window canvas
    widgetShowAll window
    let (width,height,gLength) = (200,200,70)
    imageVar <- newMVar (listArray ((1,1),(width,height)) (replicate (width*height) [0,0,0]), listArray (1,gLength) (replicate width [0]))
    let gs = GraphicsState window canvas imageVar
    timeoutAdd (displayLoop gs) 100
    forkIO mainGUI
    return gs

updateGraphics :: PopulationState -> GraphicsState -> IO ()
updateGraphics pop GraphicsState{gData = imageVar} = modifyMVar_ imageVar (\(img,graph) -> return (updateImage pop img, updateGraph pop graph))

updateImage :: PopulationState -> Array (Int,Int) [Double] -> Array (Int,Int) [Double]
updateImage pop img = accum (zipWith (+)) img (concatMap walkerPix (walkerSet pop))
    where (_,(w,h)) = bounds img
          walkerPix (Walker (Conf ps) a _) = filter (inBounds . fst) $ map (\(x,t) -> (\d -> map (*(d*a/ansatzValue (Conf ps))) (particleColour t)) <$> particlePix x) ps
          ansatzValue = case pop of PopState{ansatz=a} -> aValue a
          particlePix p = let
                  (x,y) = particlePos p 
                  x' = floor $ fromIntegral w * x/2
                  y' = floor $ fromIntegral h * y/2
              in ((x' + div w 2, y' + div h 2), fromIntegral((x'*2+1)^2+(y'*2+1)^2)**((2-fromIntegral (length p))/2))
          particlePos p = let
                  [x,y] = take 2 (p ++ [0,0])
                  r = sqrt $ sum $ map (^2) p
                  r2d = sqrt $ x*x+y*y
                  s = 0.1
              in if r == 0 then (0,0) else (x*s*r/r2d,y*s*r/r2d)
          inBounds (x,y) = 1 <= x && x <= w && 1 <= y && y <= h
          particleColour Electron = [1,1,0]
          particleColour (Nucleus _ _ _) = [0,0,1]

updateGraph :: PopulationState -> Array Int [Double] -> Array Int [Double]
updateGraph pop graph = accum (zipWith (+)) graph (concatMap walkerPix (walkerSet pop))
    where (_,l) = bounds graph
          walkerPix (Walker (Conf ps) a _) = filter (inBounds . fst) $ concatMap (\(x,t) -> if t == Electron then [(\d -> [a/d/ansatzValue (Conf ps)]) <$> particlePix x] else []) ps
          ansatzValue = case pop of PopState{ansatz=a} -> aValue a
          particlePix p = let
                  x = particlePos p
                  x' = ceiling $ fromIntegral l * x
              in (x', (volume x' - volume (x'-1))/(sum (map (^2) p)**0.75))
          particlePos p = (/3) $ (+1) $ logBase 100 $ sum $ map (^2) p
          d = fromIntegral $ confDimension $ (\(Walker c _ _:_) -> c) $ walkerSet pop
          volume x = 10**(((fromIntegral x / fromIntegral l)*2-1)*d)
          inBounds x = x > 0 && x <= l

displayLoop :: GraphicsState -> IO Bool
displayLoop gs@(GraphicsState window canvas imageVar) = do
    drawWindow <- widgetGetDrawWindow canvas
    (w,h) <- widgetGetSize canvas
    drawWindowBeginPaintRect drawWindow (Rectangle 0 0 w h)
    gc <- gcNew drawWindow
    withMVar imageVar $ drawImage w h (toDrawable drawWindow) gc
    drawWindowEndPaint drawWindow
    widgetShowAll window
    return True

drawImage :: Int -> Int -> Drawable -> GC -> (Array (Int,Int) [Double], Array Int [Double]) -> IO ()
drawImage w h drawable gc (img,graph) = mapM_ (\((x,y),b) -> drawCell x y b) (assocs img) >> setGraphColour >> mapM_ drawGraph (transpose $ map (\(i,ys) -> (i,) <$> ys) $ assocs graph)
    where (_,(iw,ih)) = bounds img
          (_,l) = bounds graph
          brightness = map (1/) $ foldr1 (zipWith (max)) (elems img)
          drawCell x y b = do
              gcSetValues gc (newGCValues{foreground = colour (zipWith (*) b brightness)})
              drawRectangle drawable gc True (div ((x-1)*w) iw) (div ((y-1)*h) ih) (div (x*w) iw - div ((x-1)*w) iw) (div (y*h) ih - div ((y-1)*h) ih)
          setGraphColour = gcSetValues gc (newGCValues{foreground = colour [1,1,1]})
          drawGraph ps = drawLines drawable gc $ map (\(i,y) -> (div (i*w) l, round $ fromIntegral h * (1-y/maximum (map snd ps)))) ps

colour :: [Double] -> Color
colour xs = Color r g b
    where [r,g,b] = map (\x -> round $ (min 1 x) * fromIntegral (maxBound :: Word16)) xs
