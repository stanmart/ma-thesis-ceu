import Data.List as L
import Data.Map  as M

type ValueGrid = [Double]
type PolicyGrid = [Double]
type Consumption = Double
type Util = Double

-- parameters

u :: Consumption -> Util
u = crra 2

r      = 0.05
beta   = 1/1.1
mu     = 100
tsigma = 100
trho   = 0

y = [50, 75, 100, 125, 150]
yProb = [0.000088417285201, 0.105561356381654, 0.788700452666289, 0.105561356381655, 0.000088417285201]
yProc = zip y yProb

-- grid

mMin = 20
mMax = 3*mu
numM = 30
mGrid = linspace mMin mMax numM
numC = 300

-- choice set

gamma k = drop 1 $ linspace 0 (2*mu) (numC+1)
pf = ((1+r)*)

-- distance in L inf norm

dist :: ValueGrid -> ValueGrid -> Double
dist v v' = maximum . L.map abs $ zipWith (\x y -> (x-y)/y) v v'

-- VFI

vfun :: ValueGrid -> Double -> Double
vfun v = interp mGrid v 

t :: ValueGrid -> ValueGrid
t v = [maximum $ L.map (\c -> (u c) + beta * (evnext k c)) (cGrid k) | k <- mGrid]
      where evnext k c = sum $ [yProb * (vfun v $ (pf $ k-c) + y) | (y, yProb) <- yProc]
            cGrid  k = drop 1 $ linspace 0 (pf k) (numC+1)
      

vList :: ValueGrid -> [ValueGrid]
vList v0 = iterate t v0

vListDef :: [ValueGrid]
vListDef = vList $ L.map u mGrid

findv :: Double -> [ValueGrid] -> ValueGrid
findv eps (v1 : v2 : vs) = if   dist v1 v2 < eps
                           then v2
                           else findv eps (v2 : vs)

-- policy function

g :: ValueGrid -> PolicyGrid
g v = [maximumC (cGrid k) $ L.map (\c -> u c + beta * evnext k c) (cGrid k) | k <- mGrid]
      where evnext k c = sum [yProb * vfun v (pf (k - c) + y) | (y, yProb) <- yProc]
            cGrid  k = drop 1 $ linspace 0 (pf k) (numC+1)
            maximumC cGrid l = snd $ maximum $ zip l cGrid

-- helper functions

linspace :: Double -> Double -> Double -> [Double]
linspace a b n = [a, a+(b-a)/(n-1) .. b]

crra :: Double -> Consumption -> Util
crra 1     c = log c
crra sigma c = (c ** (1-sigma)) / (1-sigma)

mapInterp :: Map Double Double -> Double -> Double
mapInterp m x = 
    case (M.lookupLE x m, M.lookupGE x m) of
      (Just (a, av), Just (b, bv)) -> 
        if a == b
        then av
        else interpolate (a, av) (b, bv) x
      (Nothing, Just (b, bv)) -> bv
      (Just (a, av), Nothing) -> av
      _ -> error "mapInterp"
    where interpolate (a, av) (b, bv) x = av + (x-a) * (bv-av) / (b-a)
    
interp :: [Double] -> [Double] -> Double -> Double
interp x0 y0 = mapInterp M.fromAscList al
    where al = zip x0 y0

main = do 
    eRaw <- getLine
    let v = findv (read eRaw :: Double) vListDef
    let c = g v
    print v
    putStrLn ""
    print c





