-- import Data.Vector as V
import Data.List as L
import Data.Map  as M

type ValueGrid = [Float]
type PolicyGrid = [Float]
type Consumption = Float
type Util = Float

crra :: Float -> Consumption -> Util
crra sigma c = (c ** (1-sigma)) / (1-sigma)

u :: Consumption -> Util
u = crra 3

r = 0.05
beta = 1/1.1
pf = ((1+r)*)
grid = [0.01, 0.02 .. 3]
gamma k = takeWhile (<= (pf k)) grid

dist :: ValueGrid -> ValueGrid -> Float
dist v v' = maximum . L.map abs $ zipWith (-) v v'

vfun :: ValueGrid -> Float -> Float
vfun v = interp grid v -- no need for interpolation
                       -- can only be on grid points

t :: ValueGrid -> ValueGrid
t v = [maximum $ L.map (\k' -> (f k k') + beta * (vfun v k')) (gamma k) | k <- grid]
      where f k k' = u $ (pf k) - k

vList :: ValueGrid -> [ValueGrid]
vList v0 = iterate t v0

g :: ValueGrid -> PolicyGrid
g _ = []

interp :: [Float] -> [Float] -> Float -> Float
interp x0 y0 x = mapInterp (M.fromAscList al) x
                 where al = zip x0 y0

mapInterp :: Map Float Float -> Float -> Float
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

main = print $ v112
       where v112 = take 1 . drop 111 $ vList grid





