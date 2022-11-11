import Control.Applicative ((<$>), (<*>))
import qualified Data.Map as M
import Data.Maybe (fromMaybe)


-- Type definitions

type LiquidAsset = Double
type FixedAsset = Double

type Assets = (LiquidAsset, FixedAsset)
type Income = Double
type Consumption = Double

type Utils = Double
type MU = Double

type AssetGrid = [Assets]
type ValueGrid = [Utils]

data State = H | L


-- Model parameters

ph = [0.925, 0.075]
pl = [0.5  , 0.5  ]
ye = [1.0  , 0.1  ]

beta    = 0.7
betahat = 0.8
delta   = 0.99

gamma = 3

u   = crra     gamma
u'  = crra'    gamma
iu  = invcrra  gamma
iu' = invcrra' gamma

rl = 1.00
rf = 1.03

p   af af' = (af-af') ** 2
p'  af af' = (af-af') * 2
ip' af' q  = 1/q + af'

albar = 0
afbar = 0

aGrid = [(al, af) | al <- linspace 0 10 30,
                    af <- linspace 0 10 30]

-- Main functions

gammaset :: Income -> Assets -> Assets -> Bool
gammaset y (al, af) (al', af') = and [
    al' > albar, -- borrowing constraint
    af' > afbar, -- borrowing constraint
    c > 0,       -- Inada-conditions
    c <= al + y  -- liquidity constraint
    ] where c = h y (al, af) (al', af') -- might need nat. borr. constr.!!!! - actually I don't think so

wnext :: (Assets -> Utils) -> (Assets -> (MU, MU)) -> Income -> Assets -> Assets -> Income -> Utils
wnext wfun wfun' y a@(al, af) a'@(al', af') y' =
    u c + betahat*delta/beta * (wfun .- (((1-beta)*) . ubar y a . iubar' y a . wfun') $ (al'+y', af'+y'))
    where c = h y a a'

ewnext :: (ValueGrid, ValueGrid) -> State -> Assets -> Assets -> Utils
ewnext (vl, vh) H a a' = (wnext (findW aGrid vl) (findW' aGrid vl) yh a a' yl) * (ph !! 1) + (wnext (findW aGrid vh) (findW' aGrid vh) yh a a' yh) * (ph !! 0)
ewnext (vl, vh) L a a' = (wnext (findW aGrid vl) (findW' aGrid vl) yl a a' yl) * (pl !! 1) + (wnext (findW aGrid vh) (findW' aGrid vh) yl a a' yh) * (pl !! 0)

t :: State -> (ValueGrid, ValueGrid) -> ValueGrid
t H v = [maximum $ map (ewnext v H a) (filter (gammaset yh a) aGrid) | a <- aGrid]
t L v = [maximum $ map (ewnext v L a) (filter (gammaset yl a) aGrid) | a <- aGrid]

t2 :: (ValueGrid, ValueGrid) -> (ValueGrid, ValueGrid)
t2 v = (t H v, t L v)


-- Default values and finding the solution

-- vlist :: (ValueGrid, ValueGrid) -> [(ValueGrid, ValueGrid)]
-- vlist = iterate t2

-- vListDef :: [(ValueGrid, ValueGrid)]
-- vListDef = vlist $ zip aGrid aGrid

-- findv :: Double -> [(ValueGrid, ValueGrid)] -> (ValueGrid, ValueGrid)
-- findv (v:v':vs) eps = if dist2 v v' < eps
--                       then v'
--                       else findv eps (v':vs)

-- Helper functions

dist :: ValueGrid -> ValueGrid -> Double
dist v v' = maximum . map abs $ zipWith (\x y -> (x-y)/y) v v'

dist2 :: (ValueGrid, ValueGrid) -> (ValueGrid, ValueGrid) -> Double
dist2 (vh, vl) (vh', vl') = max (dist vh vh') (dist vl vl')

crra :: Double -> Consumption -> Utils
crra 1 c = log c
crra g c = c**(1-g) / (1-g)

crra' :: Double -> Consumption -> MU
crra' g c = c ** (-g)

invcrra :: Double -> Utils -> Consumption
invcrra 1 u = exp u
invcrra g u = ((1-g) * u) ** (1/(1-g))

invcrra' :: Double -> MU -> Consumption
invcrra' g mu = mu ** (-1/g)

(.-) :: (Num b) => (a -> b) -> (a -> b) -> (a -> b)
f .- g = (-) <$> f <*> g -- function minus

h :: Income -> Assets -> Assets -> Consumption
h y (al, af) (al', af') = y + al + af - al'/rl - af'/rf - p af af'

ubar :: Income -> Assets -> Assets -> Utils
ubar y a a' = u $ h y a a'  -- ugly but works

-- iubar' :: Income -> Assets -> MU -> Assets
-- iubar' y (al', af') (m1, m2) =
--     let af = ip' af' (m1-m2) in (al, af) where
--         al = (iu' m1) + al'/rl + af'/rf - y - af

iubar' :: Income -> Assets -> (MU, MU) -> Assets
iubar' y (al', af') (m1, m2) =
    let af = ip' af' (m1-m2) in
        let al = (iu' m1) + al'/rl + af'/rf - y - af in
            (al, af)

dp :: [Double] -> [Double] -> Double
l1 `dp` l2 = sum $ zipWith (*) l1 l2

linspace :: Double -> Double -> Double -> [Double]
linspace a b n = [a, a+(b-a)/(n-1) .. b]
    
findW :: [Assets] -> [Utils] -> Assets -> Utils
findW k0 w0 = mapFind $ M.fromAscList (zip k0 w0)
    where mapFind m a = fromMaybe (error "Cannot find value in Map") (M.lookup a m) 

-- problem: W : R^2 -> R => W' : R^2 -> R^2
mapDiff :: M.Map Double Double -> Double -> Double
mapDiff m a =
    case (M.lookupLT a m, M.lookup a m, M.lookupGT a m) of 
        (Just (al, avl), _, Just (ah, avh)) -> (avh-avl) / (ah-al)
        (Just (al, avl), Just av, Nothing) -> (av-avl) / (a-al)
        (Nothing, Just av, Just (ah, avh)) -> (avh-av) / (ah-a)
        _ -> error "Error in mapDiff"

findW' :: [Assets] -> [Utils] -> Assets -> (MU, MU)
findW' k0 w0 (al, ah) =
    let ml = M.fromAscList $ zip (map fst . filter (\a -> snd a == ah) $ k0) w0
        mh = M.fromAscList $ zip (map snd . filter (\a -> fst a == al) $ k0) w0
    in  (mapDiff ml al, mapDiff mh ah)
        
yh = maximum ye
yl = minimum ye