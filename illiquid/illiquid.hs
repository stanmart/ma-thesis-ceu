import Data.Matrix (fromList)
import Control.Applicative ((<$>), (<*>)))


-- Type definitions

type LiquidAsset = Float
type FixedAsset = Float

type Assets = (LiquidAsset, FixedAsset)
type Income = Float
type Consumption = Float

type Utils = Float
type MU = Float

type AssetGrid = [Assets]
type ValueGrid = [Utils]

data State = H | L


-- Constants

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


-- Main functions

gammaset :: Income -> Assets -> Assets -> Bool
gammaset y (al, af) (al', af') = and [
    al' > albar, -- borrowing constraint
    af' > afbar, -- borrowing constraint
    c > 0,       -- Inada-conditions
    c =< al + y  -- liquidity constraint
] where c = al + af + y - al'/rl - af'/rf -- need nat. borr. constr.!!!!

wfun :: Assets -> Assets -> Utils

wfun' :: Assets -> Assets -> MU

wnext :: Income -> Assets -> Assets -> Income -> Utils
wnext y a@(al, af) a'@(al', af') y' =
    u c + betahat*delta/beta * 
    (wfun .- (((1-beta)*) . u . iu' . wfun') $ a')
    where c = al + af + y' - al'/rl - af'/rf

ewnext :: State -> Assets -> Assets -> Utils
ewnext H a a' = (map (wnext yh a a') ye) `dp` ph
ewnext L a a' = (map (wnext yl a a') ye) `dp` pl
    where l1 `dp` l2 = sum $ zipWith (*) l1 l2


t :: State -> ValueGrid -> ValueGrid
t H v = [maximum $ map (ewnext H a) (filter (gammaset yh a) aGrid) | a <- aGrid]
t L v = [maximum $ map (ewnext L a) (filter (gammaset yl a) aGrid) | a <- aGrid]

t2 :: (ValueGrid, ValueGrid) -> (ValueGrid, ValueGrid)
t2 (vh, vl) = (t H vh, t L vl)


-- Default values and finding the solution

vlist :: (ValueGrid, ValueGrid) -> [(ValueGrid, ValueGrid)]
vlist v0 = iterate t2 v0

vListDef :: [(ValueGrid, ValueGrid)]
vListDef = vlist $ zip aGrid aGrid

findv :: Float -> [(ValueGrid, ValueGrid)] -> (ValueGrid, ValueGrid)
findv v:v':vs eps = if dist2 v v' < eps
                    then v'
                    else findv eps (v':vs)

-- Helper functions

dist :: ValueGrid -> ValueGrid -> Float
dist v v' = maximum . L.map abs $ zipWith (\x y -> (x-y)/y) v v'

dist2 :: (ValueGrid, ValueGrid) -> (ValueGrid, ValueGrid) -> Float
dist2 (vh, vl) (vh', vl') = max (dist vh vh') (dist vl vl')

crra :: Float -> Consumption -> Utils
crra 1 c = log c
crra g c = c**(1-g) / (1-g)

crra' :: Float -> Consumption -> MU
crra g c = c ** (-g)

invcrra :: Float -> Utils -> Consumption
invcrra 1 u = exp u
invcrra g u = ((1-g) * u) ** (1/(1-g))

invcrra' :: Float -> MU -> Consumption
invcrra' g mu :: mu ** (-1/g)

(.-) :: (Num b) => (a -> b) -> (a -> b) -> (a -> b)
f .- g = (-) <$> f <*> g -- function minus

yh = maximum ye
yl = minimum ye