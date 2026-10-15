#import "/typst/notes.typ": *

The angle from $a$ to $b$ is $"atan2"(a times b, a dot b)$, and
$u^perp = (-u_y, u_x)$ is $u$ turned a quarter counter-clockwise.
For $u times v != 0$, $w = x u + y v$ with
$x = (w times v)\/(u times v)$ and $y = (u times w)\/(u times v)$
(Cramer): segment parameters, or which shifted copy of a line
passes through a point.

=== Polygons
#fact[Centroid][Of the area: $G = 1\/(6A) sum_i (p_i + p_(i + 1))
  (p_i times p_(i + 1))$ with the signed $A$. The mean of the
  vertices is a different point.]
#fact[Pick][A lattice polygon has area $A = I + B\/2 - 1$ with $I$
  lattice points inside and $B = sum_i gcd(|Delta x_i|, |Delta y_i|)$
  on the boundary. A lattice segment holds
  $gcd(|Delta x|, |Delta y|) + 1$ lattice points.]
#fact[Quadrilateral][A simple quadrilateral $A B C D$, convex or
  not, has signed area $(C - A) times (D - B)\/2$. The largest one
  on four given points is the best of the three ways to pair them
  into diagonals. With $A, B, C, D$ from four sets, maximise
  $u times w$ over $u in C - A$, $w in D - B$ (Minkowski sums): for
  a vertex $u$ the best $w$ is extreme in direction $u^perp$, and
  it turns monotonically with $u$.]
#fact[Convex][A closed polygon is strictly convex iff every turn
  is a strict left turn and its edge directions go around once
  (a pentagram goes around twice). Three or more lengths are the
  sides of a non-degenerate polygon iff the longest is shorter
  than the sum of the others; then a convex one in a circle
  exists.]
#fact[Erdős-Szekeres][Of 5 points, no three collinear, some 4 are
  in convex position. A convex $k$-gon is forced by 5, 9, 17
  points for $k = 4, 5, 6$, and by $binom(2k - 4, k - 2) + 1$ in
  general.]

=== Triangles and circles
Triangle $A B C$ with opposite sides $a, b, c$, angle $alpha$ at
$A$, $s = (a + b + c)\/2$, $u = B - A$, $v = C - A$
and area $S = |u times v|\/2$:
- Incenter $(a A + b B + c C)\/(a + b + c)$, inradius $S\/s$.
- Circumradius $R = a b c\/(4S)$, circumcenter
  $O = A + (|v|^2 u^perp - |u|^2 v^perp)\/(2 u times v)$,
  orthocenter $H = A + B + C - 2O$. With the centroid $G$:
  $H - G = 2(G - O)$ (Euler line; one point iff equilateral).
- Median to $a$: $sqrt(2 b^2 + 2 c^2 - a^2)\/2$; bisector of $alpha$:
  $2 b c cos(alpha\/2)\/(b + c)$.
- Circular segment between a chord and its arc of central angle
  $theta$: $r^2 (theta - sin theta)\/2$.
- Lens of two circles at distance $d$, if
  $|r_1 - r_2| < d < r_1 + r_2$: area
  $r_1^2 alpha_1 + r_2^2 alpha_2 - d r_1 sin alpha_1$ with
  $alpha_1 = arccos((d^2 + r_1^2 - r_2^2)\/(2 d r_1))$, half the
  angle the lens spans at center 1 ($alpha_2$ likewise).
#fact[Area by arcs][Twice the area is $integral.cont (x dif y -
  y dif x)$ counter-clockwise along the boundary: a segment
  $p -> q$ adds $p times q$, an arc of circle $(c, r)$ from angle
  $alpha$ to $beta > alpha$ adds $r^2 (beta - alpha) +
  r (c_x (sin beta - sin alpha) - c_y (cos beta - cos alpha))$.
  Union of circles: the arcs of each circle inside no other
  (remove duplicate circles first).]
#fact[Around a point][Of $n$ points with no two on a line through
  $O$, $binom(n, 3) - sum_i binom(c_i, 2)$ triangles contain $O$, where
  $c_i$ counts the points less than a half-turn counter-clockwise
  of $p_i$ (angular sort, two pointers).]

=== 3D
- Tetrahedron volume $|(b - a) dot ((c - a) times (d - a))|\/6$;
  the triple product (a determinant) is positive iff $a b c$
  looks counter-clockwise from $d$.
- Plane through $a, b, c$: normal $n = (b - a) times (c - a)$,
  signed distance of $p$: $n dot (p - a)\/|n|$. Line $a + t d$:
  distance of $p$ is $|d times (p - a)|\/|d|$; between two lines
  $|(a_2 - a_1) dot (d_1 times d_2)|\/|d_1 times d_2|$ (not parallel).
- Spherical cap of height $h$: $V = pi h^2 (3r - h)\/3$,
  $S = 2 pi r h$.
- Great-circle distance: $r "atan2"(|u times v|, u dot v)$ for unit
  vectors $u, v$; latitude $phi$ and longitude $lambda$ give
  $u = (cos phi cos lambda, cos phi sin lambda, sin phi)$.

=== Precision and overflow
- Stay in integers: turns, sides, containment and distance
  comparisons (squared) are exact signs. Divide or take roots
  only for the printed answer.
- Coordinates up to $C$: a cross or dot product of differences
  reaches $8 C^2$, $8 dot 10^18$ for $C = 10^9$. #cpp[That is just
  inside `ll`, and a sum of two overflows: use `__int128`.]#py[In
  PyPy, values beyond $2^63$ silently become slow big integers.]
  An exact line intersection is a degree-3 term over a degree-2
  term. A 3D triple product reaches $48 C^3$, below $2^63$ up to
  $C approx 5 dot 10^5$.
- A double keeps about 16 digits, so its error is about $10^(-16)$
  times the largest intermediate value, input rounding included.
  Pick `eps` per quantity, well above that: with coordinates near
  $10^4$, a cross product (near $10^8$) needs about $10^(-6)$,
  while a distance (near $10^4$) can use $10^(-9)$. Never compare
  computed floats with `==`.
- `acos` and `asin` near $plus.minus 1$ lose half the digits; take
  angles as `atan2(cross, dot)`. Clamp arguments to $[-1, 1]$ and
  use `sqrt(max(0, x))`.
- Sort by angle with a half-plane test plus a cross product
  (#cpp[`angle_less`]#py[`angle_cmp`]), not `atan2`, once
  coordinates exceed about $10^7$. Decide where zero vectors and equal directions go;
  subtract the center first to sort around a point.

=== Counting regions
Euler's formula (_Graph theory facts_) counts regions: every
crossing and every end of a segment or arc is a vertex, every
piece of a curve between them an edge, and $k$ lines share one
more vertex at infinity. With closed curves and full lines only,
simple crossings only (no tangency, no three curves through a
point) and a connected drawing, $I$ crossing points give
$I + k + 2 - [k > 0]$ regions: $n$ lines in general position make
$1 + n + binom(n, 2)$, $n$ pairwise crossing circles $n^2 - n + 2$.
A line crosses a circle or a convex polygon at most twice, so two
circles cross at most twice and a convex $p$-gon crosses a circle
or a convex polygon at most $2p$ times (6 for a triangle).

=== Rotating calipers
On a convex polygon in counter-clockwise order, when the answer
for edge $i$ is a vertex $j(i)$ that only moves forward as $i$
does, walk $j$ with a pointer: $O(n)$ in total, no binary search.
Uses: the diameter (`hull_diameter`), the width (for each edge
its farthest vertex), the smallest enclosing rectangle by area or
perimeter (one side on an edge, four pointers), the largest
triangle on the hull, the distance and the common tangents of two
convex polygons, and for each vertex the range of vertices beyond
a line through it.

=== Reflections
#fact[Unfolding][A path from $p$ that reflects off the lines
  $ell_1, ..., ell_k$ in this order and ends at $q$ unfolds into
  the segment from $p$ to $sigma_1(sigma_2(... sigma_k (q)))$,
  where $sigma_i$ reflects in $ell_i$. Check that the segment
  crosses the unfolded mirrors $sigma_1 ... sigma_(i - 1)(ell_i)$
  in this order, each within the mirror's extent.]

#full[
=== Cues
- A beam or ball bouncing off walls in a given order: unfold.
- "Most regions": Euler's formula with every pair crossing the
  most times at distinct points.
- Area seen from a point: angular sweep over all corners.
- Largest quadrilateral, area from four groups: the cross product
  of the diagonals, Minkowski sums.
- Sum or expectation of hull areas over subsets: per directed
  edge, count the subsets for which it is a hull edge.
- Groups with pairwise disjoint hulls: enumerate separating lines
  through two input points.
- "No convex quadrilateral" or $k$-gon: Erdős-Szekeres keeps the
  answer tiny except for collinear points.
- Order points so every turn is left: always possible (spiral).

=== More facts
#fact[Hull edges][Points in general position: $p_i -> p_j$ is a
  counter-clockwise edge of the hull of a subset $S$ with
  $|S| >= 3$ iff $p_i, p_j in S$ and no point of $S$ lies right of
  it. So twice the hull area summed over subsets is
  $sum_(i != j) (p_i times p_j) N_(i j)$, with $N_(i j) =
  2^(l_(i j)) - 1$ for $l_(i j)$ points left of $p_i -> p_j$. For a
  convex polygon, chord $i -> i + d$ is an edge iff none of the
  $d - 1$ vertices between is chosen.]
#fact[Separation][Two point sets, no three points collinear, have
  disjoint hulls iff some line through a point of each has the
  other points of each set strictly on its own side. A split into
  groups with pairwise disjoint hulls is fixed by one such line
  per pair of groups: $O(n^2)$ choices each.]
#fact[Around a disk][Points at distances $R_1, R_2 >= r$ from the
  center of a disk of radius $r$, an angle $Delta in [0, pi]$
  apart, $alpha_i = arccos(r\/R_i)$. The shortest path avoiding the
  open disk is straight if $Delta <= alpha_1 + alpha_2$, else
  tangent, arc, tangent:
  $sqrt(R_1^2 - r^2) + sqrt(R_2^2 - r^2) + r(Delta - alpha_1 -
  alpha_2)$. With $R_2 = r$ it is convex in $Delta$.]

=== More formulas
- Heron with floats loses digits on needle-like triangles; with
  $a >= b >= c$ use $S = 1/4 sqrt((a + (b + c))(c - (a - b))
  (c + (a - b))(a + (b - c)))$, brackets as written (Kahan).
- Excenter opposite $A$: $(-a A + b B + c C)\/(-a + b + c)$, radius
  $S\/(s - a)$.
- Power of $P$: $|P O|^2 - r^2$, the product of the signed
  distances from $P$ to the two points where any line through $P$
  meets the circle, and the squared tangent length.
- Polyhedron with triangle faces $a b c$ counter-clockwise seen
  from outside: $V = 1/6 sum a dot (b times c)$. Planar polygon in
  3D: area $|sum_i p_i times p_(i + 1)|\/2$.
- Rotation about the unit axis $k$ by $theta$ (Rodrigues):
  $v cos theta + (k times v) sin theta + k (k dot v)(1 - cos theta)$.
- Manhattan distance is the Chebyshev distance of
  $(x + y, x - y)$: $|Delta x| + |Delta y| =
  max(|Delta x + Delta y|, |Delta x - Delta y|)$. Rotating by
  $45 degree$ turns diamonds into squares, which split per axis.
- Inversion in the circle $(O, r)$: $P' = O + r^2 (P - O)\/|P - O|^2$.
  Circles through $O$ become lines not through $O$ and back;
  other circles stay circles. It turns tangent-circle problems
  into line problems.
- Spherical triangle with angles $alpha, beta, gamma$ on a sphere
  of radius $r$: area $r^2 (alpha + beta + gamma - pi)$.
- Ellipse with foci at distance $2c$, points with
  $|P F_1| + |P F_2| <= 2a$: semi-axes $a$ and $b = sqrt(a^2 - c^2)$,
  area $pi a b$; the shortest chord through a focus is
  perpendicular to the major axis, of length $2b^2\/a$.
- Extreme vertex of a convex polygon in direction $d$ in
  $O(log n)$: binary search for the first edge $e$ with
  $d dot e <= 0$ on the lower and on the upper chain (both by
  increasing $x$) and take the better one. Tangents from an
  outside point are the ends of the run of edges it sees.
- Duality: the line $y = k x + m$ is the point $(k, m)$, and its
  value at $x$ is $(k, m) dot (x, 1)$; the upper envelope of lines
  is the upper hull of their points (convex hull trick). The line
  $x\/p + y\/q = 1$ is the point $(1\/p, 1\/q)$: $(b, c)$ lies beyond
  it iff $(b, c) dot (1\/p, 1\/q) > 1$, an extreme dot product.
- Lattice points inside a counter-clockwise polygon, column by
  column: on $x = X$ each interior interval runs up from a
  crossing of a right-going edge to one of a left-going edge. If
  no vertex has an integer $x$ and no edge passes a lattice point,
  the count is $sum floor(y)$ over all crossings of left-going
  edges with integer $x = X$ minus the same over right-going
  ones.
- Region seen from $G$ (star-shaped): area
  $1/2 integral rho(theta)^2 dif theta$. Between consecutive event
  angles (corners, switches between boundary curves) each piece
  is a triangle at $G$, or one plus a circular segment.
- Green: $integral.double_D (Q_x - P_y) dif A =
  integral.cont P dif x + Q dif y$, counter-clockwise. A radial
  density $g(|p|)$ over a region:
  $integral.cont F(|p|^2) (x dif y - y dif x)$ with
  $F(s) = 1\/(2s) integral_0^s g(sqrt(sigma)) dif sigma$; over a disk this
  is a smooth periodic integral in the angle (trapezoid rule).
- Midpoint polygon $Q_i = (P_i + P_(i + 1))\/2$:
  $P_(i + 1) = 2 Q_i - P_i$. For odd $n$, $P_0 = sum_i (-1)^i Q_i$
  is unique; for even $n$ a preimage exists iff
  $sum_i (-1)^i Q_i = 0$, and then $P_0$ is free.
- Any points can be ordered with only left or straight turns:
  start at the lowest point heading in $+x$, and repeatedly move
  to the remaining point with the smallest counter-clockwise turn
  from the current direction, the nearer one on ties.
- Disjoint curves in a disk joining given pairs of distinct points
  exist iff no two pairs with both ends on the boundary alternate
  along it; pairs with an end inside never block.
- Shortest paths among polygonal obstacles bend only at obstacle
  corners: Dijkstra on the visibility graph. With a drift $u$ (a
  current, or the frame of an obstacle moving at $-u$) and own
  speed $v > |u|$, a move $Delta$ takes the least $t >= 0$ with
  $|Delta - u t| <= v t$; this time is convex and positively
  homogeneous in $Delta$, so the same holds.
- Rectilinear Steiner tree of terminals on the boundary of a
  rectangle: some optimal tree is the boundary minus its longest
  gap between terminals, or one full line across plus boundary
  pieces and stubs, or several parallel full lines, or an L from
  one side to the adjacent one that cuts off a corner (Cohoon,
  Richards, Salowe); each is a linear DP along the boundary.
- Pairs closer than $d$ in $k$ dimensions: hash points into
  half-open cells of side $s$. Points closer than $d$ lie in cells
  whose indices differ by at most $ceil(d\/s)$ per axis, and two
  points of one cell are closer than $s sqrt(k)$. With $s = d$ scan
  $3^k$ cells; with $s = d\/sqrt(k)$ points pairwise $>= d$ apart
  sit in distinct cells, so when close pairs are forbidden (or a
  conflict graph has no triangle) each cell holds $O(1)$ points;
  scan $(2 ceil(sqrt(k)) + 1)^k$ cells.
- Surface words: a word of $2 n$ letters using each label twice,
  as $x$ or $x^(-1)$, glues a $2 n$-gon into a closed surface. It
  is orientable iff every label appears once with each sign, and
  $chi = V - n + 1$, with $V$ the corners left after identifying
  the tails and the heads of the two sides of each label (DSU).
  Closed surfaces are equal iff orientability and $chi$ agree
  (sphere 2, torus 0, projective plane 1, Klein bottle 0).
- Pitfalls: input given clockwise (check the sign of the area),
  repeated points, all points collinear, polygons of 1 or 2
  points, open versus closed boundaries, and printing `-0.000`
  (print 0 for values that round to zero).
]
