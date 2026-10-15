#import "/typst/notes.typ": *

$a dot b = a_x b_x + a_y b_y = |a| |b| cos theta$ and
$a times b = a_x b_y - a_y b_x = |a| |b| sin theta$, where $theta$
turns $a$ into $b$ counter-clockwise; the angle itself is
$"atan2"(a times b, a dot b)$. $u^perp = (-u_y, u_x)$ is $u$ turned
a quarter counter-clockwise.

=== Polygons
#fact[Shoelace][$2A = sum_i p_i times p_(i + 1)$ (indices mod $n$),
  positive for counter-clockwise order. With floats, subtract $p_0$
  from every point first.]
#fact[Centroid][Of the area: $G = 1\/(6A) sum_i (p_i + p_(i + 1))
  (p_i times p_(i + 1))$ with the signed $A$. The mean of the
  vertices is a different point.]
#fact[Pick][A lattice polygon has area $A = I + B\/2 - 1$ with $I$
  lattice points inside and $B = sum_i gcd(|Delta x_i|, |Delta y_i|)$
  on the boundary. A lattice segment holds
  $gcd(|Delta x|, |Delta y|) + 1$ lattice points.]

=== Triangles and circles
Triangle $A B C$ with opposite sides $a, b, c$, angles
$alpha, beta, gamma$, $s = (a + b + c)\/2$, $u = B - A$, $v = C - A$
and area $S = |u times v|\/2$:
- Centroid $(A + B + C)\/3$. Incenter
  $(a A + b B + c C)\/(a + b + c)$, inradius $S\/s$.
- Circumradius $R = a b c\/(4S)$, circumcenter
  $O = A + (|v|^2 u^perp - |u|^2 v^perp)\/(2 u times v)$,
  orthocenter $H = A + B + C - 2O$.
- Heron $S = sqrt(s (s - a)(s - b)(s - c))$. Cosines
  $a^2 = b^2 + c^2 - 2 b c cos alpha$. Sines
  $a\/sin alpha = b\/sin beta = c\/sin gamma = 2R$. Median to $a$:
  $sqrt(2 b^2 + 2 c^2 - a^2)\/2$; bisector of $alpha$:
  $2 b c cos(alpha\/2)\/(b + c)$.
- Circle of radius $r$, central angle $theta$: chord
  $2r sin(theta\/2)$, arc $r theta$, sector $r^2 theta\/2$, segment
  (between chord and arc) $r^2 (theta - sin theta)\/2$.
- Lens of two circles at distance $d$, if
  $|r_1 - r_2| < d < r_1 + r_2$: area
  $r_1^2 alpha_1 + r_2^2 alpha_2 - d r_1 sin alpha_1$ with
  $alpha_1 = arccos((d^2 + r_1^2 - r_2^2)\/(2 d r_1))$, half the
  angle the lens spans at center 1 ($alpha_2$ likewise).

=== 3D
- $a times b = (a_y b_z - a_z b_y, a_z b_x - a_x b_z,
  a_x b_y - a_y b_x)$; $|a times b|$ is twice the area of the
  triangle $0, a, b$.
- Tetrahedron volume $|(b - a) dot ((c - a) times (d - a))|\/6$;
  the triple product (a determinant) is positive iff $a b c$
  looks counter-clockwise from $d$.
- Plane through $a, b, c$: normal $n = (b - a) times (c - a)$,
  signed distance of $p$: $n dot (p - a)\/|n|$. Line $a + t d$:
  distance of $p$ is $|d times (p - a)|\/|d|$; between two lines
  $|(a_2 - a_1) dot (d_1 times d_2)|\/|d_1 times d_2|$ (not parallel).
- Sphere: $V = 4/3 pi r^3$, $S = 4 pi r^2$; cap of height $h$:
  $V = pi h^2 (3r - h)\/3$, $S = 2 pi r h$.
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

=== Euler's formula
A planar graph with $V$ vertices, $E$ edges, $F$ faces (the outer
one included) and $C$ components has $V - E + F = 1 + C$, so 2 if
connected. To count regions cut out by segments or circles, make
every crossing a vertex and every piece between crossings an
edge. $n$ lines in general position make $1 + n + binom(n, 2)$
regions. A simple planar graph has $E <= 3V - 6$ (for $V >= 3$),
and $E <= 2V - 4$ without triangles.

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

#full[
=== More formulas
- Heron with floats loses digits on needle-like triangles; with
  $a >= b >= c$ use $S = 1/4 sqrt((a + (b + c))(c - (a - b))
  (c + (a - b))(a + (b - c)))$, brackets as written (Kahan).
- Excenter opposite $A$: $(-a A + b B + c C)\/(-a + b + c)$, radius
  $S\/(s - a)$. $O$, the centroid and $H$
  lie on one line (Euler line).
- Power of $P$: $|P O|^2 - r^2$, the product of the signed
  distances from $P$ to the two points where any line through $P$
  meets the circle, and the squared tangent length.
- Polyhedron with triangle faces $a b c$ counter-clockwise seen
  from outside: $V = 1/6 sum a dot (b times c)$.
- Line through $p, q$: $a x + b y = c$ with $a = q_y - p_y$,
  $b = p_x - q_x$, $c = a p_x + b p_y$; a point's distance to it is
  $|a x + b y - c|\/sqrt(a^2 + b^2)$. Two such lines meet at
  $(c_1 b_2 - c_2 b_1, a_1 c_2 - a_2 c_1)\/(a_1 b_2 - a_2 b_1)$.
- Rotation by $theta$: $(x cos theta - y sin theta,
  x sin theta + y cos theta)$. In 3D about the unit axis $k$
  (Rodrigues): $v cos theta + (k times v) sin theta +
  k (k dot v)(1 - cos theta)$.
- Reflection of $p$ in the line through $a$ with unit direction
  $d$: $2(a + d (d dot (p - a))) - p$.
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
- Planar polygon in 3D: area $|sum_i p_i times p_(i + 1)|\/2$.
- Extreme vertex of a convex polygon in direction $d$ in
  $O(log n)$: binary search for the first edge $e$ with
  $d dot e <= 0$ on the lower and on the upper chain (both by
  increasing $x$) and take the better one. Tangents from an
  outside point are the ends of the run of edges it sees.
- Duality: the line $y = k x + m$ is the point $(k, m)$, and its
  value at $x$ is $(k, m) dot (x, 1)$; the upper envelope of lines
  is the upper hull of their points (convex hull trick).
- Pitfalls: input given clockwise (check the sign of the area),
  repeated points, all points collinear, polygons of 1 or 2
  points, open versus closed boundaries, and printing `-0.000`
  (print 0 for values that round to zero).
]
