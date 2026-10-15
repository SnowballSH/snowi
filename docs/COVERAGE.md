# Coverage of the UCup Season 4-5 training set

Goal: every theorem or niche fact needed by the 474 training problems of
UCup Seasons 4-5 (from the "UCup S4-S5 Training" material: problem
catalog, full solutions, 23 study packets) is stated somewhere in the
packet, and every problem's key idea is reachable from a recognition cue.

## Method

1. Each problem was assigned to one or more areas by its topic tags (problems
   with only generic tags such as greedy or constructive went to the
   constructive area); every problem got at least one area.
2. One agent per area read every assigned solution in full, listed the
   non-basic facts it relies on, checked whether snowi already states each
   one, and otherwise added it to the area's notes after verifying it by
   brute force, exact computation, or against the standard statement.
   Facts that belonged to another area's notes were handed off as exact text
   and applied by that area (or by a mop-up pass) after re-verification.
3. A guide agent mapped all 474 problems (and the 14 problems of S5 Extra
   Stage 1) to rows of the cue table in *Recognition cues*.
4. Three read-only reviewers tried to refute every added claim; their
   confirmed findings were fixed. A visual check of every ICPC page
   followed.

The same pass also removed textbook facts and moved recognition advice out
of the ICPC packet (see CONTRIBUTING.md).

## Numbers

| Area | Problems read | Facts already stated | Facts added | Handed off |
|---|---|---|---|---|
| Algebra and probability | 70 | 82 | 26 | 8 |
| Combinatorics | 82 | 91 | 28 | 19 |
| Constructive and ad hoc | 81 | 28 | 33 | 4 |
| Dynamic programming | 90 | 72 | 24 | 5 |
| Data structures | 91 | 89 | 46 | 13 |
| Games, matroids, randomisation | 56 | 46 | 9 | 11 |
| Graphs and flows | 53 | 52 | 32 | 9 |
| Number theory | 93 | 41 | 28 | 21 |
| Strings and geometry | 58 | 50 | 39 | 3 |
| Trees | 43 | 46 | 12 | 4 |
| **Total** (problems counted per area) | 717 | 597 | 277 | 97 |

Fact counts are (problem, fact) pairs, so a fact used by several problems
counts several times. Handed-off facts were applied by their target area or
by the mop-up (15 applied there, 47 already stated by then, 1 rejected as
too elementary). Every one of the 474 problems appears in an area report,
and every one maps to a cue row.

Review: 58 confirmed findings in the notes (1 high, 17 medium, 40 low),
all fixed except one that belonged to another file and was fixed there.
The high one: the best rational approximation with a bounded denominator
is the last convergent *or a semiconvergent*, not truncation alone.

Known gap: none recorded. A fact the DP area could not verify (the
rectilinear Steiner tree on a rectangle boundary, S5-1 M) was later added to
*Geometry facts* by the strings-and-geometry area, checked against exact
Steiner trees on small grids.
