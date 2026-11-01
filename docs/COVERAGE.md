# Coverage of the training sets

One section per problem set.

## UCup Seasons 4-5

Goal: every theorem or niche fact needed by the 474 training problems of
UCup Seasons 4-5 (from the "UCup S4-S5 Training" material: problem
catalog, full solutions, 23 study packets) is stated somewhere in the
packet, and every problem's key idea is reachable from a recognition cue.

### Method

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

### Numbers

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

## KAIST 16th ICPC Mock Competition (2026)

Goal: find what the 13 problems (A to M) need that snowi does not state, add
the ideas worth keeping, and make every problem's key idea reachable from a
recognition cue. Additions went to the full packet only: the ICPC packets are
at their 25-page limit and did not change.

### Method

1. Every problem was solved and its solution checked against one independent
   brute force. Every sample was rerun separately, and the two interactive
   problems (D, I) ran against local interactors.
2. Each fact and tool the solutions use was looked up in snowi and given a
   print status. A read-only review then searched `content/` for every gap,
   under alternative wordings, and confirmed the print statuses against the
   text of all four PDFs.
3. 14 reusable ideas went into the notes, written as general facts rather than
   contest write-ups. Each statement has its own fresh brute force or exact
   check of its wording (a measurement for the relabelling speed-up), with a
   deliberately wrong variant to show the check can fail. The same holds for
   the claims in the new cue rows.
4. Three read-only reviews of the additions raised 10 findings, all fixed, on
   the wording of statements, one layout regression in the full PyPy packet
   and cue rows with no note behind them. The statements those fixes changed
   were then checked a second time by separate code.

### Numbers

| Status before the additions | Facts | Tools |
|---|---|---|
| Both ICPC packets | 1 | 12 |
| C++ ICPC packet only | 0 | 1 |
| Full packet only | 0 | 2 |
| Partly stated | 10 | 2 |
| Missing | 66 | 0 |
| Textbook, not needed | 12 | 17 |
| **Total** | 89 | 34 |

Review: all 78 missing or partly stated rows were searched and none moved to
stated. 7 rows had their notes corrected, with no change of status. 32 rows
were spot-checked, and all 68 file and line citations resolved.

Cues: before the additions the key ideas of 3 problems (F, H, M) were
reachable from a cue. All 13 now are, through 17 new rows in *Recognition
cues* and 2 extended lines in *Technique finder*.

### Ideas added

| # | Idea | Problem | Note |
|---|---|---|---|
| 1 | All minimum cuts from the residual SCCs, global cuts, their count | B | *Flow and cut modelling*, pointer in `flow/global-min-cut` |
| 2 | Exactly X ways by double-and-add, exactly X domino tilings on two rails, snake graphs | L | *Constructive problems and invariants*, *Number theory facts* |
| 3 | Commuting moves | A | *Constructive problems and invariants* |
| 4 | Wake-up lists | A | *Sqrt decomposition and small tricks* |
| 5 | Paths through a chain of rectangles, binomial kernel | E | *Counting formulas*, *Polynomial techniques* |
| 6 | Weighted tree medians | K | *Tree techniques* |
| 7 | Three-part split | M | *Constructive problems and invariants* |
| 8 | Neighbour-distinct grid colouring | I | *Constructive problems and invariants* |
| 9 | One question per level in a query, injective transcripts | D | *Constructive problems and invariants*, pointer in `misc/interactive` |
| 10 | Removing a block so a scan succeeds | H | *DP techniques* |
| 11 | Sequential overwrite on a functional graph | J | *Graph theory facts* |
| 12 | Greedy walk to the root, stable bucket split and parity in centroid decomposition, relabelling first | G | *Tree techniques* |
| 13 | Best in a group | J | *Designing a segment tree node* |
| 14 | Partition into collapsing segments | C | *DP techniques* |

The other missing facts were judged problem-specific and left out.

Known gaps:

- The snake construction for L is tested, not proven, for every X. A continued
  fraction that fits was found for every X below 12,000,000, for the top
  1,500,000 values below 2^30 and on random samples in between, so about 1.3%
  of the range is exhaustive. *Number theory facts* marks its bound as
  empirical. L itself has no gap: the two-rail construction in *Constructive
  problems and invariants* is proven from the exact transfer table of each
  block, and it was accepted by the judge.
- That is the only judge verdict. The other solutions were checked locally,
  since no judge was available.
- The cactus representation in *All minimum cuts* is a citation, not computed.
  Flow statements were tested with integer weights only.
- *Sequential overwrite* is exhaustive up to n = 4 and random above. The
  Hamiltonian path claims in the colouring bullet are exhaustive up to 36
  cells, and larger boards rest on a proof.
- The relabelling speed-ups were measured on one loaded box.
- The comment in `tree/centroid.py` says relabelling is 2x faster, which holds
  on paths only. It prints in the ICPC packet, so it was left alone.
