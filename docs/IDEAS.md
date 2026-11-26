# Ideas behind the jewels

For a Boolean assignment x, let u_j(x)=1 when clause j is violated and 0 otherwise. With strictly positive weights,

E_w(x) = sum_j w_j u_j(x)

vanishes exactly at the satisfying assignments. We can alter the penalties while retaining those points.

The experiments try fixed cosine weights, adaptive positive boosts, and diffusion across a tiling's adjacency graph. The identity

3+4 cos(a)+cos(2a) = 2(1+cos(a))^2 >= 0

supplies one positive boost field. A small extra positive term ensures a violated clause always receives a boost. The broad sweep normalizes the mean boost over violated clauses to one, so different schedules do not win merely by applying larger increments.

Hat and pinwheel orientations supply phases; their shared edges supply the averaging graph for diffusion. The assignment of logical clauses to geometric tiles remains an experimental choice, and is tested under relabelings.

Golden-ratio powers and reciprocals can encode the same or reversed rotations modulo one. They are related representations, not independent sources of evidence. The other constants provide controls and alternative schedules.

Complex projection and cut-and-project geometry make the objects visible. They must retain labels and constraints: nearby drawn points need not be adjacent in the SAT search graph. The full protocol and proof scope are in the experiment reports.
