# knapster

Simple knapsack linear problem solver using branch & bound.

```rust
// the constraints of this problem might be a little weird, but it's what i needed.
// max z = 2x1 + 3x2 + 3x3 + 5x4 + 2x5 + 4x6
// s.t. 11x1 + 8x2 + 6x3 + 14x4 + 10x5 + 10x6 <= 40
// x1, x2, x3, x4, x5, x6 <= 1

let values = vec![2.0, 3.0, 3.0, 5.0, 2.0, 4.0];
let weights = vec![11.0, 8.0, 6.0, 14.0, 10.0, 10.0];
let max_weight = 40.0;

// ...setup coefficient matrices

// output written to branches.txt

// === Processing Problem 0 ===
// Problem 0: Branching on variable 5 with value 0.20000005
// --- Creating branch 0.1 (x5 ≤ 0) ---
// <> Initial pivoting for dual problem with entering variable x1 and leaving row 8
// Problem 0.1: Found optimal solution with objective value: 15.364
// Problem 0.1: Variable values: x1 = 0.182 x2 = 1.000 x3 = 1.000 x4 = 1.000 x5 = 0.000 x6 = 1.000 s1 = 0.000 s2 = 0.818 s3 = 0.000 s4 = 0.000 s5 = 0.000 s6 = 1.000 s7 = 0.000 s8 = 0.000 
// --- Creating branch 0.2 (x5 ≥ 1) ---
// <> Initial pivoting for dual problem with entering variable s5 and leaving row 8
// Problem 0.2: Found optimal solution with objective value: 14.143
// Problem 0.2: Variable values: x1 = 0.000 x2 = 1.000 x3 = 1.000 x4 = 0.429 x5 = 1.000 x6 = 1.000 s1 = 0.000 s2 = 1.000 s3 = 0.000 s4 = 0.000 s5 = 0.571 s6 = 0.000 s7 = 0.000 e8 = 0.000 

// === Processing Problem 0.1 ===
// Problem 0.1: Branching on variable 1 with value 0.18181822
// --- Creating branch 0.1.1 (x1 ≤ 0) ---
// <> Initial pivoting for dual problem with entering variable s1 and leaving row 9
// Problem 0.1.1: Found optimal solution with objective value: 15.000
// Problem 0.1.1: Variable values: x1 = 0.000 x2 = 1.000 x3 = 1.000 x4 = 1.000 x5 = 0.000 x6 = 1.000 s1 = 2.000 s2 = 1.000 s3 = 0.000 s4 = 0.000 s5 = 0.000 s6 = 1.000 s7 = 0.000 s8 = 0.000 s9 = 0.000 
// --- Creating branch 0.1.2 (x1 ≥ 1) ---
// <> Initial pivoting for dual problem with entering variable s8 and leaving row 9
// Problem 0.1.2: Found optimal solution with objective value: 13.786
// Problem 0.1.2: Variable values: x1 = 1.000 x2 = 1.000 x3 = 1.000 x4 = 0.357 x5 = 0.000 x6 = 1.000 s1 = 0.000 s2 = 0.000 s3 = 0.000 s4 = 0.000 s5 = 0.643 s6 = 1.000 s7 = 0.000 s8 = 0.000 e9 = 0.000 

// === Processing Problem 0.2 ===
// Problem 0.2: Branching on variable 4 with value 0.42857146
// --- Creating branch 0.2.1 (x4 ≤ 0) ---
// <> Initial pivoting for dual problem with entering variable e8 and leaving row 9
// Problem 0.2.1: Found optimal solution with objective value: 13.091
// Problem 0.2.1: Variable values: x1 = 0.545 x2 = 1.000 x3 = 1.000 x4 = 0.000 x5 = 1.000 x6 = 1.000 s1 = 0.000 s2 = 0.455 s3 = 0.000 s4 = 0.000 s5 = 1.000 s6 = 0.000 s7 = 0.000 e8 = 0.000 s9 = 0.000 
// --- Creating branch 0.2.2 (x4 ≥ 1) ---
// <> Initial pivoting for dual problem with entering variable s3 and leaving row 9
// Problem 0.2.2: Found optimal solution with objective value: 14.000
// Problem 0.2.2: Variable values: x1 = 0.000 x2 = 0.000 x3 = 1.000 x4 = 1.000 x5 = 1.000 x6 = 1.000 s1 = 0.000 s2 = 1.000 s3 = 1.000 s4 = 0.000 s5 = 0.000 s6 = 0.000 s7 = 0.000 e8 = 0.000 e9 = 0.000 

// === Processing Problem 0.1.1 ===
// Problem 0.1.1: Branching on variable 1 with value 0.000000014901161
// --- Creating branch 0.1.1.1 (x1 ≤ 0) ---
// Problem 0.1.1.1: Found optimal solution with objective value: 15.000
// Problem 0.1.1.1: Variable values: x1 = 0.000 x2 = 1.000 x3 = 1.000 x4 = 1.000 x5 = 0.000 x6 = 1.000 s1 = 2.000 s2 = 1.000 s3 = 0.000 s4 = 0.000 s5 = 0.000 s6 = 1.000 s7 = 0.000 s8 = 0.000 s9 = 0.000 s10 = -0.000 
// --- Creating branch 0.1.1.2 (x1 ≥ 1) ---
// <> Infeasible dual problem, no entering variable found. Attempted with leaving row 10
// Problem 0.1.1.2: Infeasible or unbounded
// ...
```