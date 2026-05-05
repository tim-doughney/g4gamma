// Bateman.hh -- analytic activity solver for a (branched) decay chain.
//
// Given a topologically-ordered chain (parents before daughters) with mean
// lives lambda_i = 1/tau_i and edge branching ratios BR(p -> d), compute the
// activity A_i(t) of each chain member when at t=0 only the root is active
// (with A_root(0) = 1).
//
// Methods:
//   solve(t)           -- finite time t. Uses analytic generalised Bateman:
//                         A_i(t) = sum_j c_{i,j} exp(-lambda_j t)
//   solveSecularEq()   -- t -> infinity, but interpreting as "equilibrium with
//                         the root pinned at A=1". Implementation walks the
//                         topology and propagates: A_i = sum_{p parent} A_p * BR.
//                         Stable nodes get A=0 (no decays).
//
// Edge case: if two unstable chain members have lambda values that differ by
// less than `kRepeatedEigenvalueEps`, we add a tiny perturbation to one of
// them. For real ENSDF half-lives this never bites in practice.
#pragma once

#include "g4gamma/ChainBuilder.hh"
#include <vector>

namespace g4gamma {

class Bateman {
public:
    explicit Bateman(const std::vector<ChainNode>& chain);

    // A_i(t) for each i, with A_root(0)=1. t in internal units.
    std::vector<double> solve(double t) const;

    // A_i at t -> infty under "input from primary at constant activity 1"
    // interpretation, which for one-shot initial activity is the same as
    // walking the chain with running BRs (i.e. activity reaches *transient*
    // equilibrium for each daughter relative to its parent).
    //
    // When the parent has finite mean life and the daughter's mean life is
    // shorter, the daughter reaches A_p in the secular limit; if longer, the
    // daughter reaches a transient-equilibrium ratio. Here we adopt the
    // simpler convention used in the user's rdecay01 setup with
    // thresholdForVeryLongDecayTime: every unstable chain member has
    // activity equal to the sum over root->i path products of BRs (the limit
    // as the parent half-life >> all daughter half-lives, which is true for
    // the natural decay series given a long-lived primary parent).
    std::vector<double> solveSecularEq() const;

    // Lambda for each chain node (1/meanLife). 0 for stable nodes.
    const std::vector<double>& lambdas() const { return fLambda; }

private:
    const std::vector<ChainNode>& fChain;
    std::vector<double> fLambda;
    std::vector<std::vector<int>> fParents;        // node -> list of parent indices in chain
    std::vector<std::vector<double>> fParentBR;    // node -> branching ratios from each parent (matches fParents)
    std::vector<int> fTopoOrder;                   // topological order (root first)
};

} // namespace g4gamma
