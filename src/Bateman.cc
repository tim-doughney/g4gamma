#include "g4gamma/Bateman.hh"
#include <cmath>
#include <queue>
#include <unordered_set>
#include <stdexcept>
#include <algorithm>

namespace g4gamma {

namespace {
constexpr double kRepeatedEigenvalueEps = 1e-12;  // in 1/internal-time units

// Topological sort using Kahn's algorithm. The chain may have multiple paths
// to a daughter (branched decays), so we need a real toposort, not BFS depth.
std::vector<int> topoSort(const std::vector<ChainNode>& chain) {
    int n = chain.size();
    std::vector<int> indeg(n, 0);
    for (int i = 0; i < n; ++i) {
        for (const auto& e : chain[i].edges) ++indeg[e.daughterIndex];
    }
    std::queue<int> q;
    for (int i = 0; i < n; ++i) if (indeg[i] == 0) q.push(i);
    std::vector<int> order;
    while (!q.empty()) {
        int u = q.front(); q.pop();
        order.push_back(u);
        for (const auto& e : chain[u].edges) {
            if (--indeg[e.daughterIndex] == 0) q.push(e.daughterIndex);
        }
    }
    if (static_cast<int>(order.size()) != n) {
        // Cycle detected: fall back to insertion order.
        order.clear();
        for (int i = 0; i < n; ++i) order.push_back(i);
    }
    return order;
}
} // namespace

Bateman::Bateman(const std::vector<ChainNode>& chain) : fChain(chain) {
    int n = chain.size();
    fLambda.assign(n, 0.0);
    fParents.assign(n, {});
    fParentBR.assign(n, {});
    for (int i = 0; i < n; ++i) {
        if (chain[i].meanLife > 0.0) fLambda[i] = 1.0 / chain[i].meanLife;
        for (const auto& e : chain[i].edges) {
            fParents[e.daughterIndex].push_back(i);
            fParentBR[e.daughterIndex].push_back(e.branchingRatio);
        }
    }
    fTopoOrder = topoSort(chain);

    // Perturb identical lambdas slightly to avoid 0/0 in finite-time solver.
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (fLambda[i] > 0.0 && fLambda[j] > 0.0 &&
                std::abs(fLambda[i] - fLambda[j]) <
                    kRepeatedEigenvalueEps * std::max(fLambda[i], fLambda[j])) {
                fLambda[i] *= (1.0 + 1e-9);  // tiny perturbation
            }
        }
    }
}

std::vector<double> Bateman::solveSecularEq() const {
    int n = fChain.size();
    std::vector<double> A(n, 0.0);
    A[0] = 1.0;  // root has activity 1 by convention
    // Topologically propagate: A_i = sum_{p parent} A_p * BR(p->i)
    for (int idx : fTopoOrder) {
        if (idx == 0) continue;
        double sum = 0.0;
        for (size_t k = 0; k < fParents[idx].size(); ++k) {
            int p = fParents[idx][k];
            sum += A[p] * fParentBR[idx][k];
        }
        A[idx] = sum;
    }
    // Stable nodes (no decay channels) emit no gammas. We use the chain's
    // `stable` flag rather than lambda==0, because the half-life column in
    // the Geant4 RadioactiveDecay zZ.aA files is often a placeholder (0 or
    // dummy) -- the real half-life lives in ENSDFSTATE which we don't read.
    // What matters for "do gammas come from this node" is whether it has any
    // decay channels at all.
    for (int i = 0; i < n; ++i) {
        if (fChain[i].stable) A[i] = 0.0;
    }
    return A;
}

// Finite-time generalised Bateman with branching.
//
// We work with N_i(t) = number of atoms of species i, with N_root(0) = 1/lambda_root
// (so A_root(0) = 1 Bq) and N_i(0) = 0 for i != root.
//
// Solution form:  N_i(t) = sum_j n_{i,j} exp(-lambda_j t)
//                 A_i(t) = lambda_i * N_i(t)
//
// Recursion (in topological order):
//   For root r: n_{r,r} = 1/lambda_r, all other n_{r,j} = 0.
//   For each i != root:
//     For each lambda_j != lambda_i (j < i in topological order, only those
//     reachable from root):
//        n_{i,j} = (1/(lambda_i - lambda_j)) * sum_{p in parents(i)}
//                                                 BR(p->i) * lambda_p * n_{p,j}
//     n_{i,i} = -sum_{j != i} n_{i,j}    (so that N_i(0) = 0)
//
// Stable nodes (lambda = 0) absorb activity. For our purposes (gamma flux from
// decay), stable nodes contribute zero gammas, so we can leave them out of the
// activity vector entirely. We still track their N coefficients for accounting.
std::vector<double> Bateman::solve(double t) const {
    int n = fChain.size();
    if (t <= 0.0) {
        std::vector<double> A(n, 0.0);
        if (fLambda[0] > 0.0) A[0] = 1.0;
        return A;
    }

    // n_coeff[i][j] = coefficient of exp(-lambda_j t) in N_i(t)
    // Use a dense n x n matrix. n is typically <100 even for 238U chain.
    std::vector<std::vector<double>> nc(n, std::vector<double>(n, 0.0));

    // Root: N_root(0) = 1/lambda_root, single exp(-lambda_root t)
    int root = fTopoOrder.front();
    if (fLambda[root] <= 0.0) {
        // Root has no known mean life (likely the half-life column in the
        // RadioactiveDecay zZ.aA file is a dummy 0; the real value lives in
        // ENSDFSTATE which we don't read). For finite-time evolution we can't
        // do better than to fall back to secular equilibrium.
        if (!fChain[root].stable) {
            return solveSecularEq();
        }
        // Genuinely stable root (no channels) -> nothing decays.
        std::vector<double> A(n, 0.0);
        return A;
    }
    nc[root][root] = 1.0 / fLambda[root];

    for (int idx : fTopoOrder) {
        if (idx == root) continue;
        // Sum over parents
        std::vector<double> contrib(n, 0.0);
        for (size_t k = 0; k < fParents[idx].size(); ++k) {
            int p = fParents[idx][k];
            double br = fParentBR[idx][k];
            for (int j = 0; j < n; ++j) {
                if (nc[p][j] == 0.0) continue;
                // n_{i,j} += (1/(lambda_i - lambda_j)) * BR * lambda_p * n_{p,j}
                if (fLambda[idx] == fLambda[j]) continue;
                double denom = fLambda[idx] - fLambda[j];
                contrib[j] += br * fLambda[p] * nc[p][j] / denom;
            }
        }
        // Sum into nc[idx]
        for (int j = 0; j < n; ++j) nc[idx][j] += contrib[j];
        // Determine n_{idx,idx} so N_idx(0) = 0:
        double sumOther = 0.0;
        for (int j = 0; j < n; ++j) if (j != idx) sumOther += nc[idx][j];
        nc[idx][idx] = -sumOther;
    }

    // Now compute A_i(t) = lambda_i * sum_j nc[i][j] * exp(-lambda_j t)
    std::vector<double> A(n, 0.0);
    for (int i = 0; i < n; ++i) {
        if (fLambda[i] <= 0.0) continue;
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (nc[i][j] == 0.0) continue;
            sum += nc[i][j] * std::exp(-fLambda[j] * t);
        }
        A[i] = fLambda[i] * sum;
        if (A[i] < 0.0 && A[i] > -1e-12) A[i] = 0.0;  // numerical noise
    }
    return A;
}

} // namespace g4gamma
