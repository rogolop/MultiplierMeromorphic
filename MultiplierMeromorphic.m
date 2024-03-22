import "SingularitiesDim2/ProximityMatrix.m": ProximityMatrixImpl, CoefficientsVectorBranch;
import "SingularitiesDim2/LogResolution.m": ComputeLogResolutionData, ExpandWeightedCluster;
import "SingularitiesDim2/IntegralClosure.m": Unloading; 
// IntegralClosureIrreducible, ProductIdeals, ClusterFactorization;


intrinsic LogResolutionMeromorphic(f::RngMPolLocElt, g::RngMPolLocElt : Coefficients := false) -> []
{ Computes the weighted cluster of base points of a bivariate
  polynomial ideal I }
  // Generators in G & fixed part F.
  G := [f, g]; //Basis(I);
  F := Gcd(G); G := [ExactQuotient(h, F) : h in G];

  ////////////// Compute all information ////////////////
  S := PuiseuxExpansion(G: Polynomial := true);
  P, EE, CC := ProximityMatrixImpl([*<s[1], 1> : s in S*]: ExtraPoint := true);
  
  // printf "P = \n%o\n", P;
  
  // printf "LogResolutionMeromorphic EE = \n%o\n", EE;

  E := []; // Multiplicities of each generator.
  C := []; // Coefficients of BP(I).
  V := []; // Vector a values for each generator.
  v := []; // Virtual values of BP(I).
  ComputeLogResolutionData(~P, ~EE, ~CC, ~S, #G, ~E, ~C, ~V, ~v);

  // /////////////// Add new free points /////////////////////
  // lastFree := [i : i in [1..Ncols(P)] | (&+P[1..Ncols(P)])[i] eq 1];
  // points2test := #lastFree; idx := 1;
  // // For each last free point on a branch...
  // while points2test gt 0 do
  //   // Values for each gen. at p.
  //   p := lastFree[idx]; Vp := [vi[1][p] : vi in V];
  //   // Generators achieving the minimum.
  //   GG := [i : i in [1..#Vp] | Vp[i] eq Min(Vp)];
  //   // If the multiplicities of all the generators achieving the minimum
  //   // at p is > 0 add new point.
  //   if &and[E[g][1][p] ne 0 : g in GG] then
  //     // The (unique) branch of the generator 'g' where 'p' belongs.
  //     assert(#[i : i in [1..#EE] | EE[i][1, p] ne 0] eq 1);
  //     b := [i : i in [1..#EE] | EE[i][1, p] ne 0][1];
  //     ExpandWeightedCluster(~P, ~EE, ~CC, ~S, b); P[Ncols(P)][p] := -1;
  //     ComputeLogResolutionData(~P, ~EE, ~CC, ~S, #G, ~E, ~C, ~V, ~v);
  //     // We may need to add more free points after the points we added.
  //     lastFree cat:= [Ncols(P)]; points2test := points2test + 1;
  //   end if;
  // points2test := points2test - 1; idx := idx + 1;
  // end while;

  // /////////////// Add new satellite points /////////////////////
  // points2test := Ncols(P) - 1; p := 2; // Do not start at the origin.
  // while points2test gt 0 do
  //   // Values for the generators at point p.
  //   Vp := [vi[1][p] - v[1][p] : vi in V];
  //   // Points p is proximate to && Points proximate to p.
  //   p_prox := [i : i in [1..Ncols(P)] | P[p][i] eq -1];
  //   prox_p := [i : i in [1..Ncols(P)] | P[i][p] eq -1];
  //   Q := [q : q in p_prox | &+Eltseq(Submatrix(P, prox_p, [q])) eq 0];
  //   for q in Q do
  //     // Values for the generators at point q.
  //     Vq := [vi[1][q] - v[1][q] : vi in V];
  //     if &*[Vp[i] + Vq[i] : i in [1..#Vp]] ne 0 then
  //       ExpandWeightedCluster(~P, ~EE, ~CC, ~S, -1);
  //       P[Ncols(P)][p] := -1; P[Ncols(P)][q] := -1;
  //       ComputeLogResolutionData(~P, ~EE, ~CC, ~S, #G, ~E, ~C, ~V, ~v);
  //       // We may need to add more satellite points after the points we added.
  //       points2test := points2test + 1;
  //     end if;
  //   end for;
  // points2test := points2test - 1; p := p + 1;
  // end while;

  /////////////// Remove non base points ////////////////
  // Multiplicities for the cluster of base points.
  e := v * Transpose(P);
  // I := [i : i in [1..Ncols(P)] ]; // | e[1][i] ne 0
  // Remove points not in the cluster of base points.
  // P := Submatrix(P, I, I); v := Submatrix(v, [1], I); C := C[I];

  // Select 1 as affine part iff F is a unit.
  F := Evaluate(F, <0, 0>) ne 0 select Parent(F)!1 else F;
  
  A := V[1]; // values of f
  B := V[2]; // values of g
  
  // printf "e = %o\n", e;
  // printf "v = %o\n", v;
  printf "A = %o\n", A;
  // printf "eA = %o\n", A * Transpose(P); // multiplicities of f
  Excess_f := A* Transpose(P)*P;
  printf "B = %o\n", B;
  // printf "eB = %o\n", B * Transpose(P); // multiplicities of g
  
  printf "C = %o\n", C;
  
  // AmB := Matrix(IntegerRing(), 1, Ncols(A), [Max([A[1,i]-B[1,i], 0]) : i in [1..Ncols(A)]]);
  AmB := ZeroMatrix(IntegerRing(), 1, Ncols(P));
  for i in [1..Ncols(P)] do AmB[1][i] := Max([A[1,i]-B[1,i], 0]); end for;
  
  if Coefficients then return P, AmB, F, C, Excess_f, A;
  else return P, AmB, F, Excess_f, A; end if;
end intrinsic;




intrinsic MultiplierIdealsMeromorphic(f::RngMPolLocElt, g::RngMPolLocElt : MinJN:=0, MaxJN:=1, ComputeIdeals:=true) -> List
{ Computes the Multiplier Ideals and its associated Jumping Number for an
  plane curve in a smooth complex surface using the algorithm
  of Alberich-Alvarez-Dachs.
  
  Returns:
      [* dataJN1, dataJN2, ... *]
  
  "dataJNi" is the data corresponding to a jumping number:
      <
        data of the multiplier ideal (only if ComputeIdeals is true),
        jumping number
      >
  
  Data of the multiplier ideal:
      <
        < f, exponent of common f >,
        [ generator1, generator2, ... ]
      >
  
  The multiplier ideal represented by this data is:
      f^(exponent of common f) * (generator1, generator2, ...)
  }

  // With the extra point there is no confusion whether and affine component
  // has multiplicity or not.
  P, E, _, C, Excess_f, A := LogResolutionMeromorphic(f, g: Coefficients := true);
  // printf "----------------\n";
  printf "P = \n%o\n", P;
  // printf "E = \n%o\n", E;
  // printf "Excess_f = \n%o\n", Excess_f;
  
  // P, E, C := ProximityMatrix(f: Coefficients := true, ExtraPoint := true);
  QQ := Rationals();
  EQ := ChangeRing(E, QQ);
  PQ := ChangeRing(P, QQ);
  ZZ := Integers();
  PQTinv := Transpose(PQ)^-1;
  k := Parent(f);
  // k := Universe(Basis(I));
  n := Ncols(P);
  // Compute relative canonical divisor.
  K := Matrix([[QQ | 1 : i in [1..n]]]);
  K := K*PQTinv;
  // printf "K = \n%o\n", K;
  
  F := EQ;
  // F := EQ*PQTinv;
  
  // Compute the extended intersection matrix by the stict transform components.
  N := Transpose(PQ)*PQ;
  StrF := Excess_f;
  // StrF := EQ*PQ;
  nAffComp := #[1 : i in [1..n] | StrF[1][i] ne 0];
  N := DiagonalJoin(N, ZeroMatrix(QQ, nAffComp)); //-IdentityMatrix(QQ, nAffComp));
  idxAff := [i : i in [1..n] | StrF[1][i] ne 0];
  for i in [1..nAffComp] do N[n + i][idxAff[i]] := -1; end for;
  // for i in [1..nAffComp] do N[n + i][idxAff[i]] := -StrF[1][idxAff[i]]; end for;
  // printf "N = \n%o\n", N;
  
  // F := HorizontalJoin(F, Matrix(QQ, [[1 : i in [1..nAffComp]]]));
  F := HorizontalJoin(F, Matrix(QQ, [[StrF[1][i] : i in idxAff]]));
  K := HorizontalJoin(K, Matrix(QQ, [[0 : i in [1..nAffComp]]]));
  
  A := ChangeRing(A, QQ);
  A := HorizontalJoin(A, Matrix(QQ, [[StrF[1][i] : i in idxAff]]));

  // printf "F = %o\n", F;
  
  printf "\n";
  JN := MinJN; S := [**];
  while JN lt MaxJN do
    printf ".";
    D := Unloading(N, Matrix([[QQ | Floor(ei) : ei in Eltseq(JN*F - K)]]));
    // printf "JN*F - K = %o\n", Matrix([[QQ | Floor(ei) : ei in Eltseq(JN*F - K)]]);
    printf "D = %o\n", D;
    
    jnfmk := Matrix([[QQ | Floor((JN*F - K)[1][i]) - Floor(JN) * A[1][i] : i in [1..(n+nAffComp)]]]);
    printf "JN*F - K menys f... = %o\n", jnfmk;
    printf "Unload(JN*F - K menys f...) = %o\n", Unloading(N, jnfmk);
    
    if JN ne 0 then
      if ComputeIdeals then
        D2 := Matrix([[QQ | D[1][i] - Floor(JN) * A[1][i] : i in [1..(n+nAffComp)]]]);
        printf "Divisor per O(D): %o\n", D2;
        DZZ := ColumnSubmatrix(ChangeRing(D2, ZZ), n);
        gen := GeneratorsOXD(P, DZZ, C, k);
        // gen := [Factorization(h) : h in gen];
        // if gen eq [] then gen := [Parent(f)!1]; end if;
        // if gen eq [] then gen := [[<Parent(f)!1,1>]]; end if;
        gen := < <f,Floor(JN)>, gen >;
        S cat:= [*<gen, JN>*];
      else
        S cat:= [*<0, JN>*];
      end if;
    end if;
    
    lastJN := JN;
    JN, i := Min([(K[1][i] + 1 + D[1][i])/F[1][i] : i in [1..(n+nAffComp)] | F[1][i] gt 0]);
    // printf "%o\n", [(K[1][i] + 1 + D[1][i])/F[1][i] : i in [1..(n+nAffComp)] | F[1][i] ne 0];
    // printf "JN = %o, i = %o\n", JN, i;
    
    printf "\nJN = %-15o, i = %o\n", JN, [j : j in [1..(n+nAffComp)] | F[1][j] ne 0 and (K[1][j] + 1 + D[1][j])/F[1][j] eq JN];
        
  end while;
  return S;
end intrinsic;




