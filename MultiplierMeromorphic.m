import "SingularitiesDim2/ProximityMatrix.m": ProximityMatrixImpl, CoefficientsVectorBranch;
import "SingularitiesDim2/LogResolution.m": ComputeLogResolutionData, ExpandWeightedCluster;
import "SingularitiesDim2/IntegralClosure.m": Unloading; 
// IntegralClosureIrreducible, ProductIdeals, ClusterFactorization;


intrinsic LogResolutionMeromorphic(f::RngMPolLocElt, g::RngMPolLocElt) -> []
{ Computes the weighted cluster of base points of a bivariate meromorphic ideal (f/g) }
	// Generators in G & fixed part commonFactor.
	G := [f, g]; //Basis(I);
	commonFactor := Gcd(G); G := [ExactQuotient(h, commonFactor) : h in G];

	////////////// Compute all information ////////////////
	S := PuiseuxExpansion(G: Polynomial := true);
	Prox, EE, CC := ProximityMatrixImpl([*<s[1], 1> : s in S*]: ExtraPoint := true);
	
	// printf "Prox = \n%o\n", Prox;
	
	// printf "LogResolutionMeromorphic EE = \n%o\n", EE;

	E := []; // Multiplicities of each generator.
	Coeffs := []; // Coefficients of BP(I).
	V := []; // Vector a values for each generator.
	v := []; // Virtual values of BP(I).
	ComputeLogResolutionData(~Prox, ~EE, ~CC, ~S, #G, ~E, ~Coeffs, ~V, ~v);

	/////////////// Add new free points /////////////////////
	lastFree := [i : i in [1..Ncols(Prox)] | (&+Prox[1..Ncols(Prox)])[i] eq 1];
	points2test := #lastFree; idx := 1;
	// For each last free point on a branch...
	while points2test gt 0 do
		// Values for each gen. at p.
		p := lastFree[idx]; Vp := [vi[1][p] : vi in V];
		// Generators achieving the minimum.
		GG := [i : i in [1..#Vp] | Vp[i] eq Min(Vp)];
		// If the multiplicities of all the generators achieving the minimum
		// at p is > 0 add new point.
		if &and[E[g][1][p] ne 0 : g in GG] then
			// The (unique) branch of the generator 'g' where 'p' belongs.
			assert(#[i : i in [1..#EE] | EE[i][1, p] ne 0] eq 1);
			b := [i : i in [1..#EE] | EE[i][1, p] ne 0][1];
			ExpandWeightedCluster(~Prox, ~EE, ~CC, ~S, b); Prox[Ncols(Prox)][p] := -1;
			ComputeLogResolutionData(~Prox, ~EE, ~CC, ~S, #G, ~E, ~Coeffs, ~V, ~v);
			// We may need to add more free points after the points we added.
			lastFree cat:= [Ncols(Prox)]; points2test := points2test + 1;
		end if;
	points2test := points2test - 1; idx := idx + 1;
	end while;

	/////////////// Add new satellite points /////////////////////
	points2test := Ncols(Prox) - 1; p := 2; // Do not start at the origin.
	while points2test gt 0 do
		// Values for the generators at point p.
		Vp := [vi[1][p] - v[1][p] : vi in V];
		// Points p is proximate to && Points proximate to p.
		p_prox := [i : i in [1..Ncols(Prox)] | Prox[p][i] eq -1];
		prox_p := [i : i in [1..Ncols(Prox)] | Prox[i][p] eq -1];
		Q := [q : q in p_prox | &+Eltseq(Submatrix(Prox, prox_p, [q])) eq 0];
		for q in Q do
			// Values for the generators at point q.
			Vq := [vi[1][q] - v[1][q] : vi in V];
			if &*[Vp[i] + Vq[i] : i in [1..#Vp]] ne 0 then
				ExpandWeightedCluster(~Prox, ~EE, ~CC, ~S, -1);
				Prox[Ncols(Prox)][p] := -1; Prox[Ncols(Prox)][q] := -1;
				ComputeLogResolutionData(~Prox, ~EE, ~CC, ~S, #G, ~E, ~Coeffs, ~V, ~v);
				// We may need to add more satellite points after the points we added.
				points2test := points2test + 1;
			end if;
		end for;
	points2test := points2test - 1; p := p + 1;
	end while;

	/////////////// Remove non base points ////////////////
	// Multiplicities for the cluster of base points.
	e := v * Transpose(Prox);
	// I := [i : i in [1..Ncols(Prox)] ]; // | e[1][i] ne 0
	// Remove points not in the cluster of base points.
	// Prox := Submatrix(Prox, I, I); v := Submatrix(v, [1], I); Coeffs := Coeffs[I];

	// Select 1 as affine part iff commonFactor is a unit.
	commonFactor := Evaluate(commonFactor, <0, 0>) ne 0 select Parent(commonFactor)!1 else commonFactor;
	
	Nf := V[1]; // values of f
	Ng := V[2]; // values of g
	
	// printf "e = %o\n", e;
	// printf "v = %o\n", v;
	// printf "Nf = %o\n", Nf;
	// printf "eA = %o\n", Nf * Transpose(Prox); // multiplicities of f
	// Excess_f := Nf* Transpose(Prox)*Prox;
	// printf "Ng = %o\n", Ng;
	// printf "eB = %o\n", Ng * Transpose(Prox); // multiplicities of g
	
	// printf "Coeffs = %o\n", Coeffs;
	
	// N := Matrix(IntegerRing(), 1, Ncols(Nf), [Max([Nf[1,i]-Ng[1,i], 0]) : i in [1..Ncols(Nf)]]);
	N := ZeroMatrix(IntegerRing(), 1, Ncols(Prox));
	for i in [1..Ncols(Prox)] do N[1][i] := Max([Nf[1,i]-Ng[1,i], 0]); end for;
	// printf "N (>0) = %o\n", N;
	
	// if Coefficients then 
	return Nf, Ng, N, Prox, commonFactor, Coeffs; // return Prox, N, commonFactor, Coeffs, Excess_f, Nf;
	// else return Prox, N, commonFactor, Excess_f, Nf; end if;
end intrinsic;




intrinsic MultiplierIdealsMeromorphic(f::RngMPolLocElt, g::RngMPolLocElt : MinJN:=0, MaxJN:=1, ComputeIdeals:=true) -> List
{ Computes the Multiplier Ideals and their associated Jumping Numbers for a meromorphic function f/g, using the algorithm of Alberich-Alvarez-Gomez. Returns a list of tuples of the form: <jumping number, <fExponent, generators>>. The multiplier ideal represented by this data is: f^fExponent * (generators). If ComputeIdeals is set to false, the output tuples are: <jumping number, 0>.
}

	// With the extra point there is no confusion whether and affine component
	// has multiplicity or not.
	Nf, Ng, N, Prox, commonFactor, Coeffs := LogResolutionMeromorphic(f, g);
	Excess_f := Nf* Transpose(Prox)*Prox;
	// printf "----------------\n";
	// printf "Prox = \n%o\n", Prox;
	// printf "N = \n%o\n", N;
	// printf "Excess_f = \n%o\n", Excess_f;
	
	// Prox, N, Coeffs := ProximityMatrix(f: Coefficients := true, ExtraPoint := true);
	QQ := Rationals();
	F := ChangeRing(N, QQ);
	ProxQ := ChangeRing(Prox, QQ);
	ZZ := Integers();
	PQTinv := Transpose(ProxQ)^-1;
	k := Parent(f);
	// k := Universe(Basis(I));
	n := Ncols(Prox);
	// Compute relative canonical divisor.
	K := Matrix([[QQ | 1 : i in [1..n]]]);
	K := K*PQTinv;
	// printf "K = \n%o\n", K;
	
	// F := EQ;
	// F := EQ*PQTinv;
	
	// Compute the extended intersection matrix by the stict transform components.
	Intersect := Transpose(ProxQ)*ProxQ;
	// StrF := Excess_f;
	// Excess_f := EQ*ProxQ;
	nAffComp := #[1 : i in [1..n] | Excess_f[1][i] ne 0];
	Intersect := DiagonalJoin(Intersect, ZeroMatrix(QQ, nAffComp)); //-IdentityMatrix(QQ, nAffComp));
	idxAff := [i : i in [1..n] | Excess_f[1][i] ne 0];
	for i in [1..nAffComp] do Intersect[n + i][idxAff[i]] := -1; end for;
	// for i in [1..nAffComp] do Intersect[n + i][idxAff[i]] := -Excess_f[1][idxAff[i]]; end for;
	// printf "Intersect = \n%o\n", Intersect;
	
	// F := HorizontalJoin(F, Matrix(QQ, [[1 : i in [1..nAffComp]]]));
	F := HorizontalJoin(F, Matrix(QQ, [[Excess_f[1][i] : i in idxAff]]));
	K := HorizontalJoin(K, Matrix(QQ, [[0 : i in [1..nAffComp]]]));
	
	Nf := ChangeRing(Nf, QQ);
	Nf := HorizontalJoin(Nf, Matrix(QQ, [[Excess_f[1][i] : i in idxAff]]));

	// printf "F = %o\n", F;
	
	// printf "\n";
	JN := MinJN;
	Out := [**];
	while JN lt MaxJN do
		// printf "JN = %o\n", JN;
		
		D := Unloading(Intersect, Matrix([[QQ | Floor(ei) : ei in Eltseq(JN*F - K)]]));
		// printf "D = %o\n", D;
		
		// Testing if D coincides with the following (?)
		// jnfmk := Matrix([[QQ | Floor((JN*F - K)[1][i]) - Floor(JN) * Nf[1][i] : i in [1..(n+nAffComp)]]]);
		// printf "JN*F - K menys f... = %o\n", jnfmk;
		// printf "Unload(JN*F - K menys f...) = %o\n", Unloading(Intersect, jnfmk);
		
		if JN ne 0 then
			if ComputeIdeals then
				D2 := Matrix([[QQ | D[1][i] - Floor(JN) * Nf[1][i] : i in [1..(n+nAffComp)]]]);
				// printf "Divisor per O(D): %o\n", D2;
				DZZ := ColumnSubmatrix(ChangeRing(D2, ZZ), n);
				gen := GeneratorsOXD(Prox, DZZ, Coeffs, k);
				if gen eq [] then gen := [Parent(f)| 1 ]; end if;
				gen := < Floor(JN), gen >;
				Out cat:= [*<JN, gen>*];
			else
				Out cat:= [*<JN, 0>*];
			end if;
		end if;
		
		lastJN := JN;
		JN, i := Min([(K[1][i] + 1 + D[1][i])/F[1][i] : i in [1..(n+nAffComp)] | F[1][i] gt 0]);

		// printf "JN = %o, i = %o\n", JN, i;
		// printf "\nJN = %-15o, i = %o\n", JN, [j : j in [1..(n+nAffComp)] | F[1][j] ne 0 and (K[1][j] + 1 + D[1][j])/F[1][j] eq JN];
	end while;
	return Out;
end intrinsic;




