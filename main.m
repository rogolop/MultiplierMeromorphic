AttachSpec("SingularitiesDim2/IntegralClosureDim2.spec");
Attach("MultiplierMeromorphic.m");
P<x, y> := LocalPolynomialRing(RationalField(), 2, "lglex");

// Input
f := (y^2-x^3)^4 + x^8*y^5;
g := (y^2-x^3)^1;

// Log resolution
Nf, Ng, N, Prox, commonFactor, Coeffs := LogResolutionMeromorphic(f, g);
printf "Nf = %o\n\n", Nf;
printf "Ng = %o\n\n", Ng;
printf "N  = %o\n\n", N;

// Jumping numbers and multiplier ideals
S := MultiplierIdealsMeromorphic(f, g : MinJN:=9/10, MaxJN:=5/3, ComputeIdeals:=true);
// Print output data nicely
for JNData in S do
	JN := JNData[1];
	multiplierIdealData := JNData[2];
	fExponent := multiplierIdealData[1];
	generators := multiplierIdealData[2];
	printf "JN = %o\n", JN;
	printf "MI = f^%o * ( ", fExponent;
	for i->gen in generators do
		printf "%o", gen;
		printf (i lt #generators) select ", " else " )\n\n";
	end for;
end for;

