AttachSpec("SingularitiesDim2/IntegralClosureDim2.spec");
Attach("MultiplierMeromorphic.m");
P<x, y> := LocalPolynomialRing(RationalField(), 2, "lglex");

// Settings
printJN           := true;
printJN_fract     := true;

f := (y^2-x^3)^4 + x^8*y^5;
g := (y^2-x^3)^1;


Nf, Ng, N, Prox, commonFactor, Coeffs := LogResolutionMeromorphic(f, g);

S := MultiplierIdealsMeromorphic(f, g : MinJN:=0, MaxJN:=2, ComputeIdeals:=true);



candidateNames := ["h","r","s","t","u","v","w","a","b","c","d","e"];
candidateIndex := 1;
factorNames := AssociativeArray();
factorNames[P!1] := "1";
factorNames[x] := "x";
factorNames[y] := "y";
factorNames[f] := "f";
if g ne (P!1) then factorNames[g] := "g"; end if;

printf "\n";
JNPrintLength := - Max([5] cat [#Sprint(t[1]) : t in S]);
for t in S do
	JN := t[1];
	gen := t[2];
	printf "%*o | ", JNPrintLength, JN;
	powerOfF := gen[1][2];
	gen := gen[2];
	
	// printf "gen = %o\n", gen;
	gen := [Factorization(h) : h in gen];
	if gen eq [] then gen := [[<Parent(f)!1,1>]]; end if;
	
	if powerOfF gt 0 then
		printf "f";
		if powerOfF gt 1 then
			if powerOfF lt 10 then printf "^%o", powerOfF;
			else printf "^{%o}", powerOfF;
			end if;
		end if;
		printf " * ";
	end if;
	printf "( ";
	for i->generatorFactorization in gen do
		for factor in generatorFactorization do // h^m
			h := factor[1];
			m := factor[2];
			if not IsDefined(factorNames, h) then
				if candidateIndex gt # candidateNames then error "Ran out of factor names"; end if;
				factorNames[h] := candidateNames[candidateIndex];
				candidateIndex +:= 1;
			end if;
			printf "%o", factorNames[h];
			if m gt 1 then
				if m lt 10 then printf "^%o", m;
				else printf "^{%o}", m;
				end if;
			end if;
		end for;
		if i lt #gen then printf ", "; end if;
	end for;
	printf " )\n";
end for;

strs := [];
for h in Keys(factorNames) do
	if h notin {P!1, x, y, f, g} then
		strs := strs cat [ Sprintf("%o = %o\n", factorNames[h], h) ];
	end if;
end for;
if #strs gt 0 then
	Sort(~strs);
	printf "\nDefinitions:\n";
	for str in strs do printf "%o", str; end for;
end if;

printf "\n--------------------------------------------------\n";

MinJN := Min([t[1] : t in S]);
MaxJN := Max([t[1] : t in S]);

if printJN then
	printf "\nJN:\n\n";
	for m in [Floor(MinJN)..(Ceiling(MaxJN)-1)] do
		printf "%o\n\n", [t[1] : t in S | (m le t[1]) and (t[1] lt (m+1))];
	end for;
end if;

// den := LCM([Denominator(t[1]) : t in S]);
// // printf "\nNumerators:\n%o\n", [t[1]*den : t in S];
// printf "\nNumerators:\n";
// for m in [Floor(MinJN)..(Ceiling(MaxJN)-1)] do
//     printf "%o\n", [t[1]*den : t in S | (m le t[1]) and (t[1] lt (m+1))];
// end for;

if printJN_fract then
	printf "JN fract:\n\n";
	for m in [Floor(MinJN)..(Ceiling(MaxJN)-1)] do
		printf "%o + %o\n\n", m, [(t[1]-m) : t in S | (m le t[1]) and (t[1] lt (m+1))];
	end for;
end if;
