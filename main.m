AttachSpec("SingularitiesDim2/IntegralClosureDim2.spec");
Attach("MultiplierMeromorphic.m");
P<x, y> := LocalPolynomialRing(RationalField(), 2, "lglex");

// Settings
ComputeIdeals     := true;
printJN           := true;
printJN_fract     := true;

f := (y^2-x^3)^4 + x^8*y^5;
g := (y^2-x^3)^1;



Prox, N, _, _, Excess_f, Nf := LogResolutionMeromorphic(f, g);
Js := MultiplierIdealsMeromorphic(f, g : MinJN:=0, MaxJN:=4, ComputeIdeals:=true);



candidateNames := ["h","r","s","t","u","v","w","a","b","c","d","e"];
candidateIndex := 1;
factorNames := AssociativeArray();
factorNames[P!1] := "1";
factorNames[x] := "x";
factorNames[y] := "y";
factorNames[f] := "f";
if g ne (P!1) then factorNames[g] := "g"; end if;

printf "\n";
JNPrintLength := - Max([5] cat [#Sprint(t[2]) : t in Js]);
for t in Js do
    gen := t[1];
    JN := t[2];
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

if printJN then
    printf "\nJN:\n\n";
    for m in [Floor(MinJN)..(Ceiling(MaxJN)-1)] do
        printf "%o\n\n", [t[2] : t in Js | (m le t[2]) and (t[2] lt (m+1))];
    end for;
end if;

// den := LCM([Denominator(t[2]) : t in Js]);
// // printf "\nNumerators:\n%o\n", [t[2]*den : t in Js];
// printf "\nNumerators:\n";
// for m in [Floor(MinJN)..(Ceiling(MaxJN)-1)] do
//     printf "%o\n", [t[2]*den : t in Js | (m le t[2]) and (t[2] lt (m+1))];
// end for;

if printJN_fract then
    printf "JN fract:\n\n";
    for m in [Floor(MinJN)..(Ceiling(MaxJN)-1)] do
        printf "%o + %o\n\n", m, [(t[2]-m) : t in Js | (m le t[2]) and (t[2] lt (m+1))];
    end for;
end if;
