AttachSpec("SingularitiesDim2/IntegralClosureDim2.spec");
Attach("MultiplierMeromorphic.m");
P<x, y> := LocalPolynomialRing(RationalField(), 2, "lglex");

// Settings
printToFile       := false;
quitWhenFinished  := true;
ComputeIdeals     := true;
printJN           := true;
printJN_fract     := true;
MinJN             := 0;
MaxJN             := 3;
outFileNamePrefix := "examples/";
outFileNameSufix  := ".txt";

////////////////////

A := [];

// Append(~A, "y^3+x^5, x");
// Append(~A, "y^3+x^5, y");
// Append(~A, "y^3+x^5, 1");

// Append(~A, "x^5*y - y^3, x^8 + 2*x^5*y");

// Append(~A, "x^2-y^3, y^2-x^3");
// Append(~A, "x^3, y^2");
// Append(~A, "y^2, x^3");
// Append(~A, "x^3+y^2, y^2");
// Append(~A, "x^3+x*y^2, y^2");
// Append(~A, "(x^2-y^3)*(y^2-x^3), 1");
// Append(~A, "(x^2-y^3)*(y^2-x^3), x^2*y^2");

// Append(~A, "(x^2-y^3)*(x^2+y^3), y^2-x^3");
// Append(~A, "(x^2-y^3)^2, y^2-x^3");
// Append(~A, "(x^2-y^3)*x, 1");
// Append(~A, "(x^2-y^3)^2, 1");
// Append(~A, "(x^2-y^3)*y, 1");
// Append(~A, "y^2-x^3, (x^2-y^3)*(x^2+y^3)");
// Append(~A, "y*(x^2-y^3)*(x^2+y^3), y^2-x^3");
// Append(~A, "(x^9-y^4-x^3*y^3)^2 - x^16*y, y^2-x^3");

// Append(~A, "(x^2-y^3)*(x^2+y^3), 1");
// Append(~A, "(x^2-y^3)^2, 1");
// Append(~A, "(x^2-y^3), 1");

// Append(~A, "(y^2-x^3)^5 + x^18, 1");
// Append(~A, "(y^2-x^3)^5 + x^18, y^2-x^3");
// Append(~A, "(y^2-x^3)^5 + x^18, (y^2-x^3)^2");
// Append(~A, "(y^2-x^3)^5 + x^18, (y^2-x^3)^3");
// Append(~A, "(y^2-x^3)^5 + x^18, (y^2-x^3)^4");
// Append(~A, "(y^2-x^3)^5 + x^18, (y^2-x^3)^5");
// Append(~A, "(y^2-x^3)^5 + x^18, (y^2-x^3)^6");
// Append(~A, "(y^2-x^3)^5 + x^18, y^2+x^3");
// Append(~A, "(y^2-x^3)^5 + x^18, x^2+y^3");

// Append(~A, "(y^2-x^3)^4 + x^8*y^5, y^2-x^3");
// Append(~A, "(y^2-x^3)^4 + x^8*y^5, (y^2-x^3)^1");
// Append(~A, "(y^2-x^3)^4 + x^8*y^5, (y^2-x^3)^2");
// Append(~A, "(y^2-x^3)^4 + x^8*y^5, (y^2-x^3)^3");
// Append(~A, "(y^2-x^3)^4 + x^8*y^5, (y^2-x^3)^4");
// Append(~A, "(y^2-x^3)^4 + x^8*y^5, y^2+x^3");
// Append(~A, "(y^2-x^3)^4 + x^8*y^5, y^8*x^8");
// Append(~A, "(y^2-x^3)^4 + x^8*y^5, y^8*x");
// Append(~A, "(y^2-x^3)^4 + x^8*y^5, y^8");
// Append(~A, "(y^2-x^3)^4 + x^8*y^5, 1");

// Append(~A, "x^3*(x+y)^5, y^2");
// Append(~A, "x^3*(x+y)^5*(y^2-x^3)^7, y^2");
// Append(~A, "x*(x+y)*(y^2-x^3), y^2");
// Append(~A, "x^6*(x^2-y^3)^4 , (x-y^3)^9");

// Append(~A, "(x^2+y^3)^4, (y^2-x^3)^2");
// Append(~A, "(x^2+y^3)^4, (y^2-x^3)*(y^2-2*x^3)");

// Append(~A, "((y^2-x^3)^5 + x^18)*(y^6-x^11), (y^2-x^3)^3*y^3");

// Append(~A, "((y^2-x^3)^4 + x^8*y^5)*((x^2-y^3)^5 + y^18), (x+y)^17");


for s in A do
    f, g := eval s;
    f := P! f;
    g := P! g;
    
    // Setup output
    sShort := &cat Split(s, "^* ");
    outFileName := outFileNamePrefix cat sShort cat outFileNameSufix;
    if printToFile then
        SetOutputFile(outFileName : Overwrite := true);
    end if;
    
    printf "f = %o\n", Split(s, ",")[1];
    printf "g =%o\n\n", Split(s, ",")[2];
    
    // if printToFile then UnsetOutputFile(); end if;
    S := MultiplierIdealsMeromorphic(f, g : MinJN:=MinJN, MaxJN:=MaxJN, ComputeIdeals:=ComputeIdeals);
    // if printToFile then SetOutputFile(outFileName : Overwrite := false); end if;
    
    printf "\n--------------------------------------------------\n";
    
    if ComputeIdeals then
        candidateNames := ["h","r","s","t","u","v","w","a","b","c","d","e"];
        candidateIndex := 1;
        factorNames := AssociativeArray();
        factorNames[P!1] := "1";
        factorNames[x] := "x";
        factorNames[y] := "y";
        factorNames[f] := "f";
        if g ne (P!1) then factorNames[g] := "g"; end if;
        
        printf "\n";
        JNPrintLength := - Max([5] cat [#Sprint(t[2]) : t in S]);
        for t in S do
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
    end if;
    
    if printJN then
        printf "\nJN:\n\n";
        for m in [Floor(MinJN)..(Ceiling(MaxJN)-1)] do
            printf "%o\n\n", [t[2] : t in S | (m le t[2]) and (t[2] lt (m+1))];
        end for;
    end if;
    
    // den := LCM([Denominator(t[2]) : t in S]);
    // // printf "\nNumerators:\n%o\n", [t[2]*den : t in S];
    // printf "\nNumerators:\n";
    // for m in [Floor(MinJN)..(Ceiling(MaxJN)-1)] do
    //     printf "%o\n", [t[2]*den : t in S | (m le t[2]) and (t[2] lt (m+1))];
    // end for;
    
    if printJN_fract then
        printf "JN fract:\n\n";
        for m in [Floor(MinJN)..(Ceiling(MaxJN)-1)] do
            printf "%o + %o\n\n", m, [(t[2]-m) : t in S | (m le t[2]) and (t[2] lt (m+1))];
        end for;
    end if;
    
    if printToFile then
        UnsetOutputFile();
        printf "\nPrinted to: %o\n", outFileName;
    end if;
end for;

if quitWhenFinished then
	quit;
end if;

