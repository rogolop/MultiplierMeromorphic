AttachSpec("SingularitiesDim2/IntegralClosureDim2.spec");
Attach("MultiplierMeromorphic.m");
P<x, y> := LocalPolynomialRing(RationalField(), 2, "lglex");

// Settings
quitWhenFinished   := false;
printToFile        := false;
outFileNamePrefix  := "examples/";
outFileNameSufix   := ".txt";

MinJN := 0;
MaxJN := 5;

////////////////////

A := [];

// Bernstein-Sato...

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

// Append(~A, "(y^2-x^3)^5 + x^18, y^2+x^3");

// Append(~A, "(y^2-x^3)^4 + x^8*y^5, 1");
// Append(~A, "(y^2-x^3)^4 + x^8*y^5, y^2-x^3");
// Append(~A, "(y^2-x^3)^4 + x^8*y^5, y^2+x^3");




for s in A do
    
    f, g := eval s;
    f := P! f;
    g := P! g;
    // f, g := Explode(I);
    
    // Setup output
    sShort := &cat Split(s, "^* ");
    outFileName := outFileNamePrefix cat sShort cat outFileNameSufix;
    if printToFile then
        SetOutputFile(outFileName : Overwrite := true);
    end if;
    
    printf "f = %o\n", Split(s, ",")[1];
    printf "g =%o\n", Split(s, ",")[2];
    
    if printToFile then
        UnsetOutputFile();
    end if;
    
    // printf "\nIdeal I=<%o, %o>\n", Basis(I)[1], Basis(I)[2];
    
    // printf "\n### LogResolution\n";
    // LogResolution(I: Coefficients:=true);
    
    // printf "\n### LogResolutionMeromorphic\n";
    // LogResolutionMeromorphic(I: Coefficients:=true);
    
    // printf "\n### ProximityMatrix\n";
    // ProximityMatrix(f: Coefficients := true, ExtraPoint := true);
    
    // printf "\n### MultiplierIdeals I\n";
    // S := MultiplierIdealsMeromorphic_I(I : MaxJN:=1);
    // printf "%o\n", [t[2] : t in S];
    
    // printf "\n### MultiplierIdeals f\n";
    
    ////////////////////
    S := MultiplierIdealsMeromorphic(f, g : MinJN:=MinJN, MaxJN:=MaxJN, ComputeIdeals:=false);
    ////////////////////
    
    // printf "%o\n", [t[2] : t in S];
    // S;
    
    // f := (y^2-x^3)*x;
    // f := (x^2-y^3)*(y^2-x^3);
    // I := ideal<P | x^2*y^2, x^5, y^5, x*y^4, x^4*y>;
    
    // printf "\n### MultiplierIdeals\n";
    // printf "%o\n", [t[2] : t in S];
    // S;
    // S := MultiplierIdeals(I : MaxJN:=2);
    
    // S := MultiplierIdeals(f : MaxJN:=3);
    // S;
    
    if printToFile then
        SetOutputFile(outFileName : Overwrite := false);
    end if;
    
    ////////////////////
    // printf "\nJN:\n%o\n", [t[2] : t in S];
    ////////////////////
    printf "\nJN:\n";
    for m in [Floor(MinJN)..(Ceiling(MaxJN)-1)] do
        printf "%o\n", [t[2] : t in S | (m le t[2]) and (t[2] lt (m+1))];
    end for;
    den := LCM([Denominator(t[2]) : t in S]);
    // printf "\nNumerators:\n%o\n", [t[2]*den : t in S];
    printf "\nNumerators:\n";
    for m in [Floor(MinJN)..(Ceiling(MaxJN)-1)] do
        printf "%o\n", [t[2]*den : t in S | (m le t[2]) and (t[2] lt (m+1))];
    end for;
    
    printf "\nNumerators - floor:\n";
    for m in [Floor(MinJN)..(Ceiling(MaxJN)-1)] do
        printf "%o + %o\n", m, [(t[2]-m)*den : t in S | (m le t[2]) and (t[2] lt (m+1))];
    end for;
    
    if printToFile then
        UnsetOutputFile();
        printf "\nPrinted to: %o\n", outFileName;
    end if;
end for;

if quitWhenFinished then
	quit;
end if;

